use std::{
    fs::File,
    io::{self},
    path::Path,
    sync::{Arc, RwLock},
};

use clap::Parser;
use dashmap::DashMap;
use memchr;
use memmap::MmapOptions;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json;

mod utils;


#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None, arg_required_else_help = true)]
struct Args {
    /// Path to "detected_transcripts.csv" file to process
    #[arg(short = 'i', long = "input", required = true)]
    input: String,

    /// Path to save the metadata JSON file.
    /// If not provided, the output will be saved in the same directory as the input file
    /// Use special value "-" to print the metadata JSON to console without saving to a file
    #[arg(short = 'o', long = "output", verbatim_doc_comment)]
    output: Option<String>,

    /// Number of lines to process in each chunk
    #[arg(short, long, default_value = "500000")]
    chunk_size: usize,

    /// Run quitely without printing to stdout (except for errors)
    #[arg(short, long, default_value = "false")]
    quiet: bool,
}



#[derive(Debug, Serialize, Deserialize)]
struct TranscriptMetadata {
    total_transcripts: RwLock<usize>,
    transcripts_within_cells: RwLock<usize>,
    layers_transcripts_count: DashMap<i32, usize>,
    layers_transcripts_within_cells_count: DashMap<i32, usize>,
    genes_frequency: DashMap<String, usize>,
    genes_frequency_per_layer: DashMap<i32, DashMap<String, usize>>,
}

impl TranscriptMetadata {
    fn new() -> Self {
        Self {
            total_transcripts: RwLock::new(0),
            transcripts_within_cells: RwLock::new(0),
            layers_transcripts_count: DashMap::new(),
            layers_transcripts_within_cells_count: DashMap::new(),
            genes_frequency: DashMap::new(),
            genes_frequency_per_layer: DashMap::new(),
        }
    }

    fn add_data_chunk(&self, data_chunk: &[Vec<String>]) {
        data_chunk.par_iter().for_each(|row| {
            let global_z = match row[0].parse::<f32>() {
                Ok(z) => z.round() as i32,
                Err(_) => {
                    eprintln!("Error: Could not parse global_z value to int32: {}", row[0]);
                    return;
                },
            };
            let cell_id = &row[1];
            let gene = &row[2];

            // Increment total transcripts
            *self.total_transcripts.write().unwrap() += 1;

            // Increment transcripts within cells
            if !cell_id.is_empty() && cell_id != "NA" && cell_id != "N/A" && cell_id != "-1" {
                *self.transcripts_within_cells.write().unwrap() += 1;
                self.layers_transcripts_within_cells_count
                    .entry(global_z)
                    .and_modify(|e| *e += 1)
                    .or_insert(1);
            }

            // Increment layer transcript count
            self.layers_transcripts_count
                .entry(global_z)
                .and_modify(|e| *e += 1)
                .or_insert(1);

            // Increment gene frequency
            self.genes_frequency
                .entry(gene.clone())
                .and_modify(|e| *e += 1)
                .or_insert(1);

            // Increment gene frequency per layer
            self.genes_frequency_per_layer
                .entry(global_z)
                .or_insert(DashMap::new())
                .entry(gene.clone())
                .and_modify(|e| *e += 1)
                .or_insert(1);
        });
    }
}



fn main() -> io::Result<()> {
    // ----------------- Read input file ----------------- // 

    // Set headers of interest
    let col_of_interest = vec!["global_z", "cell_id", "gene"];

    let args: Args = Args::parse();

    let file_path = Path::new(&args.input);

    // Validate CSV file headers is formatted for detected_transcripts.csv
    match utils::validate_detected_transcripts_csv_file(file_path, &col_of_interest) {
        Ok(_) => (),
        Err(e) => {
            eprintln!("Error: {}\nFile is not a valid detected_transcripts.csv or is not formatted properly", e);
            return Ok(());
        },
    }
    if !args.quiet {
        println!("Using input file: {:?}", file_path);
    }


    // If the -o value is not a directory, then it is a file, use the provided filename for the metadata file
    // Otherwise, it is a directory, so make sure the directory exists and use the input filename + "_metadata.csv" for the metadata file
    let metadata_file_path = if let Some(ref output) = args.output {
        if output == "-" {
            None
        } else {
            let output_path = Path::new(output);
            if output_path.is_dir() {
                let file_name = file_path.file_stem().unwrap().to_str().unwrap();
                let metadata_file_name = format!("{}_metadata.json", file_name);
                Some(output_path.join(metadata_file_name))
            } else {
                Some(output_path.to_path_buf())
            }
        }
    } else {
        let file_name = file_path.file_stem().unwrap().to_str().unwrap();
        let metadata_file_name = format!("{}_metadata.json", file_name);
        Some(file_path.with_file_name(metadata_file_name))
    };
    
    let print_json_output = args.output.as_ref().map_or(false, |output| output == "-");


    // ----------------- Process data ----------------- //

    let col_of_interest_indices = utils::get_col_indices(&utils::get_headers(file_path).unwrap(), &col_of_interest).unwrap();

    let metadata = Arc::new(TranscriptMetadata::new());

    // Read file in chunks
    let chunk_size = args.chunk_size;
    let num_threads: usize = std::thread::available_parallelism().unwrap().into();



    let file = File::open(file_path)?;
    let mmap = unsafe { MmapOptions::new().map(&file)? };
    let mmap_len = mmap.len();


    // Without loading the entire file into memory, read the file in chunks:
    // Discover the lines upto set of chunk_size * num_threads
    // Stop the seeking, now divide the mmap into num_threads chunks of chunk_size of the current set of lines
    // Process each chunk in parallel, only use the columns of interest data for processing
    // After processing the set, seek to the next set of lines and repeat the process
    let mut start = 0;

    // Skip the first line (headers)
    if let Some(first_newline) = memchr::memchr(b'\n', &mmap) {
        start = first_newline + 1;
    }

    let mut total_lines_processed = 0;

    while start < mmap_len {
        let mmap = unsafe { MmapOptions::new().map(&file)? };

        // Discover the lines upto set of chunk_size * num_threads
        let mut lines = Vec::new();
        let mut line_count = 0;

        while line_count < chunk_size * num_threads {
            match memchr::memchr(b'\n', &mmap[start..]) {
                Some(pos) => {
                    lines.push(&mmap[start..start + pos]);
                    start += pos + 1;
                    line_count += 1;
                },
                None => {
                    if start < mmap.len() {
                        lines.push(&mmap[start..]);
                    }
                    start = mmap.len();
                    break;
                },
            }
        }

        if lines.is_empty() {
            eprintln!("Encountered empty lines, assuming end of file. Stopped scanning lines in file.");
            break;
        }

        // Divide the mmap into num_threads chunks of chunk_size of the current set of lines
        let mut chunks: Vec<Vec<&[u8]>> = vec![Vec::new(); num_threads];
        for (i, line) in lines.iter().enumerate() {
            chunks[i % num_threads].push(line);
        }

        // Process each chunk in parallel!!!
        chunks.into_par_iter().for_each(|chunk| {
            let mut data_chunk = Vec::new();
            for line in chunk {
                let line = String::from_utf8_lossy(line).trim().to_string();
                let row = utils::get_col_data(&line);
                let data: Result<Vec<String>, &str> = col_of_interest_indices
                    .iter()
                    .map(|&col| row.get(col).cloned().ok_or("Index out of bounds"))
                    .collect();
                match data {
                    Ok(data) => {
                        if data.is_empty() {
                            eprintln!("Error: Could not parse data from line: \"{}\"", line);
                            continue;
                        }
                        data_chunk.push(data);
                    }
                    Err(e) => {
                        eprintln!("Error: {} for line: \"{}\"", e, line);
                        continue;
                    }
                }
            }
            metadata.add_data_chunk(&data_chunk);
            data_chunk.clear();
        });
        total_lines_processed += lines.len();
        lines.clear();
        if !args.quiet {
            println!("Processed {} lines", total_lines_processed);
        }
    }


    // ----------------- Output results ----------------
    
    let processed_metadata = Arc::try_unwrap(metadata).unwrap();


    // If the output value is the special value "-", then do not save the metadata JSON to a file
    if print_json_output {
        let metadata_json = serde_json::to_string_pretty(&processed_metadata).unwrap();
        println!("{}", metadata_json);
        return Ok(());
    }


    // Otherwise, save metadata as JSON file
    let metadata_file_path = metadata_file_path.unwrap();
    let metadata_file = File::create(metadata_file_path).unwrap();
    serde_json::to_writer(metadata_file, &processed_metadata).unwrap();

    Ok(())
}