use std::{
    fs::File,
    io::{self},
    path::Path,
    sync::{Arc, RwLock},
};

use clap::Parser;
use dashmap::DashMap;
use memchr;
use memmap2::Mmap;
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

    /// Path to save the metadata JSON file
    /// If not provided, the output will be saved in the same directory as the input file
    /// Use special value "-" to print the metadata JSON to console without saving to a file
    #[arg(short = 'o', long = "output", verbatim_doc_comment)]
    output: Option<String>,

    /// Number of lines to process in each chunk
    #[arg(short, long, default_value = "500000")]
    chunk_size: usize,

    /// Number of processes to run in parallel
    /// Default is the number of logical CPUs in the system
    #[arg(short = 'p', long = "processes", verbatim_doc_comment)]
    num_threads: Option<usize>,

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

            // Increment fields with new chunk data
            *self.total_transcripts.write().unwrap() += 1;

            if !cell_id.is_empty() && cell_id != "NA" && cell_id != "N/A" && cell_id != "-1" {
                *self.transcripts_within_cells.write().unwrap() += 1;
                self.layers_transcripts_within_cells_count
                    .entry(global_z)
                    .and_modify(|e| *e += 1)
                    .or_insert(1);
            }

            self.layers_transcripts_count
                .entry(global_z)
                .and_modify(|e| *e += 1)
                .or_insert(1);

            self.genes_frequency
                .entry(gene.clone())
                .and_modify(|e| *e += 1)
                .or_insert(1);

            self.genes_frequency_per_layer
                .entry(global_z)
                .or_insert(DashMap::new())
                .entry(gene.clone())
                .and_modify(|e| *e += 1)
                .or_insert(1);
        });
    }
}



/// Given a Vec<TranscriptMetadata>, merge all the metadata into a single TranscriptMetadata
fn merge_metadata(metadata: Vec<TranscriptMetadata>) -> TranscriptMetadata {
    let merged_metadata = TranscriptMetadata::new();

    metadata.into_par_iter().for_each(|data| {
        *merged_metadata.total_transcripts.write().unwrap() += *data.total_transcripts.read().unwrap();
        *merged_metadata.transcripts_within_cells.write().unwrap() += *data.transcripts_within_cells.read().unwrap();

        // Extend + Combine the DashMaps:
        
        for entry in data.layers_transcripts_count.iter() {
            let (layer, count) = entry.pair();
            merged_metadata.layers_transcripts_count
                .entry(*layer)
                .and_modify(|e| *e += *count)
                .or_insert(*count);
        }

        for entry in data.layers_transcripts_within_cells_count.iter() {
            let (layer, within_cells_count) = entry.pair();
            merged_metadata.layers_transcripts_within_cells_count
                .entry(*layer)
                .and_modify(|e| *e += *within_cells_count)
                .or_insert(*within_cells_count);
        }

        for entry in data.genes_frequency.iter() {
            let (gene, count) = entry.pair();
            merged_metadata.genes_frequency
                .entry(gene.clone())
                .and_modify(|e| *e += *count)
                .or_insert(*count);
        }

        for entry in data.genes_frequency_per_layer.iter() {
            let (layer, genes_frequency) = entry.pair();
            for gene_entry in genes_frequency.iter() {
                let (gene, count) = gene_entry.pair();
                merged_metadata.genes_frequency_per_layer
                    .entry(*layer)
                    .or_insert(DashMap::new())
                    .entry(gene.clone())
                    .and_modify(|e| *e += *count)
                    .or_insert(*count);
            }
        }
    });

    merged_metadata
}



fn main() -> io::Result<()> {
    // ----------------- Read input file ----------------- // 

    // Set headers of interest. These columns are required and will be used to validate the CSV file
    let col_of_interest = vec!["global_z", "cell_id", "gene"];

    let args: Args = Args::parse();

    let file_path = Path::new(&args.input);

    // Validate CSV file is formatted for detected_transcripts.csv
    match utils::validate_detected_transcripts_csv_file(file_path, &col_of_interest) {
        Ok(_) => (),
        Err(e) => {
            eprintln!("Error: {}\nFile is not a valid detected_transcripts.csv or is not formatted properly", e);
            return Ok(());
        },
    }
    if !args.quiet {
        println!("Processing input file: {:?}", file_path);
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

    // Set number of parallel processes
    let available_cpus = std::thread::available_parallelism().unwrap().into();
    let mut parallel_processes = args.num_threads.unwrap_or_else(|| available_cpus);
    if parallel_processes == 0 {
        parallel_processes = available_cpus;
    }



    // ----------------- Process data ----------------- //

    let col_of_interest_indices = utils::get_col_indices(&utils::get_headers(file_path).unwrap(), &col_of_interest).unwrap();

    // Chunk settings
    let num_threads: usize = parallel_processes;
    let process_chunk_size = args.chunk_size;


    // Strategy: Divide the file into num_threads chunks using the length of the file
    // At the estimated start byte of each chunk, use memchr to find next newline character
    // These will be the actual start byte of each chunk for each thread -> Save the chunk start and end bytes
    // Then, each thread will then discover/scan the lines in their chunk upto the chunk_size at a time
    // Process the data, repeating the process until the end of the chunk

    let file = File::open(file_path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let mmap_len = mmap.len();

    let mut start_byte = 0;
    let mut end_byte: usize;
    let mut start_bytes = Vec::with_capacity(num_threads);
    let mut end_bytes = Vec::with_capacity(num_threads);

    // Skip the first line (headers)
    if let Some(first_newline) = memchr::memchr(b'\n', &mmap) {
        start_byte = first_newline + 1;
    }

    // Split the file into num_threads chunks
    let thread_chunk_size = mmap_len / num_threads as usize;
    for _ in 0..num_threads {
        start_bytes.push(start_byte);

        end_byte = start_byte + thread_chunk_size;
        if end_byte > mmap_len {
            end_byte = mmap_len;
        }
        if end_byte == mmap_len {
            end_bytes.push(end_byte);
            break;
        }
        // Refine the next start byte to align with the newline character + 1
        match memchr::memchr(b'\n', &mmap[end_byte..]) {
            Some(pos) => {
                end_byte += pos + 1;
                end_bytes.push(end_byte);
            },
            None => {
                eprintln!("Error: Could not find newline character in chunk. Stopped scanning lines in file.");
                break;
            },
        }
        start_byte = end_byte + 1;
    }


    // Print chunk settings
    if !args.quiet {
        println!("Number of threads: {}", num_threads);
        println!("File length: {} (bytes)", mmap_len);
        for (i, (start, end)) in start_bytes.iter().zip(end_bytes.iter()).enumerate() {
            println!("Thread {}: [{}, {}]", i, start, end);
        }
        println!();
    }

    drop(mmap);


    // Vec of Metadata for each thread to be merged at the end
    // Instead of adding to a shared metadata, this reduces lock contention
    let metadata_vec = Arc::new(RwLock::new(Vec::with_capacity(num_threads)));


    // Process chunks in parallel
    let progress_counter = Arc::new(RwLock::new(0));
    let thread_complete = Arc::new(RwLock::new(0));

    start_bytes.par_iter().zip(end_bytes.par_iter()).for_each(|(start, end)| {
        let mut local_start = *start;
        let metadata = TranscriptMetadata::new();

        while local_start < *end {
            let mmap = unsafe { Mmap::map(&file).unwrap() };

            let mut lines = Vec::new();
            let mut line_count = 0;

            while line_count < process_chunk_size && local_start < *end {
                match memchr::memchr(b'\n', &mmap[local_start..]) {
                    Some(pos) => {
                        lines.push(&mmap[local_start..local_start + pos]);
                        local_start += pos + 1;
                        line_count += 1;
                    },
                    None => {
                        if local_start < *end {
                            lines.push(&mmap[local_start..]);
                        }
                        local_start = *end;
                        break;
                    },
                }
            }

            if lines.is_empty() {
                eprintln!("Encountered empty lines, assuming end of file. Stopped scanning lines in file.");
                break;
            }

            // Parse the data chunk
            // Only collect col_of_interest data for efficiency
            let mut data_chunk = Vec::with_capacity(lines.len());
            for line in lines {
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

            // Increment progress counter with lines processed
            if !args.quiet {
                *progress_counter.write().unwrap() += data_chunk.len();
                println!("Processed {} lines", *progress_counter.read().unwrap());
            }

            // Commit the data chunk to the metadata
            metadata.add_data_chunk(&data_chunk);

            // Cleanup
            data_chunk.clear();
        }

        // Push the metadata to the shared metadata vector
        metadata_vec.write().unwrap().push(metadata);

        if !args.quiet {
            *thread_complete.write().unwrap() += 1;
            println!("[{}/{}] threads completed the task", *thread_complete.read().unwrap(), num_threads);
        }
    });



    // ----------------- Merge metadata ----------------

    let metadata = Arc::try_unwrap(metadata_vec).unwrap().into_inner().unwrap();
    let processed_metadata = merge_metadata(metadata);

    if !args.quiet {
        println!("\nMetadata compilation of {} transcripts completed.", processed_metadata.total_transcripts.read().unwrap());
    }



    // ----------------- Output results ----------------
    
    // If the output value is the special value "-", print to console instead of saving to file
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