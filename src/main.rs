use std::{
    fs::File, io::{self}, path::Path, sync::{Arc, RwLock}
};

use clap::Parser;
use memchr;
use memmap2::Mmap;
use rayon::prelude::*;

mod utils;
mod transcript_metadata;


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

    /// Run quitely without printing to stdout (except for errors)
    #[arg(short, long, default_value = "false")]
    quiet: bool,

    /// Number of lines to process in each chunk
    #[arg(short, long, default_value = "500000")]
    chunk_size: usize,

    /// Number of processes to run in parallel
    /// Default is the number of logical CPUs in the system
    #[arg(short = 'p', long = "processes", verbatim_doc_comment)]
    num_threads: Option<usize>,

    /// Resolution for the transcript mask area calculation
    /// Default is 60.0 micrometers
    #[arg(short, long, default_value = "60.0")]
    resolution: f64,

    /// Minimum number of transcripts required per resolution area to be considered for the mask area calculation
    /// Default is 1 transcript
    #[arg(short, long, default_value = "1")]
    density_threshold: usize,
}




fn main() -> io::Result<()> {
    // ----------------- Read input file ----------------- // 

    // Set headers of interest. These columns are required and will be used to validate the CSV file
    let col_of_interest = vec!["global_x", "global_y", "global_z", "cell_id", "gene"];

    let args: Args = Args::parse();

    let file_path = Path::new(&args.input);

    // Validate CSV file is formatted for detected_transcripts.csv
    match utils::validate_detected_transcripts_csv_file(file_path, &col_of_interest) {
        Ok(_) => (),
        Err(e) => {
            eprintln!("Error: {}\nFile is not a valid detected_transcripts.csv or is not formatted properly.\nThe following column headers are expected: [ {} ]", e, col_of_interest.join(", "));
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



    // ----------------- Process transcripts data ----------------- //

    let col_of_interest_indices = utils::get_col_indices(&utils::get_headers(file_path).unwrap(), &col_of_interest).unwrap();

    // Chunk settings
    let num_threads: usize = parallel_processes;
    rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global().unwrap();
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
    let thread_chunk_size = (mmap_len - start_byte) / num_threads as usize;
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
            let percentage = ((end - start) as f64 / (mmap_len - start_bytes[0]) as f64 * 100.0 * 100.0).round() / 100.0;
            println!("Thread {}: [{}, {}] (~{}%)", i, start, end, percentage);
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
        let metadata = transcript_metadata::TranscriptMetadata::new();

        while local_start < *end {
            let mmap = unsafe { Mmap::map(&file).unwrap() };
            let mut slice = &mmap[local_start..*end];
            let mut lines = Vec::new();
            let mut line_count = 0;

            // Dynamically scale the process_chunk_size by the number of threads that have completed processing (at the end of the file)
            let process_chunk_size = process_chunk_size * (1 + *thread_complete.read().unwrap() / num_threads);

            while line_count < process_chunk_size && !slice.is_empty() {
                match memchr::memchr(b'\n', slice) {
                    Some(pos) => {
                        lines.push(&slice[..pos]);
                        slice = &slice[pos + 1..];
                        local_start += pos + 1;
                        line_count += 1;
                    },
                    None => {
                        if !slice.is_empty() {
                            lines.push(slice);
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

    if !args.quiet {
        println!("\nMerging result from {} threads...", num_threads);
    }

    let metadata = Arc::try_unwrap(metadata_vec).unwrap().into_inner().unwrap();
    let processed_metadata = transcript_metadata::merge_metadata(metadata);

    if !args.quiet {
        println!("Metadata compilation of {} transcripts completed.", processed_metadata.total_transcripts.read().unwrap());
    }


    // ----------------- Calculate area covered by the transcripts ----------------

    if !args.quiet {
        println!("\nCalculating area masked by transcripts...");
    }
    processed_metadata.calculate_transcript_mask_area(args.resolution, args.density_threshold);
    if !args.quiet {
        println!("Area covered by transcripts: {} sq. micron", processed_metadata.transcripts_mask_area.read().unwrap());
    }



    // ----------------- Output results ----------------

    processed_metadata.export_json(
        metadata_file_path.unwrap().to_str().unwrap(),
        &print_json_output
    )
}