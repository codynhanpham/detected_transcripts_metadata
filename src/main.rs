use clap::Parser;
use once_cell::sync::Lazy;
use std::time::Instant;
use std::{fs::File, io, path::Path};

mod line_counter;
mod transcript_metadata;
mod utils;

static DEFAULT_THREADS: Lazy<usize> = once_cell::sync::Lazy::new(|| std::thread::available_parallelism().unwrap().into());
static DEFAULT_CHUNK_SIZE: Lazy<usize> = once_cell::sync::Lazy::new(|| 512);

#[derive(Parser, Debug)]
#[command(
    name = "detected_transcript_metadata",
    version = "2.0.0",
    about = "Extract metadata from the detected_transcripts.csv file generated on the Vizgen MERFISH platform."
)]
struct Args {
    /// Path to "detected_transcripts.csv" file to process
    #[arg(short = 'i', long = "input", required = true)]
    input: String,

    /// Count the number of lines in the input file and exit
    #[arg(
        short = 'l',
        long,
        default_value = "false",
        alias = "line",
        verbatim_doc_comment
    )]
    lines: bool,

    /// Path to save the metadata JSON file
    /// If not provided, the output will be saved in the same directory as the input file
    /// Use special value "-" to print the metadata JSON to console without saving to a file
    #[arg(short = 'o', long = "output", verbatim_doc_comment)]
    output: Option<String>,

    /// Chunk size (in KiB) for each processing thread
    #[arg(short = 'c', long, default_value_t = *DEFAULT_CHUNK_SIZE)]
    chunk_size: usize,

    /// Number of processing threads
    #[arg(short = 'p', long, alias = "thread", alias = "threads", alias = "processes", default_value_t = *DEFAULT_THREADS)]
    process: usize,

    /// Run quitely without printing to stdout (except for errors)
    #[arg(short, long, default_value = "false")]
    quiet: bool,
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    // ----------------- Validate Input File ----------------- //

    let file_path = Path::new(&args.input);

    // Set headers of interest. These columns are required and will be used to validate the CSV file
    let col_of_interest = vec!["global_x", "global_y", "global_z", "gene", "cell_id"];

    match utils::validate_detected_transcripts_csv_file(file_path, &col_of_interest) {
        Ok(_) => (),
        Err(e) => {
            eprintln!("Error: {}\nFile is not a valid detected_transcripts.csv or is not formatted properly.\nThe following column headers are expected: [ {} ]", e, col_of_interest.join(", "));
            return Ok(());
        }
    }
    if !args.quiet {
        println!("Processing input file: {:?}", file_path);
    }

    let file = File::open(&args.input).expect("Failed to open file");


    // ----------------- Count Lines Only and Exit ----------------- //
    if args.lines {
        let start = Instant::now();
        let total_lines = line_counter::count_total_lines(file, args.process, args.chunk_size, args.quiet);
        println!("Line Count: {}", total_lines);
        println!("Without header and assuming each line is a transcript");
        println!("    --> the total number of transcripts is {}", total_lines - 1);

        let duration = start.elapsed();
        println!("\nTime taken: {:?}", duration);
        return Ok(());
    }


    // ----------------- Handle What to Do with Output ----------------- //
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


    // ----------------- Metadata Processing ----------------- //
    let start = Instant::now();

    let col_of_interest_indices = utils::get_col_indices(&utils::get_headers(file_path).unwrap(), &col_of_interest).unwrap();

    let metadata = transcript_metadata::process_metadata_file(file, col_of_interest_indices, args.process, args.chunk_size, args.quiet);

    if !args.quiet {
        let duration = start.elapsed();
        println!("\nTime taken: {:?}", duration);
    }


    // ----------------- Export Metadata to JSON ----------------- //

    metadata.export_json(
        metadata_file_path.unwrap_or("".into()).to_str().unwrap(),
        &print_json_output
    )
}