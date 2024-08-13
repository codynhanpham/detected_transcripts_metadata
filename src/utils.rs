use std::{
    fs,
    io::{self, BufRead, BufReader},
    path::Path,
    sync::{Arc, RwLock},
};

use memchr;
use memmap2::Mmap;
use rayon::prelude::*;



/// Check if CSV file is valid
fn validate_csv_file(file_path: &Path) -> io::Result<()> {
    let file = fs::File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Check if the first line is a header
    let first_line = lines.next().unwrap()?;
    let headers: Vec<&str> = first_line.split(',').collect();
    if headers.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "CSV file has no headers"));
    }

    // Check if there is more than 1 row
    if lines.next().is_none() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "CSV file has no data"));
    }

    Ok(())
}

/// Check if CSV file is valid detected_transcripts.csv file
pub fn validate_detected_transcripts_csv_file(file_path: &Path, col_of_interest: &Vec<&str>) -> io::Result<()> {
    // Validate CSV file
    validate_csv_file(file_path)?;

    // Read the first line of the CSV file for headers
    let headers = get_headers(file_path)?;

    // Check if the headers of interest are in the CSV file
    match get_col_indices(&headers, col_of_interest) {
        Ok(_) => Ok(()),
        Err(e) => Err(io::Error::new(io::ErrorKind::InvalidData, e)),
    }
}

/// Split line by comma and return a Vec of columns data
pub fn get_col_data(line: &str) -> Vec<String> {
    let mut cols = Vec::new();
    let mut start = 0;
    let mut end = 0;

    while end < line.len() {
        match memchr::memchr(b',', &line.as_bytes()[end..]) {
            Some(pos) => {
                end += pos;
                cols.push(line[start..end].to_string());
                start = end + 1;
                end += 1;
            }
            None => {
                cols.push(line[start..].to_string());
                break;
            }
        }
    }

    cols
}

/// Print column data as (index, value)
pub fn _print_col_data(col_data: &Vec<String>) {
    for (i, data) in col_data.iter().enumerate() {
        println!("{}:  {}", i, data);
    }
}


/// Get the headers of a CSV file
pub fn get_headers(file_path: &Path) -> io::Result<Vec<String>> {
    let file = fs::File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Check if the first line is a header
    let first_line = lines.next().unwrap()?;
    let headers: Vec<String> = first_line.split(',').map(|s| s.to_string()).collect();
    Ok(headers)
}

/// Print header of a CSV file
pub fn _print_header(file_path: &Path) -> io::Result<()> {
    let headers = get_headers(file_path)?;
    _print_col_data(&headers);
    Ok(())
}


/// Return the index of cols from a given Vec of column values
/// If any col value is not found, error is returned
pub fn get_col_indices(headers: &Vec<String>, cols: &Vec<&str>) -> Result<Vec<usize>, String> {
    let mut indices = Vec::new();
    for col in cols {
        match headers.iter().position(|h| h == col) {
            Some(i) => indices.push(i),
            None => return Err(format!("Column '{}' not found", col)),
        }
    }
    Ok(indices)
}


/// Get a column of interest from a CSV file in parallel
/// Return as a Vec of Strings
pub fn _get_single_col_from_file(file_path: &Path, col_index: usize, off_set_rows: usize, num_threads: usize, process_chunk_size: usize, quite_mode: bool) -> io::Result<Vec<String>> {
    // Validate CSV file
    validate_csv_file(file_path)?;

    // Read the first line of the CSV file for headers
    let headers = get_headers(file_path)?;

    // Make sure the column index is within the range of the headers
    if col_index >= headers.len() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Column index out of range"));
    }


    // Read the CSV file in parallel
    let file = fs::File::open(file_path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let mmap_len = mmap.len();

    let mut start_byte = 0;
    let mut end_byte: usize;
    let mut start_bytes = Vec::with_capacity(num_threads);
    let mut end_bytes = Vec::with_capacity(num_threads);

    // Off set rows from the start of the file specified with off_set_rows
    let mut off_set_rows = off_set_rows;
    while off_set_rows > 0 {
        match memchr::memchr(b'\n', &mmap[start_byte..]) {
            Some(pos) => {
                start_byte += pos + 1;
                off_set_rows -= 1;
            }
            None => {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "Off set rows out of range"));
            }
        }
    }

    // Split the file into num_threads chunks
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
    if quite_mode == false {
        println!("Number of threads: {}", num_threads);
        println!("File length: {} (bytes)", mmap_len);
        for (i, (start, end)) in start_bytes.iter().zip(end_bytes.iter()).enumerate() {
            let percentage = ((end - start) as f64 / (mmap_len - start_bytes[0]) as f64 * 100.0 * 100.0).round() / 100.0;
            println!("Thread {}: [{}, {}] (~{}%)", i, start, end, percentage);
        }
        println!();
    }

    drop(mmap);

    
    // Process the file in parallel
    let result: Arc<RwLock<Vec<Vec<String>>>> = Arc::new(RwLock::new(Vec::with_capacity(num_threads)));

    let progress_counter = Arc::new(RwLock::new(0));
    let thread_complete = Arc::new(RwLock::new(0));

    start_bytes.par_iter().zip(end_bytes.par_iter()).enumerate().for_each(|(_, (start, end))| {
        let mut local_start = *start;
        let mut local_result = Vec::with_capacity(num_threads);

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

            let line_len = lines.len();

            // Process the lines
            for line in lines {
                let cols = get_col_data(&String::from_utf8_lossy(line));
                // Scan to get the column of interest
                if cols.len() > col_index {
                    local_result.push(cols[col_index].clone());
                }
            }

            // Update the progress counter
            if !quite_mode {
                *progress_counter.write().unwrap() += line_len;
                println!("Processed {} lines", *progress_counter.read().unwrap());
            }

            drop(mmap);
        }

        // Append the local result to the global result
        result.write().unwrap().push(local_result);

        if !quite_mode {
            *thread_complete.write().unwrap() += 1;
            println!("[{}/{}] threads completed the task", *thread_complete.read().unwrap(), num_threads);
        }
    });

    // Collect the results
    let mut final_result = Vec::new();
    for res in result.read().unwrap().iter() {
        final_result.extend(res.clone());
    }

    Ok(final_result)
}