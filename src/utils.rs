use std::{
    fs,
    io::{self, BufRead, BufReader},
    path::Path,
};

/// Check if CSV file is valid
fn validate_csv_file(file_path: &Path) -> io::Result<()> {
    let file = fs::File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Check if the first line is a header
    let first_line = lines.next().unwrap()?;
    let headers: Vec<&str> = first_line.split(',').collect();
    if headers.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "CSV file has no headers",
        ));
    }

    // Check if there is more than 1 row
    if lines.next().is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "CSV file has no data",
        ));
    }

    Ok(())
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


/// Print column data as (index, value)
pub fn _print_col_data(col_data: &Vec<String>) {
    for (i, data) in col_data.iter().enumerate() {
        println!("{}:  {}", i, data);
    }
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

/// Split line by comma and return a Vec of columns data
pub fn get_col_data(line: &&[u8]) -> Vec<Vec<u8>> {
    // let line = String::from_utf8_lossy(line).trim().to_string();
    
    let mut cols = Vec::new();
    let mut start = 0;
    let mut end = 0;

    while end < line.len() {
        match memchr::memchr(b',', &line[end..]) {
            Some(pos) => {
                end += pos;
                cols.push(line[start..end].to_vec());
                start = end + 1;
                end += 1;
            }
            None => {
                cols.push(line[start..].to_vec());
                break;
            }
        }
    }

    cols
}

/// Check if CSV file is valid detected_transcripts.csv file
pub fn validate_detected_transcripts_csv_file(
    file_path: &Path,
    col_of_interest: &Vec<&str>,
) -> io::Result<()> {
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
