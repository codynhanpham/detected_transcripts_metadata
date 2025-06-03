use fast_float;
use indicatif::{ProgressBar, ProgressState, ProgressStyle};
use memchr::memchr;
use rustc_hash::FxHashMap;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::io::Read;
use std::{fs::File, io};

use crate::utils;


#[derive(Debug, Clone)]
pub struct TranscriptMetadata {
    pub cell_map: FxHashMap<String, usize>, // Cell: Transcript count
    pub total_transcripts: usize, // Total number of transcripts (== number of rows in the file)
    pub transcripts_within_cells: usize, // Number of transcripts within cells (cell_id != "NA" or "N/A" or "-1")
    pub layers_transcripts_count: FxHashMap<i32, usize>, // Number of transcripts for each z-layer
    pub layers_transcripts_within_cells_count: FxHashMap<i32, usize>, // Number of transcripts within cells for each z-layer
    pub genes_frequency: FxHashMap<String, usize>, // Frequency of each gene
    pub genes_frequency_per_layer: FxHashMap<i32, FxHashMap<String, usize>>, // Frequency of each gene per z-layer
}

impl TranscriptMetadata {
    pub fn new() -> Self {
        Self {
            cell_map: FxHashMap::default(),
            total_transcripts: 0,
            transcripts_within_cells: 0,
            layers_transcripts_count: FxHashMap::default(),
            layers_transcripts_within_cells_count: FxHashMap::default(),
            genes_frequency: FxHashMap::default(),
            genes_frequency_per_layer: FxHashMap::default(),
        }
    }

    /// Handle exporting metadata to JSON, either to console or to a file as determined by the ``print_json_output_only`` argument
    pub fn export_json(&self, metadata_file_output_path: &str, print_json_output_only: &bool) -> io::Result<()> {
        let processed_metadata = TranscriptMetadataExport::new(&self);

        if *print_json_output_only || metadata_file_output_path == "" || metadata_file_output_path == "-" {
            let metadata_json = processed_metadata.as_json();
            println!("{}", metadata_json);
            return Ok(());
        }

        let metadata_file = File::create(metadata_file_output_path).unwrap();
        serde_json::to_writer(metadata_file, &processed_metadata).unwrap();

        Ok(())
    }
}


fn process_chunk(data: &[u8], col_of_interest_indices: &Vec<usize>, result: &mut TranscriptMetadata) -> () {
    
    let mut buffer = &data[..];
    
    loop {
        match memchr(b'\n', buffer) {
            None => break,
            Some(line_end) => {
                let line = &buffer[..line_end];
                
                let row = utils::get_col_data(&line);
                let data: Result<Vec<Vec<u8>>, &str> = col_of_interest_indices
                    .iter()
                    .map(|&col| row.get(col).cloned().ok_or("Index out of bounds"))
                    .collect();
                
                // Order of the columns in 'data' depends on order of these header set by col_of_interest in main()
                let data = data.expect("Could not parse some row data");
                
                // Process found values if we have all required fields
                if data.len() != col_of_interest_indices.len() {
                    panic!("Found a row that does not have matching number of columns");
                }

                // // If collecting the coordinates is needed. Must update the TranscriptMetadata struct accordingly
                // let global_x = fast_float::parse::<f64, _>(data[0].as_slice()).expect("Could not parse some global_x value to float64");
                // let global_y = fast_float::parse::<f64, _>(data[1].as_slice()).expect("Could not parse some global_y value to float64");
                let global_z = fast_float::parse::<f32, _>(data[2].as_slice()).expect("Could not parse some global_z value to float32");
                let global_z = global_z.round() as i32;

                let gene = String::from_utf8_lossy(data[3].as_slice()).to_string();
                let cell_id = String::from_utf8_lossy(data[4].as_slice()).to_string();


                // Update this thread result
                result.total_transcripts += 1;

                if !cell_id.is_empty() && cell_id != "NA" && cell_id != "N/A" && cell_id != "-1" {
                    result.transcripts_within_cells += 1;

                    result.cell_map
                        .entry(cell_id.clone())
                        .and_modify(|e| *e += 1)
                        .or_insert(1);

                    result.layers_transcripts_within_cells_count
                        .entry(global_z)
                        .and_modify(|e| *e += 1)
                        .or_insert(1);
                }

                result.layers_transcripts_count
                    .entry(global_z)
                    .and_modify(|e| *e += 1)
                    .or_insert(1);

                result.genes_frequency
                    .entry(gene.clone())
                    .and_modify(|e| *e += 1)
                    .or_insert(1);

                result.genes_frequency_per_layer
                    .entry(global_z)
                    .or_insert(FxHashMap::default())
                    .entry(gene.clone())
                    .and_modify(|e| *e += 1)
                    .or_insert(1);

                // result.genes_transcripts_coords
                //     .entry(gene.clone())
                //     .and_modify(|e| e.push((global_x, global_y, global_z)))
                //     .or_insert(vec![(global_x, global_y, global_z)]);
            

                buffer = &buffer[line_end + 1..];
            }
        }
    }

}


fn find_new_line_pos(bytes: &[u8]) -> Option<usize> {
    bytes.iter().rposition(|&b| b == b'\n')
}


/// Given a Vec<TranscriptMetadata>, merge all the metadata into a single TranscriptMetadata
pub fn merge_metadata_chunks(handles: Vec<std::thread::JoinHandle<TranscriptMetadata>>) -> TranscriptMetadata {
    let merged_metadata = handles.into_par_iter()
        .map(|handle| handle.join().unwrap())
        .reduce(
            || TranscriptMetadata::new(),
            |mut acc, map| {
                acc.total_transcripts += map.total_transcripts;
                acc.transcripts_within_cells += map.transcripts_within_cells;

                for (cell, count) in map.cell_map.iter() {
                    acc.cell_map
                        .entry(cell.clone())
                        .and_modify(|e| *e += *count)
                        .or_insert(*count);
                }
                
                for (layer, count) in map.layers_transcripts_count.iter() {
                    acc.layers_transcripts_count
                        .entry(*layer)
                        .and_modify(|e| *e += *count)
                        .or_insert(*count);
                }

                for (layer, within_cells_count) in map.layers_transcripts_within_cells_count.iter() {
                    acc.layers_transcripts_within_cells_count
                        .entry(*layer)
                        .and_modify(|e| *e += *within_cells_count)
                        .or_insert(*within_cells_count);
                }

                for (gene, count) in map.genes_frequency.iter() {
                    acc.genes_frequency
                        .entry(gene.clone())
                        .and_modify(|e| *e += *count)
                        .or_insert(*count);
                }

                for (layer, genes_frequency) in map.genes_frequency_per_layer.iter() {
                    for (gene, count) in genes_frequency.iter() {
                        acc.genes_frequency_per_layer
                            .entry(*layer)
                            .or_insert(FxHashMap::default())
                            .entry(gene.clone())
                            .and_modify(|e| *e += *count)
                            .or_insert(*count);
                    }
                }
                
                // for (gene, coords) in map.genes_transcripts_coords.iter() {
                //     acc.genes_transcripts_coords
                //         .entry(gene.clone())
                //         .and_modify(|e| e.extend(coords.iter()))
                //         .or_insert(coords.clone());
                // }

                acc
            }
        );

    merged_metadata
}


/// Run the full metadata extraction process
pub fn process_metadata_file(mut file: std::fs::File, col_of_interest_indices: Vec<usize>, n_threads: usize, chunk_size: usize, quiet: bool) -> TranscriptMetadata {
    // Start the processor threads
    let (sender, receiver) = crossbeam_channel::bounded::<Box<[u8]>>(2_000);
    let chunk_size = chunk_size * 1024;

    if !quiet {
        print!("Using {} threads ", n_threads);
        println!("x {} KiB chunks", chunk_size / 1024);
    }

    let mut handles = Vec::with_capacity(n_threads);
    
    // Spawn worker threads
    for _ in 0..n_threads {
        let receiver = receiver.clone();
        let indices = col_of_interest_indices.clone();

        let handle = std::thread::spawn(move || {
            let mut result: TranscriptMetadata = TranscriptMetadata::new();
            // wait until the sender sends the chunk
            for buf in receiver {
                process_chunk(&buf, &indices, &mut result);
            }
            result
        });
        handles.push(handle);
    }


    // Read the file in chunks and send the chunks to the processor threads
    let mut buf = vec![0; chunk_size];
    let mut bytes_not_processed: usize;


    // Remove the first line (header) from the initial read
    let initial_read = file.read(&mut buf).expect("Failed to read file");
    if let Some(first_newline) = memchr::memchr(b'\n', &buf[..initial_read]) {
        // Copy everything after the first newline to the start of the buffer
        buf.copy_within(first_newline + 1..initial_read, 0);
        bytes_not_processed = initial_read - first_newline - 1;
    } else {
        panic!("Could not find header line terminator. You will need to increase the chunk size by adding the -c <size-in-KiB> argument.");
    }

    let mut progress: u64 = 0;
    let mut last_printed_progress: u64 = 0;
    let mut pb: ProgressBar = ProgressBar::hidden();
    if !quiet {
        let file_size = file.metadata().unwrap().len();
        progress = 0;
        last_printed_progress = 0;

        pb = ProgressBar::new(file_size);
        pb.set_style(
            ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {bytes}/{total_bytes} ({eta})")
            .unwrap()
            .with_key("eta", |state: &ProgressState, w: &mut dyn std::fmt::Write| write!(w, "{:.1}s", state.eta().as_secs_f64()).unwrap())
            .progress_chars("#>-")
        );
    }

    loop {
        let bytes_read = file.read(&mut buf[bytes_not_processed..]).expect("Failed to read file");
        if bytes_read == 0 {
            break;
        }

        let actual_buf = &mut buf[..bytes_not_processed+bytes_read];
        let last_new_line_index = match find_new_line_pos(&actual_buf) {
            Some(index) => index,
            None => {
                println!("No new line found in the read buffer");
                bytes_not_processed += bytes_read;
                if bytes_not_processed == buf.len(){
                    panic!("No new line found in the read buffer");
                }
                continue; // try again, maybe we next read will have a newline
            }
        };

        let buf_boxed = Box::<[u8]>::from(&actual_buf[..(last_new_line_index + 1)]);
        sender.send(buf_boxed).expect("Failed to send buffer");

        actual_buf.copy_within(last_new_line_index+1.., 0);
        // You cannot use bytes_not_processed = bytes_read - last_new_line_index
        // - 1; because the buffer will contain unprocessed bytes from the
        // previous iteration and the new line index will be calculated from the
        // start of the buffer
        bytes_not_processed = actual_buf.len() - last_new_line_index - 1;

        if !quiet {
            progress += bytes_read as u64;
            if progress - last_printed_progress >= 150 * 1024 * 1024 { // 150 MB ~ 1 second on typical HDD
                pb.set_position(progress);
                last_printed_progress = progress;
            }
        }
    }
    drop(sender);
    pb.finish();

    // Combine data from all threads
    if !quiet {
        println!("\nMerging result from {} threads...", n_threads);
    }
    let result = merge_metadata_chunks(handles);

    return result;
}




/// Struct for exporting metadata to JSON
#[derive(Debug, Serialize, Deserialize)]
pub struct TranscriptMetadataExport {
    cell_count: usize, // Number of cells with transcripts
    transcripts_per_cell: Vec<usize>, // Number of transcripts per cell
    total_transcripts: usize, // Total number of transcripts (== number of rows in the file)
    transcripts_within_cells: usize, // Number of transcripts within cells (cell_id != "NA" or "N/A" or "-1")
    layers_transcripts_count: FxHashMap<i32, usize>, // Number of transcripts within cells for each z-layer
    layers_transcripts_within_cells_count: FxHashMap<i32, usize>, // Number of transcripts within cells for each z-layer
    genes_frequency: FxHashMap<String, usize>, // Frequency of each gene
    genes_frequency_per_layer: FxHashMap<i32, FxHashMap<String, usize>>, // Frequency of each gene per z-layer
}

impl TranscriptMetadataExport {
    pub fn new(metadata: &TranscriptMetadata) -> Self {
        Self {
            cell_count: metadata.cell_map.len(),
            transcripts_per_cell: metadata.cell_map.iter().map(|entry| *entry.1).collect(),
            total_transcripts: metadata.total_transcripts,
            transcripts_within_cells: metadata.transcripts_within_cells,
            layers_transcripts_count: metadata.layers_transcripts_count.clone(),
            layers_transcripts_within_cells_count: metadata.layers_transcripts_within_cells_count.clone(),
            genes_frequency: metadata.genes_frequency.clone(),
            genes_frequency_per_layer: metadata.genes_frequency_per_layer.clone(),
        }
    }

    /// Serialize the TranscriptMetadataExport struct to JSON
    pub fn as_json(&self) -> String {
        serde_json::to_string_pretty(self).unwrap()
    }
}