use std::{
    fs::File, io::{self}, sync::RwLock
};

use dashmap::DashMap;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json;
use ndarray::Array2;


#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct TranscriptHistogram2dData {
    pub resolution_microns: f64,
    pub data: Vec<(f64, f64, usize)>, // (x_start_edge, y_start_edge, count)
}

impl TranscriptHistogram2dData {
    pub fn new(resolution_microns: f64, data: Vec<(f64, f64, usize)>) -> Self {
        Self {
            resolution_microns,
            data,
        }
    }
}


#[derive(Debug, Serialize, Deserialize)]
pub struct TranscriptMetadata {
    pub cell_map: DashMap<String, usize>, // Cell: Transcript count
    pub total_transcripts: RwLock<usize>, // Total number of transcripts (== number of rows in the file)
    pub transcripts_within_cells: RwLock<usize>, // Number of transcripts within cells (cell_id != "NA" or "N/A" or "-1")
    pub layers_transcripts_count: DashMap<i32, usize>, // Number of transcripts for each z-layer
    pub layers_transcripts_within_cells_count: DashMap<i32, usize>, // Number of transcripts within cells for each z-layer
    pub genes_frequency: DashMap<String, usize>, // Frequency of each gene
    pub genes_frequency_per_layer: DashMap<i32, DashMap<String, usize>>, // Frequency of each gene per z-layer
    pub genes_transcripts_coords: DashMap<String, Vec<(f64, f64, i32)>>, // <gene: (x,y,z)> coord for each transcript
    pub transcripts_mask_area: RwLock<f64>, // Area the transcripts occupy in all layers (area behind the scatter plot)
    pub transcripts_histogram2d_data: RwLock<TranscriptHistogram2dData>, // 2D histogram data for the transcripts
}

impl TranscriptMetadata {
    pub fn new() -> Self {
        Self {
            cell_map: DashMap::new(),
            total_transcripts: RwLock::new(0),
            transcripts_within_cells: RwLock::new(0),
            layers_transcripts_count: DashMap::new(),
            layers_transcripts_within_cells_count: DashMap::new(),
            genes_frequency: DashMap::new(),
            genes_frequency_per_layer: DashMap::new(),
            genes_transcripts_coords: DashMap::new(),
            transcripts_mask_area: RwLock::new(-1.0),
            transcripts_histogram2d_data: RwLock::new(TranscriptHistogram2dData::new(-1.0, vec![])),
        }
    }

    pub fn add_data_chunk(&self, data_chunk: &[Vec<String>]) {
        // Order of the columns in the data_chunk depends on order of these header set by col_of_interest in main()
        // col_of_interest = vec!["global_y", "global_y", "global_z", "cell_id", "gene"];
        data_chunk.par_iter().for_each(|row| {
            // Parse the global_x, global_y to f64
            let global_x = match row[0].parse::<f64>() {
                Ok(x) => x,
                Err(_) => {
                    eprintln!("Error: Could not parse global_x value to float32: {}", row[0]);
                    return;
                },
            };
            let global_y = match row[1].parse::<f64>() {
                Ok(y) => y,
                Err(_) => {
                    eprintln!("Error: Could not parse global_y value to float32: {}", row[0]);
                    return;
                },
            };
            // Parse the global_z to i32
            let global_z = match row[2].parse::<f32>() {
                Ok(z) => z.round() as i32,
                Err(_) => {
                    eprintln!("Error: Could not parse global_z value to int32: {}", row[0]);
                    return;
                },
            };
            let cell_id = &row[3];
            let gene = &row[4];

            // Increment fields with new chunk data

            *self.total_transcripts.write().unwrap() += 1;

            if !cell_id.is_empty() && cell_id != "NA" && cell_id != "N/A" && cell_id != "-1" {
                *self.transcripts_within_cells.write().unwrap() += 1;

                self.cell_map
                .entry(cell_id.clone())
                .and_modify(|e| *e += 1)
                .or_insert(1);

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

            self.genes_transcripts_coords
                .entry(gene.clone())
                .and_modify(|e| e.push((global_x, global_y, global_z)))
                .or_insert(vec![(global_x, global_y, global_z)]);
        });
    }

    /// Given the `TranscriptMetadata`, return the gene coordinates matrix of all genes (concatenated)
    fn get_genes_coords_matrix(&self) -> Vec<(f64, f64, i32)> {
        let genes_coords_matrix: Vec<(f64, f64, i32)> = self.genes_transcripts_coords.iter()
            .flat_map(|entry| entry.value().clone())
            .collect();
        
        genes_coords_matrix
    }

    /// Calculate the area the transcripts occupy in all layers (area behind the scatter plot)
    /// Update the `transcripts_mask_area` field in the `TranscriptMetadata`
    pub fn calculate_transcript_mask_area(&self, resolution_microns: f64, density_threshold: usize) -> () {
        let genes_coords_matrix = self.get_genes_coords_matrix();
        let (x, y): (Vec<f64>, Vec<f64>) = genes_coords_matrix.par_iter().map(|(x, y, _)| (*x, *y)).unzip();

        let x_range = (x.iter().cloned().fold(f64::INFINITY, f64::min), x.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
        let y_range = (y.iter().cloned().fold(f64::INFINITY, f64::min), y.iter().cloned().fold(f64::NEG_INFINITY, f64::max));

        let bin_count = ((x_range.1 - x_range.0) / resolution_microns).round() as usize;

        let hist = std::sync::Mutex::new(Array2::<usize>::zeros((bin_count, bin_count)));
        let x_bin_width = (x_range.1 - x_range.0) / bin_count as f64;
        let y_bin_width = (y_range.1 - y_range.0) / bin_count as f64;

        // Parallelize the histogram filling process
        x.par_iter().zip(y.par_iter()).for_each(|(&x, &y)| {
            let x_bin = ((x - x_range.0) / x_bin_width).floor() as usize;
            let y_bin = ((y - y_range.0) / y_bin_width).floor() as usize;
            if x_bin < bin_count && y_bin < bin_count {
                let mut hist = hist.lock().unwrap();
                hist[[x_bin, y_bin]] += 1;
            }
        });

        let hist = hist.into_inner().unwrap();
        let over_threshold = hist.mapv(|x| if x >= density_threshold { 1 } else { 0 });
        let area = over_threshold.sum() as f64 * x_bin_width * y_bin_width;

        *self.transcripts_mask_area.write().unwrap() = area;

        // Update the 2D histogram data with x_start_edge, y_start_edge, count
        let mut data = vec![];
        for (i, row) in over_threshold.outer_iter().enumerate() {
            for (j, &count) in row.iter().enumerate() {
                if count > 0 {
                    let x_start_edge = x_range.0 + i as f64 * x_bin_width;
                    let y_start_edge = y_range.0 + j as f64 * y_bin_width;
                    data.push((x_start_edge, y_start_edge, count));
                }
            }
        }

        *self.transcripts_histogram2d_data.write().unwrap() = TranscriptHistogram2dData::new(resolution_microns, data);
    }

    pub fn export_json(&self, metadata_file_output_path: &str, print_json_output_only: &bool) -> io::Result<()> {
        let processed_metadata = TranscriptMetadataExport::new(&self);

        // If the output value is the special value "-", print to console instead of saving to file
        if *print_json_output_only {
            let metadata_json = serde_json::to_string_pretty(&processed_metadata).unwrap();
            println!("{}", metadata_json);
            return Ok(());
        }

        // Otherwise, save metadata as JSON file
        let metadata_file = File::create(metadata_file_output_path).unwrap();
        serde_json::to_writer(metadata_file, &processed_metadata).unwrap();

        Ok(())
    }
}



/// Struct for exporting metadata to JSON
#[derive(Debug, Serialize, Deserialize)]
pub struct TranscriptMetadataExport {
    cell_count: usize, // Number of cells with transcripts
    transcripts_per_cell: Vec<usize>, // Number of transcripts per cell
    total_transcripts: usize, // Total number of transcripts (== number of rows in the file)
    transcripts_within_cells: usize, // Number of transcripts within cells (cell_id != "NA" or "N/A" or "-1")
    layers_transcripts_count: DashMap<i32, usize>, // Number of transcripts for each z-layer
    layers_transcripts_within_cells_count: DashMap<i32, usize>, // Number of transcripts within cells for each z-layer
    genes_frequency: DashMap<String, usize>, // Frequency of each gene
    genes_frequency_per_layer: DashMap<i32, DashMap<String, usize>>, // Frequency of each gene per z-layer
    // genes_transcripts_coords: DashMap<String, Vec<(f32, f32, i32)>>, // <gene: (x,y,z)> coord for each transcript
    transcripts_mask_area: f64, // Area the transcripts occupy in all layers (area behind the scatter plot)
    transcripts_histogram2d_data: TranscriptHistogram2dData, // 2D histogram data for the transcripts
}

impl TranscriptMetadataExport {
    pub fn new(metadata: &TranscriptMetadata) -> Self {
        Self {
            cell_count: metadata.cell_map.len(),
            transcripts_per_cell: metadata.cell_map.iter().map(|entry| *entry.value()).collect(),
            total_transcripts: *metadata.total_transcripts.read().unwrap(),
            transcripts_within_cells: *metadata.transcripts_within_cells.read().unwrap(),
            layers_transcripts_count: metadata.layers_transcripts_count.clone(),
            layers_transcripts_within_cells_count: metadata.layers_transcripts_within_cells_count.clone(),
            genes_frequency: metadata.genes_frequency.clone(),
            genes_frequency_per_layer: metadata.genes_frequency_per_layer.clone(),
            // genes_transcripts_coords: metadata.genes_transcripts_coords.clone(),
            transcripts_mask_area: *metadata.transcripts_mask_area.read().unwrap(),
            transcripts_histogram2d_data: metadata.transcripts_histogram2d_data.read().unwrap().clone(),
        }
    }
}




/// Given a Vec<TranscriptMetadata>, merge all the metadata into a single TranscriptMetadata
pub fn merge_metadata(metadata: Vec<TranscriptMetadata>) -> TranscriptMetadata {
    let merged_metadata = TranscriptMetadata::new();

    metadata.into_par_iter().for_each(|data| {
        // Combine the RwLocks:
        *merged_metadata.total_transcripts.write().unwrap() += *data.total_transcripts.read().unwrap();
        *merged_metadata.transcripts_within_cells.write().unwrap() += *data.transcripts_within_cells.read().unwrap();

        // Extend + Combine the DashMaps:

        for entry in data.cell_map.iter() {
            let (cell, count) = entry.pair();
            merged_metadata.cell_map
                .entry(cell.clone())
                .and_modify(|e| *e += *count)
                .or_insert(*count);
        }
        
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

        for entry in data.genes_transcripts_coords.iter() {
            let (gene, coords) = entry.pair();
            merged_metadata.genes_transcripts_coords
                .entry(gene.clone())
                .and_modify(|e| e.extend(coords.iter()))
                .or_insert(coords.clone());
        }
    });

    merged_metadata
}



