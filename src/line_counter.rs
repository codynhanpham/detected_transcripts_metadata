use indicatif::{ProgressBar, ProgressState, ProgressStyle};
use memchr::memchr_iter;
use std::io::Read;
use rayon::prelude::*;

fn count_newlines(bytes: &[u8]) -> u64 {
    memchr_iter(b'\n', bytes).count() as u64
}


fn find_new_line_pos(bytes: &[u8]) -> Option<usize> {
    bytes.iter().rposition(|&b| b == b'\n')
}


pub fn count_total_lines(mut file: std::fs::File, n_threads: usize, chunk_size: usize, quiet: bool) -> u64 {
    let (sender, receiver) = crossbeam_channel::bounded::<Box<[u8]>>(2_000);
    let chunk_size = chunk_size * 1024;
    
    if !quiet {
        print!("Using {} threads ", n_threads);
        println!("x {} KiB chunks", chunk_size / 1024);
    }

    let mut handles = Vec::with_capacity(n_threads);
    
    for _ in 0..n_threads {
        let receiver = receiver.clone();
        let handle = std::thread::spawn(move || {
            let mut line_count = 0u64;
            for buf in receiver {
                line_count += count_newlines(&buf);
            }
            line_count
        });
        handles.push(handle);
    }

    let mut buf = vec![0; chunk_size];
    let mut bytes_not_processed = 0;

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

        let actual_buf = &mut buf[..bytes_not_processed + bytes_read];
        let last_new_line_index = match find_new_line_pos(&actual_buf) {
            Some(index) => index,
            None => {
                bytes_not_processed += bytes_read;
                if bytes_not_processed == buf.len() {
                    panic!("No new line found in the read buffer");
                }
                continue;
            }
        };
        let buf_boxed = Box::<[u8]>::from(&actual_buf[..(last_new_line_index + 1)]);
        sender.send(buf_boxed).expect("Failed to send buffer");

        actual_buf.copy_within(last_new_line_index + 1.., 0);
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

    // Combine results from all threads
    let total_lines: u64 = handles.into_par_iter()
        .map(|h| h.join().unwrap())
        .sum();

    total_lines
}