use std::env;
use std::fs;
use std::io::{BufRead, BufReader};
use std::process::{Command, Stdio};

const SPACE: usize = 1 << 24;

fn encode_base(base: u8) -> Option<usize> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn encode_barcode(sequence: &[u8]) -> Option<usize> {
    if sequence.len() < 42 {
        return None;
    }
    let positions = [0usize, 1, 2, 3, 19, 20, 21, 22, 38, 39, 40, 41];
    let mut code = 0usize;
    for position in positions {
        code = (code << 2) | encode_base(sequence[position])?;
    }
    Some(code)
}

fn encode_whitelist_barcode(sequence: &[u8]) -> Option<usize> {
    if sequence.len() != 12 {
        return None;
    }
    let mut code = 0usize;
    for base in sequence {
        code = (code << 2) | encode_base(*base)?;
    }
    Some(code)
}

fn load_whitelist(path: &str) -> Result<(Vec<u8>, Vec<u8>, usize), String> {
    let content = fs::read_to_string(path).map_err(|error| format!("{}: {}", path, error))?;
    let mut exact = vec![0u8; SPACE];
    let mut corrected = vec![0u8; SPACE];
    let mut count = 0usize;
    for (line_number, line) in content.lines().enumerate() {
        let barcode = line
            .split(',')
            .nth(2)
            .ok_or_else(|| format!("{}:{} missing barcode column", path, line_number + 1))?
            .trim();
        let code = encode_whitelist_barcode(barcode.as_bytes())
            .ok_or_else(|| format!("{}:{} invalid 12-base barcode", path, line_number + 1))?;
        if exact[code] == 0 {
            count += 1;
        }
        exact[code] = 1;
        corrected[code] = 1;
        for position in 0..12 {
            let shift = 2 * (11 - position);
            let original = (code >> shift) & 3;
            let cleared = code & !(3usize << shift);
            for replacement in 0..4 {
                if replacement != original {
                    corrected[cleared | (replacement << shift)] = 1;
                }
            }
        }
    }
    Ok((exact, corrected, count))
}

#[derive(Default)]
struct Metrics {
    exact: u64,
    corrected: u64,
    exact_poly_t: u64,
    corrected_poly_t: u64,
}

fn main() -> Result<(), String> {
    let arguments: Vec<String> = env::args().collect();
    if arguments.len() < 5 {
        return Err(
            "usage: chip_assignment_scan SAMPLE --targets=N TARGET1[=WHITELIST] [TARGET2=WHITELIST] R1.fastq.gz [R1.fastq.gz ...]"
                .to_string(),
        );
    }
    let sample = &arguments[1];
    let (target_start, target_count) = if let Some(value) = arguments[2].strip_prefix("--targets=") {
        let count = value.parse::<usize>().map_err(|_| "invalid --targets value".to_string())?;
        if !(1..=2).contains(&count) {
            return Err("--targets must be 1 or 2".to_string());
        }
        (3usize, count)
    } else {
        (2usize, 2usize)
    };
    let target_end = target_start + target_count;
    if arguments.len() <= target_end {
        return Err("at least one R1 FASTQ is required".to_string());
    }
    let mut target_names = Vec::new();
    let mut exact_tables = Vec::new();
    let mut corrected_tables = Vec::new();
    let mut whitelist_sizes = Vec::new();
    for specification in &arguments[target_start..target_end] {
        let (target, path) = specification
            .split_once('=')
            .ok_or_else(|| format!("invalid target specification: {}", specification))?;
        let (exact, corrected, size) = load_whitelist(path)?;
        target_names.push(target.to_string());
        exact_tables.push(exact);
        corrected_tables.push(corrected);
        whitelist_sizes.push(size);
    }

    let mut command = Command::new("pigz");
    let pigz_threads = env::var("CHIP_QC_PIGZ_THREADS").unwrap_or_else(|_| "4".to_string());
    command.args(["-p", &pigz_threads, "-dc"]);
    for path in &arguments[target_end..] {
        command.arg(path);
    }
    let mut child = command
        .stdout(Stdio::piped())
        .spawn()
        .map_err(|error| format!("unable to launch pigz: {}", error))?;
    let stdout = child.stdout.take().ok_or_else(|| "pigz stdout unavailable".to_string())?;
    let mut reader = BufReader::with_capacity(8 * 1024 * 1024, stdout);
    let mut buffer = Vec::with_capacity(256);
    let mut line_index = 0u64;
    let mut total_reads = 0u64;
    let mut poly_t_reads = 0u64;
    let mut metrics = (0..target_count).map(|_| Metrics::default()).collect::<Vec<_>>();

    loop {
        buffer.clear();
        let bytes = reader
            .read_until(b'\n', &mut buffer)
            .map_err(|error| format!("FASTQ read failure: {}", error))?;
        if bytes == 0 {
            break;
        }
        if line_index % 4 == 1 {
            while matches!(buffer.last(), Some(b'\n' | b'\r')) {
                buffer.pop();
            }
            total_reads += 1;
            if let Some(code) = encode_barcode(&buffer) {
                let poly_t = buffer.len() >= 70
                    && buffer[52..70].iter().filter(|base| **base == b'T').count() >= 12;
                if poly_t {
                    poly_t_reads += 1;
                }
                for index in 0..target_count {
                    if exact_tables[index][code] != 0 {
                        metrics[index].exact += 1;
                        if poly_t {
                            metrics[index].exact_poly_t += 1;
                        }
                    }
                    if corrected_tables[index][code] != 0 {
                        metrics[index].corrected += 1;
                        if poly_t {
                            metrics[index].corrected_poly_t += 1;
                        }
                    }
                }
            }
        }
        line_index += 1;
    }
    let status = child.wait().map_err(|error| format!("pigz wait failure: {}", error))?;
    if !status.success() {
        return Err(format!("pigz exited with {}", status));
    }
    if line_index % 4 != 0 {
        return Err(format!("truncated FASTQ stream: {} lines", line_index));
    }

    println!(
        "sample_id\ttested_chip_id\twhitelist_barcodes\ttotal_reads\tpoly_t_reads\texact_reads\tcorrected_reads\texact_poly_t_reads\tcorrected_poly_t_reads"
    );
    for index in 0..target_count {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            sample,
            target_names[index],
            whitelist_sizes[index],
            total_reads,
            poly_t_reads,
            metrics[index].exact,
            metrics[index].corrected,
            metrics[index].exact_poly_t,
            metrics[index].corrected_poly_t
        );
    }
    Ok(())
}
