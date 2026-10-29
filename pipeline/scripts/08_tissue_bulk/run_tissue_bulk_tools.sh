#!/usr/bin/env bash
set -euo pipefail
unset DISPLAY

usage() {
  cat >&2 <<'EOF'
Usage: run_tissue_bulk_tools.sh \
  --sample SAMPLE --stage-dir DIR --tmp-dir DIR --barcode-table FILE \
  --tissue-barcodes FILE --fastq FILE --reference-fasta FILE --gtf FILE \
  --gedi-genome NAME [--threads N] [--trim-5p N]

Runs tissue-only FASTQ filtering, slamdunk and GRAND-SLAM with ordered stage
outputs. The large filtered FASTQ is retained under --tmp-dir and recorded in
05_metadata/large_file_manifest.tsv.
EOF
}

sample= stage_dir= tmp_dir= barcode_table= tissue_barcodes= fastq=
reference_fasta= gtf= gedi_genome= threads=32 trim5=22
interpretation_mode=${INTERPRETATION_MODE:-standard}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) sample=$2; shift 2 ;;
    --stage-dir) stage_dir=$2; shift 2 ;;
    --tmp-dir) tmp_dir=$2; shift 2 ;;
    --barcode-table) barcode_table=$2; shift 2 ;;
    --tissue-barcodes) tissue_barcodes=$2; shift 2 ;;
    --fastq) fastq=$2; shift 2 ;;
    --reference-fasta) reference_fasta=$2; shift 2 ;;
    --gtf) gtf=$2; shift 2 ;;
    --gedi-genome) gedi_genome=$2; shift 2 ;;
    --threads) threads=$2; shift 2 ;;
    --trim-5p) trim5=$2; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

for value in sample stage_dir tmp_dir barcode_table tissue_barcodes fastq reference_fasta gtf gedi_genome; do
  [[ -n "${!value}" ]] || { echo "ERROR: --${value//_/-} is required" >&2; exit 2; }
done
for file in "$barcode_table" "$tissue_barcodes" "$fastq" "$reference_fasta" "$gtf"; do
  [[ -s "$file" ]] || { echo "ERROR: missing or empty input: $file" >&2; exit 2; }
done
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --threads must be a positive integer" >&2; exit 2; }
[[ "$trim5" =~ ^[0-9]+$ ]] || { echo "ERROR: --trim-5p must be a nonnegative integer" >&2; exit 2; }
[[ "$interpretation_mode" == "standard" || "$interpretation_mode" == "unlabeled_control" ]] || {
  echo "ERROR: INTERPRETATION_MODE must be standard or unlabeled_control" >&2
  exit 2
}
if [[ "$interpretation_mode" == "unlabeled_control" ]]; then
  grand_limitation="Unlabeled negative control: GRAND-SLAM NTR is apparent background or false-positive signal, not newly synthesized RNA."
else
  grand_limitation="No matched no4sU control or validated error model was supplied; NTR is exploratory."
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$stage_dir"/{01_tables,02_figures,03_logs,04_tool_output,05_metadata} "$tmp_dir"
utr_bed="$stage_dir/01_tables/02_utr_regions.bed"
tissue_fastq="$tmp_dir/${sample}.tissue.fastq.gz"
slamdunk_dir="$stage_dir/04_tool_output/slamdunk"
grand_dir="$stage_dir/04_tool_output/grand_slam"
mkdir -p "$slamdunk_dir" "$grand_dir"

input_signature_file="$stage_dir/05_metadata/bulk_input_signature.sha256"
input_signature=$(
  {
    printf 'sample=%s\ntrim5=%s\ngedi_genome=%s\n' "$sample" "$trim5" "$gedi_genome"
    for file in "$barcode_table" "$fastq" "$reference_fasta" "$gtf"; do
      stat -c '%n\t%s\t%Y' "$file"
    done
    sha256sum "$tissue_barcodes"
  } | sha256sum | cut -d' ' -f1
)
previous_signature=""
if [[ -s "$input_signature_file" ]]; then read -r previous_signature < "$input_signature_file"; fi
if [[ "$previous_signature" != "$input_signature" ]]; then
  # Generated caches are valid only for the exact tissue mask, source inputs,
  # reference and parameters above. This makes a new manual lasso propagate
  # through tissue FASTQ, slamdunk and GRAND-SLAM instead of reusing stale data.
  [[ "$tissue_fastq" == "$tmp_dir/"* ]] || { echo "ERROR: unsafe tissue FASTQ cache path: $tissue_fastq" >&2; exit 2; }
  [[ "$slamdunk_dir" == "$stage_dir/04_tool_output/"* ]] || { echo "ERROR: unsafe slamdunk cache path: $slamdunk_dir" >&2; exit 2; }
  [[ "$grand_dir" == "$stage_dir/04_tool_output/"* ]] || { echo "ERROR: unsafe GRAND-SLAM cache path: $grand_dir" >&2; exit 2; }
  rm -f -- "$tissue_fastq"
  rm -rf -- "$slamdunk_dir" "$grand_dir"
  mkdir -p "$slamdunk_dir" "$grand_dir"
fi

run_timed() {
  local name=$1
  shift
  local started ended rc
  started=$(date -u +%FT%TZ)
  set +e
  if [[ -x /usr/bin/time ]]; then
    /usr/bin/time -v -o "$stage_dir/03_logs/${name}.time.txt" \
      "$@" >"$stage_dir/03_logs/${name}.stdout.log" 2>"$stage_dir/03_logs/${name}.stderr.log"
  else
    python3 "$script_dir/../00_shared/run_with_resources.py" \
      --output "$stage_dir/03_logs/${name}.time.txt" -- "$@" \
      >"$stage_dir/03_logs/${name}.stdout.log" 2>"$stage_dir/03_logs/${name}.stderr.log"
  fi
  rc=$?
  set -e
  ended=$(date -u +%FT%TZ)
  printf 'step\tstarted_at_utc\tended_at_utc\texit_code\n%s\t%s\t%s\t%s\n' \
    "$name" "$started" "$ended" "$rc" >"$stage_dir/03_logs/${name}.status.tsv"
  return "$rc"
}

cat >"$stage_dir/05_metadata/bulk_tool_input_manifest.tsv" <<EOF
role\tpath
barcode_table\t$barcode_table
tissue_barcodes\t$tissue_barcodes
source_fastq\t$fastq
reference_fasta\t$reference_fasta
annotation_gtf\t$gtf
EOF
cat >"$stage_dir/05_metadata/large_file_manifest.tsv" <<EOF
path\tpurpose\tdeliverable
$tissue_fastq\tTissue-only R2 reads used by slamdunk\tfalse
EOF
cat >"$stage_dir/05_metadata/bulk_method.json" <<EOF
{
  "schema_version": 1,
  "sample": "$sample",
  "scope": "tissue_only_bulk",
  "interpretation_mode": "$interpretation_mode",
  "slamdunk_trim_5p": $trim5,
  "trim_5p_rationale": "Remove the 22 nt template-switch sequence AAGCAGTGGTATCAACGCAGAG observed at R2 starts.",
  "grand_slam_error_mode": "exploratory_sample_background",
  "grand_slam_limitation": "$grand_limitation",
  "threads": $threads
}
EOF

if [[ ! -s "$utr_bed" ]]; then
  run_timed build_utr_bed python3 "$script_dir/build_utr_bed.py" --gtf "$gtf" --output "$utr_bed"
fi
if [[ ! -s "$tissue_fastq" ]]; then
  run_timed filter_tissue_fastq python3 "$script_dir/filter_tissue_fastq.py" \
    --barcode-table "$barcode_table" --tissue-barcodes "$tissue_barcodes" \
    --fastq "$fastq" --output "$tissue_fastq" \
    --summary "$stage_dir/05_metadata/tissue_fastq_summary.json"
fi
if ! find "$slamdunk_dir/filter" -type f -name '*_filtered.bam' -size +0c -print -quit 2>/dev/null | grep -q .; then
  run_timed slamdunk_all slamdunk all \
    -r "$reference_fasta" -b "$utr_bed" -o "$slamdunk_dir" \
    -5 "$trim5" -t "$threads" -sn "$sample" -sy 4sU "$tissue_fastq"
fi

filtered_bam=$(find "$slamdunk_dir/filter" -type f -name '*_filtered.bam' -size +0c -print -quit)
[[ -n "$filtered_bam" ]] || { echo "ERROR: slamdunk filtered BAM not found" >&2; exit 1; }

if [[ ! -s "$HOME/.gedi/genomic/${gedi_genome}.oml" ]]; then
  run_timed prepare_gedi_genome "$script_dir/prepare_gedi_genome.sh" \
    "$reference_fasta" "$gtf" "$gedi_genome"
fi
grand_prefix="$grand_dir/$sample"
if [[ ! -s "${grand_prefix}.tsv" && ! -s "${grand_prefix}.tsv.gz" && ! -s "${grand_prefix}.ext.tsv" && ! -s "${grand_prefix}.ext.tsv.gz" ]]; then
  run_timed grand_slam "$script_dir/run_grand_slam.sh" \
    "$filtered_bam" "$gedi_genome" "$grand_prefix" \
    exploratory_sample_background "$reference_fasta"
fi

ngm_version=$(ngm --version 2>&1 | head -1 || true)
{
  printf 'tool\tversion\n'
  printf 'slamdunk\t%s\n' "$(slamdunk --version 2>&1 | head -1)"
  printf 'NextGenMap\t%s\n' "$ngm_version"
  printf 'GRAND-SLAM\t2.0.7b\n'
  printf 'Gedi\t1.0.6f\n'
  printf 'java\t%s\n' "$(java -version 2>&1 | head -1)"
} >"$stage_dir/05_metadata/software_versions.tsv"

printf '%s\n' "$input_signature" > "$input_signature_file"
manifest="$stage_dir/05_metadata/stage_output_manifest.tsv"
{
  printf 'relative_path\tbytes\n'
  find "$stage_dir" -type f ! -path "$manifest" -print0 \
    | sort -z \
    | while IFS= read -r -d '' file; do
        printf '%s\t%s\n' "${file#"$stage_dir"/}" "$(stat -c %s "$file")"
      done
} >"$manifest"
