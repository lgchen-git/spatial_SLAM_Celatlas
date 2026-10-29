#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 <bam> <genomic-name> <output-prefix> <error-mode> <reference-fasta>" >&2
  exit 2
fi

bam=$1
genomic_name=$2
prefix=$3
error_mode=$4
reference_fasta=$5

command -v gedi >/dev/null 2>&1 || { echo "ERROR: gedi/GRAND-SLAM is not on PATH" >&2; exit 127; }
mkdir -p "$(dirname "$prefix")"

# GRAND-SLAM can read BAM directly. The BAM must use the same reference as the
# prepared Gedi genome and expose MD/NH tags. no4sU controls are absent for this
# project, so an exploratory run is deliberately annotated in its metadata.
if [[ "$error_mode" != "exploratory_sample_background" ]]; then
  echo "ERROR: unsupported GRAND-SLAM error mode: $error_mode" >&2
  exit 2
fi
gedi -e Slam -genomic "$genomic_name" -prefix "$prefix" -progress -plot -D -reads "$bam"
cat > "${prefix}.run_metadata.json" <<EOF
{"error_mode":"exploratory_sample_background","reference_fasta":"$reference_fasta","interpretation":"NTR is exploratory because no no4sU control or validated errlm was supplied."}
EOF
