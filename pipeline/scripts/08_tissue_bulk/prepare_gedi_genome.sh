#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <genome.fa> <genes.gtf> <genomic-name>" >&2
  exit 2
fi

fasta=$1
gtf=$2
name=$3
command -v gedi >/dev/null 2>&1 || { echo 'ERROR: gedi is not on PATH' >&2; exit 127; }
mkdir -p "$HOME/.gedi/genomic"

# IndexGenome receives explicit files, so it cannot silently substitute an
# Ensembl release that differs from the STAR/slamdunk reference.
gedi -e IndexGenome -s "$fasta" -a "$gtf" -n "$name" -nobowtie -nostar -nokallisto
