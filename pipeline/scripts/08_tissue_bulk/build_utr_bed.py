#!/usr/bin/env python3
"""Export explicit three_prime_utr records from a reference-matched GTF."""
import argparse
import re
from pathlib import Path


def attr(value, key):
    match = re.search(rf'{re.escape(key)} "([^"]+)"', value)
    return match.group(1) if match else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    rows = []
    with open(args.gtf, encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2].lower() not in {"three_prime_utr", "3utr", "three_prime_utr_region"}:
                continue
            gene = attr(fields[8], "gene_id") or attr(fields[8], "transcript_id")
            if gene:
                rows.append((fields[0], int(fields[3]) - 1, fields[4], gene, 0, fields[6]))
    if not rows:
        raise SystemExit("No explicit three_prime_utr features found in GTF; refusing an exon fallback for slamdunk.")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8") as handle:
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
