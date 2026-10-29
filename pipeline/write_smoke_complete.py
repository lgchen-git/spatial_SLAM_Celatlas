#!/usr/bin/env python3
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("--sample", required=True)
parser.add_argument("--barcode-table", type=Path, required=True)
parser.add_argument("--filtered-r2", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
for path in (args.barcode_table, args.filtered_r2):
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"missing smoke output: {path}")
payload = {
    "schema_version": 1,
    "ok": True,
    "sample": args.sample,
    "scope": "real_fastq_subset_through_barcode_stage",
    "completed_at_utc": datetime.now(timezone.utc).isoformat(),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
