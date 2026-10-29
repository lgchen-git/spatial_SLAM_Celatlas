#!/usr/bin/env python3
"""Adapt authoritative Celatlas read-pair metrics to one conserved flow."""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--celatlas-data", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--plotly-js", required=True, type=Path)
    parser.add_argument("--bin-label", default="bin50")
    return parser.parse_args()


def load_json(path):
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f"Missing metric source: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def integer(mapping, key, source):
    value = mapping.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or int(value) != value or value < 0:
        raise ValueError(f"{source}: {key} is not a non-negative integer: {value!r}")
    return int(value)


def summary_mapping(payload, name):
    summary = payload.get(name, {})
    if isinstance(summary, dict) and "metric_list" in summary:
        return {item.get("name"): item.get("value") for item in summary["metric_list"]}
    return summary


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def node(view, node_id, stage, label, value, unit, kind, scope="sample_library"):
    return {"view_id": view, "node_id": node_id, "stage": stage, "label": label,
            "value": value, "unit": unit, "scope": scope, "node_kind": kind}


def link(view, link_id, source, target, value, reason):
    return {"view_id": view, "link_id": link_id, "source": source, "target": target,
            "value": value, "reason": reason}


def main():
    args = arguments()
    celatlas = load_json(args.celatlas_data)
    barcode = summary_mapping(celatlas, "barcode_summary")
    cutadapt = summary_mapping(celatlas, "cutadapt_summary")
    star = summary_mapping(celatlas, "star_summary")
    feature = summary_mapping(celatlas, "featureCounts_summary")

    raw = integer(barcode, "Raw Reads", "Celatlas barcode_summary")
    exact = integer(barcode, "Exactly Matched Reads", "Celatlas barcode_summary")
    corrected = integer(barcode, "Mismatched Reads", "Celatlas barcode_summary")
    valid = integer(barcode, "Valid Reads", "Celatlas barcode_summary")
    no_barcode = integer(barcode, "No Barcode Reads", "Celatlas barcode_summary") if "No Barcode Reads" in barcode else raw - valid
    if no_barcode < 0:
        raise ValueError("Celatlas barcode_summary has more valid reads than raw reads")
    if exact + corrected != valid or valid + no_barcode != raw:
        raise ValueError("Celatlas barcode_summary is not conserved")

    trimmed = integer(cutadapt, "Reads Written", "Celatlas cutadapt_summary")
    trim_short = integer(cutadapt, "Reads Too Short", "Celatlas cutadapt_summary")
    if trimmed + trim_short != valid:
        raise ValueError("Celatlas cutadapt_summary is not conserved")

    star_fates = {
        "unique_mapped": integer(star, "Uniquely Mapped Reads", "Celatlas star_summary"),
        "multi_mapped": integer(star, "Multi-Mapped Reads", "Celatlas star_summary"),
        "too_many_loci": integer(star, "Reads Mapped To Too Many Loci", "Celatlas star_summary"),
        "unmapped_too_short": integer(star, "Unmapped Reads - Too Short", "Celatlas star_summary"),
        "unmapped_other": integer(star, "Unmapped Reads - Other", "Celatlas star_summary"),
    }
    if sum(star_fates.values()) != trimmed:
        raise ValueError("Celatlas star_summary is not conserved")

    feature_fates = {
        "exonic": integer(feature, "Reads Assigned To Exonic Regions", "Celatlas featureCounts_summary"),
        "intronic": integer(feature, "Reads Assigned To Intronic Regions", "Celatlas featureCounts_summary"),
        "intergenic": integer(feature, "Reads Assigned To Intergenic Regions", "Celatlas featureCounts_summary"),
        "ambiguity": integer(feature, "Reads Unassigned Ambiguity", "Celatlas featureCounts_summary"),
    }
    if sum(feature_fates.values()) != star_fates["unique_mapped"]:
        raise ValueError("Celatlas featureCounts_summary is not conserved")

    reads = "celatlas_read_pairs"
    nodes = [
        node(reads, "raw", 0, "Raw FASTQ read pairs", raw, "read_pairs", "retained"),
        node(reads, "barcode_exact", 1, "Exact spatial barcode", exact, "read_pairs", "branch"),
        node(reads, "barcode_corrected", 1, "Corrected spatial barcode", corrected, "read_pairs", "branch"),
        node(reads, "no_barcode", 1, "No valid spatial barcode", no_barcode, "read_pairs", "loss"),
        node(reads, "barcode_valid", 2, "Celatlas barcode-valid", valid, "read_pairs", "retained"),
        node(reads, "trimmed", 3, "Cutadapt retained", trimmed, "read_pairs", "retained"),
        node(reads, "trim_too_short", 3, "Cutadapt too short", trim_short, "read_pairs", "loss"),
    ]
    read_labels = {
        "unique_mapped": "STAR uniquely mapped", "multi_mapped": "STAR multi-mapped",
        "too_many_loci": "STAR mapped to too many loci",
        "unmapped_too_short": "STAR unmapped: too short", "unmapped_other": "STAR unmapped: other",
        "exonic": "Assigned exonic", "intronic": "Assigned intronic",
        "intergenic": "Assigned intergenic", "ambiguity": "Unassigned ambiguity",
    }
    nodes += [node(reads, key, 4, read_labels[key], value, "read_pairs", "retained" if key == "unique_mapped" else "loss") for key, value in star_fates.items()]
    nodes += [node(reads, key, 5, read_labels[key], value, "read_pairs", "branch" if key in {"exonic", "intronic"} else "loss") for key, value in feature_fates.items()]
    links = [
        link(reads, "raw_to_exact", "raw", "barcode_exact", exact, "exact Celatlas spatial-barcode match"),
        link(reads, "raw_to_corrected", "raw", "barcode_corrected", corrected, "Celatlas-corrected spatial-barcode match"),
        link(reads, "raw_to_no_barcode", "raw", "no_barcode", no_barcode, "no valid Celatlas spatial barcode"),
        link(reads, "exact_to_valid", "barcode_exact", "barcode_valid", exact, "retained as barcode-valid"),
        link(reads, "corrected_to_valid", "barcode_corrected", "barcode_valid", corrected, "retained after barcode correction"),
        link(reads, "valid_to_trimmed", "barcode_valid", "trimmed", trimmed, "written by Celatlas cutadapt"),
        link(reads, "valid_to_trim_short", "barcode_valid", "trim_too_short", trim_short, "too short after adapter/quality trimming"),
    ]
    links += [link(reads, f"trimmed_to_{key}", "trimmed", key, value, read_labels[key]) for key, value in star_fates.items()]
    links += [link(reads, f"unique_to_{key}", "unique_mapped", key, value, read_labels[key]) for key, value in feature_fates.items()]
    inputs = args.output_dir / "inputs"
    write_tsv(inputs / "read_flow_nodes.input.tsv", ["view_id", "node_id", "stage", "label", "value", "unit", "scope", "node_kind"], nodes)
    write_tsv(inputs / "read_flow_links.input.tsv", ["view_id", "link_id", "source", "target", "value", "reason"], links)
    write_tsv(inputs / "read_flow_views.tsv", ["view_id", "display_name", "order", "is_default"], [
        {"view_id": reads, "display_name": "Celatlas FASTQ → feature assignment", "order": 0, "is_default": 1},
    ])
    stale_transformations = inputs / "read_flow_transformations.input.tsv"
    if stale_transformations.exists():
        stale_transformations.unlink()

    # Renderer output is view-keyed. Remove only stale module SVGs so a changed
    # view registry cannot leave extra figures that fail delivery validation.
    figure_dir = args.output_dir / "figures" / "read_flow_sankey"
    if figure_dir.is_dir():
        for old_svg in figure_dir.glob("*.svg"):
            old_svg.unlink()
    renderer = Path(__file__).with_name("build_read_flow_sankey.py")
    command = [sys.executable, str(renderer), "--nodes", str(inputs / "read_flow_nodes.input.tsv"), "--links", str(inputs / "read_flow_links.input.tsv"), "--views", str(inputs / "read_flow_views.tsv"), "--default-mode", "view", "--node-layout", "main_top_losses", "--orientation", "horizontal", "--output-dir", str(args.output_dir), "--plotly-js", str(args.plotly_js), "--title", f"{args.sample} conserved workflow flow", "--bin-label", args.bin_label]
    subprocess.run(command, check=True)
    subprocess.run([sys.executable, str(Path(__file__).with_name("validate_read_flow_outputs.py")), "--output-dir", str(args.output_dir)], check=True)


if __name__ == "__main__":
    main()
