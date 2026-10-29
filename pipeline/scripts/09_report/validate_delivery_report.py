#!/usr/bin/env python3
"""Validate the offline report and delivery manifest (Bio Report v3 contract)."""
import argparse
import base64
import csv
import json
import re
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--require-section", action="append", default=[])
    parser.add_argument("--forbid-text", action="append", default=[])
    parser.add_argument("--require-expandable-figures", action="store_true")
    parser.add_argument("--require-print-css", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    report, manifest = args.report.resolve(), args.manifest.resolve()
    errors = []
    if not report.is_file() or not report.stat().st_size:
        errors.append(f"report missing or empty: {report}")
        text = ""
    else:
        text = report.read_text(encoding="utf-8", errors="replace")
        markup = re.sub(r"<script\b[^>]*>.*?</script>", "", text, flags=re.I | re.S)
        hrefs = re.findall(r'\bhref="([^"]+)"', markup)
        srcs = re.findall(r'\bsrc="([^"]+)"', markup)
        if any(not value.startswith("#") for value in hrefs):
            errors.append("report contains external or file href values")
        if any(not value.startswith("data:") for value in srcs):
            errors.append("report contains non-data src values")
        ids = set(re.findall(r'\bid="([^"]+)"', markup))
        stale = [value for value in hrefs if value.startswith("#") and value[1:] not in ids]
        if stale:
            errors.append(f"report contains stale anchors: {', '.join(stale[:5])}")
        sidebar = re.search(r'<nav\b([^>]*)class="[^"]*\bsidebar\b[^"]*"[^>]*>(.*?)</nav>', text, flags=re.I | re.S)
        if not sidebar:
            errors.append("report sidebar is missing")
        else:
            section_ids = [value for attrs, value in re.findall(
                r'<section\b([^>]*)id="([^"]+)"[^>]*>', text, flags=re.I
            ) if re.search(r'\bclass="[^"]*\breport-section\b', attrs, flags=re.I)]
            nav_ids = re.findall(r'href="#([^"]+)"', sidebar.group(2), flags=re.I)
            if section_ids and nav_ids != section_ids:
                errors.append("sidebar anchors must match report sections exactly and in order")
        if 'class="footer"' not in text:
            errors.append("report footer is missing")
        if args.require_expandable_figures:
            if ('class="expandable"' not in text or "figure-card" not in text or
                    "openFigure" not in text or "Escape" not in text):
                errors.append("expandable figure controls are incomplete")
        if args.require_print_css and "@media print" not in text:
            errors.append("report print CSS is missing")
        for uri in re.findall(r'\bsrc="(data:[^"]+)"', text):
            if ";base64," not in uri:
                errors.append("report contains a non-base64 data URI")
                continue
            try:
                if not base64.b64decode(uri.split(";base64,", 1)[1], validate=True):
                    errors.append("report contains an empty data URI")
            except ValueError:
                errors.append("report contains an invalid data URI")
        for index, (attrs, body) in enumerate(
            re.findall(r'<table\b([^>]*)>(.*?)</table>', text, flags=re.I | re.S), start=1
        ):
            bodies = re.findall(r'<tbody\b[^>]*>(.*?)</tbody>', body, flags=re.I | re.S)
            row_count = sum(len(re.findall(r'<tr\b', item, flags=re.I)) for item in bodies)
            if row_count > 20 and 'data-page-size' not in attrs.lower() and 'paginated-table' not in attrs.lower():
                errors.append(f"table {index} has {row_count} rows without pagination")
            if 'data-page-size' in attrs.lower() or 'paginated-table' in attrs.lower():
                table_id = re.search(r'\bid="([^"]+)"', attrs, flags=re.I)
                pager = (re.search(
                    rf'<div\b(?=[^>]*data-pager-for="{re.escape(table_id.group(1))}")[^>]*>(.*?)</div>',
                    text, flags=re.I | re.S
                ) if table_id else None)
                if not table_id or not pager or any(control not in pager.group(1) for control in ("data-page-prev", "data-page-next", "data-page-status")):
                    errors.append(f"paginated table {index} has incomplete pager")
        for section in args.require_section:
            if section not in text:
                errors.append(f"required section is missing: {section}")
        for term in args.forbid_text:
            if term and term in text:
                errors.append(f"forbidden text found in report: {term}")
    rows = []
    if not manifest.is_file() or not manifest.stat().st_size:
        errors.append(f"manifest missing or empty: {manifest}")
    else:
        with manifest.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        for row in rows:
            if row.get("status") in {"copied", "up_to_date"}:
                target = Path(row.get("target", ""))
                if not target.is_absolute():
                    target = manifest.parent / target
                target = target.resolve()
                if not target.is_file() or not target.stat().st_size:
                    errors.append(f"manifest target missing: {target}")
                elif row.get("bytes") and int(row["bytes"]) != target.stat().st_size:
                    errors.append(f"manifest byte size mismatch: {target}")

        def resolved_target(row):
            target = Path(row.get("target", ""))
            if not target.is_absolute():
                target = manifest.parent / target
            return target.resolve()

        if sum(resolved_target(row) == report for row in rows) != 1:
            errors.append("manifest must contain exactly one final-report target")
        elif not rows or resolved_target(rows[-1]) != report:
            errors.append("final-report target must be the last manifest row")
    payload = {"report": str(report), "manifest_rows": len(rows), "data_uri_count": len(re.findall(r'\bsrc="data:', text)), "errors": errors, "ok": not errors}
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        print("OK" if not errors else "FAILED")
        for error in errors:
            print(f"ERROR: {error}")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
