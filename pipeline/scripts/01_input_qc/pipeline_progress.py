#!/usr/bin/env python3
"""Framework-agnostic stage progress and resource recorder."""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def runtime_dir(result_dir: Path) -> Path:
    path = result_dir / "01_input_qc" / "05_metadata" / "progress_runtime"
    path.mkdir(parents=True, exist_ok=True)
    return path


def disk_bytes(path: Path) -> int:
    total = 0
    for root, _, files in os.walk(path):
        for filename in files:
            try:
                total += (Path(root) / filename).stat().st_size
            except FileNotFoundError:
                pass
    return total


def append_event(result_dir: Path, event: dict) -> None:
    event["timestamp"] = now()
    with (runtime_dir(result_dir) / "stage_metrics.jsonl").open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(event, ensure_ascii=False, sort_keys=True) + "\n")


def parse_time(path: Path) -> dict:
    values: dict[str, float | int | str] = {}
    if not path.is_file():
        return values
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.lstrip().startswith("Elapsed (wall clock)"):
            values["wall_clock"] = line.split("):", 1)[-1].strip()
            continue
        key, _, value = line.partition(":")
        key = key.strip()
        value = value.strip()
        if key.startswith("User time"):
            values["user_cpu_seconds"] = float(value)
        elif key.startswith("System time"):
            values["system_cpu_seconds"] = float(value)
        elif key.startswith("Maximum resident set size"):
            values["max_rss_kb"] = int(value)
    if "user_cpu_seconds" in values or "system_cpu_seconds" in values:
        values["cpu_seconds"] = round(float(values.get("user_cpu_seconds", 0)) + float(values.get("system_cpu_seconds", 0)), 3)
    return values


def wall_clock_seconds(value: str | None) -> float | None:
    """Convert GNU time's h:mm:ss or m:ss wall-clock value to seconds."""
    if not value:
        return None
    try:
        parts = [float(part) for part in value.split(":")]
    except ValueError:
        return None
    if len(parts) == 3:
        return parts[0] * 3600 + parts[1] * 60 + parts[2]
    if len(parts) == 2:
        return parts[0] * 60 + parts[1]
    return parts[0] if len(parts) == 1 else None


def load_manifest(path: Path) -> dict:
    if not path.is_file():
        return {}
    try:
        import yaml  # type: ignore
        return yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except Exception:
        result: dict[str, object] = {"stages": []}
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith("pipeline:"):
                result["pipeline"] = line.split(":", 1)[1].strip().strip('"')
            elif "{" in line and "id:" in line and "name:" in line:
                stage_id = re.search(r"id:\s*[\"']?([^,}\"']+)", line)
                name = re.search(r"name:\s*[\"']?([^}\"']+)", line)
                if stage_id and name:
                    result["stages"].append({"id": stage_id.group(1).strip(), "name": name.group(1).strip()})
        return result


def load_events(result_dir: Path) -> list[dict]:
    path = runtime_dir(result_dir) / "stage_metrics.jsonl"
    if not path.is_file():
        return []
    events = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        try:
            events.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return events


def fmt_bytes(value: int | None) -> str:
    if value is None:
        return "未采集"
    units = ["B", "KB", "MB", "GB", "TB"]
    amount = float(value)
    for unit in units:
        if amount < 1024 or unit == units[-1]:
            return f"{amount:.1f} {unit}"
        amount /= 1024
    return "未采集"


def fmt_seconds(value: float | int | None) -> str:
    if value is None:
        return "未采集"
    seconds = int(round(float(value)))
    hours, seconds = divmod(seconds, 3600)
    minutes, seconds = divmod(seconds, 60)
    return f"{hours}h {minutes}m {seconds}s" if hours else f"{minutes}m {seconds}s"


def command_init(args: argparse.Namespace) -> int:
    result_dir = Path(args.result_dir).resolve()
    info = {"pipeline": args.pipeline, "run_id": args.run_id or result_dir.name, "manifest": args.manifest or "", "started_at": now()}
    target = runtime_dir(result_dir) / "run.json"
    if not target.exists():
        target.write_text(json.dumps(info, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return 0


def command_event(args: argparse.Namespace) -> int:
    result_dir = Path(args.result_dir).resolve()
    event = {"event": args.event, "stage_id": args.stage_id, "stage_name": args.stage_name, "threads": args.threads, "log": args.log or ""}
    if args.event == "start":
        event["disk_bytes_before"] = disk_bytes(result_dir)
    else:
        event.update(parse_time(Path(args.time_file)) if args.time_file else {})
        before = args.disk_before
        after = disk_bytes(result_dir)
        event["disk_bytes_before"] = before
        event["disk_bytes_after"] = after
        event["disk_bytes_delta"] = after - before if before is not None else None
        event["exit_code"] = args.exit_code
    append_event(result_dir, event)
    return 0


def command_status(args: argparse.Namespace) -> int:
    result_dir = Path(args.result_dir).resolve()
    manifest = load_manifest(Path(args.manifest)) if args.manifest else load_manifest(result_dir / "progress_manifest.yaml")
    stages = manifest.get("stages", []) if isinstance(manifest, dict) else []
    events = load_events(result_dir)
    by_stage: dict[str, dict] = {}
    for event in events:
        by_stage[str(event.get("stage_id"))] = event
    print(f"## 流程进展：{manifest.get('pipeline', result_dir.name) if isinstance(manifest, dict) else result_dir.name}")
    print(f"结果目录：`{result_dir}`\n")
    if events:
        latest_event = events[-1]
        state = {"start": "运行中", "finish": "完成", "fail": "失败"}.get(latest_event.get("event"), "未知")
        print(f"当前阶段：{latest_event.get('stage_id', '?')} {latest_event.get('stage_name', '')}（{state}）\n")
    else:
        print("当前阶段：未采集运行事件；请检查历史日志或在下次运行时接入监控。\n")
    print("| 阶段 | 状态 | 耗时 | CPU 核时 | 峰值内存 | 磁盘增量 |")
    print("|---|---|---:|---:|---:|---:|")
    for stage in stages:
        stage_id = str(stage.get("id"))
        latest = by_stage.get(stage_id, {})
        status = {"start": "运行中", "finish": "完成", "fail": "失败"}.get(latest.get("event"), "等待")
        if latest.get("event") == "finish" and latest.get("exit_code") not in (None, 0):
            status = "失败"
        cpu = latest.get("cpu_seconds")
        rss = latest.get("max_rss_kb")
        disk = latest.get("disk_bytes_delta")
        print(f"| {stage_id} {stage.get('name', '')} | {status} | {latest.get('wall_clock', '未采集')} | {fmt_seconds(cpu)} | {fmt_bytes(int(rss) * 1024) if rss is not None else '未采集'} | {fmt_bytes(disk) if isinstance(disk, int) else '未采集'} |")
    if not stages:
        print("\n未找到 manifest；可读取 `01_input_qc/05_metadata/progress_runtime/stage_metrics.jsonl`，但无法生成预定义阶段表。")
    return 0


def command_resources(args: argparse.Namespace) -> int:
    """Collect detailed GNU-time and sampled resource records into one TSV."""
    result_dir = Path(args.result_dir).resolve()
    rows: list[dict[str, object]] = []
    recorded_steps: set[tuple[str, str]] = set()
    for path in sorted(result_dir.rglob("*.time.txt")):
        values = parse_time(path)
        if not values:
            continue
        component = path.relative_to(result_dir).parts[0]
        step = path.name.removesuffix(".time.txt")
        rows.append({
            "source_type": "gnu_time",
            "component": component,
            "step": step,
            "wall_clock": values.get("wall_clock", ""),
            "wall_seconds": wall_clock_seconds(str(values.get("wall_clock", ""))) or "",
            "cpu_seconds": values.get("cpu_seconds", ""),
            "peak_cpu_percent": "",
            "max_rss_gb": round(int(values.get("max_rss_kb", 0)) / 1024 / 1024, 3) if values.get("max_rss_kb") else "",
            "source": str(path),
        })
        recorded_steps.add((component, step))
    for source in args.sampled_resource_tsv:
        path = Path(source).resolve()
        if not path.is_file():
            continue
        with path.open(encoding="utf-8", newline="") as handle:
            for sampled in csv.DictReader(handle, delimiter="\t"):
                component = path.relative_to(result_dir).parts[0] if path.is_relative_to(result_dir) else "external"
                step = sampled.get("substage", "") or sampled.get("step", "")
                if (component, step) in recorded_steps:
                    continue
                wall_clock = sampled.get("elapsed", "")
                wall_seconds = sampled.get("observed_wall_seconds", "") or wall_clock_seconds(wall_clock) or ""
                cpu_seconds: float | str = ""
                if sampled.get("user_seconds") or sampled.get("system_seconds"):
                    cpu_seconds = round(float(sampled.get("user_seconds", 0) or 0) + float(sampled.get("system_seconds", 0) or 0), 3)
                max_rss_gb: float | str = sampled.get("peak_rss_gb", "")
                if not max_rss_gb and sampled.get("max_rss_kb"):
                    max_rss_gb = round(float(sampled["max_rss_kb"]) / 1024 / 1024, 3)
                rows.append({
                    "source_type": "process_sampling" if sampled.get("substage") else "resource_table",
                    "component": component,
                    "step": step,
                    "wall_clock": wall_clock,
                    "wall_seconds": wall_seconds,
                    "cpu_seconds": cpu_seconds,
                    "peak_cpu_percent": sampled.get("peak_cpu_percent", "") or sampled.get("cpu_percent", "").rstrip("%"),
                    "max_rss_gb": max_rss_gb,
                    "source": str(path),
                })
                recorded_steps.add((component, step))
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["source_type", "component", "step", "wall_clock", "wall_seconds", "cpu_seconds", "peak_cpu_percent", "max_rss_gb", "source"]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(output)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    init = sub.add_parser("init")
    init.add_argument("--result-dir", required=True); init.add_argument("--pipeline", required=True)
    init.add_argument("--run-id"); init.add_argument("--manifest"); init.set_defaults(func=command_init)
    event = sub.add_parser("event")
    event.add_argument("--result-dir", required=True); event.add_argument("--event", choices=("start", "finish", "fail"), required=True)
    event.add_argument("--stage-id", required=True); event.add_argument("--stage-name", required=True); event.add_argument("--threads", type=int, default=1)
    event.add_argument("--log"); event.add_argument("--time-file"); event.add_argument("--disk-before", type=int); event.add_argument("--exit-code", type=int)
    event.set_defaults(func=command_event)
    status = sub.add_parser("status")
    status.add_argument("--result-dir", required=True); status.add_argument("--manifest"); status.set_defaults(func=command_status)
    resources = sub.add_parser("resources", help="Consolidate detailed resource logs into one TSV")
    resources.add_argument("--result-dir", required=True); resources.add_argument("--output", required=True)
    resources.add_argument("--sampled-resource-tsv", action="append", default=[]); resources.set_defaults(func=command_resources)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
