#!/usr/bin/env python3
"""Portable subset of GNU time -v used when /usr/bin/time is unavailable."""
from __future__ import annotations

import argparse
import resource
import subprocess
import time
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("missing command")
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    completed = subprocess.run(command)
    elapsed = time.monotonic() - started
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    user = after.ru_utime - before.ru_utime
    system = after.ru_stime - before.ru_stime
    cpu = 100.0 * (user + system) / elapsed if elapsed else 0.0
    minutes, seconds = divmod(elapsed, 60)
    hours, minutes = divmod(int(minutes), 60)
    elapsed_text = f"{hours}:{minutes:02d}:{seconds:05.2f}" if hours else f"{minutes}:{seconds:05.2f}"
    Path(args.output).write_text(
        f"User time (seconds): {user:.6f}\n"
        f"System time (seconds): {system:.6f}\n"
        f"Percent of CPU this job got: {cpu:.1f}%\n"
        f"Elapsed (wall clock) time (h:mm:ss or m:ss): {elapsed_text}\n"
        f"Maximum resident set size (kbytes): {after.ru_maxrss}\n"
        f"File system inputs: {after.ru_inblock - before.ru_inblock}\n"
        f"File system outputs: {after.ru_oublock - before.ru_oublock}\n"
        f"Exit status: {completed.returncode}\n",
        encoding="utf-8",
    )
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
