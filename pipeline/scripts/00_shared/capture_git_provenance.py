#!/usr/bin/env python3
"""Record the exact pipeline Git state used by one result stage."""
from __future__ import annotations

import argparse
import hashlib
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def git(repo: Path, *args: str, binary: bool = False):
    result = subprocess.run(
        ["git", "-C", str(repo), *args], check=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=not binary,
    )
    return result.stdout


def digest(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    commit = git(repo, "rev-parse", "--verify", "HEAD").strip()
    try:
        branch = git(repo, "symbolic-ref", "--quiet", "--short", "HEAD").strip()
    except subprocess.CalledProcessError:
        branch = "DETACHED"
    status = git(repo, "status", "--porcelain=v1", "--untracked-files=all")
    diff = git(repo, "diff", "--binary", "HEAD", "--", binary=True)
    describe = git(repo, "describe", "--always", "--dirty", "--tags").strip()
    try:
        remote = git(repo, "remote", "get-url", "origin").strip().split("@")[-1]
    except subprocess.CalledProcessError:
        remote = ""
    tags = git(repo, "tag", "--points-at", "HEAD").splitlines()
    release = next((x for x in sorted(tags) if x.startswith("pipeline-v")), "")
    rows = (
        ("schema_version", "1"),
        ("captured_at_utc", datetime.now(timezone.utc).isoformat()),
        ("repo_name", repo.name), ("commit", commit), ("branch", branch),
        ("describe", describe), ("dirty", str(bool(status)).lower()),
        ("tracked_diff_sha256", digest(diff)),
        ("status_sha256", digest(status.encode())),
        ("changed_path_count", str(len(status.splitlines()))),
        ("remote_url", remote), ("pipeline_release_tag", release),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("key\tvalue\n" + "".join(f"{k}\t{v}\n" for k, v in rows), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
