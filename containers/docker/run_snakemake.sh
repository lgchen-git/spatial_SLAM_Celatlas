#!/usr/bin/env bash
set -euo pipefail

source_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
project_root=$source_root
image=celatlas-spatial-slam:20261023-v1.8.0
config=config/run_example_run.yaml
cores=32
memory_gib=96
outputs_host=${OUTPUTS_HOST:-}
reuse_outputs_host=${REUSE_OUTPUTS_HOST:-}
dry_run=0
detach=0
target=all
force_target=""
rerun_triggers=""

usage() {
  cat <<'EOF'
Usage: run_snakemake.sh [--project-root DIR] [--image TAG] [--config FILE]
                        [--outputs-host DIR] [--reuse-outputs-host DIR]
                        [--cores N] [--memory-gib N]
                        [--target TARGET] [--force-target RULE] [--dry-run] [--detach]
                        [--rerun-triggers MODE]
EOF
}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root) project_root=$2; shift 2 ;;
    --image) image=$2; shift 2 ;;
    --config) config=$2; shift 2 ;;
    --outputs-host) outputs_host=$2; shift 2 ;;
    --reuse-outputs-host) reuse_outputs_host=$2; shift 2 ;;
    --cores) cores=$2; shift 2 ;;
    --memory-gib) memory_gib=$2; shift 2 ;;
    --target) target=$2; shift 2 ;;
    --force-target) force_target=$2; shift 2 ;;
    --rerun-triggers) rerun_triggers=$2; shift 2 ;;
    --dry-run) dry_run=1; shift ;;
    --detach) detach=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

project_root=$(realpath "$project_root")
for name in config inputs log outputs pipeline reference tmp; do
  [[ -d "$project_root/$name" ]] || { echo "ERROR: missing project directory: $project_root/$name" >&2; exit 2; }
done
[[ -f "$project_root/$config" ]] || { echo "ERROR: missing run config: $project_root/$config" >&2; exit 2; }
[[ -f "$project_root/pipeline/Snakefile" ]] || { echo "ERROR: missing pipeline/Snakefile" >&2; exit 2; }
[[ "$cores" =~ ^[1-9][0-9]*$ && "$memory_gib" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: cores and memory must be positive integers" >&2; exit 2; }

inputs=$(realpath "$project_root/inputs")
reference=$(realpath "$project_root/reference")
config_host=$(realpath "$project_root/config")
pipeline_host=$(realpath "$project_root/pipeline")
[[ -z "$outputs_host" ]] || mkdir -p "$outputs_host"
mkdir -p "$project_root/log" "$project_root/tmp/runtime"/{home,tmp,numba,matplotlib,pycache,snakemake}
outputs=$(realpath "${outputs_host:-$project_root/outputs}")
reuse_mount=()
if [[ -n "$reuse_outputs_host" ]]; then
  reuse_outputs=$(realpath "$reuse_outputs_host")
  [[ -d "$reuse_outputs" ]] || { echo "ERROR: reuse outputs directory not found: $reuse_outputs" >&2; exit 2; }
  reuse_mount=(--mount "type=bind,src=$reuse_outputs,dst=/workspace/reuse,readonly")
fi
logs=$(realpath "$project_root/log")
runtime_tmp=$(realpath "$project_root/tmp")
git_commit=$(git -C "$source_root" rev-parse HEAD 2>/dev/null || printf unknown)
git_dirty=false
[[ -z "$(git -C "$source_root" status --porcelain 2>/dev/null)" ]] || git_dirty=true

args=(--snakefile pipeline/Snakefile --configfile "$config" --cores "$cores"
      --config "config_path=$config"
      --resources "mem_mb=$((memory_gib * 1024))" --rerun-incomplete --keep-incomplete --printshellcmds "$target")
[[ "$dry_run" -eq 0 ]] || args+=(--dry-run)
[[ -z "$force_target" ]] || args+=(--forcerun "$force_target")
[[ -z "$rerun_triggers" ]] || args+=(--rerun-triggers "$rerun_triggers")

docker_args=()
if [[ "$detach" -eq 1 ]]; then
  docker_args+=(--detach)
else
  docker_args+=(--rm)
fi

docker run "${docker_args[@]}" --network none \
  --user "$(id -u):$(id -g)" --cpus "$cores" --memory "${memory_gib}g" \
  --ulimit nofile=65536:65536 \
  --workdir /workspace \
  --env USER="$(id -un)" --env LOGNAME="$(id -un)" \
  --env HOME=/workspace/tmp/runtime/home --env TMPDIR=/workspace/tmp/runtime/tmp \
  --env NUMBA_CACHE_DIR=/workspace/tmp/runtime/numba --env MPLCONFIGDIR=/workspace/tmp/runtime/matplotlib \
  --env PYTHONPYCACHEPREFIX=/workspace/tmp/runtime/pycache \
  --env CELATLAS_MODEL=/opt/celatlas-spatial-slam/celatlas/src/swin_tiny.pth \
  --env PIPELINE_GIT_COMMIT="$git_commit" --env PIPELINE_GIT_DIRTY="$git_dirty" --env CELATLAS_SLAM_IMAGE="$image" \
  --mount "type=bind,src=$config_host,dst=/workspace/config,readonly" \
  --mount "type=bind,src=$pipeline_host,dst=/workspace/pipeline,readonly" \
  --mount "type=bind,src=$inputs,dst=/workspace/inputs,readonly" \
  --mount "type=bind,src=$reference,dst=/workspace/reference,readonly" \
  --mount "type=bind,src=$logs,dst=/workspace/log" \
  --mount "type=bind,src=$outputs,dst=/workspace/outputs" \
  --mount "type=bind,src=$runtime_tmp,dst=/workspace/tmp" \
  --mount "type=bind,src=$runtime_tmp/runtime/snakemake,dst=/workspace/.snakemake" \
  "${reuse_mount[@]}" \
  "$image" "${args[@]}"
