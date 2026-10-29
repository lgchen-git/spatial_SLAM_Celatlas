#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
python=${PYTHON:-python3}
smoke_tmp=${SMOKE_TMPDIR:-"$root/tmp/smoke"}
export TMPDIR=${TMPDIR:-"$smoke_tmp/runtime"}
export PYTHONPYCACHEPREFIX=${PYTHONPYCACHEPREFIX:-"$smoke_tmp/pycache"}
export NUMBA_CACHE_DIR=${NUMBA_CACHE_DIR:-"$smoke_tmp/numba"}
export MPLCONFIGDIR=${MPLCONFIGDIR:-"$smoke_tmp/matplotlib"}
mkdir -p "$TMPDIR" "$PYTHONPYCACHEPREFIX" "$NUMBA_CACHE_DIR" "$MPLCONFIGDIR"

"$python" "$root/pipeline/scripts/01_input_qc/check_samples.py" \
  --samples "$root/config/samples_example_run.tsv" \
  --libraries "$root/config/libraries_example_run.tsv" --sample example_sample
bash -n "$root/containers/docker/build_export_image.sh" "$root/containers/docker/run_snakemake.sh" \
  "$root/pipeline/scripts/07_spatial_dynamics/run_bin50_showcase.sh" \
  "$root/pipeline/scripts/08_tissue_bulk/run_grand_slam.sh" \
  "$root/pipeline/scripts/08_tissue_bulk/prepare_gedi_genome.sh"
mapfile -d '' python_sources < <(find "$root/pipeline" -type f -name '*.py' -print0)
"$python" -m py_compile "${python_sources[@]}"
"$python" "$root/tests/test_core_logic.py"
"$python" "$root/tests/test_pipeline_contract.py"
"$python" "$root/tests/test_build_read_flow_sankey.py"
"$python" "$root/tests/test_snakemake_delivery.py"
layout_auditor=/opt/project/.codex/skills/bio-pipeline-layout/scripts/audit_analysis_layout.py
if [[ -f "$layout_auditor" ]]; then
  "$python" "$layout_auditor" --config "$root/config/analysis_layout.yaml" \
    --project-root "$root" --outputs-root "$root/outputs" --json
fi
