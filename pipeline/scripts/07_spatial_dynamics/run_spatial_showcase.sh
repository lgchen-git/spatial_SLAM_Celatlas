#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
SCRIPTS="$ROOT/pipeline/scripts"
export NUMBA_DISABLE_JIT=${NUMBA_DISABLE_JIT:-0}
export NUMBA_DISABLE_CACHING=${NUMBA_DISABLE_CACHING:-1}
export NUMBA_CACHE_DIR=${NUMBA_CACHE_DIR:-$ROOT/tmp/celatlas-numba-cache}
export MPLCONFIGDIR=${MPLCONFIGDIR:-$ROOT/tmp/celatlas-matplotlib}
mkdir -p "$NUMBA_CACHE_DIR" "$MPLCONFIGDIR"

if [[ $# -lt 6 ]]; then
  echo "Usage: $0 SAMPLE SAMPLE_RESULT_DIR COORDINATE_UM OFFICIAL_BIN_DIR SOURCE_MATRIX_DIR SOURCE_MOLECULES [--force]" >&2
  exit 2
fi
SAMPLE=$1
OUT=$2
COORDINATE_UM=$3
CELATLAS_BIN_DIR=$4
SOURCE_MATRIX_DIR=$5
SOURCE_MOLECULES=$6
FORCE=${7:-}
[[ -d "$CELATLAS_BIN_DIR" ]] || { echo "ERROR: official Celatlas bin directory not found: $CELATLAS_BIN_DIR" >&2; exit 1; }
[[ -d "$SOURCE_MATRIX_DIR" ]] || { echo "ERROR: source matrix directory not found: $SOURCE_MATRIX_DIR" >&2; exit 1; }
[[ -s "$SOURCE_MOLECULES" ]] || { echo "ERROR: source molecule table not found: $SOURCE_MOLECULES" >&2; exit 1; }
BIN_SIZE_RAW=${BIN_SIZE_RAW:-50}
BIN_SIZE_UM=${BIN_SIZE_UM:-50}
BIN_SIZE_PIXELS=${BIN_SIZE_PIXELS:-100}
BIN_LABEL=bin${BIN_SIZE_RAW}
SPATIAL_HVG=${SPATIAL_HVG:-3000}
SPATIAL_HVG_FLAVOR=${SPATIAL_HVG_FLAVOR:-seurat}
SPATIAL_PCA=${SPATIAL_PCA:-25}
SPATIAL_NEIGHBORS=${SPATIAL_NEIGHBORS:-15}
LEIDEN_RES=${LEIDEN_RES:-1.0}
LEIDEN_SEED=${LEIDEN_SEED:-0}
LEIDEN_ITER=${LEIDEN_ITER:-2}
SPATIAL_POINT=${SPATIAL_POINT:-14}
UMAP_POINT=${UMAP_POINT:-10}
REP_BG_POINT=${REP_BG_POINT:-8}
REP_POINT=${REP_POINT:-12}
SMOOTHING_NEIGHBORS=${SMOOTHING_NEIGHBORS:-24}
SMOOTHING_MAX_RADIUS=${SMOOTHING_MAX_RADIUS:-}
INTERPRETATION_MODE=${INTERPRETATION_MODE:-standard}
REVIEWED_MODE=${REVIEWED_MODE:-0}
TC_COLUMN=${TC_COLUMN:-tc_events}
T_COLUMN=${T_COLUMN:-t_covered}
DEFAULT_BASE_QUALITY=${DEFAULT_BASE_QUALITY:-20}
BASE_QUALITY_THRESHOLDS=${BASE_QUALITY_THRESHOLDS:-15,20,25,30}
GRAND_GENE_NTR=${GRAND_GENE_NTR:-}
[[ "$INTERPRETATION_MODE" == "standard" || "$INTERPRETATION_MODE" == "unlabeled_control" ]] || {
  echo "ERROR: INTERPRETATION_MODE must be standard or unlabeled_control" >&2
  exit 2
}
SHOW="$OUT/07_spatial_dynamics"
TABLES="$SHOW/01_tables"
FIGURES="$SHOW/02_figures"
LOGS="$SHOW/03_logs"
TOOLS="$SHOW/04_tool_output"
META="$SHOW/05_metadata"
mkdir -p "$TABLES" "$FIGURES" "$LOGS" "$TOOLS/scanpy" "$TOOLS/correction" "$TOOLS/kinetics" "$META"

run_timed() {
  local output=$1
  shift
  if [[ -x /usr/bin/time ]]; then
    /usr/bin/time -v -o "$output" "$@"
  else
    python3 "$SCRIPTS/00_shared/run_with_resources.py" --output "$output" -- "$@"
  fi
}

run_scanpy() {
  run_timed "$LOGS/scanpy_total_${BIN_LABEL}.time.txt" \
    python3 "$SCRIPTS/07_spatial_dynamics/spatial_slam_scanpy_analysis.py" \
      --sample "$SAMPLE" --matrix-dir "$SOURCE_MATRIX_DIR" --outdir "$TOOLS/scanpy" \
      --celatlas-bin-dir "$CELATLAS_BIN_DIR" --bin-size "$BIN_SIZE_RAW" --bin-size-um "$BIN_SIZE_UM" --bin-size-pixels "$BIN_SIZE_PIXELS" --hvg "$SPATIAL_HVG" --pca-components "$SPATIAL_PCA" \
      --neighbors "$SPATIAL_NEIGHBORS" --leiden-resolution "$LEIDEN_RES" --leiden-random-state "$LEIDEN_SEED" --leiden-n-iterations "$LEIDEN_ITER" \
      --hvg-flavor "$SPATIAL_HVG_FLAVOR" --umap-point-size "$UMAP_POINT" --formats png,pdf --assignments total \
      >"$LOGS/scanpy_total_${BIN_LABEL}.log" 2>&1
}

run_correction() {
  run_timed "$LOGS/correction_${BIN_LABEL}.time.txt" \
    python3 "$SCRIPTS/07_spatial_dynamics/spatial_slam_correct_new_fraction.py" \
      --molecule-calls "$SOURCE_MOLECULES" \
      --sample "$SAMPLE" --outdir "$TOOLS/correction" --bin-size "$BIN_SIZE_RAW" --bin-size-pixels "$BIN_SIZE_PIXELS" \
      --tc-column "$TC_COLUMN" --t-column "$T_COLUMN" \
      --official-bin-positions "$CELATLAS_BIN_DIR/spatial/tissue_positions_list.csv" \
      >"$LOGS/correction_${BIN_LABEL}.log" 2>&1
}

pids=()
if [[ "$FORCE" == "--force" || ! -s "$TOOLS/scanpy/06_spatial_total_domains_${BIN_LABEL}.tsv" ]]; then
  run_scanpy & pids+=("$!")
fi
if [[ "$FORCE" == "--force" || ! -s "$TOOLS/correction/${SAMPLE}.corrected_bin_gene_counts.tsv.gz" ]]; then
  run_correction & pids+=("$!")
fi
fail=0
for pid in "${pids[@]}"; do
  wait "$pid" || fail=1
done
[[ "$fail" -eq 0 ]] || { echo "ERROR: ${BIN_LABEL} parallel prerequisite failed; inspect $LOGS" >&2; exit 1; }

if [[ "$REVIEWED_MODE" == "1" && ( "$FORCE" == "--force" || ! -s "$TOOLS/model_validation/model_validation_summary.json" ) ]]; then
  run_timed "$LOGS/model_validation_${BIN_LABEL}.time.txt" \
    python3 "$SCRIPTS/07_spatial_dynamics/validate_mixture_model.py" \
      --molecule-calls "$SOURCE_MOLECULES" \
      --substitution-spectrum "$OUT/05_molecule_call/05_metadata/spatial_substitution_spectrum.tsv" \
      --outdir "$TOOLS/model_validation" --base-quality-thresholds "$BASE_QUALITY_THRESHOLDS" \
      --default-base-quality "$DEFAULT_BASE_QUALITY" --grand-slam-p-old 0.0013034 --bootstrap-replicates 200 \
      >"$LOGS/model_validation_${BIN_LABEL}.log" 2>&1
fi

if [[ "$FORCE" == "--force" || ! -s "$TOOLS/kinetics/domain_gene_kinetics.tsv" ]]; then
  radius_args=()
  [[ -z "$SMOOTHING_MAX_RADIUS" ]] || radius_args=(--smoothing-max-radius "$SMOOTHING_MAX_RADIUS")
  run_timed "$LOGS/kinetics_${BIN_LABEL}.time.txt" \
    python3 "$SCRIPTS/07_spatial_dynamics/spatial_slam_kinetics.py" \
      --sample "$SAMPLE" \
      --corrected-bin-gene "$TOOLS/correction/${SAMPLE}.corrected_bin_gene_counts.tsv.gz" \
      --domains "$TOOLS/scanpy/06_spatial_total_domains_${BIN_LABEL}.tsv" \
      --outdir "$TOOLS/kinetics" --min-gene-total-umis 100 --top-genes-per-domain 8 \
      --smoothing-neighbors "$SMOOTHING_NEIGHBORS" --official-bin-positions "$CELATLAS_BIN_DIR/spatial/tissue_positions_list.csv" \
      "${radius_args[@]}" --formats png,pdf --dpi 450 \
      >"$LOGS/kinetics_${BIN_LABEL}.log" 2>&1
fi

if [[ "$REVIEWED_MODE" == "1" && ( "$FORCE" == "--force" || ! -s "$TOOLS/domain_validation/domain_stability_marker_summary.json" ) ]]; then
  run_timed "$LOGS/domain_validation_${BIN_LABEL}.time.txt" \
    python3 "$SCRIPTS/07_spatial_dynamics/analyze_domain_stability_markers.py" \
      --celatlas-bin-dir "$CELATLAS_BIN_DIR" --domains "$TOOLS/scanpy/06_spatial_total_domains_${BIN_LABEL}.tsv" \
      --outdir "$TOOLS/domain_validation" --resolutions 0.4,0.6,0.8,1.0,1.2,1.4 \
      --seeds 0,1,2,3,4 --final-resolution "$LEIDEN_RES" --hvg "$SPATIAL_HVG" \
      --pca-components "$SPATIAL_PCA" --neighbors "$SPATIAL_NEIGHBORS" \
      >"$LOGS/domain_validation_${BIN_LABEL}.log" 2>&1
fi

if [[ "$REVIEWED_MODE" == "1" && ( "$FORCE" == "--force" || ! -s "$META/reviewed_biology_summary.json" ) ]]; then
  [[ -s "$GRAND_GENE_NTR" ]] || { echo "ERROR: GRAND-SLAM gene table missing: $GRAND_GENE_NTR" >&2; exit 1; }
  run_timed "$LOGS/reviewed_biology_${BIN_LABEL}.time.txt" \
    python3 "$SCRIPTS/07_spatial_dynamics/analyze_reviewed_biology.py" \
      --domains "$TOOLS/scanpy/06_spatial_total_domains_${BIN_LABEL}.tsv" \
      --bin-summary "$TOOLS/kinetics/bin_kinetics_summary.tsv.gz" \
      --domain-gene "$TOOLS/kinetics/domain_gene_kinetics.tsv" \
      --global-gene "$TOOLS/kinetics/global_gene_kinetics.tsv" \
      --markers "$TOOLS/domain_validation/domain_markers.tsv" \
      --marker-universe "$TOOLS/domain_validation/domain_marker_gene_universe.tsv" \
      --raw-bin-gene "$TOOLS/kinetics/bin_gene_kinetics_raw.tsv.gz" \
      --smoothed-bin-gene "$TOOLS/kinetics/bin_gene_kinetics_smoothed.tsv.gz" \
      --positions "$CELATLAS_BIN_DIR/spatial/tissue_positions_list.csv" \
      --gene-info "$ROOT/pipeline/resources/Mus_musculus.gene_info.gz" \
      --reactome-gmt "$ROOT/pipeline/resources/mouse.reactome.gmt" --grand-gene-ntr "$GRAND_GENE_NTR" \
      --output-dir "$SHOW" >"$LOGS/reviewed_biology_${BIN_LABEL}.log" 2>&1
fi

if [[ "$REVIEWED_MODE" == "1" && ( "$FORCE" == "--force" || ! -s "$TOOLS/model_validation/method_spatial_sensitivity.tsv" ) ]]; then
  run_timed "$LOGS/method_spatial_sensitivity_${BIN_LABEL}.time.txt" \
    python3 "$SCRIPTS/07_spatial_dynamics/analyze_method_spatial_sensitivity.py" \
      --molecule-calls "$SOURCE_MOLECULES" --model-summary "$TOOLS/model_validation/model_validation_summary.json" \
      --domains "$TOOLS/scanpy/06_spatial_total_domains_${BIN_LABEL}.tsv" \
      --positions "$CELATLAS_BIN_DIR/spatial/tissue_positions_list.csv" \
      --domain-gene "$TOOLS/kinetics/domain_gene_kinetics.tsv" --outdir "$TOOLS/model_validation" \
      --bin-size-pixels "$BIN_SIZE_PIXELS" >"$LOGS/method_spatial_sensitivity_${BIN_LABEL}.log" 2>&1
fi

if [[ "$REVIEWED_MODE" == "1" ]]; then
  cp "$TOOLS/model_validation"/*.tsv "$TABLES/"
  cp "$TOOLS/model_validation"/*.tsv.gz "$TABLES/"
  cp "$TOOLS/model_validation"/*.png "$FIGURES/"
  cp "$TOOLS/domain_validation"/*.tsv "$TABLES/"
  cp "$TOOLS/domain_validation"/*.png "$FIGURES/"
fi

run_timed "$LOGS/showcase_figures.time.txt" \
  python3 "$SCRIPTS/07_spatial_dynamics/make_bin50_showcase.py" \
    --sample "$SAMPLE" --domains "$TOOLS/scanpy/06_spatial_total_domains_${BIN_LABEL}.tsv" \
    --corrected-bin-gene "$TOOLS/correction/${SAMPLE}.corrected_bin_gene_counts.tsv.gz" \
    --bin-summary "$TOOLS/kinetics/bin_kinetics_summary.tsv.gz" \
    --domain-gene-kinetics "$TOOLS/kinetics/domain_gene_kinetics.tsv" \
    --smoothed-bin-gene "$TOOLS/kinetics/bin_gene_kinetics_smoothed.tsv.gz" \
    --output-dir "$SHOW" --bin-size-raw "$BIN_SIZE_RAW" --bin-size-um "$BIN_SIZE_UM" --coordinate-um "$COORDINATE_UM" \
    --spatial-point-size "$SPATIAL_POINT" --umap-point-size "$UMAP_POINT" --interpretation-mode "$INTERPRETATION_MODE" --dpi 450 \
    >"$LOGS/showcase_figures.log" 2>&1

run_timed "$LOGS/representative_gene_old_new_ntr.time.txt" \
  python3 "$SCRIPTS/07_spatial_dynamics/plot_representative_gene_old_new_ntr.py" \
    --smoothed-bin-gene "$TOOLS/kinetics/bin_gene_kinetics_smoothed.tsv.gz" \
    --domain-gene-kinetics "$TOOLS/kinetics/domain_gene_kinetics.tsv" \
    --domains "$TOOLS/scanpy/06_spatial_total_domains_${BIN_LABEL}.tsv" \
    --output-dir "$FIGURES" --metadata-dir "$META" --background-point-size "$REP_BG_POINT" --point-size "$REP_POINT" --n-genes 4 --interpretation-mode "$INTERPRETATION_MODE" --bin-label "$BIN_LABEL" --dpi 450 \
    >"$LOGS/representative_gene_old_new_ntr.log" 2>&1

python3 - "$FIGURES" "$META/output_manifest.tsv" <<'PY'
import hashlib,sys
from pathlib import Path
figures=Path(sys.argv[1]); output=Path(sys.argv[2])
rows=[]
for path in sorted(figures.iterdir()):
    if path.suffix.lower() not in {'.png','.pdf'}: continue
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    rows.append((path.name,path.stat().st_size,digest))
output.write_text('file\tbytes\tsha256\n'+''.join(f'{name}\t{size}\t{digest}\n' for name,size,digest in rows),encoding='utf-8')
PY
python3 - "$LOGS" "$META/resource_usage.tsv" <<'PY'
import re,sys
from pathlib import Path
logs=Path(sys.argv[1]); output=Path(sys.argv[2])
patterns={
    'user_seconds': r'^\s*User time \(seconds\):\s*(.+)$',
    'system_seconds': r'^\s*System time \(seconds\):\s*(.+)$',
    'cpu_percent': r'^\s*Percent of CPU this job got:\s*(.+)$',
    'elapsed': r'^\s*Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(.+)$',
    'max_rss_kb': r'^\s*Maximum resident set size \(kbytes\):\s*(.+)$',
    'fs_inputs': r'^\s*File system inputs:\s*(.+)$',
    'fs_outputs': r'^\s*File system outputs:\s*(.+)$',
}
rows=[]
for path in sorted(logs.glob('*.time.txt')):
    text=path.read_text(encoding='utf-8',errors='replace')
    values={key:(match.group(1).strip() if (match:=re.search(pattern,text,re.M)) else '') for key,pattern in patterns.items()}
    rows.append((path.name.removesuffix('.time.txt'),*[values[key] for key in patterns]))
output.write_text('step\t'+'\t'.join(patterns)+'\n'+''.join('\t'.join(row)+'\n' for row in rows),encoding='utf-8')
PY
echo "Completed ${BIN_LABEL} showcase: $SHOW"
