#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
image=${IMAGE_TAG:-celatlas-spatial-slam:20261023-v1.8.0}
archive=${IMAGE_ARCHIVE:-$root/results/docker_delivery/celatlas-spatial-slam_20261023-v1.8.0.tar}
source_root=${CELATLAS_SOURCE:-/opt/project/software/celatlas_spatial_v1.8_DV/celatlas_spatial-1.8.0/celatlas_spatial-1.8.0}
runtime_env=${CELATLAS18_ENV:-/opt/project/software/micromamba-root/envs/celatlas18-stardist-py311}
model=${CELATLAS_MODEL_FILE:-/opt/project/software/celatlas_spatial_v1.8_DV/src/swin_tiny.pth}
logo=${FOOTER_LOGO:-/opt/project/projects/epi-logo.png}
context="$root/tmp/docker-build-context.$$"
archive_partial="$archive.partial.$$"

command -v docker >/dev/null || { echo "ERROR: docker is unavailable" >&2; exit 127; }
command -v rsync >/dev/null || { echo "ERROR: rsync is unavailable" >&2; exit 127; }
[[ -f "$source_root/celatlas_spatial/__init__.py" && -d "$runtime_env" && -f "$model" ]] || { echo "ERROR: Celatlas v1.8 source/runtime/model missing" >&2; exit 2; }
[[ -f "$logo" ]] || { echo "ERROR: footer logo missing: $logo" >&2; exit 2; }
[[ ! -e "$archive" ]] || { echo "ERROR: archive already exists: $archive" >&2; exit 2; }
if docker image inspect "$image" >/dev/null 2>&1; then
  echo "ERROR: image tag already exists: $image" >&2
  exit 2
fi

mkdir -p "$context/gedi-bin" "$context/runtime_celatlas18" "$context/system_libs" "$(dirname "$archive")" "$root/results/docker_delivery"
trap 'rm -rf -- "$context"; rm -f -- "$archive_partial"' EXIT
rsync -a --exclude='__pycache__/' --exclude='*.pyc' "$runtime_env/" "$context/runtime_celatlas18/"
cp "$model" "$context/swin_tiny.pth"
for lib in libGL.so.1.7.0 libGLdispatch.so.0.0.0 libGLX.so.0.0.0 \
           libX11.so.6.4.0 libxcb.so.1.1.0 libXau.so.6.0.0 libXdmcp.so.6.0.0 \
           libbsd.so.0.11.5 libmd.so.0.0.5; do
  [[ -e "/lib/x86_64-linux-gnu/$lib" ]] || { echo "ERROR: required host library is missing: $lib" >&2; exit 2; }
  cp -L "/lib/x86_64-linux-gnu/$lib" "$context/system_libs/$lib"
done
rsync -a --exclude=/.git/ "$root/vendor/gedi-bin/" "$context/gedi-bin/"
mkdir -p "$context/subread/bin" "$context/subread/lib"
cp "$root/vendor/subread/bin/featureCounts" "$context/subread/bin/featureCounts"
cp "$root/vendor/subread/lib/libz.so.1" "$context/subread/lib/libz.so.1"
cp "$root/containers/docker/Dockerfile" "$context/Dockerfile"
cp "$root/pipeline/envs/snakemake.yaml" "$context/snakemake.yaml"
cp "$root/pipeline/envs/runtime.yaml" "$context/runtime.yaml"
cp "$logo" "$context/epi-logo.png"

docker build --pull=false --tag "$image" "$context"
docker run --rm "$image" --help >/dev/null
docker run --rm --entrypoint /opt/celatlas18/bin/python "$image" -c \
  'import anndata,scanpy,pysam,yaml,celatlas_spatial; assert celatlas_spatial.__version__ == "1.8.0"; print("Celatlas", celatlas_spatial.__version__, "runtime imports: OK")'
docker run --rm --network none --entrypoint /bin/bash "$image" -c \
  'set -e; celatlas count --help >/dev/null; celatlas_spatial rna analysis --help >/dev/null; celatlas_spatial_report --help >/dev/null'
docker run --rm --network none --entrypoint /bin/bash "$image" -c \
  'set -e; test -s "$CELATLAS_MODEL"; for tool in celatlas_spatial STAR samtools bedtools featureCounts slamdunk ngm gedi java Rscript rustc pigz; do command -v "$tool" >/dev/null; done; featureCounts -v; gedi -e Slam -h >/tmp/gedi-help.txt 2>&1; test -s /tmp/gedi-help.txt; Rscript -e "library(cowplot);library(ggplot2);library(reshape2)"'
docker run --rm --network none --entrypoint /opt/celatlas18/bin/python "$image" -c \
  'import numpy as np, anndata as ad, scanpy as sc; a=ad.AnnData(np.array([[1,0],[0,1],[1,1],[2,1],[1,2]],dtype=float)); sc.pp.neighbors(a,n_neighbors=2,use_rep="X"); sc.tl.leiden(a,resolution=0.5,random_state=0,n_iterations=2); assert "leiden" in a.obs'
docker save --output "$archive_partial" "$image"
mv -- "$archive_partial" "$archive"
printf 'image\t%s\narchive\t%s\n' "$image" "$archive"
