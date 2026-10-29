# Celatlas Spatial SLAM

Reusable Docker/Snakemake workflow for spatial SLAM analysis. The pipeline
separates input QC, Celatlas integration, barcode processing, alignment,
molecule calls, matrices, spatial dynamics, tissue-bulk analysis, and report
generation into explicit stages.

This public repository contains reusable workflow code, container definitions,
tests, and layout contracts. Raw data, references, result trees, and
server-specific run profiles are excluded.

~~~bash
make lint
make smoke
bash containers/docker/run_snakemake.sh --help
~~~

Provide a local run configuration and mount external inputs and references
read-only. Runtime logs, outputs, and temporary files stay outside Git.
