#!/usr/bin/env python3
"""Execute exactly one Celatlas Spatial SLAM stage for Snakemake."""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
from pathlib import Path

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = PROJECT_ROOT / "pipeline" / "scripts"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(command: list[str], *, env: dict[str, str] | None = None, stdout=None) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True, env=env, stdout=stdout)


class StageContext:
    def __init__(self, config_path: Path):
        self.config_path = config_path.resolve()
        self.config = yaml.safe_load(self.config_path.read_text(encoding="utf-8"))
        self.sample = self.config["sample_id"]
        self.project = self.config["project_id"]
        self.samples_path = (PROJECT_ROOT / self.config["samples"]).resolve()
        self.libraries_path = (PROJECT_ROOT / self.config["libraries"]).resolve()
        self.sample_row = next(r for r in rows(self.samples_path) if r["sample_id"] == self.sample)
        self.libraries = [r for r in rows(self.libraries_path) if r["sample_id"] == self.sample]
        self.chip = self.sample_row["chip_id"]
        self.species = self.sample_row["species"]
        self.chemistry = self.sample_row["chemistry"]
        self.coordinate_um = self.sample_row["coordinate_um_per_unit"]
        self.mask_dir = (PROJECT_ROOT / self.sample_row["mask_dir"]).resolve()
        self.reference_dir = (PROJECT_ROOT / self.sample_row["reference_dir"]).resolve()
        self.outputs = (PROJECT_ROOT / self.config.get("outputs_root", "outputs")).resolve()
        self.source_outputs = (PROJECT_ROOT / self.config.get("source_outputs_root", self.outputs)).resolve()
        self.tmp = (PROJECT_ROOT / self.config.get("tmp_root", "tmp")).resolve()
        self.out = self.outputs / self.sample
        self.source_out = self.source_outputs / self.sample
        celatlas_source = (self.sample_row.get("celatlas_result_dir") or "").strip()
        self.celatlas_reuse = bool(celatlas_source)
        if celatlas_source:
            self.celatlas = (PROJECT_ROOT / celatlas_source).resolve()
        else:
            self.celatlas = (
                self.out / "02_celatlas/04_tool_output/celatlas_full/results"
                / self.project / self.chip
            )
        self.run_scope = self.config.get("run_scope", "full")
        self.threads = int(self.config.get("threads", 4))
        self.interpretation = self.config.get("interpretation_mode", "standard")
        self.fastq_r1 = [str((PROJECT_ROOT / r["fastq_r1"]).resolve()) for r in self.libraries]
        self.fastq_r2 = [str((PROJECT_ROOT / r["fastq_r2"]).resolve()) for r in self.libraries]
        self.upstream_out = self.source_out if self.run_scope == "reviewed" else self.out
        self.filtered_r2 = self.tmp / self.sample / "03_barcode" / f"{self.sample}.filtered_r2.fastq.gz"
        self.tagged_r2 = self.tmp / self.sample / "04_trim_align" / f"{self.sample}.tagged_r2.fastq.gz"
        self.barcode_table = self.upstream_out / "03_barcode" / "01_tables" / f"{self.sample}.read_barcode_umi.tsv.gz"
        self.star_bam = self.upstream_out / "04_trim_align" / "04_tool_output" / "celatlas" / f"{self.sample}_Aligned.sortedByCoord.out.bam"
        self.molecules = self.out / "05_molecule_call" / "01_tables" / f"{self.sample}.molecule_calls.tsv.gz"
        self.matrix_dir = self.out / "06_matrix" / "04_tool_output" / "matrix"
        self.spatial_source_out = self.source_out if self.run_scope == "spatial_report" else self.out
        self.source_molecules = self.spatial_source_out / "05_molecule_call" / "01_tables" / f"{self.sample}.molecule_calls.tsv.gz"
        self.source_matrix_dir = self.spatial_source_out / "06_matrix" / "04_tool_output" / "matrix"
        self.bin_size = int(self.config.get("spatial", {}).get("bin_size_raw", 50))
        self.bin_label = f"bin{self.bin_size}"
        self.bin_dir = self.celatlas / "06.segment" / "01.binsegment" / "square_bin" / f"{self.chip}_{self.bin_label}"
        self.env = os.environ.copy()
        runtime = self.tmp / "runtime"
        for path in (runtime / "home", runtime / "tmp", runtime / "numba", runtime / "matplotlib"):
            path.mkdir(parents=True, exist_ok=True)
        self.env.update({
            "HOME": str(runtime / "home"), "TMPDIR": str(runtime / "tmp"),
            "TMP": str(runtime / "tmp"), "TEMP": str(runtime / "tmp"),
            "NUMBA_CACHE_DIR": str(runtime / "numba"),
            "MPLCONFIGDIR": str(runtime / "matplotlib"), "MPLBACKEND": "Agg",
            "CELATLAS_MODEL": os.environ.get("CELATLAS_MODEL", self.config["celatlas_model"]),
            "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
            "INTERPRETATION_MODE": self.interpretation,
        })

    def stage_dir(self, number: str, name: str) -> Path:
        stage = self.out / f"{number}_{name}"
        for sub in ("01_tables", "02_figures", "03_logs", "04_tool_output", "05_metadata"):
            (stage / sub).mkdir(parents=True, exist_ok=True)
        return stage

    def finish(self, stage: Path, number: str, inputs: list[Path]) -> None:
        metadata = stage / "05_metadata"
        provenance = metadata / "git_provenance.tsv"
        provenance.write_text(
            "key\tvalue\n"
            f"pipeline_git_commit\t{os.environ.get('PIPELINE_GIT_COMMIT', 'unavailable_in_container')}\n"
            f"pipeline_git_dirty\t{os.environ.get('PIPELINE_GIT_DIRTY', 'unknown')}\n"
            f"container_image\t{os.environ.get('CELATLAS_SLAM_IMAGE', 'unknown')}\n",
            encoding="utf-8",
        )
        command = ["python3", str(SCRIPTS / "00_shared" / "stage_contract.py"),
                   "--stage-dir", str(stage), "--stage-id", number]
        for path in inputs:
            command.extend(["--input", str(path)])
        run(command, env=self.env)


def stage_01(c: StageContext) -> None:
    stage = c.stage_dir("01", "input_qc")
    run(["python3", str(SCRIPTS / "01_input_qc" / "check_samples.py"),
         "--samples", str(c.samples_path), "--libraries", str(c.libraries_path), "--sample", c.sample], env=c.env)
    project_qc = c.outputs / "_project" / "01_input_qc"
    run(["python3", str(SCRIPTS / "01_input_qc" / "run_chip_assignment_qc.py"),
         "--samples", str(c.samples_path), "--libraries", str(c.libraries_path),
         "--project-id", c.project, "--output-dir", str(project_qc),
         "--tmp-dir", str(c.tmp / c.project / "01_input_qc" / "chip_assignment"),
         "--minimum-fold", "2.0", "--threads", str(c.threads)], env=c.env)
    shutil.copy2(project_qc / "01_tables" / "chip_assignment_crosscheck.tsv", stage / "01_tables")
    shutil.copy2(project_qc / "05_metadata" / "chip_assignment_qc.json", stage / "05_metadata")
    inputs = [c.samples_path, c.libraries_path, *map(Path, c.fastq_r1),
              c.mask_dir / f"{c.chip}_FilterBarcodes.csv"]
    c.finish(stage, "01", inputs)


def stage_02(c: StageContext) -> None:
    stage = c.stage_dir("02", "celatlas")
    if not c.celatlas_reuse and not _celatlas_full_outputs_ready(c):
        workspace = stage / "04_tool_output/celatlas_full"
        fastq_dir = workspace / "input_fastq"
        image_dir = workspace / "images"
        model_source = Path(c.env["CELATLAS_MODEL"]).resolve()
        model_dir = workspace / "src"
        model_dir.mkdir(parents=True, exist_ok=True)
        model_link = model_dir / "swin_tiny.pth"
        if model_link.is_symlink():
            if model_link.resolve() != model_source:
                raise RuntimeError(f"Celatlas model link points to a different model: {model_link}")
        elif model_link.exists():
            raise RuntimeError(f"Celatlas model path already exists and is not the expected symlink: {model_link}")
        else:
            model_link.symlink_to(model_source)
        runner_config = workspace / "celatlas_runner.yaml"
        runner_config.write_text(
            yaml.safe_dump({"workspace": str(workspace), "src_dir": str(model_dir),
                            "default_insert_r2": 150}, sort_keys=True),
            encoding="utf-8",
        )
        fastq_dir.mkdir(parents=True, exist_ok=True)
        image_dir.mkdir(parents=True, exist_ok=True)

        for index, library in enumerate(c.libraries, start=1):
            lane = f"{index:03d}"
            for mate, key in (("R1", "fastq_r1"), ("R2", "fastq_r2")):
                source = (PROJECT_ROOT / library[key]).resolve()
                staged = fastq_dir / f"{c.chip}_S1_L{lane}_{mate}_001.fastq.gz"
                if staged.is_symlink():
                    if staged.resolve() != source:
                        raise RuntimeError(f"Celatlas FASTQ link points to a different input: {staged}")
                elif staged.exists():
                    raise RuntimeError(f"Celatlas FASTQ staging path already exists and is not a symlink: {staged}")
                else:
                    staged.symlink_to(source)

        reference_root = c.tmp / "runtime/celatlas_reference_root"
        reference_root.mkdir(parents=True, exist_ok=True)
        species_reference = reference_root / c.species
        if species_reference.is_symlink():
            if species_reference.resolve() != c.reference_dir:
                raise RuntimeError(f"Celatlas reference link points to a different reference: {species_reference}")
        elif species_reference.exists():
            raise RuntimeError(f"Celatlas reference path already exists and is not a symlink: {species_reference}")
        else:
            species_reference.symlink_to(c.reference_dir, target_is_directory=True)

        celatlas_env = c.env.copy()
        celatlas_env.update({
            "CELATLAS_WORKSPACE": str(workspace),
            "CELATLAS_SRC_DIR": str(model_dir),
            "CELATLAS_PYTHON": os.environ.get(
                "CELATLAS_PYTHON", "/opt/conda/envs/runtime/bin/python"
            ),
            "CLEAN_INTERMEDIATE": "0",
        })
        celatlas_cli = c.config.get("celatlas_cli", "/opt/conda/envs/runtime/bin/celatlas")
        command = [
            celatlas_cli, "count",
            "--id", c.chip,
            "--sample-name", c.sample,
            "--tissue", c.config.get("tissue", "unspecified"),
            "--targetdir", str(workspace / "results" / c.project / c.chip),
            "--project", c.project,
            "--chemistry", c.chemistry,
            "--species", c.species,
            "--image", "gene_expr",
            "--workflow", "ST",
            "--mode", "strna",
            "--thread", str(c.threads),
            "--featurecounts-thread", str(c.threads),
            "--bin", "50",
            "--pixel-size", str(c.coordinate_um),
            "--reference-dir", str(reference_root),
            "--mask-dir", str(c.mask_dir),
            "--image-dir", str(image_dir),
            "--fastqs", str(fastq_dir),
            "--fastq-name", c.chip,
            "--cell-num", "50000",
            "--config", str(runner_config),
            "--log-dir", str(workspace / "logs"),
            "--resume-existing",
            "--allow-existing-step-outputs",
        ]
        run(command, env=celatlas_env)
    elif not c.celatlas_reuse:
        print(f"[Celatlas stage] {c.sample}: complete v1.8 outputs already exist; skip rerunning Celatlas")

    validation = stage / "05_metadata" / "celatlas_reuse_validation.json"
    run(["python3", str(SCRIPTS / "02_celatlas" / "validate_reused_celatlas.py"),
         "--source", str(c.celatlas), "--sample", c.sample, "--chip", c.chip,
         "--chemistry", c.chemistry, "--species", c.species, "--output", str(validation)], env=c.env)
    (stage / "05_metadata" / "celatlas_result_dir.txt").write_text(str(c.celatlas) + "\n", encoding="utf-8")
    execution = "validated_existing_result" if c.celatlas_reuse else "full_official_run"
    (stage / "05_metadata" / "celatlas_full_method.tsv").write_text(
        "key\tvalue\n"
        f"tool\tCelatlas Spatial\nversion\t1.8.0\nexecution\t{execution}\n"
        f"source_result_dir\t{c.celatlas}\nbin_um\t50\npixel_size_um\t{c.coordinate_um}\n",
        encoding="utf-8",
    )
    run(["python3", str(SCRIPTS / "07_spatial_dynamics" / "record_celatlas_image_provenance.py"),
         "--positions", str(c.celatlas / "06.segment/01.binsegment/square_bin" / f"{c.chip}_Raw/spatial/tissue_positions_list.csv"),
         "--registration-image", str(c.celatlas / "06.segment/01.binsegment/images" / f"{c.chip}_regist.tif"),
         "--hires-image", str(c.celatlas / "06.segment/01.binsegment/images/tissue_hires_image.png"),
         "--lowres-image", str(c.celatlas / "06.segment/01.binsegment/images/tissue_lowres_image.png"),
         "--scalefactors", str(c.celatlas / "06.segment/01.binsegment/square_bin" / f"{c.chip}_Raw/spatial/scalefactors_json.json"),
         "--count-detail", str(c.celatlas / "05.count" / f"{c.chip}_count_detail.txt"),
         "--output", str(stage / "05_metadata/celatlas_image_provenance.json")], env=c.env)
    c.finish(stage, "02", [c.out / "01_input_qc/05_metadata/stage_complete.json", c.celatlas / ".metrics.json"])


def _celatlas_full_outputs_ready(c: StageContext) -> bool:
    required = (
        c.celatlas / ".metrics.json",
        c.celatlas / "05.count" / f"{c.chip}_count_detail.txt",
        c.bin_dir / "filtered_feature_bc_matrix/matrix.mtx.gz",
        c.bin_dir / "filtered_feature_bc_matrix/features.tsv.gz",
        c.bin_dir / "filtered_feature_bc_matrix/barcodes.tsv.gz",
        c.bin_dir / "spatial/tissue_positions_list.csv",
        c.bin_dir.parent / f"{c.chip}_Raw/spatial/tissue_positions_list.csv",
        c.bin_dir.parent / f"{c.chip}_Raw/spatial/scalefactors_json.json",
        c.celatlas / "07.outs/binned_outputs/h5ad" / f"{c.chip}_bin50.h5ad",
        c.celatlas / f"{c.chip}_spatial_analysis_report.html",
    )
    return all(path.is_file() and path.stat().st_size > 0 for path in required)


def stage_03(c: StageContext) -> None:
    stage = c.stage_dir("03", "barcode")
    c.filtered_r2.parent.mkdir(parents=True, exist_ok=True)
    run(["python3", str(SCRIPTS / "03_barcode" / "run_parallel_barcode.py"),
         "--r1", ",".join(c.fastq_r1), "--r2", ",".join(c.fastq_r2),
         "--whitelist", str(c.mask_dir / f"{c.chip}_FilterBarcodes.csv"),
         "--sample", c.sample, "--outdir", str(stage), "--filtered-r2-output", str(c.filtered_r2),
         "--tmp-dir", str(c.tmp / c.sample / "03_barcode/lane_parts"),
         "--log-dir", str(stage / "03_logs"), "--workers", str(c.threads),
         "--pattern", c.config["barcode_pattern"]], env=c.env)
    (stage / "05_metadata/large_file_manifest.tsv").write_text(
        f"path\tpurpose\tdeliverable\n{c.filtered_r2}\tBarcode-filtered R2 for downstream alignment\tno\n",
        encoding="utf-8")
    c.finish(stage, "03", [*map(Path, c.fastq_r1 + c.fastq_r2), c.mask_dir / f"{c.chip}_FilterBarcodes.csv"])


def stage_04(c: StageContext) -> None:
    stage = c.stage_dir("04", "trim_align")
    cutadapt = c.tmp / c.sample / "04_trim_align/cutadapt"
    cutadapt.mkdir(parents=True, exist_ok=True)
    run(["python3", str(SCRIPTS / "04_trim_align/tag_fastq_headers.py"),
         "--fastq", str(c.filtered_r2), "--barcode-table", str(c.barcode_table),
         "--output", str(c.tagged_r2), "--stats", str(stage / "05_metadata/qname_tagging.json"),
         "--threads", str(c.threads)], env=c.env)
    run(["celatlas_spatial", "rna", "cutadapt", "--outdir", str(cutadapt), "--sample", c.sample,
         "--thread", str(c.threads), "--minimum_length", "20", "--nextseq_trim", "20", "--gzip",
         "--overlap", "10", "--insert", "150", "--fq", str(c.tagged_r2)], env=c.env)
    star_out = stage / "04_tool_output/celatlas"
    star_tmp = star_out / f"{c.sample}__STARtmp"
    if star_tmp.exists():
        failed_root = c.tmp / c.sample / "04_trim_align/failed_star_attempts"
        failed_root.mkdir(parents=True, exist_ok=True)
        attempt = 1
        while (failed_root / f"{star_tmp.name}.attempt_{attempt:02d}.txt").exists():
            attempt += 1
        (failed_root / f"{star_tmp.name}.attempt_{attempt:02d}.txt").write_text(
            f"Removed incomplete STAR temporary directory before retry: {star_tmp}\n",
            encoding="utf-8",
        )
        shutil.rmtree(star_tmp)
    run(["celatlas_spatial", "rna", "star", "--outdir", str(star_out),
         "--sample", c.sample, "--thread", str(c.threads), "--genomeDir", str(c.reference_dir),
         "--outFilterMatchNmin", "30", "--outFilterMultimapNmax", "1",
         "--STAR_param", "--outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --outSAMattributes MD NH HI AS nM",
         "--starMem", "30", "--fq", str(cutadapt / f"{c.sample}_clean_2.fq.gz")], env=c.env)
    c.finish(stage, "04", [c.filtered_r2, c.barcode_table, c.reference_dir / "Genome"])


def stage_05(c: StageContext) -> None:
    stage = c.stage_dir("05", "molecule_call")
    annotation = stage / "04_tool_output/annotation"
    run(["python3", str(SCRIPTS / "05_molecule_call/build_exon_intron_annotation.py"),
         "--gtf", str(c.reference_dir / "genes.gtf"), "--outdir", str(annotation)], env=c.env)
    if not Path(f"{c.star_bam}.bai").is_file():
        run(["samtools", "index", "-@", str(c.threads), str(c.star_bam)], env=c.env)
    reviewed = c.run_scope == "reviewed"
    caller = "call_new_old_molecules_reviewed.py" if reviewed else "call_new_old_molecules.py"
    command = ["python3", str(SCRIPTS / "05_molecule_call" / caller),
         "--bam", str(c.star_bam), "--fasta", str(c.reference_dir / "genome.fa"),
         "--annotation-dir", str(annotation), "--barcode-table", str(c.barcode_table),
         "--sample", c.sample, "--outdir", str(stage),
         "--tmp-dir", str(c.tmp / c.sample / "05_molecule_call/chromosome_parts"),
         "--threads", str(c.config.get("molecule_call_workers", min(c.threads, 8))),
         "--new-threshold", str(c.config["new_threshold"])]
    if reviewed:
        molecule = c.config["molecule_call"]
        command.extend([
            "--min-base-quality", str(molecule["min_base_quality"]),
            "--sensitivity-base-qualities", ",".join(map(str, molecule["sensitivity_base_qualities"])),
        ])
        if molecule.get("snp_vcf"):
            command.extend(["--snp-vcf", str((PROJECT_ROOT / molecule["snp_vcf"]).resolve())])
    run(command, env=c.env)
    c.finish(stage, "05", [c.star_bam, c.barcode_table, c.reference_dir / "genome.fa", c.reference_dir / "genes.gtf"])


def stage_06(c: StageContext) -> None:
    stage = c.stage_dir("06", "matrix")
    run(["python3", str(SCRIPTS / "06_matrix/build_old_new_matrices.py"),
         "--molecule-calls", str(c.molecules), "--sample", c.sample, "--outdir", str(stage)], env=c.env)
    c.finish(stage, "06", [c.molecules])


def stage_07(c: StageContext) -> None:
    stage = c.stage_dir("07", "spatial_dynamics")
    spatial = c.config["spatial"]
    plotting = spatial["plotting"]
    env = c.env.copy()
    env.update({
        "BIN_SIZE_RAW": str(spatial["bin_size_raw"]), "BIN_SIZE_UM": str(spatial["bin_size_um"]),
        "BIN_SIZE_PIXELS": str(spatial["bin_size_pixels"]), "SPATIAL_HVG": str(spatial["hvg"]),
        "SPATIAL_HVG_FLAVOR": spatial["hvg_flavor"], "SPATIAL_PCA": str(spatial["pca_components"]),
        "SPATIAL_NEIGHBORS": str(spatial["neighbors"]), "LEIDEN_RES": str(spatial["leiden_resolution"]),
        "LEIDEN_SEED": str(spatial["leiden_random_state"]), "LEIDEN_ITER": str(spatial["leiden_n_iterations"]),
        "SPATIAL_POINT": str(plotting["spatial_point_size"]), "UMAP_POINT": str(plotting["umap_point_size"]),
        "REP_BG_POINT": str(plotting["representative_background_point_size"]),
        "REP_POINT": str(plotting["representative_point_size"]),
        "SMOOTHING_NEIGHBORS": str(spatial.get("smoothing_neighbors", 24)),
        "SMOOTHING_MAX_RADIUS": str(spatial.get("smoothing_max_radius_pixels", "")),
        "REVIEWED_MODE": "1" if c.run_scope == "reviewed" else "0",
    })
    if c.run_scope == "reviewed":
        molecule = c.config["molecule_call"]
        threshold = int(molecule["min_base_quality"])
        env.update({
            "TC_COLUMN": f"tc_events_bq{threshold}",
            "T_COLUMN": f"t_covered_bq{threshold}",
            "DEFAULT_BASE_QUALITY": str(threshold),
            "BASE_QUALITY_THRESHOLDS": ",".join(map(str, molecule["sensitivity_base_qualities"])),
            "GRAND_GENE_NTR": str(c.source_out / "08_tissue_bulk/04_tool_output/grand_slam" / f"{c.sample}.tsv.gz"),
        })
    run(["bash", str(SCRIPTS / "07_spatial_dynamics/run_spatial_showcase.sh"), c.sample,
         str(c.out), c.coordinate_um, str(c.bin_dir), str(c.source_matrix_dir),
         str(c.source_molecules)], env=env)
    c.finish(stage, "07", [c.source_matrix_dir / "old_umi.mtx", c.source_matrix_dir / "new_umi.mtx",
                            c.source_molecules, c.bin_dir / "spatial/tissue_positions_list.csv"])


def stage_08(c: StageContext) -> None:
    stage = c.stage_dir("08", "tissue_bulk")
    domains = c.out / "07_spatial_dynamics/04_tool_output/scanpy" / f"06_spatial_total_domains_{c.bin_label}.tsv"
    corrected = c.out / "07_spatial_dynamics/04_tool_output/correction" / f"{c.sample}.corrected_bin_gene_counts.tsv.gz"
    run(["python3", str(SCRIPTS / "08_tissue_bulk/export_tissuumaps.py"),
         "--spatial-domains", str(domains), "--corrected-bin-gene", str(corrected),
         "--output-dir", str(stage / "04_tool_output/tissuumaps"), "--coordinate-um", c.coordinate_um], env=c.env)
    converted = stage / "04_tool_output/celatlas_tissue_mask"
    run(["python3", str(SCRIPTS / "08_tissue_bulk/build_tissue_mask.py"),
         "--molecule-calls", str(c.molecules), "--celatlas-dir", str(c.celatlas),
         "--outdir", str(converted), "--roi-coordinate-um", c.coordinate_um], env=c.env)
    shutil.copy2(converted / "01_tissue_barcodes.tsv", stage / "01_tables/01_tissue_barcodes.tsv")
    shutil.copy2(converted / "tissue_mask_qc.json", stage / "05_metadata/tissue_mask_qc.json")
    run(["bash", str(SCRIPTS / "08_tissue_bulk/run_tissue_bulk_tools.sh"),
         "--sample", c.sample, "--stage-dir", str(stage),
         "--tmp-dir", str(c.tmp / c.sample / "08_tissue_bulk"),
         "--barcode-table", str(c.barcode_table),
         "--tissue-barcodes", str(stage / "01_tables/01_tissue_barcodes.tsv"),
         "--fastq", str(c.filtered_r2), "--reference-fasta", str(c.reference_dir / "genome.fa"),
         "--gtf", str(c.reference_dir / "genes.gtf"), "--gedi-genome", c.config["gedi_genomic"],
         "--threads", str(c.threads), "--trim-5p", "22"], env=c.env)
    c.finish(stage, "08", [c.molecules, c.barcode_table, c.filtered_r2, domains, corrected])


def stage_09(c: StageContext) -> None:
    stage = c.stage_dir("09", "report")
    bulk_source = (c.source_out if c.run_scope in {"spatial_report", "reviewed"} else c.out) / "08_tissue_bulk"
    if c.run_scope == "reviewed":
        runtime_logs = (PROJECT_ROOT / c.config.get("log_root", "log") / c.sample).resolve()
        delivered_logs = stage / "03_logs/snakemake"
        delivered_logs.mkdir(parents=True, exist_ok=True)
        if runtime_logs.is_dir():
            for source in sorted(runtime_logs.glob("*")):
                if source.is_file() and source.suffix in {".log", ".tsv"}:
                    shutil.copy2(source, delivered_logs / source.name)
        report = stage / c.config.get("report_filename", f"{c.sample}_report_reviewed.html")
        manifest = stage / "05_metadata/delivery_manifest.tsv"
        run([
            "python3", str(SCRIPTS / "09_report/build_reviewed_delivery.py"),
            "--sample-root", str(c.out),
            "--baseline-root", str(c.config.get("baseline_result_root", "baseline_bin20_result")),
            "--bulk-root", str(bulk_source), "--output", str(report), "--manifest", str(manifest),
        ], env=c.env)
        validator = [
            "python3", str(SCRIPTS / "09_report/validate_delivery_report.py"),
            "--report", str(report), "--manifest", str(manifest),
            "--require-expandable-figures", "--require-print-css",
        ]
        for section in (
            "01 数据概览", "02 Barcode 与基础 QC", "03 Total RNA 空间图谱",
            "04 Domain markers 与 putative annotation", "05 Spatial SLAM 模型与总体 RNA 年龄",
            "06 Total / New / NTR 空间对照", "07 Domain RNA dynamics", "08 Representative genes",
            "09 Pathway-level dynamics", "10 Bulk method cross-check", "11 Limitations",
        ):
            validator.extend(["--require-section", section])
        for term in ("orthogonal validation", "显著动态基因", "绝对更新速率"):
            validator.extend(["--forbid-text", term])
        validator.append("--json")
        validation = stage / "05_metadata/report_validation.json"
        with validation.open("w", encoding="utf-8") as handle:
            run(validator, env=c.env, stdout=handle)
        c.finish(stage, "09", [bulk_source / "05_metadata/stage_complete.json",
                                c.out / "07_spatial_dynamics/05_metadata/stage_complete.json"])
        return
    bulk = c.out / "08_tissue_bulk"
    if c.run_scope == "spatial_report":
        for sub in ("01_tables", "02_figures", "05_metadata"):
            (bulk / sub).mkdir(parents=True, exist_ok=True)
        shutil.copy2(bulk_source / "05_metadata/tissue_mask_qc.json", bulk / "05_metadata/tissue_mask_qc.json")
        shutil.copy2(bulk_source / "05_metadata/software_versions.tsv", bulk / "05_metadata/software_versions.tsv")
    spatial = c.out / "07_spatial_dynamics"
    bulk_summary = stage / "01_tables/bulk_slam_summary.json"
    run(["python3", str(SCRIPTS / "09_report/summarize_bulk_slam.py"),
         "--slamdunk-dir", str(bulk_source / "04_tool_output/slamdunk"),
         "--grand-prefix", str(bulk_source / "04_tool_output/grand_slam" / c.sample),
         "--mask-qc", str(bulk_source / "05_metadata/tissue_mask_qc.json"), "--output", str(bulk_summary),
         "--spectrum-table", str(bulk / "01_tables/05_complete_12_channel_mismatch_spectrum.tsv"),
         "--figure-dir", str(bulk / "02_figures"), "--provenance", str(bulk / "05_metadata/bulk_figure_provenance.json"),
         "--slamdunk-min-reads", "100", "--grand-min-reads", "100", "--top-genes", "15",
         "--interpretation-mode", c.interpretation, "--dpi", "450"], env=c.env)
    resource_summary = stage / "05_metadata/pipeline_resource_summary.tsv"
    resource_command = ["python3", str(SCRIPTS / "01_input_qc/pipeline_progress.py"), "resources",
                        "--result-dir", str(c.out), "--output", str(resource_summary)]
    for candidate in (spatial / "05_metadata/resource_usage.tsv",):
        if candidate.is_file():
            resource_command.extend(["--sampled-resource-tsv", str(candidate)])
    run(resource_command, env=c.env)
    read_flow = stage / "04_tool_output/read_flow"
    plotly = SCRIPTS / "09_report/read_flow_assets/plotly-3.5.1.min.js"
    run(["python3", str(SCRIPTS / "09_report/build_spatial_read_flow.py"),
         "--sample", c.sample, "--celatlas-data", str(c.celatlas / ".data.json"),
         "--plotly-js", str(plotly), "--output-dir", str(read_flow),
         "--bin-label", c.bin_label], env=c.env)
    report = stage / c.config.get("report_filename", f"{c.sample}_report.html")
    manifest = stage / "05_metadata/delivery_manifest.tsv"
    run(["python3", str(SCRIPTS / "09_report/build_spatial_slam_combined_report.py"),
         "--sample", c.sample, "--sample-metadata", str(c.samples_path),
         "--bin-label", c.bin_label, "--bin-size-um", str(c.config["spatial"]["bin_size_um"]),
         "--figure-dir", str(spatial / "02_figures"),
         "--spatial-metrics", str(spatial / "05_metadata/coordinate_and_figure_provenance.json"),
         "--barcode-qc", str(c.upstream_out / "03_barcode/05_metadata" / f"{c.sample}.barcode_qc.json"),
         "--call-qc", str(c.spatial_source_out / "05_molecule_call/05_metadata" / f"{c.sample}.qc.json"),
         "--celatlas-metrics", str(c.celatlas / ".data.json"),
         "--chip-assignment-qc", str(c.upstream_out / "01_input_qc/05_metadata/chip_assignment_qc.json"),
         "--domain-kinetics", str(spatial / "01_tables" / f"{c.bin_label}_domain_gene_kinetics.tsv.gz"),
         "--resource-summary", str(resource_summary), "--software-versions", str(bulk / "05_metadata/software_versions.tsv"),
         "--read-flow-section", str(read_flow / "reports/read_flow_sankey_section.html"),
         "--read-flow-plotly", str(plotly), "--read-flow-manifest", str(read_flow / "reports/sankey_manifest.tsv"),
         "--footer-logo", c.config["footer_logo"], "--interpretation-mode", c.interpretation,
         "--bulk-summary", str(bulk_summary), "--bulk-figure-dir", str(bulk / "02_figures"),
         "--bulk-spectrum-table", str(bulk / "01_tables/05_complete_12_channel_mismatch_spectrum.tsv"),
         "--output", str(report), "--manifest", str(manifest)], env=c.env)
    validator = ["python3", str(SCRIPTS / "09_report/validate_delivery_report.py"),
                 "--report", str(report), "--manifest", str(manifest),
                 "--require-section", "1 项目概览与核心结论", "--require-section", "2 守恒数据流与关键损失",
                 "--require-section", "3 数据质量与 Celatlas 基线", "--require-section", "4 Total RNA 空间结构",
                 "--require-section", "5 RNA 年龄分解与校正 NTR", "--require-section", "6 空间动力学与代表基因",
                 "--require-section", "7 Tissue-only bulk SLAM", "--require-section", "8 方法、限制与运行信息",
                 "--forbid-text", "old/new 辅助聚类", "--forbid-text", "single_sample_report",
                 "--forbid-text", "用于定位组织内 RNA 更新活性",
                 "--forbid-text", "备选芯片", "--json"]
    validation = stage / "05_metadata/report_validation.json"
    with validation.open("w", encoding="utf-8") as handle:
        run(validator, env=c.env, stdout=handle)
    c.finish(stage, "09", [bulk_source / "05_metadata/stage_complete.json",
                            spatial / "05_metadata/stage_complete.json"])


STAGES = {"01": stage_01, "02": stage_02, "03": stage_03, "04": stage_04,
          "05": stage_05, "06": stage_06, "07": stage_07, "08": stage_08, "09": stage_09}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--stage", choices=sorted(STAGES), required=True)
    args = parser.parse_args()
    context = StageContext(args.config)
    context.out.mkdir(parents=True, exist_ok=True)
    STAGES[args.stage](context)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
