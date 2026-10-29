#!/usr/bin/env python3
import gzip
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


class SnakemakeDockerContractTests(unittest.TestCase):
    def test_required_layout_and_entrypoint(self):
        for name in ("config", "inputs", "log", "outputs", "pipeline", "reference", "tmp"):
            self.assertTrue((ROOT / name).is_dir(), name)
        self.assertTrue((ROOT / "pipeline/Snakefile").is_file())

    def test_layout_has_nine_ordered_rules_and_one_report(self):
        text = (ROOT / "config/analysis_layout.yaml").read_text(encoding="utf-8")
        numbers = re.findall(r'^\s+- number: ["\'](\d{2})["\']$', text, re.M)
        self.assertEqual(numbers, [f"{i:02d}" for i in range(1, 10)])
        self.assertEqual(text.count("final_report: true"), 1)
        snakefile = (ROOT / "pipeline/Snakefile").read_text(encoding="utf-8")
        for rule in ("input_qc", "celatlas", "barcode", "trim_align", "molecule_call", "matrix", "spatial_dynamics", "tissue_bulk", "report"):
            self.assertRegex(snakefile, rf"(?m)^rule {rule}:")

    def test_docker_mount_permissions_and_single_snakemake_entrypoint(self):
        runner = (ROOT / "containers/docker/run_snakemake.sh").read_text(encoding="utf-8")
        for name in ("config", "pipeline", "inputs", "reference"):
            self.assertRegex(runner, rf"dst=/workspace/{re.escape(name)}[^\n]*readonly")
        for name in ("log", "outputs", "tmp"):
            self.assertRegex(runner, rf"dst=/workspace/{name}")
        self.assertNotIn("workflow/run.sh", runner)
        self.assertIn("--snakefile pipeline/Snakefile", runner)
        self.assertIn("--ulimit nofile=65536:65536", runner)
        self.assertIn("dst=/workspace/reuse,readonly", runner)
        self.assertIn("--reuse-outputs-host", runner)
        self.assertIn('--config "config_path=$config"', runner)

    def test_image_excludes_data_and_historical_celatlas_results(self):
        builder = (ROOT / "containers/docker/build_export_image.sh").read_text(encoding="utf-8")
        self.assertIn("--exclude=/results/", builder)
        self.assertIn("--exclude=/.git/", builder)
        dockerfile = (ROOT / "containers/docker/Dockerfile").read_text(encoding="utf-8")
        for forbidden in ("inputs/", "reference/", "outputs/", "results/"):
            self.assertNotIn(f"COPY {forbidden}", dockerfile)

    def test_run_config_locks_reuse_and_bin50(self):
        text = (ROOT / "config/run_example_run.yaml").read_text(encoding="utf-8")
        self.assertIn("sample_id: example_sample", text)
        self.assertIn("bin_size_raw: 50", text)
        self.assertIn("official_celatlas_bin50: true", text)
        self.assertIn("skip_bulk: false", text)

    def test_bin20_config_reuses_upstream_and_keeps_an_isolated_scale(self):
        config = yaml.safe_load((ROOT / "config/run_example_run-bin20.yaml").read_text(encoding="utf-8"))
        self.assertEqual(config["run_scope"], "spatial_report")
        self.assertEqual(config["source_outputs_root"], "/workspace/reuse")
        self.assertEqual(config["spatial"]["bin_size_raw"], 20)
        self.assertEqual(config["spatial"]["bin_size_um"], 20)
        self.assertEqual(config["spatial"]["bin_size_pixels"], 40)
        self.assertEqual(config["spatial"]["leiden_resolution"], 1.0)
        self.assertTrue(config["reuse_tissue_bulk"])

    def test_spatial_report_scope_uses_source_outputs_without_full_rerun(self):
        snakefile = (ROOT / "pipeline/Snakefile").read_text(encoding="utf-8")
        self.assertIn('RUN_SCOPE == "spatial_report"', snakefile)
        self.assertIn("SOURCE_SAMPLE_ROOT", snakefile)
        self.assertIn("SPATIAL_BIN_LABEL", snakefile)
        self.assertIn("BULK_SOURCE_ROOT", snakefile)

    def test_reviewed_bin20_config_locks_unique_site_and_sensitivity_controls(self):
        config = yaml.safe_load((ROOT / "config/run_example_run-bin20-reviewed.yaml").read_text(encoding="utf-8"))
        self.assertEqual(config["run_scope"], "reviewed")
        self.assertEqual(config["spatial"]["bin_size_raw"], 20)
        self.assertEqual(config["molecule_call"]["counting_method"], "unique_genomic_site_highest_base_quality")
        self.assertEqual(config["molecule_call"]["min_base_quality"], 20)
        self.assertEqual(config["molecule_call"]["sensitivity_base_qualities"], [15, 20, 25, 30])
        self.assertEqual(config["spatial"]["smoothing_max_radius_pixels"], 18)
        self.assertEqual(config["report_filename"], "example_sample_report_reviewed.html")

    def test_fastq_qname_tagging_is_synchronized(self):
        with tempfile.TemporaryDirectory() as tmp_name:
            tmp = Path(tmp_name)
            fastq = tmp / "reads.fq.gz"
            table = tmp / "barcodes.tsv.gz"
            output = tmp / "tagged.fq.gz"
            stats = tmp / "stats.json"
            with gzip.open(fastq, "wt", encoding="utf-8") as handle:
                handle.write("@read1 2:N:0:1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nIIII\n")
            with gzip.open(table, "wt", encoding="utf-8") as handle:
                handle.write(
                    "read_name\tbarcode\tumi\tx\ty\tmean_r1_qual\tmean_r2_qual\tpolyt_len\n"
                    "read1\tAACCGGTT\tAAAA\t10\t20\t35\t36\t18\n"
                    "read2\tTTGGCCAA\tCCCC\t30\t40\t35\t36\t18\n"
                )
            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "pipeline/scripts/04_trim_align/tag_fastq_headers.py"),
                    "--fastq", str(fastq), "--barcode-table", str(table),
                    "--output", str(output), "--stats", str(stats), "--threads", "1",
                ],
                check=True,
            )
            with gzip.open(output, "rt", encoding="utf-8") as handle:
                tagged = handle.read()
            self.assertIn("@read1|CB=AACCGGTT|UB=AAAA|X=10|Y=20 2:N:0:1", tagged)
            self.assertIn("@read2|CB=TTGGCCAA|UB=CCCC|X=30|Y=40", tagged)


if __name__ == "__main__":
    unittest.main()
