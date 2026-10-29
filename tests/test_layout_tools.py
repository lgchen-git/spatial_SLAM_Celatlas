#!/usr/bin/env python3
"""Contract tests for the project-local layout audit and reversible migration."""
from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "pipeline/scripts/00_shared/audit_analysis_layout.py"
MIGRATE = ROOT / "pipeline/scripts/00_shared/migrate_layout_v2.py"
STAGE_CONTRACT = ROOT / "pipeline/scripts/00_shared/stage_contract.py"
SUBDIRS = ("01_tables", "02_figures", "03_logs", "04_tool_output", "05_metadata")


def write_layout(root: Path) -> tuple[Path, Path, Path]:
    results = root / "results"
    tmp_root = root / "tmp"
    results.mkdir()
    tmp_root.mkdir()
    config = root / "layout.yaml"
    config.write_text(
        f"""project_id: test
results_root: {results}
tmp_root: {tmp_root}
stage_subdirectories:
  - 01_tables
  - 02_figures
  - 03_logs
  - 04_tool_output
  - 05_metadata
stages:
  - {{number: \"01\", name: one, script_dir: 01_one, result_dir: 01_one, command: run, final_report: false}}
  - {{number: \"08\", name: report, script_dir: 08_report, result_dir: 08_report, command: run, final_report: true}}
""",
        encoding="utf-8",
    )
    for name in ("01_one", "08_report"):
        (root / "pipeline" / "scripts" / name).mkdir(parents=True)
    return config, results, tmp_root


class LayoutAuditTests(unittest.TestCase):
    def test_clean_layout_and_vendor_html_pass_strict(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            config, results, _ = write_layout(root)
            sample = results / "S1"
            for stage_name in ("01_one", "08_report"):
                for subdir in SUBDIRS:
                    (sample / stage_name / subdir).mkdir(parents=True, exist_ok=True)
            (sample / "01_one/04_tool_output/vendor/report.html").parent.mkdir(parents=True)
            (sample / "01_one/04_tool_output/vendor/report.html").write_text("vendor")
            (sample / "08_report/S1_report.html").write_text("customer")
            result = subprocess.run(
                [sys.executable, str(AUDIT), "--config", str(config), "--repo-root", str(root), "--sample", "S1", "--strict", "--json"],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertTrue(json.loads(result.stdout)["ok"])

    def test_audit_detects_nested_layout_misplaced_figure_and_root_log(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            config, results, _ = write_layout(root)
            sample = results / "S1"
            for stage_name in ("01_one", "08_report"):
                for subdir in SUBDIRS:
                    (sample / stage_name / subdir).mkdir(parents=True, exist_ok=True)
            nested = sample / "01_one/bin50_showcase"
            for subdir in SUBDIRS:
                (nested / subdir).mkdir(parents=True)
            (sample / "01_one/bad.png").write_bytes(b"png")
            (sample / "pipeline.log").write_text("log")
            (sample / "08_report/S1_report.html").write_text("customer")
            result = subprocess.run(
                [sys.executable, str(AUDIT), "--config", str(config), "--repo-root", str(root), "--sample", "S1", "--strict", "--json"],
                capture_output=True, text=True,
            )
            payload = json.loads(result.stdout)
            codes = {row["code"] for row in payload["issues"]}
            self.assertNotEqual(result.returncode, 0)
            self.assertTrue({"nested_stage_layout", "misplaced_figure", "root_log"}.issubset(codes))


class LayoutMigrationTests(unittest.TestCase):
    def make_legacy_sample(self, root: Path) -> tuple[Path, Path, Path]:
        config, results, tmp_root = write_layout(root)
        # The migration mapping is intentionally fixed to the production 01-08 contract.
        sample = results / "S1"
        for stage_name in ("02_barcode", "03_trim_align", "04_molecule_call", "05_matrix", "06_spatial_dynamics", "07_tissue_bulk", "08_report"):
            for subdir in SUBDIRS:
                (sample / stage_name / subdir).mkdir(parents=True, exist_ok=True)
        (sample / "02_barcode/S1.barcode_qc.json").write_text("{}")
        (sample / "02_barcode/S1.read_barcode_umi.tsv.gz").write_bytes(b"table")
        (sample / "06_spatial_dynamics/bin50_showcase/02_figures").mkdir(parents=True)
        (sample / "06_spatial_dynamics/bin50_showcase/02_figures/01_total.png").write_bytes(b"figure")
        (sample / "06_spatial_dynamics/bin50_showcase/02_figures/07_spatial_slam_seq_capability_overview_bin50.png").write_bytes(b"old")
        (sample / "07_tissue_bulk/02_figures/01_slamdunk_conversion_qc.png").write_bytes(b"slam")
        (sample / "pipeline_spatial.log").write_text("log")
        return config, results, tmp_root

    def test_default_is_dry_run_and_apply_is_reversible(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            config, _, tmp_root = self.make_legacy_sample(root)
            base = [sys.executable, str(MIGRATE), "--config", str(config), "--repo-root", str(root), "--sample", "S1", "--timestamp", "STAMP"]
            dry = subprocess.run(base, check=True, capture_output=True, text=True)
            self.assertEqual(json.loads(dry.stdout)["mode"], "dry-run")
            self.assertTrue((root / "results/S1/02_barcode/S1.barcode_qc.json").exists())

            applied = subprocess.run(base + ["--apply"], check=True, capture_output=True, text=True)
            manifest = Path(json.loads(applied.stdout)["manifest"])
            self.assertTrue((root / "results/S1/02_barcode/05_metadata/S1.barcode_qc.json").is_file())
            self.assertTrue((root / "results/S1/06_spatial_dynamics/02_figures/01_total.png").is_file())
            self.assertTrue((tmp_root / "results_legacy/S1/STAMP/06_spatial_dynamics/bin50_showcase/02_figures/07_spatial_slam_seq_capability_overview_bin50.png").is_file())
            with manifest.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(rows and all(row["sha256"] and row["status"] == "moved" for row in rows))

            subprocess.run([sys.executable, str(MIGRATE), "--rollback-manifest", str(manifest)], check=True, capture_output=True, text=True)
            self.assertTrue((root / "results/S1/02_barcode/S1.barcode_qc.json").is_file())
            self.assertTrue((root / "results/S1/06_spatial_dynamics/bin50_showcase/02_figures/01_total.png").is_file())
            self.assertTrue((root / "results/S1/pipeline_spatial.log").is_file())

    def test_different_destination_is_a_hard_conflict(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            config, _, _ = self.make_legacy_sample(root)
            destination = root / "results/S1/02_barcode/05_metadata/S1.barcode_qc.json"
            destination.write_text("different")
            result = subprocess.run(
                [sys.executable, str(MIGRATE), "--config", str(config), "--repo-root", str(root), "--sample", "S1", "--timestamp", "STAMP"],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("different destination already exists", result.stderr)


class StageContractTests(unittest.TestCase):
    def test_completion_requires_same_inputs_and_unchanged_outputs(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            stage = root / "06_spatial_dynamics"
            (stage / "01_tables").mkdir(parents=True)
            input_path = root / "input.tsv"
            input_path.write_text("input\n", encoding="utf-8")
            output = stage / "01_tables/result.tsv"
            output.write_text("result\n", encoding="utf-8")
            base = [sys.executable, str(STAGE_CONTRACT), "--stage-dir", str(stage), "--stage-id", "06", "--input", str(input_path)]
            subprocess.run(base, check=True)
            marker = json.loads((stage / "05_metadata/stage_complete.json").read_text(encoding="utf-8"))
            self.assertEqual(marker["completion_mode"], "executed")
            self.assertEqual(subprocess.run(base + ["--check"]).returncode, 0)
            output.write_text("changed\n", encoding="utf-8")
            self.assertNotEqual(subprocess.run(base + ["--check"]).returncode, 0)


if __name__ == "__main__":
    unittest.main()
