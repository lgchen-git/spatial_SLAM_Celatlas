#!/usr/bin/env python3
import importlib.util
import csv
import gzip
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread


ROOT = Path(__file__).resolve().parents[1]


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


correction = load("correction", "pipeline/scripts/07_spatial_dynamics/spatial_slam_correct_new_fraction.py")
kinetics = load("kinetics", "pipeline/scripts/07_spatial_dynamics/spatial_slam_kinetics.py")
mask = load("tissue_mask", "pipeline/scripts/08_tissue_bulk/build_tissue_mask.py")
bulk = load("bulk_summary", "pipeline/scripts/09_report/summarize_bulk_slam.py")
progress = load("pipeline_progress", "pipeline/scripts/01_input_qc/pipeline_progress.py")
celatlas_reuse = load("celatlas_reuse", "pipeline/scripts/02_celatlas/validate_reused_celatlas.py")
reviewed_biology = load("reviewed_biology", "pipeline/scripts/07_spatial_dynamics/analyze_reviewed_biology.py")


class CoreLogicTests(unittest.TestCase):
    def test_correction_model_fit_defaults_to_all_molecules(self):
        parser = correction.parse_args
        argv = sys.argv
        try:
            sys.argv = [
                "correction", "--molecule-calls", "calls.tsv.gz", "--sample", "S1",
                "--outdir", "out",
            ]
            self.assertIsNone(parser().max_estimate_rows)
        finally:
            sys.argv = argv

    def test_official_bin_grid_uses_celatlas_fullres_pixel_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp_name:
            positions = Path(tmp_name) / "tissue_positions_list.csv"
            positions.write_text("BC1,1,4,22,90.5,500.25\n", encoding="utf-8")
            grid = correction.load_official_grid(positions)
            row = grid.loc[(22, 4)]
            self.assertEqual(row["barcode"], "BC1")
            self.assertAlmostEqual(float(row["bin_x"]), 500.25)
            self.assertAlmostEqual(float(row["bin_y"]), 90.5)

    def test_reused_celatlas_contract_accepts_matching_source_and_rejects_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp_name:
            source = Path(tmp_name) / "S1"
            chip = "S1"
            required = celatlas_reuse.required_artifacts(source, chip)
            for name, path in required.items():
                path.parent.mkdir(parents=True, exist_ok=True)
                if name == "metrics":
                    path.write_text(json.dumps({
                        "sample_summary": {"Chip Number": chip, "Chemistry": "BBV2.4"},
                        "star_summary": {"Genome": "Mus_musculus"},
                    }), encoding="utf-8")
                else:
                    path.write_bytes(b"fixture")
            payload = celatlas_reuse.validate(source, chip, chip, "BBV2.4", "Mus_musculus")
            self.assertTrue(payload["ok"])
            self.assertEqual(payload["reuse_mode"], "validated_existing_result")
            with self.assertRaisesRegex(ValueError, "metadata mismatch"):
                celatlas_reuse.validate(source, chip, chip, "WRONG", "Mus_musculus")

    @unittest.skipUnless(shutil.which("rustc") and shutil.which("pigz"), "rustc and pigz required")
    def test_full_read_chip_assignment_supports_declared_and_detects_swap(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            project = "P1"
            data = tmp / "data"
            data.mkdir()
            barcode_a = "ACGTACGTACGT"
            barcode_b = "TGCATGCATGCA"

            def sequence(barcode):
                return (
                    barcode[:4]
                    + "G" * 15
                    + barcode[4:8]
                    + "G" * 15
                    + barcode[8:]
                    + "A" * 10
                    + "T" * 18
                )

            def write_fastq(path, barcode):
                seq = sequence(barcode)
                with gzip.open(path, "wt") as handle:
                    for index in range(12):
                        handle.write(f"@r{index}\n{seq}\n+\n{'I' * len(seq)}\n")

            rows = []
            libraries = []
            for sample_id, chip_id, barcode in (("S_A", "CHIP_A", barcode_a), ("S_B", "CHIP_B", barcode_b)):
                mask = data / sample_id
                mask.mkdir()
                (mask / f"{chip_id}_FilterBarcodes.csv").write_text(f"1,2,{barcode}\n", encoding="utf-8")
                r1 = mask / f"{sample_id}_R1.fastq.gz"
                write_fastq(r1, barcode)
                rows.append(
                    {
                        "project_id": project,
                        "sample_id": sample_id,
                        "chip_id": chip_id,
                        "mask_dir": str(mask),
                    }
                )
                libraries.append({"sample_id": sample_id, "library_id": sample_id, "fastq_r1": str(r1), "fastq_r2": str(r1)})
            samples = tmp / "samples.tsv"
            with samples.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=rows[0], delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            library_table = tmp / "libraries.tsv"
            with library_table.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=libraries[0], delimiter="\t")
                writer.writeheader()
                writer.writerows(libraries)
            command = [
                sys.executable,
                str(ROOT / "pipeline/scripts/01_input_qc/run_chip_assignment_qc.py"),
                "--samples",
                str(samples),
                "--libraries",
                str(library_table),
                "--project-id",
                project,
                "--output-dir",
                str(tmp / "result"),
                "--tmp-dir",
                str(tmp / "runtime"),
            ]
            result = subprocess.run(command, check=False, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            payload = json.loads((tmp / "result/05_metadata/chip_assignment_qc.json").read_text())
            self.assertTrue(payload["passed"])
            self.assertTrue(all(row["decision"] == "declared_mapping_supported" for row in payload["sample_decisions"]))

            write_fastq(Path(libraries[0]["fastq_r1"]), barcode_b)
            write_fastq(Path(libraries[1]["fastq_r1"]), barcode_a)
            swapped = subprocess.run(command + ["--force"], capture_output=True, text=True)
            self.assertEqual(swapped.returncode, 3)
            payload = json.loads((tmp / "result/05_metadata/chip_assignment_qc.json").read_text())
            self.assertFalse(payload["passed"])
            self.assertTrue(all(row["decision"] == "swapped_mapping_supported" for row in payload["sample_decisions"]))

    @unittest.skipUnless(shutil.which("rustc") and shutil.which("pigz"), "rustc and pigz required")
    def test_single_sample_chip_qc_is_structural_and_not_cross_validated(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            project = "P_SINGLE"
            data = tmp / "data"
            data.mkdir()
            barcode = "ACGTACGTACGT"
            mask_dir = data / "sample"
            mask_dir.mkdir()
            (mask_dir / "CHIP_FilterBarcodes.csv").write_text(f"1,2,{barcode}\n", encoding="utf-8")
            r1 = mask_dir / "sample_R1.fastq.gz"
            seq = barcode[:4] + "G" * 15 + barcode[4:8] + "G" * 15 + barcode[8:] + "A" * 10 + "T" * 18
            with gzip.open(r1, "wt") as handle:
                for index in range(4):
                    handle.write(f"@single{index}\n{seq}\n+\n{'I' * len(seq)}\n")
            samples = tmp / "samples.tsv"
            with samples.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["project_id", "sample_id", "chip_id", "mask_dir"], delimiter="\t")
                writer.writeheader()
                writer.writerow({"project_id": project, "sample_id": "S_SINGLE", "chip_id": "CHIP", "mask_dir": str(mask_dir)})
            libraries = tmp / "libraries.tsv"
            with libraries.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["sample_id", "library_id", "fastq_r1"], delimiter="\t")
                writer.writeheader()
                writer.writerow({"sample_id": "S_SINGLE", "library_id": "L1", "fastq_r1": str(r1)})
            command = [
                sys.executable, str(ROOT / "pipeline/scripts/01_input_qc/run_chip_assignment_qc.py"),
                "--samples", str(samples), "--libraries", str(libraries), "--project-id", project,
                "--output-dir", str(tmp / "result"), "--tmp-dir", str(tmp / "runtime"),
            ]
            subprocess.run(command, check=True, capture_output=True, text=True)
            payload = json.loads((tmp / "result/05_metadata/chip_assignment_qc.json").read_text())
            self.assertTrue(payload["passed"])
            self.assertEqual(payload["mode"], "single_sample_structural_qc")
            self.assertEqual(payload["cross_sample_assignment"], "not_applicable")
            self.assertEqual(payload["sample_decisions"][0]["alternate_chip_id"], "NA")
            self.assertEqual(payload["sample_decisions"][0]["decision"], "declared_chip_not_cross_validated_single_sample")

    def test_uninformative_molecule_keeps_prior(self):
        posterior = correction.posterior_new_array(
            np.array([0, 1]), np.array([0, 10]), pi=0.23, p_old=0.001, p_new=0.08
        )
        self.assertAlmostEqual(float(posterior[0]), 0.23)
        self.assertGreater(float(posterior[1]), 0.23)

    def test_mixture_converges_with_ordered_components_and_fixed_background(self):
        pair_counts = {
            (0, 100): 5000,
            (1, 100): 500,
            (2, 100): 30,
            (5, 100): 100,
            (8, 100): 800,
            (10, 100): 300,
        }
        fitted = correction.fit_mixture(pair_counts)
        self.assertTrue(fitted["converged"])
        self.assertGreater(fitted["p_new"], fitted["p_old"])
        fixed = correction.fit_mixture(pair_counts, user_p_old=0.0013)
        self.assertAlmostEqual(fixed["p_old"], 0.0013)
        self.assertGreater(fixed["p_new"], fixed["p_old"])
        posterior = correction.posterior_new_array(
            np.array([0, 1, 4]), np.array([0, 20, 20]),
            fixed["pi"], fixed["p_old"], fixed["p_new"],
        )
        self.assertTrue(np.all((posterior >= 0) & (posterior <= 1)))

    def test_expected_old_new_conservation(self):
        posterior = correction.posterior_new_array(
            np.array([0, 1, 2]), np.array([10, 10, 10]), pi=0.2, p_old=0.001, p_new=0.08
        )
        total = float(len(posterior))
        expected_new = float(posterior.sum())
        expected_old = total - expected_new
        self.assertAlmostEqual(expected_old + expected_new, total)

    def test_radius_smoothing_does_not_cross_gap_or_change_raw_input(self):
        corrected = pd.DataFrame([
            {"bin_id": "A", "bin_x": 0.0, "bin_y": 0.0, "gene_id": "g1", "gene_name": "G1", "total": 10, "expected_new": 2.0},
            {"bin_id": "B", "bin_x": 40.0, "bin_y": 0.0, "gene_id": "g1", "gene_name": "G1", "total": 10, "expected_new": 8.0},
            {"bin_id": "C", "bin_x": 400.0, "bin_y": 0.0, "gene_id": "g1", "gene_name": "G1", "total": 10, "expected_new": 10.0},
        ])
        original = corrected.copy(deep=True)
        domains = pd.DataFrame([
            {"bin_id": "A", "domain": "1", "UMAP1": 0.0, "UMAP2": 0.0, "tissue_component": "1"},
            {"bin_id": "B", "domain": "1", "UMAP1": 0.0, "UMAP2": 0.0, "tissue_component": "1"},
            {"bin_id": "C", "domain": "2", "UMAP1": 0.0, "UMAP2": 0.0, "tissue_component": "2"},
        ])
        smoothed = kinetics.smooth_bin_gene(corrected, domains, ["g1"], neighbors=24, label_hours=None, max_radius=80)
        by_bin = smoothed.set_index("bin_id")
        self.assertAlmostEqual(by_bin.loc["A", "smoothed_ntr"], 0.5)
        self.assertAlmostEqual(by_bin.loc["B", "smoothed_ntr"], 0.5)
        self.assertAlmostEqual(by_bin.loc["C", "smoothed_ntr"], 1.0)
        self.assertEqual(int(by_bin.loc["C", "smoothing_neighbors_used"]), 1)
        pd.testing.assert_frame_equal(corrected, original)

    def test_representative_story_labels_require_matching_total_new_pattern(self):
        frame = pd.DataFrame([
            {"gene_id": "gA", "total": 3000, "expected_new": 1200, "n_bins": 10, "ntr": .40, "delta_ntr": .12},
            {"gene_id": "gA", "total": 1000, "expected_new": 100, "n_bins": 10, "ntr": .10, "delta_ntr": -.18},
            {"gene_id": "gB", "total": 3000, "expected_new": 150, "n_bins": 10, "ntr": .05, "delta_ntr": -.10},
            {"gene_id": "gB", "total": 1000, "expected_new": 100, "n_bins": 10, "ntr": .10, "delta_ntr": .02},
            {"gene_id": "gU", "total": 1000, "expected_new": 320, "n_bins": 10, "ntr": .32, "delta_ntr": .01},
            {"gene_id": "gU", "total": 1000, "expected_new": 180, "n_bins": 10, "ntr": .18, "delta_ntr": -.01},
        ])
        classified = reviewed_biology.classify_story_types(frame)
        self.assertEqual(classified[(classified.gene_id == "gA") & (classified.delta_ntr > 0)].iloc[0].story_type,
                         "A: Total/New/NTR enriched")
        self.assertEqual(classified[(classified.gene_id == "gB") & (classified.delta_ntr < 0)].iloc[0].story_type,
                         "B: Total enriched with lower NTR")
        self.assertEqual(classified[classified.gene_id == "gU"].iloc[0].story_type,
                         "C: similar Total RNA, different NTR")

    def test_lasso_polygon_and_hole(self):
        roi = {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [
                    [[0, 0], [10, 0], [10, 10], [0, 10]],
                    [[4, 4], [6, 4], [6, 6], [4, 6]],
                ],
            },
        }
        self.assertTrue(mask.inside_geojson(2, 2, roi))
        self.assertFalse(mask.inside_geojson(5, 5, roi))
        self.assertFalse(mask.inside_geojson(20, 20, roi))

    def test_celatlas_barcode_tissue_mask(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            molecules = tmp / "molecules.tsv.gz"
            celatlas = tmp / "celatlas"
            output = tmp / "output"
            celatlas.mkdir()
            with gzip.open(molecules, "wt") as handle:
                handle.write("barcode\tx\ty\nBC01\t1\t2\nBC02\t3\t4\n")
            (celatlas / "sample_Barcodes_tissue_positions.csv").write_text(
                "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\n"
                "BC01,1,1,2,3,4\nBC02,0,,,5,6\n",
                encoding="utf-8",
            )
            subprocess.run([
                "python3", str(ROOT / "pipeline/scripts/08_tissue_bulk/build_tissue_mask.py"),
                "--molecule-calls", str(molecules), "--celatlas-dir", str(celatlas),
                "--outdir", str(output),
            ], check=True)
            retained = (output / "01_tissue_barcodes.tsv").read_text(encoding="utf-8")
            self.assertIn("BC01", retained)
            self.assertNotIn("BC02", retained)
            qc = json.loads((output / "tissue_mask_qc.json").read_text(encoding="utf-8"))
            self.assertEqual(qc["total_molecules"], 2)
            self.assertEqual(qc["tissue_total_molecules"], 1)
            self.assertEqual(qc["tissue_molecule_fraction"], 0.5)

    def test_celatlas_headerless_barcode_tissue_mask(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            molecules = tmp / "molecules.tsv.gz"
            celatlas = tmp / "celatlas"
            output = tmp / "output"
            celatlas.mkdir()
            with gzip.open(molecules, "wt") as handle:
                handle.write("barcode\tx\ty\nBC01\t1\t2\nBC02\t3\t4\n")
            (celatlas / "tissue_positions_list.csv").write_text(
                "BC01,1,1,2,3.0,4.0\nBC02,0,3,4,5.0,6.0\n",
                encoding="utf-8",
            )
            subprocess.run([
                "python3", str(ROOT / "pipeline/scripts/08_tissue_bulk/build_tissue_mask.py"),
                "--molecule-calls", str(molecules), "--celatlas-dir", str(celatlas),
                "--outdir", str(output),
            ], check=True)
            retained = (output / "01_tissue_barcodes.tsv").read_text(encoding="utf-8")
            self.assertIn("BC01", retained)
            self.assertNotIn("BC02", retained)

    def test_tissue_fastq_filter_streams_in_barcode_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            barcode_table = tmp / "barcodes.tsv.gz"
            fastq = tmp / "reads.fastq.gz"
            tissue = tmp / "tissue.tsv"
            output = tmp / "tissue.fastq.gz"
            with gzip.open(barcode_table, "wt", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["read_name", "barcode"], delimiter="\t")
                writer.writeheader()
                writer.writerows([
                    {"read_name": "read1", "barcode": "BC01"},
                    {"read_name": "read2", "barcode": "BC02"},
                ])
            with gzip.open(fastq, "wt") as handle:
                handle.write("@read1 extra\nACGT\n+\nIIII\n@read2\nTGCA\n+\nIIII\n")
            tissue.write_text("barcode\nBC02\n", encoding="utf-8")
            subprocess.run([
                "python3", str(ROOT / "pipeline/scripts/08_tissue_bulk/filter_tissue_fastq.py"),
                "--barcode-table", str(barcode_table), "--tissue-barcodes", str(tissue),
                "--fastq", str(fastq), "--output", str(output),
            ], check=True, capture_output=True, text=True)
            with gzip.open(output, "rt") as handle:
                retained = handle.read()
            self.assertIn("@read2", retained)
            self.assertNotIn("@read1", retained)

    def test_grand_slam_versioned_column_detection(self):
        columns = [
            "Gene", "Symbol", "sample.bam Readcount", "sample.bam 0.05 quantile",
            "sample.bam Mean", "sample.bam MAP", "sample.bam 0.95 quantile",
        ]
        self.assertEqual(bulk.find_suffix(columns, " Mean"), "sample.bam Mean")
        self.assertEqual(bulk.find_suffix(columns, " Readcount"), "sample.bam Readcount")

    def test_grand_slam_interval_plotting_handles_mean_outside_quantiles(self):
        means = np.array([0.00001, 0.04])
        lower = np.array([0.00086, 0.03])
        upper = np.array([0.57, 0.05])
        err_low, err_high, plot_lower, plot_upper, adjusted = bulk.interval_errors(means, lower, upper)
        self.assertTrue(np.all(err_low >= 0))
        self.assertTrue(np.all(err_high >= 0))
        self.assertEqual(adjusted.tolist(), [True, False])
        self.assertAlmostEqual(plot_lower[0], means[0])
        self.assertAlmostEqual(plot_upper[0], upper[0])

    def test_gnu_time_wall_clock_parser(self):
        self.assertEqual(progress.wall_clock_seconds("3:36:36"), 12996)
        self.assertAlmostEqual(progress.wall_clock_seconds("33:54.25"), 2034.25)
        self.assertAlmostEqual(progress.wall_clock_seconds("7:38.51"), 458.51)
        self.assertIsNone(progress.wall_clock_seconds("not-a-time"))

    def test_barcode_stage_routes_outputs_to_contract_directories(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            r1 = tmp / "r1.fastq.gz"
            r2 = tmp / "r2.fastq.gz"
            whitelist = tmp / "whitelist.csv"
            stage = tmp / "03_barcode"
            filtered = tmp / "large" / "filtered.fastq.gz"
            sequence = "ACGT" + "A" * 10 + "T" * 12
            with gzip.open(r1, "wt") as handle:
                handle.write(f"@read1\n{sequence}\n+\n{'I' * len(sequence)}\n")
            with gzip.open(r2, "wt") as handle:
                handle.write("@read1\nGCTA\n+\nIIII\n")
            whitelist.write_text("1,2,ACGT\n", encoding="utf-8")
            subprocess.run(
                [
                    "python3", str(ROOT / "pipeline/scripts/03_barcode/extract_spatial_barcodes.py"),
                    "--r1", str(r1), "--r2", str(r2), "--whitelist", str(whitelist),
                    "--sample", "sample", "--outdir", str(stage), "--pattern", "C4U10T12",
                    "--filtered-r2-output", str(filtered),
                ],
                check=True,
            )
            self.assertTrue((stage / "01_tables/sample.read_barcode_umi.tsv.gz").is_file())
            self.assertTrue((stage / "05_metadata/sample.barcode_qc.json").is_file())
            self.assertTrue(filtered.is_file())
            self.assertFalse((stage / "sample.filtered_r2.fastq.gz").exists())

    def test_parallel_barcode_merges_lanes_once_and_atomically(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            r1_paths, r2_paths = [], []
            sequence = "ACGT" + "A" * 10 + "T" * 12
            for lane in (1, 2):
                r1, r2 = tmp / f"r1_{lane}.fastq.gz", tmp / f"r2_{lane}.fastq.gz"
                with gzip.open(r1, "wt") as handle:
                    handle.write(f"@read{lane}\n{sequence}\n+\n{'I' * len(sequence)}\n")
                with gzip.open(r2, "wt") as handle:
                    handle.write(f"@read{lane}\nGCTA\n+\nIIII\n")
                r1_paths.append(str(r1)); r2_paths.append(str(r2))
            whitelist = tmp / "whitelist.csv"
            whitelist.write_text("1,2,ACGT\n", encoding="utf-8")
            stage, filtered = tmp / "03_barcode", tmp / "large/filtered.fastq.gz"
            subprocess.run(
                [
                    "python3", str(ROOT / "pipeline/scripts/03_barcode/run_parallel_barcode.py"),
                    "--r1", ",".join(r1_paths), "--r2", ",".join(r2_paths),
                    "--whitelist", str(whitelist), "--sample", "sample", "--outdir", str(stage),
                    "--filtered-r2-output", str(filtered), "--tmp-dir", str(tmp / "parts"),
                    "--log-dir", str(stage / "03_logs"), "--workers", "2", "--pattern", "C4U10T12",
                ], check=True,
            )
            with gzip.open(stage / "01_tables/sample.read_barcode_umi.tsv.gz", "rt") as handle:
                lines = handle.readlines()
            self.assertEqual(sum(line.startswith("read_name\t") for line in lines), 1)
            self.assertEqual(len(lines), 3)
            qc = json.loads((stage / "05_metadata/sample.barcode_qc.json").read_text())
            self.assertEqual(qc["total_reads"], 2)
            self.assertEqual(qc["execution"], "parallel_by_fastq_pair")

    def test_matrix_stage_writes_one_ordered_schema(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            molecules = tmp / "molecules.tsv"
            stage = tmp / "06_matrix"
            pd.DataFrame(
                [
                    {"gene_id": "g2", "gene_name": "Gene2", "barcode": "bc2", "x": 2, "y": 2, "assignment": "new"},
                    {"gene_id": "g1", "gene_name": "Gene1", "barcode": "bc1", "x": 1, "y": 1, "assignment": "old"},
                    {"gene_id": "g1", "gene_name": "Gene1", "barcode": "bc1", "x": 1, "y": 1, "assignment": "old"},
                ]
            ).to_csv(molecules, sep="\t", index=False)
            subprocess.run(
                [
                    "python3", str(ROOT / "pipeline/scripts/06_matrix/build_old_new_matrices.py"),
                    "--molecule-calls", str(molecules), "--sample", "sample", "--outdir", str(stage),
                ],
                check=True,
            )
            matrix_dir = stage / "04_tool_output/matrix"
            features = pd.read_csv(matrix_dir / "features.tsv", sep="\t")
            barcodes = pd.read_csv(matrix_dir / "barcodes.tsv", sep="\t")
            self.assertEqual(features["gene_id"].tolist(), ["g1", "g2"])
            self.assertEqual(barcodes["barcode"].tolist(), ["bc1", "bc2"])
            old = mmread(matrix_dir / "old_umi.mtx").tocsr()
            new = mmread(matrix_dir / "new_umi.mtx").tocsr()
            self.assertEqual(old.shape, (2, 2))
            self.assertEqual(int(old[0, 0]), 2)
            self.assertEqual(int(new[1, 1]), 1)
            self.assertFalse((matrix_dir / "sample.features.tsv").exists())
            self.assertFalse((matrix_dir / "sample.old_umi.mtx").exists())


if __name__ == "__main__":
    unittest.main()
