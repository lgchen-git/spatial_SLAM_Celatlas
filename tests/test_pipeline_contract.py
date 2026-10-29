#!/usr/bin/env python3
"""Regression checks for reviewed pipeline defaults and report visibility."""
import ast
import base64
import csv
import gzip
import importlib.util
import json
import re
import subprocess
import sys
import tempfile
import unittest
from html.parser import HTMLParser
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
REMOVED_OLD_NEW_CLUSTER_FIGURES = {
    "07_cluster_new_fraction_box.png",
    "08_cluster_new_fraction_bar.png",
    "09_gene_new_fraction_rank.png",
    "10_cluster_dynamic_gene_heatmap.png",
    "11_old_new_cluster_overlap.png",
}
REMOVED_BULK_FIGURES = {
    "01_bulk_slam_gene_ntr_summary.png",
    "02_grand_slam_mismatch_spectrum.png",
    "04_slamdunk_grand_slam_tc_rate_comparison.png",
    "01_slamdunk_conversion_qc.png",
    "02_grand_slam_gene_ntr.png",
    "03_grand_slam_model_qc.png",
    "04_tc_mismatch_enrichment_over_background.png",
}
BULK_FIGURES = {
    "02_slamdunk_conversion_qc.png",
    "03_grand_slam_gene_ntr.png",
    "04_grand_slam_model_qc.png",
    "05_tc_mismatch_enrichment_over_background.png",
}
REMOVED_SPATIAL_REPORT_FIGURES = {
    "06_representative_dynamic_genes_spatial_ntr_bin50.png",
    "07_spatial_slam_seq_capability_overview_bin50.png",
    "08_representative_dynamic_genes_old_new_ntr_bin50.png",
}
REPORT_CHAPTERS = [
    "项目概览与核心结论",
    "守恒数据流与关键损失",
    "数据质量与Celatlas基线",
    "Total RNA空间结构",
    "RNA年龄分解与校正NTR",
    "空间动力学与代表基因",
    "Tissue-only bulk SLAM",
    "方法、限制与运行信息",
]
REPORT_METADATA_OPTIONS = {
    "--chip-assignment-qc",
    "--sample-metadata",
    "--domain-kinetics",
    "--resource-summary",
    "--software-versions",
    "--interpretation-mode",
    "--analysis-markdown",
    "--read-flow-section",
    "--read-flow-plotly",
    "--read-flow-manifest",
    "--footer-logo",
}
REPORT_FIGURES = [
    "01_total_rna_spatial_domains_bin50.png",
    "02_total_rna_umap_bin50.png",
    "03_total_old_new_spatial_decomposition_bin50.png",
    "04_corrected_ntr_spatial_bin50.png",
    "05_domain_ntr_and_dynamic_gene_heatmap_bin50.png",
    "06_representative_dynamic_genes_old_new_ntr_bin50.png",
    *sorted(BULK_FIGURES),
]


class HeadingParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.in_h2 = False
        self.current = []
        self.h2 = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() == "h2":
            self.in_h2 = True
            self.current = []

    def handle_data(self, data):
        if self.in_h2:
            self.current.append(data)

    def handle_endtag(self, tag):
        if tag.lower() == "h2" and self.in_h2:
            self.h2.append("".join(self.current).strip())
            self.in_h2 = False


def normalized_chapter(title):
    return re.sub(r"^\s*\d+[.、:]?\s*", "", title).replace(" ", "")


def load_report_builder():
    path = ROOT / "pipeline/scripts/09_report/build_spatial_slam_combined_report.py"
    spec = importlib.util.spec_from_file_location("report_builder_contract", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def scanpy_assignment_default():
    path = ROOT / "pipeline/scripts/07_spatial_dynamics/spatial_slam_scanpy_analysis.py"
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
            continue
        if node.func.attr != "add_argument" or not node.args:
            continue
        if isinstance(node.args[0], ast.Constant) and node.args[0].value == "--assignments":
            defaults = [kw.value for kw in node.keywords if kw.arg == "default"]
            if len(defaults) != 1 or not isinstance(defaults[0], ast.Constant):
                raise AssertionError("--assignments must have one literal default")
            return defaults[0].value
    raise AssertionError("--assignments argument not found")


class PipelineContractTests(unittest.TestCase):
    def test_test_iaa_config_locks_current_roots_and_samples(self):
        config = yaml.safe_load((ROOT / "config/analysis_test_iaa_20260821.yaml").read_text(encoding="utf-8"))
        layout = yaml.safe_load((ROOT / "config/analysis_layout_test_iaa_20260821.yaml").read_text(encoding="utf-8"))
        with (ROOT / "config/samples_test_iaa_20260821.tsv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(config["interpretation_mode"], "standard")
        self.assertEqual(layout["results_root"], "/opt/project/data/Celatlas_spatialSLAM/Test-IAA-20260821/results")
        self.assertEqual(layout["tmp_root"], "/opt/project/data/Celatlas_spatialSLAM/Test-IAA-20260821/tmp")
        self.assertEqual({row["sample_id"] for row in rows}, {"example_sample_a", "example_sample_b"})
        self.assertTrue(all(row["labeling_status"] == "labeled" for row in rows))

    def test_unlabeled_control_report_uses_background_not_dynamics_claims(self):
        report = load_report_builder()
        report.TITLES.update(report.CONTROL_TITLES)
        report.CAPTIONS.update(report.CONTROL_CAPTIONS)
        rendered = report.build_report(
            "S1", {"project_id": "P", "species": "Mus_musculus", "chemistry": "BBV2.4"},
            {"assignment_totals": {"old": 90, "new": 10}, "total_molecules": 100, "n_bins": 4, "domains": ["1"]},
            {}, {"molecules_total": 100, "molecules_new": 10}, {}, {}, {}, {}, [], [], [],
            "unlabeled_control",
        )
        self.assertIn("<title>Spatial SLAM 分析报告</title>", rendered)
        self.assertIn("Unlabeled spatial-RNA negative control", rendered)
        self.assertIn("无标记背景分解与表观 NTR", rendered)
        self.assertIn("不代表新生 RNA", rendered)
        self.assertNotIn("RNA 年龄增加新的信息轴", rendered)
        self.assertNotIn("用于定位组织内 RNA 更新活性", rendered)

    def test_reviewed_spatial_defaults_are_configuration_driven(self):
        config = yaml.safe_load((ROOT / "config/analysis.yaml").read_text(encoding="utf-8"))
        spatial = config["spatial"]
        self.assertEqual(spatial["bin_size_raw"], 50)
        self.assertEqual(spatial["bin_size_um"], 50)
        self.assertEqual(spatial["bin_size_pixels"], 100)
        self.assertTrue(spatial["official_celatlas_bin50"])
        self.assertEqual(spatial["hvg"], 3000)
        self.assertEqual(spatial["hvg_flavor"], "seurat")
        self.assertEqual(spatial["pca_components"], 25)
        self.assertEqual(spatial["neighbors"], 15)
        self.assertEqual(spatial["leiden_resolution"], 1.0)
        self.assertEqual(spatial["plotting"]["spatial_point_size"], 14)

    def test_scanpy_defaults_to_total_only(self):
        self.assertEqual(scanpy_assignment_default(), "total")

    def test_removed_old_new_cluster_figures_cannot_enter_report(self):
        report = load_report_builder()
        for mapping in (report.FIGURE_ORDER, report.TITLES, report.CAPTIONS):
            self.assertTrue(REMOVED_OLD_NEW_CLUSTER_FIGURES.isdisjoint(mapping))
            self.assertTrue(REMOVED_SPATIAL_REPORT_FIGURES.isdisjoint(mapping))

    def test_report_uses_four_current_bulk_figures_only(self):
        report = load_report_builder()
        for mapping in (report.FIGURE_ORDER, report.TITLES, report.CAPTIONS):
            self.assertTrue(REMOVED_BULK_FIGURES.isdisjoint(mapping))
            self.assertTrue(BULK_FIGURES.issubset(mapping))

    def test_bulk_summary_cli_exposes_reviewed_thresholds(self):
        result = subprocess.run(
            [sys.executable, str(ROOT / "pipeline/scripts/09_report/summarize_bulk_slam.py"), "--help"],
            check=True,
            capture_output=True,
            text=True,
        )
        for option in ("--slamdunk-min-reads", "--grand-min-reads", "--top-genes", "--spectrum-table"):
            self.assertIn(option, result.stdout)

    def test_report_cli_exposes_customer_report_metadata_inputs(self):
        result = subprocess.run(
            [sys.executable, str(ROOT / "pipeline/scripts/09_report/build_spatial_slam_combined_report.py"), "--help"],
            check=True,
            capture_output=True,
            text=True,
        )
        for option in REPORT_METADATA_OPTIONS:
            self.assertIn(option, result.stdout)
        self.assertIn("--bin-label", result.stdout)
        self.assertIn("--bin-size-um", result.stdout)

    def test_report_builder_registers_bin20_assets_without_bin50_residue(self):
        report = load_report_builder()
        report.configure_bin("bin20", 20)
        spatial_names = [name for name in report.FIGURE_ORDER if "bin20" in name]
        self.assertEqual(len(spatial_names), 6)
        self.assertTrue(all("bin20" in name for name in spatial_names))
        self.assertTrue(all("bin50" not in name for name in spatial_names))

    def test_customer_report_has_seven_chapters_pagination_and_ten_figures(self):
        one_pixel_png = base64.b64decode(
            "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII="
        )
        with tempfile.TemporaryDirectory() as tmp_name:
            tmp = Path(tmp_name)
            spatial_figures = tmp / "spatial_figures"
            bulk_figures = tmp / "bulk_figures"
            spatial_figures.mkdir()
            bulk_figures.mkdir()
            for name in REPORT_FIGURES:
                target = bulk_figures / name if name in BULK_FIGURES else spatial_figures / name
                target.write_bytes(one_pixel_png)
                target.with_suffix(".pdf").write_bytes(b"%PDF-1.4\n%%EOF\n")

            def write_json(name, payload):
                path = tmp / name
                path.write_text(json.dumps(payload), encoding="utf-8")
                return path

            spatial_metrics = write_json(
                "spatial_metrics.json",
                {
                    "assignment_totals": {"old": 80, "new": 20},
                    "total_molecules": 100,
                    "bin_size_raw": 50,
                    "bin_size_um": 25,
                    "n_bins": 12,
                },
            )
            barcode_qc = write_json("barcode_qc.json", {"total_reads": 200, "kept_reads": 150})
            call_qc = write_json("call_qc.json", {"tc_rate": 0.02})
            bulk_summary = write_json(
                "bulk_summary.json",
                {
                    "tissue_mask": {"method": "manual_lasso", "n_tissue_barcodes": 10, "tissue_total_molecules": 77},
                    "slamdunk": {"n_tcount_rows": 5},
                    "grand_slam": {"available": True, "n_genes": 4, "warning": "exploratory"},
                    "tc_mismatch_enrichment": {
                        "complete_12_channel_spectrum": [
                            {"region": region, "channel": channel, "coverage": 1000, "mismatches": 10, "mismatch_rate": 0.01}
                            for region in ("Exonic", "Intronic")
                            for channel in ("A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G")
                        ]
                    },
                },
            )
            celatlas_metrics = write_json("celatlas_metrics.json", {})
            chip_assignment_qc = write_json("chip_assignment_qc.json", {
                "sample_decisions": [{
                    "sample_id": "S1", "declared_chip_id": "S1", "alternate_chip_id": "S2",
                    "declared_exact_poly_t_rate": 0.42, "alternate_exact_poly_t_rate": 0.08,
                    "declared_corrected_poly_t_rate": 0.71, "alternate_corrected_poly_t_rate": 0.20,
                    "declared_over_alternate_exact_fold": 5.25,
                }]
            })
            sample_metadata = tmp / "sample_metadata.tsv"
            sample_metadata.write_text(
                "sample_id\tspecies\tchemistry\tcoordinate_um_per_unit\nS1\tMus_musculus\tBBV2.4\t0.5\n",
                encoding="utf-8",
            )
            domain_kinetics = tmp / "domain_kinetics.tsv"
            domain_kinetics.write_text(
                "domain\tgene_name\ttotal\tntr\tdelta_ntr\tdynamic_score\n"
                + "".join(f"D{i % 3 + 1}\tGene{i:02d}\t{100 + i}\t0.{i % 9}\t0.0{i % 7}\t{i + 1}\n" for i in range(25)),
                encoding="utf-8",
            )
            resource_summary = tmp / "resource_summary.tsv"
            resource_summary.write_text(
                "stage_id\tstage_name\telapsed_seconds\tmax_rss_kb\tcpu_percent\n06\tspatial\t120\t2048\t400\n",
                encoding="utf-8",
            )
            software_versions = tmp / "software_versions.tsv"
            software_versions.write_text(
                "tool\tversion\nCelatlas\t1.8.0\nslamdunk\t0.4.3\nGRAND-SLAM\t2.0.7b\n",
                encoding="utf-8",
            )
            spectrum_table = tmp / "04_complete_12_channel_mismatch_spectrum.tsv"
            spectrum_table.write_text("region\tchannel\tcoverage\tmismatches\tmismatch_rate\nExonic\tT>C\t1000\t10\t0.01\n", encoding="utf-8")
            output = tmp / "report.html"
            footer_logo = tmp / "epi-logo.png"
            footer_logo.write_bytes(base64.b64decode("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNk+M/wHwAEAQH/69iP5wAAAABJRU5ErkJggg=="))
            manifest = tmp / "delivery_manifest.tsv"
            command = [
                sys.executable,
                str(ROOT / "pipeline/scripts/09_report/build_spatial_slam_combined_report.py"),
                "--sample", "S1",
                "--figure-dir", str(spatial_figures),
                "--bulk-figure-dir", str(bulk_figures),
                "--spatial-metrics", str(spatial_metrics),
                "--barcode-qc", str(barcode_qc),
                "--call-qc", str(call_qc),
                "--bulk-summary", str(bulk_summary),
                "--bulk-spectrum-table", str(spectrum_table),
                "--celatlas-metrics", str(celatlas_metrics),
                "--chip-assignment-qc", str(chip_assignment_qc),
                "--sample-metadata", str(sample_metadata),
                "--domain-kinetics", str(domain_kinetics),
                "--resource-summary", str(resource_summary),
                "--software-versions", str(software_versions),
                "--footer-logo", str(footer_logo),
                "--output", str(output),
                "--manifest", str(manifest),
            ]
            result = subprocess.run(command, check=False, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            html = output.read_text(encoding="utf-8")
            headings = HeadingParser()
            headings.feed(html)
            self.assertEqual(
                [normalized_chapter(title) for title in headings.h2],
                [normalized_chapter(title) for title in REPORT_CHAPTERS],
            )
            self.assertIn("<table", html.lower())
            self.assertRegex(html.lower(), r"pagination")
            self.assertRegex(html, r"(?i)(paginate|renderPage|changePage)")
            self.assertEqual(html.count("data:image/png;base64,"), 11)
            self.assertIn("<title>Spatial SLAM 分析报告</title>", html)
            self.assertIn("<h1>Spatial SLAM 分析报告</h1>", html)
            self.assertEqual(html.count('figure-card paired'), 4)
            self.assertNotIn("full portrait", html)
            self.assertIn("<summary>参考文献</summary>", html)
            self.assertIn("10.1093/bioinformatics/bty256", html)
            self.assertIn('class="footer-logo"', html)
            self.assertIn("<div class=\"metric-label\">Tissue molecules</div><div class=\"metric-value\">77</div>", html)
            self.assertIn("<td>当前芯片</td><td>S1</td>", html)
            self.assertIn("当前芯片 exact+PolyT 命中率", html)
            self.assertNotIn("备选芯片", html)
            self.assertNotIn("<td>S2</td>", html)
            self.assertRegex(html, r"(?i)no-?4su")
            self.assertIn("仅作为探索性相对信号解释", html)
            self.assertNotRegex(html, r"(?i)\bfigure\s+[0-9x]+\b")
            self.assertNotIn("figure-kicker", html)
            with manifest.open(encoding="utf-8", newline="") as handle:
                manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(manifest_rows), 22)
            self.assertEqual([row["category"] for row in manifest_rows], ["01_tables"] + ["02_figures"] * 20 + ["09_report"])
            expected_names = [spectrum_table.name]
            for name in REPORT_FIGURES:
                expected_names.extend((name, str(Path(name).with_suffix(".pdf"))))
            expected_names.append(output.name)
            self.assertEqual([Path(row["target"]).name for row in manifest_rows], expected_names)
            self.assertTrue(all(not Path(row["target"]).is_absolute() for row in manifest_rows))
            self.assertTrue(all(int(row["bytes"]) > 0 for row in manifest_rows))

    def test_spatial_read_flow_uses_celatlas_read_pair_view_only(self):
        with tempfile.TemporaryDirectory() as tmp_name:
            tmp = Path(tmp_name)
            metrics = tmp / "metrics.json"
            metrics.write_text(json.dumps({
                "barcode_summary": {"Raw Reads": 100, "Exactly Matched Reads": 40, "Mismatched Reads": 30, "No Barcode Reads": 30, "Valid Reads": 70},
                "cutadapt_summary": {"Reads Written": 68, "Reads Too Short": 2},
                "star_summary": {"Uniquely Mapped Reads": 50, "Multi-Mapped Reads": 8, "Unmapped Reads - Too Short": 7, "Unmapped Reads - Other": 3},
                "featureCounts_summary": {"Reads Assigned To Exonic Regions": 35, "Reads Assigned To Intronic Regions": 5, "Reads Assigned To Intergenic Regions": 8, "Reads Unassigned Ambiguity": 2},
            }), encoding="utf-8")
            output = tmp / "flow"
            stale_dir = output / "figures/read_flow_sankey"
            stale_dir.mkdir(parents=True)
            (stale_dir / "obsolete_view.svg").write_text("<svg/>", encoding="utf-8")
            result = subprocess.run([
                sys.executable, str(ROOT / "pipeline/scripts/09_report/build_spatial_read_flow.py"),
                "--sample", "S1", "--celatlas-metrics", str(metrics),
                "--output-dir", str(output),
                "--plotly-js", str(ROOT / "pipeline/scripts/09_report/read_flow_assets/plotly-3.5.1.min.js"),
            ], check=False, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            payload = json.loads((output / "tables/read_flow_payload.json").read_text(encoding="utf-8"))
            self.assertEqual(list(payload["views"]), ["celatlas_read_pairs"])
            self.assertEqual(payload["views"]["celatlas_read_pairs"]["root_value"], 100)
            self.assertEqual(payload["transformations"], [])
            self.assertFalse((stale_dir / "obsolete_view.svg").exists())
            self.assertEqual(len(list((output / "figures/read_flow_sankey").glob("*.svg"))), 1)


if __name__ == "__main__":
    unittest.main()
