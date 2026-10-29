#!/usr/bin/env python3
"""Self-tests for the conserved flow renderer."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path


REPORT_SCRIPTS = Path(__file__).resolve().parents[1] / "pipeline" / "scripts" / "09_report"
SCRIPT = REPORT_SCRIPTS / "build_read_flow_sankey.py"
VALIDATOR = REPORT_SCRIPTS / "validate_read_flow_outputs.py"


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class SankeyTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory(prefix="bio-flow-test-")
        self.root = Path(self.temp.name)
        self.nodes = [
            {"view_id": "sample-a", "node_id": "raw", "stage": 0, "label": "Raw reads", "value": 100, "unit": "reads", "scope": "sample", "node_kind": "retained"},
            {"view_id": "sample-a", "node_id": "usable", "stage": 1, "label": "Usable reads", "value": 80, "unit": "reads", "scope": "sample", "node_kind": "retained"},
            {"view_id": "sample-a", "node_id": "structure_loss", "stage": 1, "label": "No complete structure", "value": 20, "unit": "reads", "scope": "sample", "node_kind": "loss"},
            {"view_id": "sample-a", "node_id": "mapped", "stage": 2, "label": "Mapped reads", "value": 60, "unit": "reads", "scope": "sample", "node_kind": "retained"},
            {"view_id": "sample-a", "node_id": "map_loss", "stage": 2, "label": "Unmapped reads", "value": 15, "unit": "reads", "scope": "sample", "node_kind": "loss"},
            {"view_id": "sample-a", "node_id": "map_qc_loss", "stage": 2, "label": "Non-primary alignment reads", "value": 5, "unit": "reads", "scope": "sample", "node_kind": "loss"},
            {"view_id": "molecules", "node_id": "called", "stage": 0, "label": "Called molecules", "value": 10, "unit": "molecules", "scope": "sample", "node_kind": "retained"},
            {"view_id": "molecules", "node_id": "passing", "stage": 1, "label": "Passing molecules", "value": 7, "unit": "molecules", "scope": "sample", "node_kind": "retained"},
            {"view_id": "molecules", "node_id": "filtered", "stage": 1, "label": "Filtered molecules", "value": 3, "unit": "molecules", "scope": "sample", "node_kind": "loss"},
        ]
        self.links = [
            {"view_id": "sample-a", "link_id": "raw_usable", "source": "raw", "target": "usable", "value": 80, "reason": "retained"},
            {"view_id": "sample-a", "link_id": "raw_loss", "source": "raw", "target": "structure_loss", "value": 20, "reason": "structure unavailable"},
            {"view_id": "sample-a", "link_id": "usable_mapped", "source": "usable", "target": "mapped", "value": 60, "reason": "primary mapped"},
            {"view_id": "sample-a", "link_id": "usable_loss", "source": "usable", "target": "map_loss", "value": 15, "reason": "unmapped"},
            {"view_id": "sample-a", "link_id": "usable_qc_loss", "source": "usable", "target": "map_qc_loss", "value": 5, "reason": "non-primary"},
            {"view_id": "molecules", "link_id": "called_passing", "source": "called", "target": "passing", "value": 7, "reason": "passes filters"},
            {"view_id": "molecules", "link_id": "called_filtered", "source": "called", "target": "filtered", "value": 3, "reason": "filtered"},
        ]
        self.views = [
            {"view_id": "sample-a", "display_name": "Sample A reads", "order": 0, "is_default": 1},
            {"view_id": "molecules", "display_name": "Sample A molecules", "order": 1, "is_default": 0},
        ]
        self.transformations = [{
            "transform_id": "reads_to_molecules", "source_view_id": "sample-a", "source_node_id": "mapped",
            "target_view_id": "molecules", "target_node_id": "called", "relation": "aggregate",
            "label": "Reads collapse to molecules", "reason": "UMI deduplication",
        }]

    def tearDown(self) -> None:
        self.temp.cleanup()

    def paths(self) -> tuple[Path, Path, Path, Path]:
        nodes = self.root / "nodes.tsv"
        links = self.root / "links.tsv"
        views = self.root / "views.tsv"
        transformations = self.root / "transformations.tsv"
        write_tsv(nodes, ["view_id", "node_id", "stage", "label", "value", "unit", "scope", "node_kind"], self.nodes)
        write_tsv(links, ["view_id", "link_id", "source", "target", "value", "reason"], self.links)
        write_tsv(views, ["view_id", "display_name", "order", "is_default"], self.views)
        write_tsv(transformations, ["transform_id", "source_view_id", "source_node_id", "target_view_id", "target_node_id", "relation", "label", "reason"], self.transformations)
        return nodes, links, views, transformations

    def run_tool(self, *extra: str, node_layout: str = "main_top_losses", orientation: str = "horizontal") -> subprocess.CompletedProcess[str]:
        nodes, links, views, transformations = self.paths()
        return subprocess.run(
            [sys.executable, str(SCRIPT), "--nodes", str(nodes), "--links", str(links), "--views", str(views), "--transformations", str(transformations), "--default-mode", "view", "--node-layout", node_layout, "--orientation", orientation, "--output-dir", str(self.root / "out"), *extra],
            text=True, capture_output=True, check=False,
        )

    def test_builds_all_outputs(self) -> None:
        result = self.run_tool("--title", "Toy conserved flow")
        self.assertEqual(result.returncode, 0, result.stderr)
        out = self.root / "out"
        expected = [
            out / "tables/read_flow_nodes.tsv", out / "tables/read_flow_links.tsv",
            out / "tables/read_flow_payload.json", out / "figures/read_flow_sankey/sample-a.svg",
            out / "tables/read_flow_transformations.tsv",
            out / "figures/read_flow_sankey/molecules.svg", out / "reports/read_flow_sankey.html",
            out / "reports/read_flow_sankey_section.html", out / "reports/sankey_manifest.tsv",
            out / "reports/sankey_validation.json",
        ]
        self.assertTrue(all(path.is_file() and path.stat().st_size > 0 for path in expected))
        svg_root = ET.parse(out / "figures/read_flow_sankey/sample-a.svg").getroot()
        self.assertEqual(svg_root.attrib["data-node-layout"], "main_top_losses")
        self.assertEqual(svg_root.attrib["data-orientation"], "horizontal")
        svg_ns = {"svg": "http://www.w3.org/2000/svg"}
        svg_style = "".join(item.text or "" for item in svg_root.findall(".//svg:style", svg_ns))
        self.assertIn("paint-order:stroke fill", svg_style)
        self.assertIn("stroke:#fff", svg_style)
        self.assertIn("stroke-width:3px", svg_style)
        svg_children = list(svg_root)
        last_band = max(
            index for index, item in enumerate(svg_children)
            if item.tag.endswith("path") and "data-link-id" in item.attrib
        )
        first_label = min(
            index for index, item in enumerate(svg_children)
            if item.tag.endswith("text") and item.attrib.get("class") == "label"
        )
        self.assertLess(last_band, first_label)
        svg_links = {
            item.attrib["data-link-id"]: item
            for item in svg_root.findall(".//svg:path", svg_ns)
            if "data-link-id" in item.attrib
        }
        svg_nodes = {
            item.attrib["data-node-id"]: item
            for item in svg_root.findall(".//svg:rect", svg_ns)
            if "data-node-id" in item.attrib
        }
        widths = {key: float(item.attrib["data-width"]) for key, item in svg_links.items()}
        self.assertAlmostEqual(widths["raw_usable"], widths["usable_mapped"] + widths["usable_loss"] + widths["usable_qc_loss"], places=5)
        self.assertAlmostEqual(widths["raw_loss"], widths["usable_loss"] + widths["usable_qc_loss"], places=5)
        self.assertAlmostEqual(float(svg_nodes["raw"].attrib["data-height"]), widths["raw_usable"] + widths["raw_loss"], places=5)
        mapped_bottom = float(svg_nodes["mapped"].attrib["y"]) + float(svg_nodes["mapped"].attrib["height"])
        self.assertAlmostEqual(float(svg_nodes["map_loss"].attrib["y"]) - mapped_bottom, 12.0, places=2)
        payload = json.loads((out / "tables/read_flow_payload.json").read_text(encoding="utf-8"))
        self.assertEqual(payload["default_view"], "sample-a")
        self.assertEqual(payload["default_mode"], "view")
        self.assertEqual(payload["node_layout"], "main_top_losses")
        self.assertEqual(payload["orientation"], "horizontal")
        self.assertEqual(payload["transformations"][0]["source_value"], 60)
        read_nodes = {row["node_id"]: row for row in payload["views"]["sample-a"]["nodes"]}
        self.assertGreater(read_nodes["raw"]["y"], 0.0)
        self.assertLess(read_nodes["raw"]["y"], read_nodes["structure_loss"]["y"])
        self.assertLess(read_nodes["raw"]["x"], read_nodes["structure_loss"]["x"])
        self.assertLess(read_nodes["structure_loss"]["x"], read_nodes["usable"]["x"])
        stage_two_losses = sorted(
            (row for row in read_nodes.values() if row["stage"] == 2 and row["node_kind"] == "loss"),
            key=lambda row: row["layout_order"],
        )
        self.assertEqual([row["value"] for row in stage_two_losses], [15, 5])
        mapped = next(row for row in payload["views"]["sample-a"]["links"] if row["link_id"] == "usable_mapped")
        self.assertAlmostEqual(mapped["pct_parent"], 75.0)
        self.assertAlmostEqual(mapped["pct_root"], 60.0)
        report = (out / "reports/read_flow_sankey.html").read_text(encoding="utf-8")
        self.assertIn("Plotly.react", report)
        self.assertIn("data.node_layout==='main_top_losses'?'fixed':'snap'", report)
        self.assertIn("nodeSpec.y=nodes.map(n=>n.y)", report)
        self.assertIn("@media print", report)
        self.assertIn("<select class='bio-flow-view-select' id='bio-flow-view'", report)
        self.assertIn("raiseSankeyLabels(el)", report)
        self.assertIn("group.querySelector('.sankey-node-set')", report)
        self.assertIn("nodes.parentNode.appendChild(nodes)", report)
        self.assertIn("coordinateRoot.appendChild(layer)", report)
        self.assertIn("bio-flow-vertical-label-layer", report)
        self.assertIn("plotly_afterplot", report)
        self.assertIn("plotView(el,view,forcedWidth)", report)
        self.assertIn("el.innerHTML=view.vertical_svg", report)
        self.assertIn("Plotly.purge(el)", report)
        self.assertIn("Plotly.newPlot", report)
        self.assertIn("installVerticalInteractions(el,view)", report)
        self.assertIn("bio-flow-svg-tooltip", report)
        self.assertIn("max-height:78vh", report)
        self.assertIn("overflow-y:auto", report)
        self.assertIn("overscroll-behavior:contain", report)
        self.assertIn("Math.min(780,Math.floor(innerHeight*.72))", report)
        self.assertNotIn("id='bio-flow-orientation'", report)
        self.assertIn("const orientation=data.orientation||'horizontal'", report)
        self.assertIn(".bio-flow-label-layer text{paint-order:stroke fill;stroke:#fff;stroke-width:4px", report)
        self.assertNotIn("svg.main-svg", report)
        self.assertNotIn("id='bio-flow-connected'", report)
        self.assertNotIn("class='bio-flow-tab'", report)
        self.assertNotIn("class='bio-flow-tabs'", report)
        self.assertNotRegex(report, r"<(?:script|link)[^>]+(?:src|href)=")
        validation = json.loads((out / "reports/sankey_validation.json").read_text(encoding="utf-8"))
        self.assertTrue(validation["ok"])
        with (out / "reports/sankey_manifest.tsv").open(encoding="utf-8", newline="") as handle:
            manifest = list(csv.DictReader(handle, delimiter="\t"))
        self.assertTrue(all(not Path(row["target"]).is_absolute() for row in manifest))
        checked = subprocess.run(
            [sys.executable, str(VALIDATOR), "--output-dir", str(out)],
            text=True, capture_output=True, check=False,
        )
        self.assertEqual(checked.returncode, 0, checked.stderr)
        self.assertTrue(json.loads(checked.stdout)["ok"])

    def test_snap_layout_remains_available(self) -> None:
        result = self.run_tool(node_layout="snap")
        self.assertEqual(result.returncode, 0, result.stderr)
        out = self.root / "out"
        payload = json.loads((out / "tables/read_flow_payload.json").read_text(encoding="utf-8"))
        self.assertEqual(payload["node_layout"], "snap")
        checked = subprocess.run(
            [sys.executable, str(VALIDATOR), "--output-dir", str(out)],
            text=True, capture_output=True, check=False,
        )
        self.assertEqual(checked.returncode, 0, checked.stderr)

    def test_vertical_orientation_builds_and_validates(self) -> None:
        result = self.run_tool(orientation="vertical")
        self.assertEqual(result.returncode, 0, result.stderr)
        out = self.root / "out"
        payload = json.loads((out / "tables/read_flow_payload.json").read_text(encoding="utf-8"))
        validation = json.loads((out / "reports/sankey_validation.json").read_text(encoding="utf-8"))
        self.assertEqual(payload["orientation"], "vertical")
        self.assertEqual(validation["orientation"], "vertical")
        svg_root = ET.parse(out / "figures/read_flow_sankey/sample-a.svg").getroot()
        self.assertEqual(svg_root.attrib["data-orientation"], "vertical")
        svg_ns = {"svg": "http://www.w3.org/2000/svg"}
        raw = next(item for item in svg_root.findall(".//svg:rect", svg_ns) if item.attrib.get("data-node-id") == "raw")
        usable = next(item for item in svg_root.findall(".//svg:rect", svg_ns) if item.attrib.get("data-node-id") == "usable")
        self.assertGreater(float(usable.attrib["y"]), float(raw.attrib["y"]))
        self.assertGreater(float(raw.attrib["width"]), float(raw.attrib["height"]))
        checked = subprocess.run(
            [sys.executable, str(VALIDATOR), "--output-dir", str(out)],
            text=True, capture_output=True, check=False,
        )
        self.assertEqual(checked.returncode, 0, checked.stderr)
        self.assertEqual(json.loads(checked.stdout)["orientation"], "vertical")

    def test_validate_only_writes_nothing(self) -> None:
        result = self.run_tool("--validate-only")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(json.loads(result.stdout)["ok"])
        self.assertFalse((self.root / "out").exists())

    def test_validate_only_needs_no_output_or_plotly(self) -> None:
        nodes, links, views, transformations = self.paths()
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--nodes", str(nodes), "--links", str(links), "--views", str(views), "--transformations", str(transformations), "--default-mode", "view", "--validate-only"],
            text=True, capture_output=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        validation = json.loads(result.stdout)
        self.assertFalse(validation["presentation_checked"])
        self.assertIsNone(validation["inputs"]["plotly"])

    def test_rejects_nonconservation(self) -> None:
        self.links[0]["value"] = 79
        result = self.run_tool("--validate-only")
        self.assertEqual(result.returncode, 2)
        self.assertIn("outgoing", result.stderr)

    def test_rejects_mixed_units(self) -> None:
        self.nodes[1]["unit"] = "cells"
        result = self.run_tool("--validate-only")
        self.assertEqual(result.returncode, 2)
        self.assertIn("mixed units", result.stderr)

    def test_rejects_noninteger(self) -> None:
        self.nodes[0]["value"] = "100.0"
        result = self.run_tool("--validate-only")
        self.assertEqual(result.returncode, 2)
        self.assertIn("non-negative integer", result.stderr)

    def test_rejects_stage_skip(self) -> None:
        self.nodes[3]["stage"] = 3
        result = self.run_tool("--validate-only")
        self.assertEqual(result.returncode, 2)
        self.assertIn("exactly one stage", result.stderr)

    def test_rejects_cycle_or_multiple_roots(self) -> None:
        self.links.append({"view_id": "sample-a", "link_id": "cycle", "source": "mapped", "target": "raw", "value": 60, "reason": "invalid cycle"})
        result = self.run_tool("--validate-only")
        self.assertEqual(result.returncode, 2)

    def test_rejects_unknown_transformation_endpoint(self) -> None:
        self.transformations[0]["target_node_id"] = "missing"
        result = self.run_tool("--validate-only")
        self.assertEqual(result.returncode, 2)
        self.assertIn("endpoint does not exist", result.stderr)

    def test_rejects_one_to_one_count_mismatch(self) -> None:
        self.transformations[0]["relation"] = "one_to_one"
        result = self.run_tool("--validate-only")
        self.assertEqual(result.returncode, 2)
        self.assertIn("equal values", result.stderr)

    def test_rejects_same_unit_transformation(self) -> None:
        self.transformations[0]["target_view_id"] = "sample-a"
        self.transformations[0]["target_node_id"] = "mapped"
        result = self.run_tool("--validate-only")
        self.assertEqual(result.returncode, 2)
        self.assertIn("different units", result.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
