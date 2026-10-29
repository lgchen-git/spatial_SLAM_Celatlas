#!/usr/bin/env python3
"""Validate and record Celatlas expression-derived image coordinates."""
import argparse
import csv
import json
from itertools import chain
from pathlib import Path

from PIL import Image


def image_meta(path):
    with Image.open(path) as image:
        return {"path": str(Path(path).resolve()), "width": image.width, "height": image.height, "mode": image.mode}


def update_extrema(extrema, x, y):
    candidates = {"min_x": x, "max_x": x, "min_y": y, "max_y": y}
    for key, value in candidates.items():
        if extrema[key] is None or (key.startswith("min") and value < extrema[key]) or (key.startswith("max") and value > extrema[key]):
            extrema[key] = value


def position_summary(path, hires_scale):
    pixel_extrema = {"min_x": None, "max_x": None, "min_y": None, "max_y": None}
    array_extrema = {"min_x": None, "max_x": None, "min_y": None, "max_y": None}
    n_rows = n_tissue = 0
    pxl_is_scaled_array = None
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        first = next(reader, None)
        if first is None:
            raise ValueError(f"Celatlas position table is empty: {path}")
        expected = {
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
        }
        if expected.issubset(first):
            indexes = {name: first.index(name) for name in expected}
            rows = reader
        elif len(first) >= 6:
            indexes = {
                "barcode": 0,
                "in_tissue": 1,
                "array_row": 2,
                "array_col": 3,
                "pxl_row_in_fullres": 4,
                "pxl_col_in_fullres": 5,
            }
            rows = chain([first], reader)
        else:
            raise ValueError(f"Unrecognized Celatlas position-table schema: {path}")
        for row in rows:
            pixel_x = float(row[indexes["pxl_col_in_fullres"]])
            pixel_y = float(row[indexes["pxl_row_in_fullres"]])
            array_x = float(row[indexes["array_col"]])
            array_y = float(row[indexes["array_row"]])
            if pxl_is_scaled_array is None:
                pxl_is_scaled_array = (
                    abs(pixel_x - array_x * hires_scale) <= 2
                    and abs(pixel_y - array_y * hires_scale) <= 2
                )
            n_rows += 1
            n_tissue += str(row[indexes["in_tissue"]]) == "1"
            update_extrema(pixel_extrema, pixel_x, pixel_y)
            update_extrema(array_extrema, array_x, array_y)
    fullres = array_extrema if pxl_is_scaled_array else pixel_extrema
    return {
        "n_positions": n_rows,
        "n_tissue_positions": n_tissue,
        "fullres_coordinate_source": "array_col,array_row" if pxl_is_scaled_array else "pxl_col_in_fullres,pxl_row_in_fullres",
        "rendered_pixel_coordinates_are_hires_scaled": bool(pxl_is_scaled_array),
        "rendered_pixel_extrema": pixel_extrema,
        "array_coordinate_extrema": array_extrema,
        **{key: int(value) for key, value in fullres.items()},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--positions", required=True)
    parser.add_argument("--registration-image", required=True)
    parser.add_argument("--hires-image", required=True)
    parser.add_argument("--lowres-image", required=True)
    parser.add_argument("--scalefactors", required=True)
    parser.add_argument("--count-detail", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    registration = image_meta(args.registration_image)
    hires = image_meta(args.hires_image)
    lowres = image_meta(args.lowres_image)
    scales = json.loads(Path(args.scalefactors).read_text(encoding="utf-8"))
    positions = position_summary(args.positions, scales["tissue_hires_scalef"])
    within_bounds = (
        0 <= positions["min_x"] <= positions["max_x"] < registration["width"]
        and 0 <= positions["min_y"] <= positions["max_y"] < registration["height"]
    )
    expected_hires = [registration["width"] * scales["tissue_hires_scalef"], registration["height"] * scales["tissue_hires_scalef"]]
    expected_lowres = [registration["width"] * scales["tissue_lowres_scalef"], registration["height"] * scales["tissue_lowres_scalef"]]
    provenance = {
        "classification": "expression-derived",
        "classification_note": "Generated from gene-expression count density and Celatlas segmentation; no histology, staining, or microscopy image was supplied.",
        "canonical_orientation": {
            "origin": "upper-left",
            "x": "pixel column; increases right",
            "y": "pixel row; increases downward",
            "physical_overlay_convention": "X increases right; Y increases downward",
        },
        "transforms": ["gene-expression rendering", "registration to raw chip coordinates", "tissue segmentation", "hires/lowres resize"],
        "registration_image": registration,
        "hires_image": {**hires, "scale_from_fullres": scales["tissue_hires_scalef"], "expected_dimensions": expected_hires},
        "lowres_image": {**lowres, "scale_from_fullres": scales["tissue_lowres_scalef"], "expected_dimensions": expected_lowres},
        "position_table": {"path": str(Path(args.positions).resolve()), **positions, "all_landmarks_within_fullres": within_bounds},
        "landmark_checks": [
            {"name": "left_extent", "fullres_x": positions["min_x"]},
            {"name": "right_extent", "fullres_x": positions["max_x"]},
            {"name": "top_extent", "fullres_y": positions["min_y"]},
            {"name": "bottom_extent", "fullres_y": positions["max_y"]},
        ],
        "source_count_detail": str(Path(args.count_detail).resolve()),
        "interpretive_limit": "This background supports expression/tissue coverage overlays but is not morphological or histological evidence.",
    }
    if not within_bounds:
        raise ValueError("Celatlas full-resolution positions fall outside the registration image")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(provenance, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
