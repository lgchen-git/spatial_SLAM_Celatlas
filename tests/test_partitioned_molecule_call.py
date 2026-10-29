#!/usr/bin/env python3
import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

try:
    import pysam
except ImportError:
    pysam = None


ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(pysam is not None, "pysam is required")
class PartitionedMoleculeCallTests(unittest.TestCase):
    def test_qname_metadata_and_cross_read_molecule_aggregation(self):
        with tempfile.TemporaryDirectory() as tmp_name:
            tmp = Path(tmp_name)
            annotation = tmp / "annotation"
            annotation.mkdir()
            header = "chrom\tstart\tend\tgene_id\tgene_name\tstrand\tfeature_class\n"
            (annotation / "exons.tsv").write_text(
                header + "chr1\t0\t500\tgene1\tGene1\t+\texon\n",
                encoding="utf-8",
            )
            (annotation / "introns.tsv").write_text(header, encoding="utf-8")
            fasta = tmp / "genome.fa"
            fasta.write_text(">chr1\n" + "T" * 1000 + "\n", encoding="utf-8")
            pysam.faidx(str(fasta))
            bam = tmp / "reads.bam"
            bam_header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
            with pysam.AlignmentFile(bam, "wb", header=bam_header) as target:
                for index, (umi, sequence) in enumerate((("AAAA", "C" * 20), ("AAAA", "T" * 20), ("CCCC", "T" * 20))):
                    read = pysam.AlignedSegment()
                    read.query_name = f"read{index}|CB=AACCGGTT|UB={umi}|X=10|Y=20"
                    read.query_sequence = sequence
                    read.flag = 0
                    read.reference_id = 0
                    read.reference_start = 10 + index * 30
                    read.mapping_quality = 60
                    read.cigar = ((0, 20),)
                    read.query_qualities = pysam.qualitystring_to_array("I" * 20)
                    target.write(read)
            pysam.index(str(bam))
            barcode_table = tmp / "barcodes.tsv.gz"
            with gzip.open(barcode_table, "wt", encoding="utf-8") as handle:
                handle.write("read_name\tbarcode\tumi\tx\ty\n")
            outdir = tmp / "result"
            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "pipeline/scripts/05_molecule_call/call_new_old_molecules.py"),
                    "--bam", str(bam), "--fasta", str(fasta),
                    "--annotation-dir", str(annotation), "--barcode-table", str(barcode_table),
                    "--sample", "sample", "--outdir", str(outdir),
                    "--tmp-dir", str(tmp / "parts"), "--threads", "2", "--progress-every", "0",
                ],
                check=True,
            )
            with gzip.open(outdir / "01_tables/sample.molecule_calls.tsv.gz", "rt", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 2)
            by_umi = {row["umi"]: row for row in rows}
            self.assertEqual(by_umi["AAAA"]["supporting_reads"], "2")
            self.assertEqual(by_umi["AAAA"]["assignment"], "new")
            self.assertEqual(by_umi["CCCC"]["assignment"], "old")
            qc = json.loads((outdir / "05_metadata/sample.qc.json").read_text(encoding="utf-8"))
            self.assertEqual(qc["total_bam_reads"], 3)
            self.assertEqual(qc["molecules_total"], 2)
            self.assertEqual(qc.get("missing_barcode_metadata", 0), 0)


if __name__ == "__main__":
    unittest.main()
