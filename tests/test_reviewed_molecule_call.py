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
class ReviewedMoleculeCallTests(unittest.TestCase):
    def test_unique_sites_quality_conflicts_snp_mask_and_negative_strand(self):
        with tempfile.TemporaryDirectory() as tmp_name:
            tmp = Path(tmp_name)
            annotation = tmp / "annotation"
            annotation.mkdir()
            header = "chrom\tstart\tend\tgene_id\tgene_name\tstrand\tfeature_class\n"
            (annotation / "exons.tsv").write_text(
                header
                + "chr1\t0\t90\tgene_plus\tPlus\t+\texon\n"
                + "chr1\t90\t200\tgene_minus\tMinus\t-\texon\n",
                encoding="utf-8",
            )
            (annotation / "introns.tsv").write_text(header, encoding="utf-8")
            fasta = tmp / "genome.fa"
            fasta.write_text(">chr1\n" + "T" * 90 + "A" * 110 + "\n", encoding="utf-8")
            pysam.faidx(str(fasta))

            bam = tmp / "reads.bam"
            bam_header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 200}]}

            def add(target, name, umi, start, sequence, quality, reverse=False):
                read = pysam.AlignedSegment()
                read.query_name = f"INST:1:FC:L01:{name}|CB=BARCODE|UB={umi}|X=10|Y=20"
                read.query_sequence = sequence
                read.flag = 16 if reverse else 0
                read.reference_id = 0
                read.reference_start = start
                read.mapping_quality = 60
                read.cigar = ((0, len(sequence)),)
                read.query_qualities = pysam.qualitystring_to_array(chr(quality + 33) * len(sequence))
                target.write(read)

            with pysam.AlignmentFile(bam, "wb", header=bam_header) as target:
                # Same molecule and same genomic sites: the Q30 C observation wins over Q20 T.
                add(target, "r1", "AAAA", 10, "CCCCC", 30)
                add(target, "r2", "AAAA", 10, "TTTTT", 20)
                # Equal-quality conflicting bases make this site ambiguous and excluded.
                add(target, "r3", "CCCC", 30, "C", 30)
                add(target, "r4", "CCCC", 30, "T", 30)
                # Same molecule, different positions: all sites remain distinct.
                add(target, "r5", "GGGG", 50, "TTT", 30)
                add(target, "r6", "GGGG", 53, "TTT", 30)
                # Negative-strand A>G is transcript-oriented T>C.
                add(target, "r7", "TTTT", 100, "GGGGG", 30, reverse=True)
            pysam.index(str(bam))

            vcf = tmp / "mask.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.1\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n"
                "chr1\t11\t.\tT\tC\t.\tPASS\t.\tGT:DP:RD:AD:FREQ:ADF:ADR\t1/1:20:1:19:95%:10:9\n",
                encoding="utf-8",
            )
            barcode_table = tmp / "barcodes.tsv.gz"
            with gzip.open(barcode_table, "wt", encoding="utf-8") as handle:
                handle.write("read_name\tbarcode\tumi\tx\ty\n")
            outdir = tmp / "result"
            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "pipeline/scripts/05_molecule_call/call_new_old_molecules_reviewed.py"),
                    "--bam", str(bam), "--fasta", str(fasta),
                    "--annotation-dir", str(annotation), "--barcode-table", str(barcode_table),
                    "--sample", "sample", "--outdir", str(outdir),
                    "--tmp-dir", str(tmp / "parts"), "--threads", "1", "--progress-every", "0",
                    "--min-base-quality", "20", "--sensitivity-base-qualities", "15,20,25,30",
                    "--snp-vcf", str(vcf),
                ],
                check=True,
            )
            with gzip.open(outdir / "01_tables/sample.molecule_calls.tsv.gz", "rt", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 4)
            by_umi = {row["umi"]: row for row in rows}
            repeated = by_umi["AAAA"]
            self.assertEqual(repeated["supporting_reads"], "2")
            self.assertEqual(repeated["raw_t_covered"], "10")
            self.assertEqual(repeated["raw_tc_events"], "5")
            self.assertEqual(repeated["unique_t_covered_bq0"], "5")
            self.assertEqual(repeated["duplicate_t_evaluations"], "5")
            self.assertEqual(repeated["masked_snp_sites"], "1")
            self.assertEqual(repeated["t_covered"], "4")
            self.assertEqual(repeated["tc_events"], "4")

            conflict = by_umi["CCCC"]
            self.assertEqual(conflict["ambiguous_t_sites"], "1")
            self.assertEqual(conflict["t_covered"], "0")
            self.assertEqual(conflict["tc_events"], "0")

            distinct = by_umi["GGGG"]
            self.assertEqual(distinct["duplicate_t_evaluations"], "0")
            self.assertEqual(distinct["t_covered"], "6")
            self.assertEqual(distinct["tc_events"], "0")

            negative = by_umi["TTTT"]
            self.assertEqual(negative["t_covered"], "5")
            self.assertEqual(negative["tc_events"], "5")
            self.assertEqual(negative["assignment"], "new")

            summary = json.loads((outdir / "05_metadata/molecule_site_duplication_summary.json").read_text())
            self.assertEqual(summary["molecules"], 4)
            self.assertEqual(summary["molecules_with_duplicated_sites"], 2)
            self.assertEqual(summary["snp_mask_records"], 1)
            for name in (
                "spatial_substitution_spectrum.tsv", "tc_by_read_position.tsv",
                "tc_by_base_quality.tsv", "lane_tc_background.tsv",
            ):
                self.assertTrue((outdir / "05_metadata" / name).is_file(), name)


if __name__ == "__main__":
    unittest.main()
