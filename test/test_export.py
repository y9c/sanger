#!/usr/bin/env python3
"""Unit tests for cfutils export module."""

import tempfile
import unittest

from cfutils.align import call_mutations
from cfutils.export import batch_summary, to_fasta, to_json, to_vcf, write_batch
from cfutils.parser import parse_abi, parse_fasta


class TestExport(unittest.TestCase):
    """Test export formats."""

    @classmethod
    def setUpClass(cls):
        cls.query = parse_abi("./data/B5-M13R_B07.ab1")
        cls.ref = parse_fasta("./data/ref.fa")

    def test_to_fasta(self):
        s = to_fasta(self.query, name="read1", width=50)
        lines = s.strip().split("\n")
        self.assertTrue(lines[0].startswith(">read1"))
        self.assertEqual("".join(lines[1:]), self.query.seq)

    def test_to_vcf(self):
        sites = call_mutations(self.query, self.ref, report_all_sites=True)
        vcf = to_vcf(sites, sample_id="S", reference_name="3k")
        self.assertTrue(vcf.startswith("##fileformat=VCFv4.2"))
        has_variant_line = any(
            "\t.\t" not in line and not line.startswith("#")
            for line in vcf.split("\n")
        )
        self.assertTrue(has_variant_line)

    def test_to_json(self):
        js = to_json(self.query)
        self.assertIn('"length"', js)
        self.assertIn('"qc"', js)

    def test_batch_summary(self):
        rows = batch_summary([self.query, self.query])
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["bases"], len(self.query.seq))

    def test_write_batch_csv(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = write_batch([self.query], tmp, fmt="csv")
            self.assertTrue(p.endswith(".csv"))
            with open(p) as fh:
                self.assertIn("name", fh.readline())


if __name__ == "__main__":
    unittest.main()
