#!/usr/bin/env python3
"""Unit tests for the high-level Chromatogram object."""

import unittest

from cfutils.chromatogram import Chromatogram
from cfutils.parser import parse_fasta


class TestChromatogram(unittest.TestCase):
    """Test the Chromatogram domain object."""

    @classmethod
    def setUpClass(cls):
        cls.cg = Chromatogram.from_abi("./data/B5-M13R_B07.ab1")
        cls.ref = parse_fasta("./data/ref.fa")

    def test_accessors(self):
        self.assertEqual(self.cg.length, 1141)
        self.assertEqual(len(self.cg.quality), 1141)
        self.assertEqual(self.cg.traces.shape, (4, self.cg.traces.shape[1]))
        self.assertEqual(len(self.cg.peaks), 1141)
        self.assertGreater(self.cg.mean_quality, 20)
        self.assertEqual(self.cg.channels, "GATC")

    def test_qc(self):
        m = self.cg.qc()
        self.assertIn("crl", m)
        self.assertIn("snr", m)

    def test_basecall(self):
        res = self.cg.basecall()
        self.assertEqual(res.n_calls, 1141)

    def test_call_mutations(self):
        sites = self.cg.call_mutations(self.ref)
        self.assertGreater(len(sites), 0)

    def test_slice_and_trim(self):
        seg = self.cg.slice(10, 20)
        self.assertEqual(seg.length, 11)
        trimmed = self.cg.trim()
        self.assertLessEqual(trimmed.length, self.cg.length)

    def test_analyze(self):
        self.assertTrue(isinstance(self.cg.analyze("translate"), str))
        self.assertIn("EcoRI", self.cg.analyze("restriction"))

    def test_fasta_and_vcf(self):
        self.assertTrue(self.cg.to_fasta().startswith(">"))
        vcf = self.cg.to_vcf(self.ref)
        self.assertTrue(vcf.startswith("##fileformat=VCFv4.2"))

    def test_repr(self):
        self.assertIn("B5-M13R_B07", repr(self.cg))


if __name__ == "__main__":
    unittest.main()
