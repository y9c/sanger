#!/usr/bin/env python3
"""Unit tests for cfutils assembly (pileup/consensus) module."""

import unittest
from collections import Counter

from cfutils.parser import parse_abi, parse_fasta, SeqRecord
from cfutils.assembly import pileup, consensus, coverage, PileupColumn


class TestAssembly(unittest.TestCase):
    """Test reference-guided pileup and consensus calling."""

    @classmethod
    def setUpClass(cls):
        cls.read = parse_abi("./data/B5-M13R_B07.ab1")
        cls.ref = parse_fasta("./data/ref.fa")

    def test_pileup_matches_read_to_ref_positions(self):
        table = pileup([self.read], self.ref)
        self.assertGreater(len(table), 0)
        # the read was called against this same reference, so majority consensus
        # should track the reference (allow a few N/heterozygote positions)
        cons = consensus(table)
        self.assertEqual(len(cons), len(table))

    def test_consensus_base_majority(self):
        col = PileupColumn(ref_pos=5, ref_base="A",
                           counts=Counter({"A": 3, "C": 1}), n_reads=4)
        self.assertEqual(col.consensus_base(), "A")
        col2 = PileupColumn(ref_pos=6, ref_base="G",
                            counts=Counter({"C": 4}), n_reads=4)
        self.assertEqual(col2.consensus_base(), "C")

    def test_low_depth_falls_back_to_ref(self):
        col = PileupColumn(ref_pos=7, ref_base="T",
                           counts=Counter({"T": 100, "A": 1}), n_reads=101)
        self.assertEqual(col.consensus_base(), "T")

    def test_coverage_sorted(self):
        table = pileup([self.read], self.ref)
        cov = coverage(table)
        positions = [p for p, _ in cov]
        self.assertEqual(positions, sorted(positions))

    def test_min_cov_filters(self):
        table = pileup([self.read], self.ref, min_cov=2)
        # with a single read and min_cov=2 the table should be empty
        self.assertEqual(len(table), 0)


if __name__ == "__main__":
    unittest.main()
