#!/usr/bin/env python3
"""Unit tests for cfutils analysis (sequence biology) module."""

import unittest

from cfutils.analysis import (
    translate, find_motifs, restriction_sites, reverse_complement, gc_windows,
)

class TestAnalysis(unittest.TestCase):
    """Test sequence-level helpers."""

    def test_translate_frame1(self):
        # ATG=Met, TTT=Phe, TAA=stop
        self.assertEqual(translate("ATGTTTTAA", frame=1), "MFX")
        self.assertEqual(translate("ATGTTT", frame=1), "MF")

    def test_translate_frame2(self):
        self.assertEqual(translate("AATGTTT", frame=2), "MF")

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("ATGC"), "GCAT")
        self.assertEqual(reverse_complement("ACGTN"), "NACGT")
        # IUPAC
        self.assertEqual(reverse_complement("R"), "Y")
        self.assertEqual(reverse_complement("N"), "N")

    def test_find_motifs(self):
        self.assertEqual(find_motifs("GCATGC", "GCA"), [1])
        self.assertEqual(find_motifs("GCAATGCA", "GCA"), [1, 6])

    def test_find_motifs_both_strands(self):
        # 'AATT' forward (pos1) and reverse-complement 'AATT' is palindrome
        self.assertEqual(find_motifs("AATTGGAA", "AATT", both_strands=False), [1])
        self.assertGreaterEqual(len(find_motifs("AATTGGAA", "AATT", both_strands=True)), 1)

    def test_restriction_sites(self):
        # EcoRI is GAATTC; "GAATTC" appears at pos 1
        sites = restriction_sites("GAATTCGAATTC", enzymes={"EcoRI": ("GAATTC", 1)})
        self.assertEqual(sites["EcoRI"], [1, 7])

    def test_gc_windows(self):
        w = gc_windows("GGGGTTTT", window=4, step=4)
        self.assertEqual(len(w), len("GGGGTTTT") // 4)
        self.assertEqual(w[0][1], 100.0)
        self.assertEqual(w[1][1], 0.0)

    def test_re_basic_available(self):
        from cfutils.analysis import RE_BASIC
        self.assertIn("EcoRI", RE_BASIC)


if __name__ == "__main__":
    unittest.main()
