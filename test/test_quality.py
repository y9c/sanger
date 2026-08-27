#!/usr/bin/env python3
"""Unit tests for cfutils quality module."""

import unittest

from cfutils.parser import parse_abi
from cfutils.align import SitePair
from cfutils.quality import QualityFilter, site_qualities, passed_filter


class TestQuality(unittest.TestCase):
    """Test quality filtering and site-quality computation."""

    @classmethod
    def setUpClass(cls):
        cls.query = parse_abi("./data/B5-M13R_B07.ab1")

    def test_site_qualities(self):
        site_q, local_q = site_qualities(self.query, 10, flank_base_num=2)
        qual = self.query.letter_annotations["phred_quality"]
        self.assertEqual(site_q, qual[9])
        self.assertIsInstance(local_q, int)

    def test_default_filter(self):
        qf = QualityFilter()
        good = SitePair(ref_pos=1, ref_base="A", cf_pos=2, cf_base="T",
                        qual_site=30, qual_local=30)
        bad = SitePair(ref_pos=1, ref_base="A", cf_pos=2, cf_base="T",
                       qual_site=5, qual_local=5)
        self.assertTrue(qf.passed(good))
        self.assertFalse(qf.passed(bad))

    def test_filter_list(self):
        qf = QualityFilter(min_base_qual=20, min_local_qual=20)
        sites = [
            SitePair(1, "A", 2, "T", 5, 5),
            SitePair(2, "G", 5, "C", 40, 40),
        ]
        kept = qf.filter(sites)
        self.assertEqual(len(kept), 1)
        self.assertEqual(kept[0].cf_pos, 5)

    def test_passed_filter_helper(self):
        self.assertTrue(passed_filter(SitePair(1, "A", 2, "T", 30, 30)))
        self.assertFalse(passed_filter(SitePair(1, "A", 2, "T", 5, 5)))


if __name__ == "__main__":
    unittest.main()
