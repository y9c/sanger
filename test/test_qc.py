#!/usr/bin/env python3
"""Unit tests for cfutils qc (per-read QC metrics) module."""

import unittest

from cfutils.parser import parse_abi
from cfutils.qc import read_metrics, trimmed_bounds, summarize

class TestQC(unittest.TestCase):
    """Test the per-read QC metrics."""

    @classmethod
    def setUpClass(cls):
        cls.query = parse_abi("./data/B5-M13R_B07.ab1")

    def test_read_metrics_shape(self):
        m = read_metrics(self.query)
        for key in ("n_bases", "mean_qual", "min_qual", "n_fraction",
                    "trim_start", "trim_end", "trimmed_len"):
            self.assertIn(key, m)
        self.assertEqual(m["n_bases"], len(self.query.seq))
        self.assertGreaterEqual(m["mean_qual"], 0)

    def test_trimmed_bounds_within_read(self):
        s, e = trimmed_bounds(self.query)
        self.assertGreaterEqual(s, 1)
        self.assertLessEqual(e, len(self.query))
        self.assertLessEqual(s, e)

    def test_summarize_batch(self):
        table = summarize([self.query, self.query])
        self.assertEqual(len(table), 2)
        self.assertAlmostEqual(table[0]["gc_percent"], table[1]["gc_percent"])


if __name__ == "__main__":
    unittest.main()
