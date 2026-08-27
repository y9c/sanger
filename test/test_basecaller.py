#!/usr/bin/env python3
"""Unit tests for cfutils basecaller module."""

import unittest

from cfutils.parser import parse_abi
from cfutils.basecaller import call_bases, detect_peaks, basecaller_score

class TestBasecaller(unittest.TestCase):
    """Test the peak-driven base caller (raw traces mode)."""

    @classmethod
    def setUpClass(cls):
        cls.query = parse_abi("./data/B5-M13R_B07.ab1", rescale=False)

    def test_call_bases_length_matches(self):
        res = call_bases(self.query)
        self.assertEqual(res.n_calls, len(self.query.seq))
        self.assertEqual(len(res.sequence), len(self.query.seq))
        self.assertEqual(len(res.qualities), len(self.query.seq))

    def test_accuracy_against_vendor_call(self):
        res = call_bases(self.query)
        score = basecaller_score(res, self.query.seq)
        # a sane re-call should agree with the vendor call most of the time
        self.assertGreater(score["accuracy"], 0.9)
        self.assertGreater(score["mean_quality"], 20)

    def test_detect_peaks(self):
        peaks = detect_peaks(self.query)
        self.assertGreater(len(peaks), 0)


if __name__ == "__main__":
    unittest.main()
