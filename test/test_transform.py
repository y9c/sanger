#!/usr/bin/env python3
"""Unit tests for cfutils transform (trim / reverse-complement) module."""

import unittest
import numpy as np

from cfutils.parser import parse_abi
from cfutils.transform import trim, reverse_complement_record
from cfutils.utils import normalize_ambiguity

class TestTransform(unittest.TestCase):
    """Test whole-record transformations."""

    @classmethod
    def setUpClass(cls):
        # raw/rescaled both should keep axes consistent; use rescaled (plotting)
        cls.query = parse_abi("./data/B5-M13R_B07.ab1")

    def test_trim_preserve_length_consistency(self):
        trimmed = trim(self.query)
        self.assertLessEqual(len(trimmed), len(self.query))
        # length of channels must equal trace_x, and peaks aligned to seq
        self.assertEqual(len(trimmed.annotations["channel 1"]),
                         len(trimmed.annotations["trace_x"]))
        self.assertEqual(len(trimmed.annotations["peak positions"]),
                         len(trimmed.seq))
        self.assertEqual(len(trimmed.letter_annotations["phred_quality"]),
                         len(trimmed.seq))

    def test_trim_keeps_middle(self):
        trimmed = trim(self.query)
        # the original read is high quality in the middle, so trimming should
        # keep a substantial portion of it
        self.assertGreater(len(trimmed.seq), 0.5 * len(self.query.seq))

    def test_reverse_complement_record(self):
        rc = reverse_complement_record(self.query)
        self.assertEqual(rc.seq, self.query.seq.translate(
            str.maketrans("ACGTN", "TGCAN"))[::-1])
        self.assertEqual(len(rc.seq), len(self.query.seq))
        self.assertEqual(len(rc.letter_annotations["phred_quality"]),
                         len(self.query.seq))
        # trace arrays preserved
        self.assertEqual(len(rc.annotations["channel 1"]),
                         len(self.query.annotations["channel 1"]))

    def test_double_rc_is_identity_seq(self):
        rc2 = reverse_complement_record(reverse_complement_record(self.query))
        self.assertEqual(rc2.seq, self.query.seq)

    def test_normalize_ambiguity(self):
        self.assertEqual(normalize_ambiguity("ACGT"), "ACGT")
        self.assertEqual(normalize_ambiguity("ACGN"), "ACGN")
        self.assertEqual(normalize_ambiguity("ARYN"), "ANNN")
        self.assertEqual(normalize_ambiguity("acgu"), "ACGT")


if __name__ == "__main__":
    unittest.main()
