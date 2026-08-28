#!/usr/bin/env python3
"""Unit tests for sanger transform (trim / reverse-complement) module."""

import unittest

import numpy as np

from sanger.parser import parse_abi
from sanger.transform import reverse_complement_record, strip_primers, trim, trim_ends
from sanger.utils import normalize_ambiguity


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
        self.assertEqual(
            len(trimmed.annotations["channel 1"]), len(trimmed.annotations["trace_x"])
        )
        self.assertEqual(len(trimmed.annotations["peak positions"]), len(trimmed.seq))
        self.assertEqual(
            len(trimmed.letter_annotations["phred_quality"]), len(trimmed.seq)
        )

    def test_trim_keeps_middle(self):
        trimmed = trim(self.query)
        # the original read is high quality in the middle, so trimming should
        # keep a substantial portion of it
        self.assertGreater(len(trimmed.seq), 0.5 * len(self.query.seq))

    def test_reverse_complement_record(self):
        rc = reverse_complement_record(self.query)
        self.assertEqual(
            rc.seq, self.query.seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
        )
        self.assertEqual(len(rc.seq), len(self.query.seq))
        self.assertEqual(
            len(rc.letter_annotations["phred_quality"]), len(self.query.seq)
        )
        # trace arrays preserved
        self.assertEqual(
            len(rc.annotations["channel 1"]), len(self.query.annotations["channel 1"])
        )
        # axis alignment: rc.seq[0] is the complement of the ORIGINAL LAST base,
        # so it must sit at the LEFT (smallest peak x) of the flipped axis
        rp = rc.annotations["peak positions"]
        rt = rc.annotations["trace_x"]
        self.assertTrue(
            all(rp[i] <= rp[i + 1] for i in range(len(rp) - 1)),
            "RC peak axis must be increasing (mirrored then reversed)",
        )
        self.assertEqual(len(rt), len(rc.annotations["channel 1"]))
        # rc.seq[0] should sit where the original last peak maps to after flip:
        # xmax - orig_peaks[-1]
        xmax = max(self.query.annotations["trace_x"])
        self.assertAlmostEqual(
            rp[0], xmax - self.query.annotations["peak positions"][-1], delta=1e-6
        )

    def test_double_rc_is_identity_seq(self):
        rc2 = reverse_complement_record(reverse_complement_record(self.query))
        self.assertEqual(rc2.seq, self.query.seq)
        # double-RC must also restore the original peak axis (mirror+reverse twice)
        np.testing.assert_allclose(
            rc2.annotations["peak positions"],
            self.query.annotations["peak positions"],
            rtol=1e-6,
            atol=1e-6,
        )

    def test_normalize_ambiguity(self):
        self.assertEqual(normalize_ambiguity("ACGT"), "ACGT")
        self.assertEqual(normalize_ambiguity("ACGN"), "ACGN")
        self.assertEqual(normalize_ambiguity("ARYN"), "ANNN")
        self.assertEqual(normalize_ambiguity("acgu"), "ACGT")

    def test_trim_ends(self):
        # build a tiny record with low-quality 5' and 3' ends
        from sanger.parser import SeqRecord

        rec = SeqRecord("ACGTACGTAC", name="t")
        rec.annotations["channel 1"] = list(range(100))
        rec.annotations["channel 2"] = list(range(100))[::-1]
        rec.annotations["channel 3"] = list(range(100))
        rec.annotations["channel 4"] = list(range(100))[::-1]
        rec.annotations["peak positions"] = [float(i) for i in range(10)]
        rec.annotations["trace_x"] = [float(i) for i in range(100)]
        rec.letter_annotations["phred_quality"] = [5] * 3 + [40] * 4 + [5] * 3
        out = trim_ends(rec, min_qual=20)
        self.assertEqual(out.seq, "TACG")

    def test_strip_primers(self):
        from sanger.parser import SeqRecord

        rec = SeqRecord("AAAAATGCGTACGTAAA", name="t")
        rec.annotations["channel 1"] = list(range(20))
        rec.annotations["channel 2"] = list(range(20))
        rec.annotations["channel 3"] = list(range(20))
        rec.annotations["channel 4"] = list(range(20))
        rec.annotations["peak positions"] = [float(i) for i in range(18)]
        rec.annotations["trace_x"] = [float(i) for i in range(20)]
        rec.letter_annotations["phred_quality"] = [40] * 18
        out = strip_primers(rec, forward="AAAAA", reverse="AAA")
        self.assertEqual(out.seq, "TGCGTACGT")


if __name__ == "__main__":
    unittest.main()
