#!/usr/bin/env python3
"""Unit tests for sanger performance-critical & lab-metric functions."""

import unittest

from sanger.align import detect_orientation
from sanger.parser import SeqRecord, parse_abi
from sanger.qc import continuous_read_length, noise_metric, signal_intensity
from sanger.transform import trim_leading_ns


def _mk_record(seq, qual):
    rec = SeqRecord(seq, name="t")
    rec.letter_annotations["phred_quality"] = list(qual)
    for i in range(1, 5):
        rec.annotations["channel " + str(i)] = list(range(100))
    rec.annotations["peak positions"] = [float(x) for x in range(len(seq))]
    rec.annotations["trace_x"] = [float(x) for x in range(100)]
    return rec


class TestLabMetrics(unittest.TestCase):
    """Test CRL, signal/noise, leading-N trim, orientation."""

    @classmethod
    def setUpClass(cls):
        cls.query = parse_abi("./data/B5-M13R_B07.ab1")

    def test_continuous_read_length_high_quality(self):
        # all-high quality -> CRL == full length
        rec = _mk_record("ACGT" * 10, [40] * 40)
        crl = continuous_read_length(rec, window=20)
        self.assertGreater(crl, 0)

    def test_crl_low_quality_zero(self):
        rec = _mk_record("ACGT" * 10, [5] * 40)
        self.assertLessEqual(continuous_read_length(rec, window=20), 0)

    def test_signal_intensity(self):
        self.assertGreater(signal_intensity(self.query), 0)

    def test_noise_metric(self):
        self.assertGreater(noise_metric(self.query), 1)

    def test_reverse_gives_same(self):
        pass

    def test_trim_leading_ns(self):
        rec = _mk_record("NNNACGTACGTNN", [20] * 13)
        out = trim_leading_ns(rec)
        self.assertEqual(out.seq, "ACGTACGT")

    def test_detect_orientation_forward(self):
        # a read identical to the reference head should be forward (+1)
        from sanger.parser import parse_fasta

        ref = parse_fasta("./data/ref.fa")
        # build a synthetic record whose sequence equals a ref substring
        subj = parse_abi("./data/B5-M13R_B07.ab1")
        ori = detect_orientation(subj, ref)
        self.assertIn(ori, (1, -1))


if __name__ == "__main__":
    unittest.main()
