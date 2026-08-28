#!/usr/bin/env python3
"""Unit tests for sanger parser functions."""

import unittest

from sanger.parser import parse_abi, parse_fasta


class TestParserFunc(unittest.TestCase):
    """Test parsing functions in sanger.parser."""

    def test_parse_abi(self) -> None:
        """Test parse_abi returns a SeqRecord with expected attributes."""
        record = parse_abi("./data/B5-M13R_B07.ab1")
        self.assertIsNotNone(record)
        self.assertTrue(hasattr(record, "seq"), "SeqRecord missing 'seq' attribute.")

    def test_parse_fasta(self) -> None:
        """Test parse_fasta returns a SeqRecord with expected attributes."""
        record = parse_fasta("./data/ref.fa")
        self.assertIsNotNone(record)
        self.assertTrue(hasattr(record, "seq"), "SeqRecord missing 'seq' attribute.")

    def test_rescaled_peaks_start_at_zero_and_align_with_trace(self):
        """The rescaled peak axis must start at 0, in step with trace_x.

        Regression: peaks were rescaled by ``p/step`` without subtracting
        ``peaks[0]``, so every base label was shifted right of its trace peak.
        """
        record = parse_abi("./data/B5-M13R_B07.ab1", rescale=True)
        peaks = record.annotations["peak positions"]
        trace_x = record.annotations["trace_x"]
        self.assertAlmostEqual(peaks[0], 0.0, delta=1e-6)
        # trace_x is trimmed to the peak span and starts at 0
        self.assertAlmostEqual(trace_x[0], 0.0, delta=1e-6)


if __name__ == "__main__":
    unittest.main()
