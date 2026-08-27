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


if __name__ == "__main__":
    unittest.main()
