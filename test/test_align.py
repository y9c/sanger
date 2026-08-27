#!/usr/bin/env python3
"""Unit tests for sanger alignment functions."""

import unittest

from sanger.align import align_chromatograph, call_mutations
from sanger.parser import parse_abi, parse_fasta


class TestAlignFunc(unittest.TestCase):
    """Test alignment and mutation calling in sanger.align."""

    def setUp(self) -> None:
        """Load test data for alignment tests."""
        self.query_record = parse_abi("./data/B5-M13R_B07.ab1")
        self.subject_record = parse_fasta("./data/ref.fa")

    def test_align_chromatograph(self) -> None:
        """Test align_chromatograph returns a list of site pairs."""
        sitepairs = align_chromatograph(self.query_record, self.subject_record)
        self.assertIsInstance(sitepairs, list)
        self.assertGreater(len(sitepairs), 0, "No site pairs found.")

    def test_call_mutations(self) -> None:
        """Test call_mutations returns a list of mutations."""
        mutations = call_mutations(self.query_record, self.subject_record)
        self.assertIsInstance(mutations, list)
        # Optionally check mutation structure if known


if __name__ == "__main__":
    unittest.main()
