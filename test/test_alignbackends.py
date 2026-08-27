#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Alignment backend equivalence tests.

Verifies the Cython Smith-Waterman accelerator produces the same alignment as
the pure-Python/NumPy fallback, and that mutation calling works through either.
"""

import unittest

from sanger.align import _HAVE_CY, _cy_swalign, _sw_align, run_align
from sanger.parser import parse_abi, parse_fasta


class TestAlignBackends(unittest.TestCase):
    """Test alignment backends (Cython vs NumPy) agree."""

    @classmethod
    def setUpClass(cls):
        cls.query = parse_abi("./data/B5-M13R_B07.ab1")
        cls.ref = parse_fasta("./data/ref.fa")

    @unittest.skipUnless(_HAVE_CY, "cython accelerator not built")
    def test_cython_matches_numpy(self):
        ref = self.ref.seq
        qry = self.query.seq[100:250]
        cy = _cy_swalign(ref, qry)
        np_ = _sw_align(ref, qry)
        # both should find the matching substring with no gaps in this region
        self.assertEqual(cy[2].replace("-", ""), np_.alignment[0].replace("-", ""))
        self.assertGreater(len(cy[2].replace("-", "")), 100)

    def test_run_align_works(self):
        sites = run_align(self.ref.seq, self.query.seq)
        self.assertGreater(len(sites), 0)

    def test_coordinates_are_one_based(self):
        # query == ref[2:] -> ref_pos must start at 3, cf_pos at 1 (1-based)
        ref = "ACGTACGTACGT"
        query = "GTACGTAC"
        sites = run_align(ref, query)
        self.assertEqual(sites[0].ref_pos, 3)
        self.assertEqual(sites[0].cf_pos, 1)
        self.assertEqual(sites[0].ref_base, "G")
        self.assertEqual(sites[0].cf_base, "G")

    def test_call_mutations_via_backend(self):
        from sanger.align import call_mutations

        sites = call_mutations(self.query, self.ref, report_all_sites=True)
        # the sample should align most of its length
        self.assertGreater(len(sites), 1000)


if __name__ == "__main__":
    unittest.main()
