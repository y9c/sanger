#!/usr/bin/env python3
"""Unit tests for sanger dnalink (dna-features-viewer integration) module."""

import unittest

try:
    import dna_features_viewer  # noqa: F401

    HAVE_DFV = True
except ImportError:
    HAVE_DFV = False

from sanger.dnalink import to_graphic_features, to_graphic_record
from sanger.features import ChromatogramFeature
from sanger.parser import parse_abi


@unittest.skipUnless(HAVE_DFV, "dna_features_viewer not installed")
class TestDnaLink(unittest.TestCase):
    """Test the DNA-features-viewer bridge."""

    @classmethod
    def setUpClass(cls):
        cls.query = parse_abi("./data/B5-M13R_B07.ab1")

    def test_to_graphic_features_coords(self):
        feats = [ChromatogramFeature(start=20, end=30, strand=+1, label="a")]
        gfs = to_graphic_features(feats, len(self.query))
        self.assertEqual(len(gfs), 1)
        self.assertEqual(gfs[0].start, 19)
        self.assertEqual(gfs[0].end, 30)
        self.assertEqual(gfs[0].label, "a")

    def test_to_graphic_record(self):
        rec = to_graphic_record(self.query, [])
        self.assertEqual(rec.sequence_length, len(self.query.seq))
        self.assertEqual(len(rec.features), 0)

    def test_plot_combined_runs(self):
        import matplotlib.pyplot as plt

        from sanger.dnalink import plot_combined

        feats = [ChromatogramFeature(start=10, end=15, strand=+1, label="p")]
        fig, (ax_feat, ax_chrom) = plot_combined(
            self.query, features=feats, region=(5, 25)
        )
        self.assertIsNotNone(fig)
        plt.close(fig)


if __name__ == "__main__":
    unittest.main()
