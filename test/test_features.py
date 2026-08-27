#!/usr/bin/env python3
"""Unit tests for sanger features (annotation overlay) module."""

import unittest

import matplotlib.pyplot as plt

from sanger.align import SitePair
from sanger.features import (
    ChromatogramFeature,
    add_feature,
    iter_features,
    peak_to_x,
    plot_features,
)
from sanger.parser import parse_abi
from sanger.show import plot_chromatograph


class TestFeatures(unittest.TestCase):
    """Test the feature annotation / overlay API."""

    @classmethod
    def setUpClass(cls):
        cls.query = parse_abi("./data/B5-M13R_B07.ab1")

    def test_add_and_iter(self):
        feat = ChromatogramFeature(start=10, end=20, label="primer")
        add_feature(self.query, feat)
        self.assertIn(feat, list(iter_features(self.query)))

    def test_peak_to_x_maps_positions(self):
        xs = peak_to_x(self.query, [1, 10, 100])
        self.assertEqual(len(xs), 3)
        peaks = self.query.annotations["peak positions"]
        self.assertAlmostEqual(xs[0], peaks[0])
        self.assertAlmostEqual(xs[1], peaks[9])

    def test_strand_validation(self):
        with self.assertRaises(ValueError):
            ChromatogramFeature(start=5, end=3)
        with self.assertRaises(ValueError):
            ChromatogramFeature(start=5, end=10, strand=7)

    def test_from_sitepair(self):
        site = SitePair(ref_pos=3, ref_base="A", cf_pos=14, cf_base="T")
        feat = ChromatogramFeature.from_sitepair(site)
        self.assertEqual(feat.start, 14)
        self.assertEqual(feat.end, 14)
        self.assertIn("A3T", feat.label)

    def test_plot_features_runs(self):
        fig, ax = plt.subplots(1, 1, figsize=(15, 6))
        feats = [
            ChromatogramFeature(
                start=10, end=20, strand=+1, label="fwd", color="#ff8888"
            ),
            ChromatogramFeature(
                start=25, end=30, strand=-1, label="rev", color="#8888ff"
            ),
            ChromatogramFeature(start=40, end=40, strand=0, label="snv"),
        ]
        plot_features(self.query, ax, features=feats)
        plt.close(fig)

    def test_plot_features_clips_out_of_range(self):
        """Features outside the visible x-range must not expand the axes."""
        fig, ax = plt.subplots(1, 1, figsize=(15, 6))
        plot_chromatograph(self.query, region=(10, 20), ax=ax)
        xlim_before = ax.get_xlim()
        # a feature far outside the plotted region
        far = ChromatogramFeature(start=800, end=810, label="faraway")
        plot_features(self.query, ax, features=[far])
        # drawing an out-of-range feature must not widen the x-axis
        self.assertAlmostEqual(ax.get_xlim()[0], xlim_before[0], delta=2)
        self.assertAlmostEqual(ax.get_xlim()[1], xlim_before[1], delta=2)
        plt.close(fig)


if __name__ == "__main__":
    unittest.main()
