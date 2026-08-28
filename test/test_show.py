#!/usr/bin/env python3
"""Unit tests for sanger show module."""

import sysconfig
import unittest

import matplotlib.pyplot as plt

from sanger.align import SitePair
from sanger.parser import parse_abi, parse_fasta
from sanger.show import annotate_mutation, highlight_base, plot_chromatograph

# matplotlib's deepcopy of tick properties recurses infinitely on
# free-threaded (Py_GIL_DISABLED) builds, so plotting tests cannot run there.
FREE_THREADED = bool(sysconfig.get_config_var("Py_GIL_DISABLED"))


class TestShowFunc(unittest.TestCase):
    """Test visualization functions in sanger.show."""

    def setUp(self) -> None:
        """Set up test data and figure for plotting tests."""
        self.query_record = parse_abi("./data/B5-M13R_B07.ab1")
        self.subject_record = parse_fasta("./data/ref.fa")
        self.fig, self.ax = plt.subplots(1, 1, figsize=(15, 6))

    @unittest.skipIf(
        FREE_THREADED, "matplotlib deepcopy recurses on free-threaded builds"
    )
    def test_plot_chromatograph(self) -> None:
        """Test plot_chromatograph function runs without error."""
        plot_chromatograph(self.query_record, region=(10, 30), ax=self.ax)
        self.assertTrue(True)

    @unittest.skipIf(
        FREE_THREADED, "matplotlib deepcopy recurses on free-threaded builds"
    )
    def test_highlight_base(self) -> None:
        """Test highlight_base overlays highlight on chromatograph."""
        plot_chromatograph(self.query_record, region=(10, 20), ax=self.ax)
        highlight_base(14, self.query_record, self.ax)
        self.assertTrue(True)

    @unittest.skipIf(
        FREE_THREADED, "matplotlib deepcopy recurses on free-threaded builds"
    )
    def test_region_selects_1based_bases(self) -> None:
        """region=(start, end) is 1-based inclusive: render bases start..end.

        Regression: the region was applied as 0-based indices, shifting the
        window one base to the right (bases 11..21 instead of 10..20).
        """
        ax = plot_chromatograph(self.query_record, region=(10, 20), ax=self.ax)
        labels = [t.get_text() for t in ax.get_xticklabels() if t.get_text()]
        self.assertEqual(labels, [str(i) for i in range(10, 21)])

    def test_annotate_mutation(self) -> None:
        """Test annotate_mutation overlays mutation annotation."""
        mutation = SitePair(ref_pos=10, ref_base="A", cf_pos=14, cf_base="T")
        annotate_mutation(mutation, self.query_record, self.ax)
        self.assertTrue(True)

    def tearDown(self) -> None:
        """Close the matplotlib figure after each test."""
        plt.close(self.fig)


if __name__ == "__main__":
    unittest.main()
