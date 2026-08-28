#!/usr/bin/env python3
"""Unit tests for sanger composite (side-by-side plotting) module."""

import sysconfig
import unittest

import matplotlib.pyplot as plt

from sanger.composite import add_panel, side_by_side
from sanger.features import ChromatogramFeature
from sanger.parser import parse_abi

# matplotlib's deepcopy of tick properties recurses infinitely on
# free-threaded (Py_GIL_DISABLED) builds, so plotting tests cannot run there.
FREE_THREADED = bool(sysconfig.get_config_var("Py_GIL_DISABLED"))


def gc_panel(ax, trace_x, peaks, seq, record, start_pos=None):
    """A minimal user panel function drawing GC content over windows."""
    n = len(seq)
    if n < 4:
        ax.set_yticks([])
        return
    start_pos = start_pos or 1
    xs = [start_pos + i for i in range(n)]
    gc = [
        100.0
        * (
            seq[max(0, i - 1) : i + 2].count("G")
            + seq[max(0, i - 1) : i + 2].count("C")
        )
        / 3
        for i in range(n)
    ]
    ax.fill_between(xs, gc, alpha=0.3)
    ax.set_ylim(0, 100)


class TestComposite(unittest.TestCase):
    """Test side-by-side composite plotting."""

    @classmethod
    def setUpClass(cls):
        cls.query = parse_abi("./data/B5-M13R_B07.ab1")

    @unittest.skipIf(
        FREE_THREADED, "matplotlib deepcopy recurses on free-threaded builds"
    )
    def test_side_by_side_runs(self):
        fig, (ax_chrom, ax_panel) = side_by_side(
            self.query,
            gc_panel,
            region=(10, 40),
            features=[
                ChromatogramFeature(start=20, end=25, label="amplicon", color="#ffcc88")
            ],
        )
        self.assertIsNotNone(fig)
        plt.close(fig)

    @unittest.skipIf(
        FREE_THREADED, "matplotlib deepcopy recurses on free-threaded builds"
    )
    def test_add_panel_runs(self):
        fig, ax = plt.subplots(1, 1, figsize=(16, 5))
        from sanger.show import plot_chromatograph

        plot_chromatograph(self.query, region=(10, 40), ax=ax)
        ax_panel = add_panel(fig, ax, gc_panel, self.query, region=(10, 40))
        self.assertIsNotNone(ax_panel)
        plt.close(fig)


if __name__ == "__main__":
    unittest.main()
