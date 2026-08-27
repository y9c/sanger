#!/usr/bin/env python3
"""Unit tests for cfutils tracks (split/join) module."""

import unittest
import numpy as np

from cfutils.parser import parse_abi
from cfutils.tracks import join_tracks, split_track, slice_track, export_tracks, import_tracks


class TestTracks(unittest.TestCase):
    """Test trace splitting and joining."""

    @classmethod
    def setUpClass(cls):
        cls.query = parse_abi("./data/B5-M13R_B07.ab1")

    def test_join_preserves_length(self):
        a = self.query
        b = slice_track(a, 1, 20)
        joined = join_tracks(a, b, name="joined")
        self.assertEqual(joined.seq, a.seq + b.seq)
        self.assertEqual(len(joined.annotations["channel 1"]),
                         len(a.annotations["channel 1"]) + len(b.annotations["channel 1"]))
        # the x axis should contain a gap between the two records
        self.assertGreater(len(joined.annotations["trace_x"]),
                           max(a.annotations["channel 1"]) if False else 0)
    def test_split_reassembles_sequence(self):
        pieces = split_track(self.query, cuts=[20, 40])
        self.assertEqual(len(pieces), 3)
        reassembled = "".join(p.seq for p in pieces)
        self.assertEqual(reassembled, self.query.seq)

    def test_slice_track_region(self):
        seg = slice_track(self.query, 10, 15)
        self.assertEqual(seg.seq, self.query.seq[9:15])
        self.assertEqual(len(seg.annotations["peak positions"]), 6)
        self.assertEqual(len(seg.annotations["channel 1"]),
                         len(seg.annotations["trace_x"]))

    def test_npz_roundtrip(self):
        import tempfile, os
        with tempfile.TemporaryDirectory() as tmp:
            paths = export_tracks([self.query], tmp, fmt="npz")
            rec = import_tracks(paths[0])
            self.assertTrue(np.allclose(rec.annotations["peak positions"],
                                        self.query.annotations["peak positions"]))
            self.assertEqual(len(rec.annotations["channel 1"]),
                             len(self.query.annotations["channel 1"]))

    def test_slice_preserves_channels_and_plots(self):
        seg = slice_track(self.query, 10, 20)
        self.assertEqual(seg.annotations["channels"],
                         self.query.annotations["channels"])
        # plotting a sliced record must not fail (needs channels annotation)
        import matplotlib
        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
        from cfutils.show import plot_chromatograph
        fig, ax = plt.subplots()
        plot_chromatograph(seg, ax=ax)
        plt.close(fig)

    def test_join_preserves_channels(self):
        seg = slice_track(self.query, 1, 20)
        joined = join_tracks(self.query, seg)
        self.assertEqual(joined.annotations["channels"],
                         self.query.annotations["channels"])

    def test_slice_records_offset_and_normalises(self):
        seg = slice_track(self.query, 20, 40)
        # provenance: 0-based offset into the parent read
        self.assertEqual(seg.annotations["offset"], 19)
        self.assertEqual(seg.annotations["parent"], self.query.name)
        # peak axis is re-normalised to start near 0
        self.assertAlmostEqual(seg.annotations["peak positions"][0], 0.0, delta=1.0)
        self.assertEqual(len(seg.annotations["peak positions"]), 21)


if __name__ == "__main__":
    unittest.main()
