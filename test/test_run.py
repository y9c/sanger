#!/usr/bin/env python3
"""Unit tests for sanger run module."""

import os
import tempfile
import unittest

from sanger.run import report_mutation


class TestRunFunc(unittest.TestCase):
    """Test report_mutation in sanger.run."""

    def test_report_mutation(self) -> None:
        """Test report_mutation creates output files without error."""
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                report_mutation(
                    query_ab1_file="./data/B5-M13R_B07.ab1",
                    subject_fasta_file="./data/ref.fa",
                    output_dir=tmpdir,
                    file_basename="test",
                    report_all_sites=True,
                    report_mut_plot=False,
                )
                tsv_path = os.path.join(tmpdir, "test.tsv")
                if not os.path.exists(tsv_path):
                    print(f"Temp dir contents: {os.listdir(tmpdir)}")
                self.assertTrue(os.path.exists(tsv_path), f"Missing TSV: {tsv_path}")
            except Exception as e:
                self.fail(f"report_mutation raised an exception: {e}")


if __name__ == "__main__":
    unittest.main()
