#!/usr/bin/env python3
"""Basic environment and sanity tests for sanger."""

import sys
import unittest


class TestBasicFunc(unittest.TestCase):
    """Test basic environment requirements for sanger."""

    def test_python_version(self) -> None:
        """Ensure Python version is >= 3.6."""
        self.assertGreaterEqual(sys.version_info, (3, 6), "Python 3.6+ is required.")


if __name__ == "__main__":
    unittest.main()
