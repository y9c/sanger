#!/usr/bin/env python3
"""Unit tests for the sanger MCP server."""

import unittest

try:
    from sanger import mcp_server

    HAVE_MCP = mcp_server._HAVE_MCP
except Exception:
    HAVE_MCP = False


@unittest.skipUnless(HAVE_MCP, "mcp not installed")
class TestMCPServer(unittest.TestCase):
    """Test the MCP server tool functions (direct calls)."""

    @classmethod
    def setUpClass(cls):
        cls._names = None

    def test_tools_registered(self):
        import asyncio

        names = asyncio.run(mcp_server.mcp.list_tools())
        n = {t.name for t in names}
        for expected in (
            "read_chromatogram",
            "qc_metrics",
            "call_mutations",
            "re_call_bases",
            "analyze_sequence",
            "trim_read",
            "export_sequence",
            "plot_chromatogram",
        ):
            self.assertIn(expected, n)

    def test_read_chromatogram(self):
        info = mcp_server.read_chromatogram("./data/B5-M13R_B07.ab1")
        self.assertEqual(info["length"], 1141)
        self.assertGreater(info["mean_qual"], 20)

    def test_qc_metrics(self):
        m = mcp_server.qc_metrics("./data/B5-M13R_B07.ab1")
        self.assertIn("crl", m)
        self.assertIn("snr", m)

    def test_call_mutations_returns_variants(self):
        r = mcp_server.call_mutations("./data/B5-M13R_B07.ab1", "./data/ref.fa")
        self.assertGreater(r["n_variants"], 0)
        self.assertIn("ref_pos", r["variants"][0])

    def test_re_call_bases(self):
        out = mcp_server.re_call_bases("./data/B5-M13R_B07.ab1")
        self.assertEqual(out["n_calls"], 1141)
        self.assertGreater(len(out["sequence"]), 1000)

    def test_analyze_sequence(self):
        r = mcp_server.analyze_sequence("./data/B5-M13R_B07.ab1", "translate")
        self.assertTrue(r["result"])


if __name__ == "__main__":
    unittest.main()
