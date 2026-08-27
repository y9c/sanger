#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Optional, lazy matplotlib loader.

cfutils' parsing, alignment, QC and export work without matplotlib.  Plotting
is an optional capability: importing ``cfutils`` never fails when matplotlib
is absent, but calling a plotting function raises a clear, actionable error.

Use::

    from ._mpl import mpl, plt, Axes, HAVE_MPL, require_matplotlib
    ...
    require_matplotlib()          # raises a helpful ImportError if missing
"""

from __future__ import annotations

try:  # pragma: no cover - exercised when matplotlib is installed
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

    HAVE_MPL = True
except Exception:  # pragma: no cover - no matplotlib installed
    from typing import Any as Axes
    from typing import Any as Figure

    mpl = None
    plt = None
    HAVE_MPL = False


def require_matplotlib():
    """Return the matplotlib root module, or raise a helpful error."""
    if not HAVE_MPL:
        raise ImportError(
            "matplotlib is required for cfutils plotting. "
            "Install it with `pip install cfutils[plot]` (or `pip install matplotlib`)."
        )
    return mpl


__all__ = ["mpl", "plt", "Axes", "Figure", "HAVE_MPL", "require_matplotlib"]
