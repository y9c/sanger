"""Build script for sanger.

The Cython Smith-Waterman accelerator (``sanger._swalign``) is optional: it is
compiled only when a working C compiler is available, otherwise the package
builds as pure Python (and uses the NumPy Smith-Waterman fallback at runtime).
This keeps ``pip install sanger`` working even without a compiler.
"""

import os
import subprocess
import sys
import tempfile

from setuptools import setup


def _can_compile() -> bool:
    """Return True if a C compiler is available and works."""
    cc = os.environ.get("CC", "cc")
    if not cc or os.environ.get("SANGER_SKIP_EXT", "").lower() in ("1", "true"):
        return False
    code = "int main(void){return 0;}\n"
    try:
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "x.c")
            with open(src, "w") as fh:
                fh.write(code)
            subprocess.run(
                [cc, "-c", src, "-o", os.path.join(d, "x.o")],
                check=True,
                capture_output=True,
                timeout=60,
            )
        return True
    except Exception:
        return False


def main():
    ext_modules = []
    if _can_compile() and sys.platform != "win32":
        from setuptools import Extension

        ext_modules = [
            Extension(
                "sanger._swalign",
                sources=["sanger/_swalign.c"],
                extra_compile_args=["-O3"],
            )
        ]
    setup(ext_modules=ext_modules)


if __name__ == "__main__":
    main()
