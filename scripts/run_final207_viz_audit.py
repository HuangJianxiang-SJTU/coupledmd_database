#!/usr/bin/env python3
"""Run the exact final-207 visualization audit through the legacy implementation."""
from __future__ import annotations

import runpy
from pathlib import Path


LEGACY_IMPLEMENTATION = Path(__file__).with_name("run_final208_viz_audit.py")


if __name__ == "__main__":
    runpy.run_path(str(LEGACY_IMPLEMENTATION), run_name="__main__")
