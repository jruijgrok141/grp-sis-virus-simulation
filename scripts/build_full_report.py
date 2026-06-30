#!/usr/bin/env python3
"""
Refresh quantitative macros, export all report figures, compile report.pdf (pdflatex x2).

Does not run NetLogo; use after `py -3 scripts/run_pipeline.py --all`.

Usage (from project root):
  py -3 scripts/build_full_report.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "report"


def main() -> int:
    subprocess.run([sys.executable, str(ROOT / "scripts" / "write_report_macros.py")], check=True, cwd=str(ROOT))
    subprocess.run([sys.executable, str(ROOT / "scripts" / "export_report_figures.py")], check=True, cwd=str(ROOT))
    for _ in range(2):
        r = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "report.tex"],
            cwd=str(REPORT),
            check=False,
        )
        if r.returncode != 0:
            print("pdflatex exit", r.returncode)
            return r.returncode
    print("Wrote", REPORT / "report.pdf")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
