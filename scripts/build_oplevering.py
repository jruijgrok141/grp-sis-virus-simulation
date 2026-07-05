#!/usr/bin/env python3
"""
Assemble the IM1312 research-project delivery folder ``oplevering/`` at the repo root.

Copies the final report, NetLogo model, analysis pipeline, simulation outputs, proposal,
and references. Excludes course prep (leerdoelen, tentamenvoorbereiding), the early
benchmark prototype, and internal review notes.

Usage (from project root):
  py -3 scripts/build_oplevering.py
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "oplevering"

REPORT_FIGURES = (
    "er_network_example.png",
    "fig_grp_sis_process_schematics.png",
    "fig_er_survival_vs_beta.png",
    "fig_er_late_prevalence_vs_beta.png",
    "fig_er_threshold_ratio_by_seed.png",
    "fig_er_extinct_median_vs_beta_by_seed.png",
    "fig_er_extinct_tick_violin_by_beta.png",
    "fig_er_prevalence_trajectory_bands.png",
)

def _submission_date() -> str:
    now = datetime.now()
    return f"{now.day} {now.strftime('%B %Y')}"


def _readme_text() -> str:
    submission_date = _submission_date()
    return f"""# Oplevering — Virus simulation on networks (grp-SIS)

**Student:** Jan Ruijgrok (852796035)  
**Course:** IM1312 Research Methods for AI, Open University of the Netherlands  
**Project:** Testing a spectral epidemic threshold for SIS on networks under heavy-tailed recovery  
**GitHub:** https://github.com/jruijgrok141/grp-sis-virus-simulation  
**Submission bundle assembled:** {submission_date}

This folder contains the materials submitted for the research project: the written report,
the NetLogo agent-based model, Python analysis scripts, simulation outputs, the original
research proposal, and key references.

## Contents

| Path | Description |
|------|-------------|
| `report/report.pdf` | **Final report** (main deliverable). |
| `report/report.tex` | LaTeX source; figures listed below; macros in `generated_quantities.tex`. |
| `proposal/` | Research proposal (`research_proposal.pdf` / `.tex`). |
| `netlogo/` | NetLogo 7 model (`virus_simulation.nlogox`) and BehaviorSpace experiment XML. |
| `scripts/` | Pipeline: headless BehaviorSpace runs, threshold aggregation, figure export, PDF build. |
| `output/` | Edge lists, BehaviorSpace CSV exports (`raw/`), and aggregated threshold tables. |
| `references/` | Core papers (PDF) and a short Tang et al. summary (Markdown). |
| `requirements.txt` | Python dependencies for the analysis scripts. |

## Software requirements

- **NetLogo 7.0.3** (or compatible 7.x) with `NetLogo_Console.exe` for headless BehaviorSpace runs.
- **Python 3.10+** with packages from `requirements.txt`.
- **pdfLaTeX** (e.g. MiKTeX) only if you want to rebuild `report/report.pdf` from source.

## Quick start

### Read the results

Open `report/report.pdf`.

### Reproduce analysis from existing simulation output

From this folder (`oplevering/`), with Python dependencies installed:

```text
py -3 scripts/build_full_report.py
```

This refreshes `report/generated_quantities.tex`, exports report figures, and compiles `report/report.pdf`. It does **not** re-run NetLogo.

### Re-run simulations (optional, time-consuming)

Set `NETLOGO_HOME` to your NetLogo installation directory, then from this folder:

```text
py -3 scripts/run_pipeline.py --all
py -3 scripts/build_full_report.py
```

`--all` runs experiments `01` (network export), `05`–`07` (threshold sweeps, three ER seeds), and `08` (prevalence trajectories). Lattice/ring experiments `02`–`04` are optional (`--with-lattice-ring`).

## Main experiments (BehaviorSpace)

| ID | Name | Role |
|----|------|------|
| `01` | export-networks | Export ER edge lists (`output/edges/edges_ER_<seed>.csv`). |
| `05` | baseline-empirical-threshold-ER | Exponential recovery; β sweep. |
| `06` | baseline-ER-power-law-Tang | Tang power-law recovery; same sweep. |
| `07` | baseline-ER-lognormal-Tang | Tang lognormal recovery; same sweep. |
| `08` | dynamics-trajectory-ER | Per-tick prevalence vs. mean-field (seed 10001). |

## Regenerating this delivery folder

From the full development repository (parent of `oplevering/`):

```text
py -3 scripts/build_oplevering.py
```

## Citation and license

Academic research project (Open University). When reusing the NetLogo model or report figures,
cite this report and the underlying literature in `references/`.

The report title page date (`\\RepReportDate`) and git revision in the PDF are refreshed automatically when you run `scripts/build_full_report.py` (via `write_report_macros.py`).
"""


def _copy_file(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def _copy_tree(src: Path, dst: Path, *, ignore=None) -> None:
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst, ignore=ignore)


def _ignore_pycache(_dir: str, names: list[str]) -> set[str]:
    skip = {n for n in names if n == "__pycache__" or n.endswith(".pyc")}
    return skip


def main() -> int:
    print("Refreshing report PDF (macros, figures, date, git hash)...")
    rebuild = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "build_full_report.py")],
        cwd=str(ROOT),
    )
    if rebuild.returncode != 0:
        return rebuild.returncode
    if not (ROOT / "report" / "report.pdf").is_file():
        print("Missing report/report.pdf after build_full_report.", file=sys.stderr)
        return 1

    if DEST.exists():
        shutil.rmtree(DEST)
    DEST.mkdir(parents=True)

    # Report
    report_src = ROOT / "report"
    report_dst = DEST / "report"
    report_dst.mkdir()
    for name in ("report.pdf", "report.tex", "generated_quantities.tex"):
        _copy_file(report_src / name, report_dst / name)
    for fig in REPORT_FIGURES:
        src = report_src / fig
        if not src.is_file():
            print(f"Warning: missing report figure {fig}", file=sys.stderr)
        else:
            _copy_file(src, report_dst / fig)

    # Proposal (source + PDF + figure; skip LaTeX aux)
    proposal_src = ROOT / "proposal"
    proposal_dst = DEST / "proposal"
    proposal_dst.mkdir()
    for name in ("research_proposal.tex", "research_proposal.pdf", "er_network_example.png"):
        src = proposal_src / name
        if src.is_file():
            _copy_file(src, proposal_dst / name)
        elif name == "research_proposal.pdf":
            alt = ROOT / "benchmark" / "research_proposal.pdf"
            if alt.is_file():
                _copy_file(alt, proposal_dst / name)
            else:
                print("Warning: no research_proposal.pdf found", file=sys.stderr)

    # NetLogo, scripts, output, references
    _copy_tree(ROOT / "netlogo", DEST / "netlogo")
    _copy_tree(ROOT / "scripts", DEST / "scripts", ignore=_ignore_pycache)
    _copy_tree(ROOT / "output", DEST / "output")
    _copy_tree(ROOT / "references", DEST / "references")
    _copy_file(ROOT / "requirements.txt", DEST / "requirements.txt")

    (DEST / "README.md").write_text(_readme_text(), encoding="utf-8")

    n_files = sum(1 for _ in DEST.rglob("*") if _.is_file())
    size_mb = sum(f.stat().st_size for f in DEST.rglob("*") if f.is_file()) / (1024 * 1024)
    print(f"Wrote {DEST} ({n_files} files, {size_mb:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
