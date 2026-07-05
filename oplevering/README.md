# Oplevering — Virus simulation on networks (grp-SIS)

**Student:** Jan Ruijgrok (852796035)  
**Course:** IM1312 Research Methods for AI, Open University of the Netherlands  
**Project:** Testing a spectral epidemic threshold for SIS on networks under heavy-tailed recovery  
**GitHub:** https://github.com/jruijgrok141/grp-sis-virus-simulation  
**Submission bundle assembled:** 5 July 2026

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

The report title page date (`\RepReportDate`) and git revision in the PDF are refreshed automatically when you run `scripts/build_full_report.py` (via `write_report_macros.py`).
