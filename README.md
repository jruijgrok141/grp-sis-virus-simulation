# Virus simulation on networks (grp-SIS)

**Author:** Jan Ruijgrok  
**Course:** IM1312 Research Methods for AI, Open University of the Netherlands

Agent-based **susceptible–infected–susceptible (SIS)** epidemic model on explicit contact networks, implemented in **NetLogo 7**, to study how well **spectral / mean-field epidemic thresholds** match simulation outcomes when infected individuals have **non-exponential (including heavy-tailed) recovery times**, while keeping the **mean** infectious period fixed.

The work builds on the grp-SIS perspective (general recovery-time distributions in SIS dynamics), following ideas from Tang et al. and related network epidemiology literature. See [`proposal/research_proposal.pdf`](proposal/research_proposal.pdf) for the formal problem statement, research questions (RQ1–RQ3), and methods narrative.

## Contents

| Path | Description |
|------|-------------|
| [`report/report.pdf`](report/report.pdf) | **Final report** (main deliverable). |
| `report/report.tex` | LaTeX source; figures; macros in `generated_quantities.tex`. |
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

Open [`report/report.pdf`](report/report.pdf).

### Reproduce analysis from existing simulation output

From the repository root, with Python dependencies installed:

```text
py -3 scripts/build_full_report.py
```

This refreshes `report/generated_quantities.tex`, exports report figures, and compiles `report/report.pdf`. It does **not** re-run NetLogo.

### Re-run simulations (optional, time-consuming)

Set `NETLOGO_HOME` to your NetLogo installation directory, then from the repository root:

```text
py -3 scripts/run_pipeline.py --all
py -3 scripts/build_full_report.py
```

`--all` runs experiments `01` (network export), `05`–`07` (threshold sweeps, three ER seeds), and `08` (prevalence trajectories). Lattice/ring experiments `02`–`04` are optional (`--with-lattice-ring`).

### Run interactively in NetLogo

1. Install **[NetLogo 7](https://ccl.northwestern.edu/netlogo/)** (the model declares version **7.0.3**).
2. Open `netlogo/virus_simulation.nlogox` in NetLogo.
3. Use the Interface tab to set population size, topology, infection probability, recovery law, limits, and seeds; run **setup**, then **go**.
4. For systematic sweeps, use **BehaviorSpace** (Tools → BehaviorSpace).

## Main experiments (BehaviorSpace)

| ID | Name | Role |
|----|------|------|
| `01` | export-networks | Export ER edge lists (`output/edges/edges_ER_<seed>.csv`). |
| `05` | baseline-empirical-threshold-ER | Exponential recovery; β sweep. |
| `06` | baseline-ER-power-law-Tang | Tang power-law recovery; same sweep. |
| `07` | baseline-ER-lognormal-Tang | Tang lognormal recovery; same sweep. |
| `08` | dynamics-trajectory-ER | Per-tick prevalence vs. mean-field (seed 10001). |

## License and citation

This repository is released under [CC BY 4.0](LICENSE). When reusing the NetLogo model, code, or report figures, cite this report and the underlying literature in `references/`.

Third-party papers in `references/` remain the property of their respective authors and publishers.
