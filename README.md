# Virus simulation on networks (grp-SIS)

Master’s research project (**RMAI**, Open University **IM1312**): an agent-based **susceptible–infected–susceptible (SIS)** epidemic model on explicit contact networks, implemented in **NetLogo 7**, to study how well **spectral / mean-field epidemic thresholds** match simulation outcomes when infected individuals have **non-exponential (including heavy-tailed) recovery times**, while keeping the **mean** infectious period fixed.

The work builds on the grp-SIS perspective (general recovery-time distributions in SIS dynamics), following ideas from Tang et al. and related network epidemiology literature. See `proposal/research_proposal.pdf` (or `research_proposal.tex`) for the formal problem statement, research questions (RQ1–RQ3), and methods narrative.

## Repository layout

| Path | Contents |
|------|----------|
| `netlogo/virus_simulation.nlogox` | Main NetLogo 7 model: grp-SIS dynamics, multiple network generators (e.g. Erdős–Rényi, lattice, ring), immunization options, BehaviorSpace hooks for batch experiments and edge export. |
| `netlogo/behaviorspace_experiments.xml` | BehaviorSpace experiment definitions (reference copy). Comments in that file note NetLogo 7 XML expectations and how experiments relate to headless workflows. |
| `proposal/` | Research proposal (LaTeX/PDF) and supporting figure(s). |
| `report/` | Written results: `report.tex` / `report.pdf`, figures, `Further_Exploration_Report.md` (extended agenda and workspace notes), and related artefacts. |
| `oplevering/` | **Submission bundle** for IM1312: self-contained copy of the final report, proposal, NetLogo model, analysis scripts, simulation outputs, and references. Start here if you want the deliverable as one folder. |
| `references/` | Key papers (PDFs) and a short summary (`Tang_Yao_Xie_Feng_grp-SIS_samenvatting.md`). |

## Running the simulation

1. Install **[NetLogo 7](https://ccl.northwestern.edu/netlogo/)** (the model declares version **7.0.3**).
2. Open `netlogo/virus_simulation.nlogox` in NetLogo.
3. Use the Interface tab to set population size, topology, infection probability, recovery law (exponential vs heavy-tailed variants), limits, and seeds; run **setup**, then **go**.
4. For systematic sweeps, use **BehaviorSpace** (Tools → BehaviorSpace). Experiment XML can be maintained in `behaviorspace_experiments.xml` and merged into the model file as needed for your NetLogo workflow (see comments at the top of that XML file).

Outputs such as exported edge lists are intended to support comparisons against spectral quantities (e.g. largest adjacency eigenvalue) as described in the proposal and in `report/Further_Exploration_Report.md`.

## Documentation to read first

- **`oplevering/report/report.pdf`** — Final report in the submission bundle (same content as `report/report.pdf` when regenerated).
- **`proposal/research_proposal.pdf`** — Abstract, background, RQ1–RQ3, and simulation setup (e.g. fixed ER-style specification, mean recovery time).
- **`report/report.pdf`** — Results, notation, and discussion aligned with the simulation (development copy under `report/`).
- **`report/Further_Exploration_Report.md`** — Consolidated framing, workspace map, and follow-up ideas (some paths there refer to optional tooling or outputs that may live outside this clone).

## Submission (IM1312)

The **deliverable folder** is `oplevering/` (local / OU portal; **not** on GitHub — contains student number). The public `report/report.pdf` omits the student number. Regenerate the report PDF and submission bundle from the repo root:

```text
py -3 scripts/build_full_report.py
py -3 scripts/build_oplevering.py
```

`build_oplevering.py` runs `build_full_report.py` first, so a single `py -3 scripts/build_oplevering.py` is enough before upload.

The report **title-page date** (`\RepReportDate`, e.g. `5 July 2026`) and **git revision** in the data-availability section are written by `scripts/write_report_macros.py` at build time (not edited by hand).

**GitHub:** [github.com/jruijgrok141/grp-sis-virus-simulation](https://github.com/jruijgrok141/grp-sis-virus-simulation)

## License / citation

This repository supports an academic research project. If you reuse the NetLogo model or figures, cite the original thesis/report and the underlying papers listed in `references/` as appropriate.
