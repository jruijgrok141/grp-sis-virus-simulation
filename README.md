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
| `references/` | Key papers (PDFs) and a short summary (`Tang_Yao_Xie_Feng_grp-SIS_samenvatting.md`). |

## Running the simulation

1. Install **[NetLogo 7](https://ccl.northwestern.edu/netlogo/)** (the model declares version **7.0.3**).
2. Open `netlogo/virus_simulation.nlogox` in NetLogo.
3. Use the Interface tab to set population size, topology, infection probability, recovery law (exponential vs heavy-tailed variants), limits, and seeds; run **setup**, then **go**.
4. For systematic sweeps, use **BehaviorSpace** (Tools → BehaviorSpace). Experiment XML can be maintained in `behaviorspace_experiments.xml` and merged into the model file as needed for your NetLogo workflow (see comments at the top of that XML file).

Outputs such as exported edge lists are intended to support comparisons against spectral quantities (e.g. largest adjacency eigenvalue) as described in the proposal and in `report/Further_Exploration_Report.md`.

## Documentation to read first

- **`proposal/research_proposal.pdf`** — Abstract, background, RQ1–RQ3, and simulation setup (e.g. fixed ER-style specification, mean recovery time).
- **`report/report.pdf`** — Results, notation, and discussion aligned with the simulation.
- **`report/Further_Exploration_Report.md`** — Consolidated framing, workspace map, and follow-up ideas (some paths there refer to optional tooling or outputs that may live outside this clone).

## License / citation

This repository supports an academic research project. If you reuse the NetLogo model or figures, cite the original thesis/report and the underlying papers listed in `references/` as appropriate.
