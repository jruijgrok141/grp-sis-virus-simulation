# Further exploration report: spectral SIS thresholds under heavy-tailed recovery

**Scope:** This document consolidates the proposed research question, the artefacts already present in the workspace, interim evidence, and a **prioritised agenda** for extending the study. It is written to support the next phase of the RMAI research project (Open University, IM1312) and to align simulation design, analysis, and reporting.

**Date:** April 2026

---

## 1. Canonical research frame

### 1.1 Core question

**How robust is the mean-field spectral benchmark** \(\tau_c \approx c/\lambda_{\max}\) **with** \(c \approx 1\) **when an agent-based susceptible–infected–susceptible (SIS) model on a fixed contact network uses non-exponential (heavy-tailed) recovery times, holding the mean infectious period fixed?**

The practical stake is whether simple “connectivity × mean infectious period” rules used in network epidemiology remain safe guides for control when real recovery-time variability is long-tailed.

### 1.2 Hypothesis (as in the proposal)

Heavy-tailed recovery **may** make persistence easier than the spectral line suggests at the **same** mean recovery duration \(\mathbb{E}[W]\)—i.e. endemic-like behaviour at lower effective transmission than \(\tau_c \approx 1/\lambda_{\max}\) would imply.

### 1.3 Research questions (RQ1–RQ3)

| ID | Question |
|----|----------|
| **RQ1** | Under **exponential** (Markovian) recovery, how closely does an **operational** simulation threshold track \(\beta_{\mathrm{pred}} = 1/(\lambda_{\max}\mathbb{E}[W])\) (equivalently \(\tau_{\mathrm{pred}} = 1/\lambda_{\max}\))? |
| **RQ2** | When recovery is replaced by **Tang-style** power-law or lognormal laws at fixed \(\mathbb{E}[W]\), how does the **empirical** threshold **shift** relative to exponential and to \(\beta_{\mathrm{pred}}\)? |
| **RQ3** | Are differences between regimes clearer in **threshold location** or in **extinction times**, **variability**, and **late-window prevalence**? |

Operational definitions already implemented in code: survival rate = fraction of runs with extinction flag `< 0.5` at the time limit; \(\hat{\beta}_{\mathrm{surv}\,50}\) = smallest infection probability on the sweep grid with survival rate \(\geq 0.5\) (with fallback to the largest grid point if the curve never crosses 50% below the top of the grid—**right-censoring**). See `scripts/threshold_estimators.py`.

---

## 2. Workspace map (what exists where)

| Area | Role | Key paths |
|------|------|-----------|
| **NetLogo ABM** | grp-SIS on exported graphs; BehaviorSpace sweeps | `netlogo/virus_simulation.nlogox`, `netlogo/behaviorspace_experiments.xml` |
| **Pipeline** | Headless runs, \(\lambda_{\max}\) from edge CSVs, aggregation | `scripts/run_pipeline.py`, `scripts/analysis_utils.py`, `scripts/threshold_estimators.py` |
| **Outputs** | Edges, \(\lambda_{\max}\), empirical thresholds | `output/edges/`, `output/lambda_max.csv`, `output/empirical_thresholds.csv` |
| **Benchmark (Python)** | Parallel numerical engine / checks (faster exploration than full ABM for some questions) | `benchmark/sis_engine.py`, `benchmark/sis_numba.py`, `benchmark/spectral.py`, `benchmark/output/threshold_summary.csv` |
| **Written report (results)** | Methods, notation, figures, discussion | `report/report.tex` |
| **Proposal & peer review** | Original aims and external critique | `peer_review/Research_Proposal_Jan Ruijgrok.tex`, `peer_review/peer_review_jr.md` |

**Experiments (XML names):** `01-export-networks`; `02`–`04` threshold sweeps on non-ER topologies; `05`–`07` ER baselines (exponential, power-law, lognormal); `08-dynamics-trajectory-ER` for microscopic vs mean-field prevalence trajectories.

---

## 3. State of play (evidence snapshot)

This section summarises **current** pipeline outputs so exploration can build on numbers that are actually in the repo.

### 3.1 Spectral reference (ER seed 10001)

From `output/lambda_max.csv`, the exported ER graph has \(\lambda_{\max} \approx 7.181\) (N = 2000, mean degree \(\approx 5.98\)). With \(\mathbb{E}[W]=5\) ticks (as in the proposal narrative), \(\beta_{\mathrm{pred}} = 1/(\lambda_{\max}\mathbb{E}[W]) \approx 0.0279\) and \(\tau_{\mathrm{pred}} = 1/\lambda_{\max} \approx 0.139\).

The same file records **heterogeneous** mean-field predictions (`tau_pred_heterogeneous_mf`) that can be compared to the homogeneous line—useful if you extend the discussion beyond \(c=1\).

### 3.2 Empirical survival thresholds (`output/empirical_thresholds.csv`)

On the **same** ER instance, \(\hat{\beta}_{\mathrm{surv}\,50}\) is **above** \(\beta_{\mathrm{pred}}\) for all three recovery laws (ratios \(\approx 1.58\)–\(1.72\) in \(\beta\)-space). Exponential and lognormal estimates **hit the upper end of the \(\beta\) sweep** (\(0.048\)), so they are **censored** on the grid: the true 50% crossing may lie beyond the largest \(\beta\) tested. Power-law Tang recovery crosses at \(0.044\) **inside** the grid.

**Interpretation sketch:** The mean-field spectral line behaves as a **conservative** (optimistic-for-control) benchmark for *this* discrete-time ABM: common survival occurs at higher \(\beta\) than \(\beta_{\mathrm{pred}}\). The ordering between heavy-tailed laws is **not** a single “heavy tail \(\Rightarrow\) lower threshold” story on this grid—lognormal matches exponential at the cap—while **extinction-time summaries** differ (e.g. median final tick among extinct runs higher under power-law in the exported table). That supports treating **RQ3** (threshold vs dynamics) as central.

*Note:* `benchmark/output/threshold_summary.csv` may reflect a different sweep or aggregation; treat it as a **cross-check** and reconcile filenames and run IDs when citing numbers in the final thesis.

### 3.3 Alignment with `report/report.tex`

The LaTeX report already documents notation, survival curves, threshold table, and limitations (censoring, single graph). **Further exploration** should either extend those analyses or replace coarse estimators—see below—not duplicate the same plots without new content.

---

## 4. Gaps that limit decisive answers

Derived from the proposal’s risk section, `peer_review/peer_review_jr.md`, and the current outputs:

1. **Grid censoring** — When survival is already high at the largest \(\beta\), \(\hat{\beta}_{\mathrm{surv}\,50}\) is not identified; extend \(\beta\) upward or use **model-based** interpolation (`logistic_ed50_linearized` in `threshold_estimators.py`).
2. **Single network realisation** — Results are for one ER draw (seed 10001). Process noise is averaged via BehaviorSpace repetitions, but **topological** uncertainty is not.
3. **Discrete time vs continuous-time theory** — The benchmark \(\tau_c \approx 1/\lambda_{\max}\) is continuous-time folklore; grp-SIS theory does not guarantee the same line for non-Markovian recovery. The study correctly uses the line as a **benchmark**; exploration should separate (i) finite-size/discrete-time effects visible in the **exponential** case from (ii) **shape-of-\(W\)** effects.
4. **Stopping time** — Runs end at `max-ticks` or absorption; “survival” mixes **true persistence** with **slow extinction**. Quasi-stationary or longer horizons may be needed for some parameters.
5. **Underused topologies** — `output/lambda_max.csv` includes **Lat4** and **Ring** edge lists; experiments `02`–`04` target threshold behaviour on those structures but are not yet the focus of the main written report.

---

## 5. Prioritised exploration agenda

### Tier A — High value, relatively low cost

| Action | Rationale | RQ |
|--------|-----------|-----|
| **Extend \(\beta\) grid** upward (and optionally refine near \(\beta_{\mathrm{pred}}\)) | Reduces censoring for exponential and lognormal | RQ1–RQ2 |
| **Report \(\hat{\beta}\) from logistic ED50** alongside grid-based \(\hat{\beta}_{\mathrm{surv}\,50}\) | Smoother threshold when points are sparse | RQ1–RQ2 |
| **Stratify extinction times by \(\beta\)** (not only global medians) | Sharpens RQ3 | RQ3 |
| **Reconcile** `output/empirical_thresholds.csv` vs `benchmark/output/` | One authoritative table for the thesis | All |
| **Short methods paragraph** on bootstrap CIs already in pipeline columns | Addresses peer-review ask for uncertainty | All |

### Tier B — Structural replication

| Action | Rationale | RQ |
|--------|-----------|-----|
| **Repeat experiments 05–07** for **several ER seeds** (new `expt-seed`s with matching edge export) | Separates graph noise from recovery-law effects | RQ2–RQ3 |
| **Report distribution** of \(\hat{\tau}_{50}/\tau_{\mathrm{pred}}\) across seeds | Quantifies “how universal” the \(c>1\) gap is | RQ1 |

### Tier C — Topology and mechanisms

| Action | Rationale | RQ |
|--------|-----------|-----|
| Run **02–04** systematically and compare to ER using the **same** recovery implementations | Tests whether spectral benchmark shifts are topology-specific | RQ2–RQ3 |
| Use **heterogeneous** \(\tau\) prediction from `lambda_max.csv` in discussion | Nuances the \(c=1\) story | RQ1 |
| **Experiment 08** trajectories at multiple \(\beta\) near threshold | Links microscopic–mean-field agreement to survival outcomes | RQ3 |

### Tier D — Theory and fast numerics

| Action | Rationale | RQ |
|--------|-----------|-----|
| Use `benchmark/` **SIS engine** to explore parameter planes **without** full ABM cost | Hypothesis generation and sanity checks | All |
| Read **Tang et al. (grp-SIS)** and **Cator et al.** for predictions testable in simulation (e.g. conditions for heavy tails to increase \(\mathbb{E}[I]\) at fixed \(\mathbb{E}[W]\)) | Grounds the hypothesis in formal results | Background |
| **Quasi-stationary** estimators (if time): condition on non-extinction | Cleaner “endemic” threshold than survival-at-fixed-horizon | RQ3 |

---

## 6. What would “settle” the hypothesis?

The original hypothesis is **not** refuted or confirmed by a single grid on one graph. A stronger test would show **at least one** of:

- A **consistent** ordering of \(\hat{\beta}_{\mathrm{surv}\,50}\) (or ED50) across **replicated** networks: heavy-tailed laws **lower** than exponential relative to \(\beta_{\mathrm{pred}}\), **or**
- A **mechanism-linked** explanation if not: e.g. censoring, time horizon, or discrete-time effects dominate **small** shifts between laws.

Pre-register (for the write-up): **primary estimator** (grid 50% vs logistic ED50), **number of seeds**, **repetitions per \((\beta,\) regime\()\), and **max-ticks**.

---

## 7. Reproducibility checklist

- Set `NETLOGO_HOME` and run `python scripts/run_pipeline.py --all` from the project root (see docstring in `run_pipeline.py`).
- After editing `behaviorspace_experiments.xml`, sync experiments into `virus_simulation.nlogox` (NetLogo comment constraints noted in the script).
- Version-control `output/*.csv` used in the thesis, or archive a run bundle with seed list and NetLogo version.

---

## 8. Suggested reading (short list)

- Tang et al. (grp-SIS / general recovery) — motivates power-law and lognormal implementations.
- Cator et al. (general infection and cure times on networks).
- Van Mieghem (*N*-intertwined SIS) and Pastor-Satorras et al. (review) — spectral thresholds and limits of mean field.
- Chakrabarti et al. — spectral threshold in real networks.

Full bibliographic entries appear in `report/report.tex`.

---

## 9. Summary

The workspace already implements a **coherent** pipeline from **network export** → **\(\lambda_{\max}\)** → **replicated BehaviorSpace sweeps** → **threshold estimators** and a **draft results report**. The **most informative** next steps are: **remove grid censoring** and/or use **logistic ED50**, **replicate across ER seeds**, and **analyse extinction and prevalence stratified by \(\beta\)** to answer **RQ3** with the same rigour as threshold location (**RQ1–RQ2**). The **benchmark** folder supports faster numerical exploration to complement—not replace—the NetLogo ABM where the research question is explicitly about **agent-level** heavy-tailed recovery on an **explicit** graph.

This report is the single place that ties **proposal intent**, **peer-review gaps**, **code**, and **numeric outputs** into one exploration roadmap for the remainder of the project.
