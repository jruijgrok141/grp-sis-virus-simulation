# Recommendations for highest possible score

Based on feedback in `review_report/result_marking.tex` (compiled as `result_marking.pdf`). Overall score noted: **8.6/10** (unweighted mean of sections roughly **8.5–9.0**).

## Cross-cutting (biggest leverage)

The written summary says the main gap is not *acknowledging* limitations, but **weakening them in the core analysis**:

1. **Graph variability (single ER instance)**  
   - Run the same pipeline on **multiple independent ER graphs** (same \(N\), same expected degree / edge count rule as now), with **reported uncertainty** (e.g. median/IQR or CI across graphs for \(\hat{\beta}^{\mathrm{surv}}_{50}\) and ratios vs. \(\beta_{\mathrm{pred}}\)).  
   - In the main text, **one figure or table** should show “across-graph” spread, not only appendix.

2. **Censoring and horizon**  
   - **Finer \(\beta\) grid near threshold crossings** so survival curves and 50% rules are not overly sensitive to censoring.  
   - Where censoring bites, add **sensitivity analysis**: e.g. alternative horizons or a clear rule for “inconclusive” runs; optionally **time-to-extinction** or **hazard-style** summary as a complement to survival-at-\(T\).

3. **Tight RQ ↔ evidence loop**  
   - End each results subsection with **one sentence**: which RQ, which figure/table, **quantitative** takeaway (effect size / ordering), and **one caveat**. That pushes “depth” and “conclusions” scores without fluff.

---

## By rubric section (how to push 8.5–8.8 → 9+)

### 1 Research question(s) (8.8)

- Add **one short paragraph** upfront: *decision the reader should make* with your results (when spectral bound is safe vs. when not), tied to your delimitations.  
- If programme weighting favours “debate”, add **one explicit contrast** with an alternative hypothesis (e.g. strong heavy-tail effect vs. your conservative finding) and what would falsify each.

### 2 Literature (8.5)

- **Synthesize** 2–3 papers in one place (not only cite): e.g. what exactly is proven vs. heuristic for non-Markovian / heterogeneous recovery in relation to \(\lambda_{\max}\) scaling.  
- Add **1–2 sentences** on **finite-size / stochastic** literature if not already central—bridges your ABM–spectral gap story.

### 3 Method (9.0)

- Already strong; preserve clarity. Small upgrades: **pseudo-code or flow diagram** for the pipeline (ABM → export → \(\lambda_{\max}\) → bootstrap), and **explicit stopping rules** for simulations (max steps, failure modes).

### 4 Results (8.5)

- Move **ensemble / multi-graph** and **near-threshold** analyses from “future work” into **primary results** where possible.  
- Ensure every key plot has **uncertainty** (bands or intervals) that is interpretable—not only noted as “degenerate” in captions.

### 5 Conclusions / discussion (8.8)

- **Quantify** main claims (e.g. typical ratio ranges across laws and across graphs).  
- Separate **“what we showed”** vs **“what we infer about real networks”** in two short paragraphs to sharpen generalisability.

### 6 Structure (9.0)

- Keep as is; only ensure the **new ensemble/censoring analyses** have a clear home in Methods + Results so the story stays linear.

### 7 Form and presentation (8.5)

- **Thin the notation** in one pass for the main line: move secondary definitions to an appendix; keep one “reader’s cheat sheet” table.  
- Check **consistent symbols** and **caption self-containment** (each figure readable without the main text).

### 8 Responsible research (8.5)

- Add a **short reproducibility box**: exact commit hash or archive, **one command** to reproduce main figures, and **runtime / hardware** order-of-magnitude.  
- If any parameter was tuned, state **pre-registration-style** choices (fixed before runs) vs. exploratory analyses.

---

## Priority order if time is limited

1. Multi-graph ER ensemble + uncertainty in **main** results.  
2. Finer \(\beta\) sweep near thresholds + censoring sensitivity.  
3. Literature synthesis paragraph + sharper quantitative conclusions.

That sequence matches the assessor’s closing remark and targets the sections at **8.5** (Literature, Results, Form, Responsible research) while preserving strongest areas (Method, Structure).
