# 2026/08/19 Lab Meeting — Speaker Notes
Aging Adipose Atlas: 3-Month vs. 16-Month WT Mouse Cell-Type Annotation

(Same text is also embedded in the "Notes" pane of each slide in `atlas_aging_report.pptx` — this file is just a convenient standalone copy for printing/rehearsing.)

---

## 1. Title
Hi everyone. Today I'll walk through the aging atlas we built from 3-month and 16-month WT mouse adipose tissue — the focus is our own cell-type annotation, and at the end I'll validate it against PanSci's original labels. This should take about 5 to 10 minutes.

## 2. Motivation
Starting with motivation. Aging is a central variable in adipose biology, but most single-cell atlases only annotate a single timepoint. We processed 3-month and 16-month WT mice together in one pipeline, 518,332 cells total, so cell populations can be compared directly across age under the same standard. The core question: does our own marker-gene plus differential-expression annotation pipeline agree with PanSci's original labels, cluster by cluster?

## 3. Experimental Design
The dataset covers three adipose tissues — brown, subcutaneous, and visceral — WT genotype only, 3-month plus 16-month, 518,332 cells total. One key design choice: Harmony integration only uses tissue as the batch key, not age, because age is the biological signal we want to observe, not a technical batch effect to remove. Annotation starts with marker-gene scoring, cross-checked with differential expression, and any ambiguous cluster gets further sub-clustering. Finally we cross-tabulate our results against PanSci's own original labels.

## 4. Annotation Pipeline
All six steps run on the same combined 3-month + 16-month WT dataset — one consistent script, not separate runs per age group.

**Step 1, QC filter**: we remove low-quality cells and genes (WT and age filtering is applied here too, before doublet detection — that way we never waste time running Scrublet on cells we'd discard anyway).

**Step 2, doublet removal**: Scrublet runs per mouse (batched by animal ID, not pooled, since doublet rates can differ mouse to mouse).

**Step 3, normalize + scale**: normalization, HVG selection, and scaling all use tissue as the batch key (not age — age is the biological signal we want to keep, not remove).

**Step 4, PCA + Harmony**: PCA down to 50 components, then Harmony integration (same batch key, tissue only, for the same reason — age has to survive this correction so it can still be compared afterward).

**Step 5, Leiden clustering**: we cluster the integrated, batch-corrected embedding.

**Step 6, cell-type annotation**: marker scoring plus differential-expression validation — any cluster without one clean matching signal gets sub-clustered further, which is exactly how the rare types on the next slide were found.

## 5. Results: Final Cell-Type Annotation
This is the final cell-type annotation — 15 cell types across 518,332 cells. Five of them are rare populations that the cluster-level pass couldn't resolve on its own; they only emerged through targeted sub-clustering: Rdh16+ epithelial-like cells, neural cells, mesothelial cells, epididymal cells, and lymphatic endothelial cells — the last two are still marked tentative since their signal is weaker. This UMAP also carries an extra check: instead of using the original HVG space, the embedding was recomputed using only each type's top 50 differentially-expressed marker genes — 608 genes total — and it still separates all 15 types cleanly, showing those markers alone are enough to support this clustering result.

## 6. Results: Validation Against PanSci's Labels
This is the similarity matrix — the x-axis is PanSci's original labels, the y-axis is our final annotation, and brighter yellow means a higher proportion. The diagonal is very clean, which means our independently-derived annotation largely agrees with PanSci's own labels. 11 out of 15 types exceed 85% concordance with PanSci. The remaining 4 — DCs, skeletal muscle cells, epididymal cells, and lymphatic endothelial cells — are genuinely rare populations that only came out through sub-clustering; the last two are still marked tentative since their signal is weaker. We honestly document that their signal isn't as clean, rather than treating it as an unfinished analysis.

## 7. Key Takeaways
To summarize in three tiers. Tier one: 7 types match PanSci at 95% or above — T cells, Macrophages, Pericytes, very clean signal. Tier two is 85 to 95%, covering ASPCs, neural cells, mesothelial cells, and adipocytes. Tier three is the genuinely rare, noisier group of 4 — DCs, skeletal muscle cells, and two tentative populations: epididymal cells and lymphatic endothelial cells — all of which only surfaced through extra sub-clustering. They sit at 56 to 80% concordance, but their destination is clearly documented.

## 8. Next Steps
This talk covered building the aging atlas and validating our cell-type annotation. The next step is to turn this into a poster, a 10-minute presentation combining this report with the July 31 preprocessing-benchmark talk.

## 9. Thank You
Thank you — that's the report. Happy to take questions.