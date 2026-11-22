# Glioblastoma cell states under temozolomide

A re-analysis of lineage-barcoded single-cell data (Larsson et al., *Molecular
Systems Biology*, 2021) with the full rate equations, solved and verified by SIRPY.

**Preprint:** *Temozolomide and glioblastoma cell states: a re-analysis of
lineage-barcoded single-cell data with the full rate equations* — bioRxiv [DOI to be
added once posted].

## What the analysis asks

The original study barcoded patient-derived glioblastoma cells (line U3065MG),
assigned each sequenced cell to one of six transcriptional states, and fitted a
linear model of growth and state switching, untreated and under temozolomide (TMZ)
or BMP4. That fit replaced the matrix exponential over one week with `I + tA`, added
a sparsity budget on the switching rates, and tied the treated networks to the
untreated one by a prior.

We asked what the same public data say without those three choices, in three steps.

| Step | Assumes | Question |
|---|---|---|
| 1 | no dynamical model at all | what do a lineage's descendants look like after 8 days? |
| 2 | the authors' rate equations, solved exactly | which switches do the data establish? |
| 3 | rates constant within each 8-day block | does the order of two drugs matter? |

## What we found

**1. Temozolomide removes cells rather than creating resistant ones.**
Over 8 days one state-2 cell leaves 12.0 descendants untreated (95% interval 10.7 to
13.4) and 5.1 under TMZ (4.2 to 6.1) — but the same number in the resistant state 5
(2.54 and 2.52; difference −0.02, −0.96 to 0.92). The entire loss is in the other
states (−6.9, −8.3 to −5.5). The state-5 share rises because the rest are removed.
The one sign of a push into state 5 comes from state-1 cells (+1.5, 0.4 to 2.5).

**2. The data establish few switches.**
Of the 30 possible switches, the untreated time course establishes five at the 5%
level, three of them only just; the TMZ data establish none. The state-2 to state-5
switch, highlighted in the original network, is established in no condition. This
does not mean those switches are absent — it means 8-day lineage counts at these cell
numbers cannot distinguish them from zero once growth rates are free.

**3. Drug order changes the outcome by about 30%.**

| Plan (8 days + 8 days) | Cells per starting cell | State-5 share |
|---|---|---|
| Control, control | 53.8 [47.2, 61.2] | 22% |
| TMZ, TMZ | 14.4 [11.7, 17.8] | 46% |
| BMP4, BMP4 | 20.8 [17.6, 24.6] | 22% |
| **TMZ, then BMP4** | **13.3 [11.1, 16.0]** | 22% |
| **BMP4, then TMZ** | **19.1 [16.7, 21.8]** | 44% |

Ratio TMZ-first to BMP4-first: **0.696 [0.585, 0.827]**. Switch days other than day 8
are model extrapolations with wide intervals; these data pin down only 8-day blocks.

## How we know the numbers

- **390 SIRPY solves of the rate model, every one verified.** A solve whose
  verdict is not VERIFIED is discarded and never used. VERIFIED means the returned
  solution, substituted back into the differential equation, leaves a residual below
  a stated tolerance — a check on the returned function, not a proven bound on its
  error.
- **An independent cross-check.** The model is linear, so the matrix exponential also
  has a closed form. It was used for the global search only; every reported value was
  recomputed with SIRPY, and the two agree at the reported fits to between 2e−14 and
  4e−10.
- **Reproduced on two computers**, matching to every printed digit.
- **Multimodality found and reported.** A first fit from 5 starting points settled on
  a worse local optimum of the BMP4 likelihood (−log L 2521.0894 against 2521.0424).
  The best optimum is reached from roughly one start in five, which is why the final
  analysis uses 40.

## Running it

```
python3 gbm_states.py
```

Needs `numpy`, `scipy`, `matplotlib` and `sirpy`; takes 5 to 10 minutes. It writes
`out_gbm/` with the three figures and `results.json`, which holds every number quoted
above.

| File | What it is |
|---|---|
| `gbm_states.py` | the analysis |
| `extract_data.py` | rebuilds the two cell tables from the authors' `.mat` file (needs `pip install mat-io`) |
| `results.json` | every number from our run (SIRPY 0.67.0) |
| `out_gbm/` | figures from that run |

## Data

The data are the authors' own and stay with them: barcode, day, treatment and state
per cell, their sampling fractions, and the lineages seen at day 8. Get them from
**<https://github.com/idalarsson/STAG>**; `extract_data.py` turns their `.mat` file
into the two CSV tables this analysis reads. We do not redistribute their files.

## Limitations

The authors' state labels and sampling fractions are taken as given; only lineages
present at day 8 contribute; each lineage's day-8 composition is itself a sample
shared by the three arms; the order predictions assume a drug's rates do not depend
on the previous block; and these are cells in culture, not a tumor.

## Citation

Larsson I, Dalmo E, Elgendy R, et al. *Modeling glioblastoma heterogeneity as a
dynamic network of cell states.* Molecular Systems Biology 17:e10105 (2021).
doi:10.15252/msb.202010105

## About

Analysis by Hamidreza Semiyari, Luminary Quantum Institute LLC. SIRPY is LQI's
differential-equation solver; it checks every solution it returns. Prepared with AI
assistance (analysis code and text); all numbers were computed and checked by the
author on his own computer.
