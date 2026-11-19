# Experiments overview

Each experiment answers one research question, compares against a stated baseline,
uses a defined metric, and limits its claim to the conditions tested. All claims
apply to synthetic, one-dimensional, rat-derived conditions and are not clinical.

## Experiments

| # | Research question | Baseline | Main metrics | Owner |
|---|---|---|---|---|
| E1 | Does the PDE solver behave correctly? | Exact solution, a finer grid, another solver | Convergence order, speed, FWHM, rise and decay | M3 |
| E2 | Does the calibrated wave resemble the rat data? | S1 features (FWHM 42.76 s, rise 35.19 s, decay 40.68 s) | Difference per feature, by animal | M2 + M3 |
| E3 | What can the forward operator resolve? | Known source shifts and widths | Rank, singular values, pattern similarity | M4 |
| E4 | Can the wave be recovered from noisy sensors? | Zero guess, free-field Tikhonov | Front error, speed error, relative L2, whitened residual | M5 |
| E5 | Can a hypothetical control suppress the wave? | No control, fixed-time control | Suppression, energy, rebound, violations | M6 + M3 |
| E6 | Does a clean checkout reproduce everything? | none | Pass or fail per stage | M6 + M1 |

The S1 values in E2 are the preliminary target values from three rats
(results/E2_rat_comparison/tables/E2_s1_per_animal.csv). FWHM is the calibration
target. Rise and decay are independent checks.

## Reporting rules used in all experiments

- The number of runs or seeds, the aggregation (median and IQR, or mean and SD) and the meaning of error bars are stated for every result.
- Small samples are shown as individual points.
- Methods are compared on paired seeds, so every method sees the same noise.
- Data used to fit parameters is kept apart from data used to evaluate them.
- Results state what was measured. Interpretation is in the Discussion.
- Every figure has a script, an input data file and a seed, listed in results/figure_index.csv.

## Test environment

| Item | Value |
|---|---|
| Operating system | Windows 11 Pro, build 26300 |
| Memory | 16 GB RAM |
| CPU | [FILL] |
| Python | 3.13 [FILL: exact patch version] |
| E2 analysis (this machine) | pandas 2.3.3, numpy 2.2.6, matplotlib 3.10.7 |
| Solver, forward model, inversion | [FILL: versions from each member's requirements file or the E6 clean-checkout run] |
| Repository commit | [FILL: hash of the commit that produced the final results] |
| Random seeds | [FILL: per experiment, from figure_index.csv] |

The commit hash is recorded when the final results are generated, not before.
