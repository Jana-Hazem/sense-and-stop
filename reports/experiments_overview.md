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
| Hardware | 13th Gen Intel Core i7-1365U, 16 GB RAM |
| E2 analysis (Member 2) | Python 3.13.12, pandas 2.3.3, numpy 2.2.6, matplotlib 3.10.7 |
| Solver (Member 3) | src/model/solver.py: excitable (FitzHugh-Nagumo type) PDE, finite differences, reflecting boundaries, RK4 time stepping; configs/model_default.yaml, dt 0.1 s, dx 0.1 mm. Package versions: [FILL: M3] |
| Forward model (Member 4) | src/forward/: A_matrix (73 x 60501, sparse) and H, stored in A_matrix.npz and H_bip.npz. Package versions: [FILL: M4] |
| Inversion (Member 5) | src/inverse/run_member5.py: Method A (Tikhonov free-field inversion) and Method B (travelling-pulse least-squares fit). Python 3.13.9, NumPy 2.3.5, SciPy 1.16.3, pandas 2.3.3, recorded in the Member 5 run metadata.json. The pulse template uses solve_ivp RK45 (rtol 1e-6, atol 1e-8) at twice the data grid spacing. |
| Code version | [FILL: repo commit hash from the E6 clean-checkout run] |
| Random seeds | E2 dot plot: 20261119 (figure_index.csv). E4: 9 evaluation seeds, values not recorded in the run metadata. E1, E3, E5: [FILL: owners] |

Package versions differ between members' machines. The versions above are those each stage was run with.