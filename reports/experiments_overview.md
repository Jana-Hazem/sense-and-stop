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

### Computing environment

| Item | Value |
|---|---|
| Operating system | Windows 11 Pro, build 26300 |
| Hardware | 13th Gen Intel Core i7-1365U, 16 GB RAM |
| Python | 3.13.9 (shared project environment) |
| NumPy | 2.3.5 |
| SciPy | 1.16.3 |
| pandas | 2.3.3 |
| Matplotlib | 3.10.7 |

### Software and computational implementation

| Component | Implementation and configuration |
|---|---|
| E1 — PDE solver (Member 3) | `src/model/solver.py`: excitable FitzHugh–Nagumo-type reaction–diffusion PDE, second-order finite differences, Neumann no-flux boundaries implemented using reflecting ghost nodes, and RK4 time integration. Configuration: `configs/model_default.yaml`. Spatial grid spacing = 0.1 mm; actual internal time step = 0.1 s. The solver validation experiments are deterministic and do not require random seeds. |
| E2 — Experimental data analysis (Member 2) | Analysis of rat cortical spreading depolarization recordings using extracted DC and laser Doppler event features. Per-animal summaries and event-level measurements are used for experimental comparisons. The FWHM figure uses plotting jitter seed = `20261119`. |
| E3 — Forward model (Member 4) | `src/forward/`: sparse area-averaging operator `A_matrix.npz` (73 × 60,501) and monopolar average-referenced lead-field matrix `H_matrix.npz` (19 × 73). The bipolar channels are constructed by subtracting rows of the monopolar H matrix. The source patch contains 73 sources on a 301 × 201 grid with 0.1 mm spacing over a 30 × 20 mm domain. The synthetic sensor record in `sensors_default.npz` has dimensions (271, 21) and was generated using seed = `0`. |
| E4 — Inversion (Member 5) | `src/inverse/run_member5.py`: Method A, Tikhonov free-field inversion, and Method B, travelling-pulse least-squares fitting. Uses `H_matrix.npz`, with bipolar channels formed by subtracting monopolar H rows. The pulse template uses `solve_ivp` with RK45 (`rtol = 1e-6`, `atol = 1e-8`) at twice the data grid spacing. Both inversion methods use deterministic computation without random initialization. |
| E5 — Control (Members 6 and 3) | The control implementation uses the PDE solver and differential-evolution optimization with Sobol initialization, followed by Nelder–Mead polishing. The search uses a base seed of 0 and a scenario-specific seed derived from the CRC32 hash of the scenario name, so seed assignment does not depend on scenario execution order. Controlled and uncontrolled simulations use matched initial conditions, parameters, grid and time axis. |

### Background-noise data and evaluation windows (E4)

| Item | Value |
|---|---|
| Dataset | CHB-MIT scalp EEG |
| Selected recordings | `chb01_01`, `chb01_02`, `chb01_08` |
| Subject | One subject (chb01) |
| Stored noise array | `noise_2s.npz` |
| Array dimensions | (3, 21, 1740) |
| Original sampling rate | 256 Hz |
| Processed sampling interval | 2 s (0.5 Hz) |
| Processing | Anti-aliased decimation from 256 Hz to one sample every 2 seconds; edge-trimmed recordings |
| Simulation window length | 271 samples per window |
| Paired windows | 3 non-overlapping pairs per recording |
| Total evaluation windows | 9 paired windows per scale and case |

Each evaluation pair consists of a data window and a separate baseline window used for noise estimation and bad-channel detection. Each recording provides three non-overlapping pairs because 1740 samples allow three pairs of 271-sample windows.

The nine pairs are deterministic window selections, not nine independent subjects or independently generated random seeds. The false-alarm denominator is 18 because both windows in each pair are evaluated as data windows.

### Sensor and operator conventions

- The noise data contain 21 bipolar channels.
- The monopolar lead-field matrix contains 19 electrodes.
- Three bipolar channels cannot be constructed because FT9 and FT10 are absent from the lead-field electrode set: `T7-FT9`, `FT9-FT10`, and `FT10-T8`.
- Therefore, 18 bipolar channels initially have corresponding lead-field rows.
- Bad channels are detected separately for each E4 noise window.
- In the E4 main case, the number of retained channels was 18 in five runs, 16 in three runs, and 15 in one run.
- `H_matrix.npz` is the monopolar operator used in the E4 run. Bipolar channels are calculated by subtracting the corresponding monopolar rows. The `H_bip.npz` label in the metadata hash record is a fallback filename and does not establish the existence of a separately generated bipolar operator.

### Reproducibility and randomization

| Experiment | Randomness and evaluation design |
|---|---|
| E1 — Solver validation | Deterministic numerical experiments; no random seed required. |
| E2 — Rat comparison | FWHM figure jitter seed = `20261119`. |
| E3 — Forward model | Synthetic sensor noise seed = `0` in `sensors_default.npz`. |
| E4 — Inversion | Nine deterministic paired noise-window selections per scale and case. No random initialization in Tikhonov or pulse fitting. |
| E5 — Control | Base search seed = `0`; scenario-specific seeds derived from the scenario-name CRC32 hash. |

For supplementary E4 analyses, the drift test uses seeded drift amplitudes, with seed = `1000 + evaluation-window index`. The E4 plotting jitter uses seed = `0`. These supplementary seeds are separate from the main inversion evaluation.

**Scope of reproducibility:** The nine E4 evaluation pairs all originate from three recordings belonging to one subject. They quantify repeated window-level performance and must not be interpreted as nine independent subjects. E5's seed assignment is documented, but its real-data control results must be reported separately once the experiment is executed.
