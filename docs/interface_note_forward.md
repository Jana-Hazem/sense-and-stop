# Forward sensing interface and hand-off to Member 5

Member 4, draft v0.1. Computational study, not clinical. Rat-derived model scale and pediatric epilepsy EEG (used only as background noise) are not CSD or migraine recordings. EEG is never labelled as CSD.

## 1. Observation model

```
sensors[t, ch] = scale * ( H_bip @ ( A @ u[t].ravel() ) )[ch]  +  noise_gain * eps[t, ch]
```

| Symbol | Meaning | Unit | Shape |
|---|---|---|---|
| `u` | hidden wave from Member 3, **not clipped** | normalised | `[271, 201, 301]` = `[t, y, x]` |
| `A` | area-averaging, PDE nodes to sources | none (rows sum to 1) | sparse `[73, 60501]` |
| `H_bip` | bipolar lead field | V per A*m | `[21, 73]` |
| `scale` | dipole moment per source at u = 1 | A*m | scalar, **illustrative** |
| `eps` | empirical EEG background, decimated to the PDE time grid | V | `[271, 21]` |
| `sensors` | synthetic observation | V | `[271, 21]` |

Node index in `A`: `iy * 301 + ix` (x varies fastest), i.e. `u[t].ravel()` order.
Time grid: 271 samples, 2 s apart, 0 to 540 s. Space grid: dx = 0.1 mm, patch 30 mm (x) by 20 mm (y).
Function: `make_sensors(u, H_bip, A, scale_am, noise_V, rng, noise_gain=1.0)` in `sensors.py` returns `(sensors, clean)`.
Scales run so far: 1e-9, 1e-8, 1e-7 A*m per source (seed 0). Report results at all three; none is calibrated to physiology.

## 2. Files (kept outside the repo, never committed)

| File | Contents |
|---|---|
| `H_bip.npz` | `H_bip` (21x73), `M` (21x21 derivation matrix), `electrodes` (21), `channels` (21), `vertno`, `units` |
| `A_matrix.npz` | sparse `A` (73 x 60501) |
| `flat_patch.npz` | `vertno`, `xy_mm` (73x2, patch coordinates), `pos_mm`, `normals`, `hemi`=0, `label`, `dx_mm`, `nx`, `ny` |
| `noise_2s.npz` | `noise_V` `[3 files, 21 ch, time]`, `channels`, `dt_s`=2 |
| `noise_cov.npz` | `C` (21x21, V^2), `R` (correlation), `channels` |
| `sensors_default.npz` | default synthetic run, including `clean_V` `[271, 21]` |
| `u_2d_default.npy` | default 2D wave, float32, about 65 MB, shared by other means |

## 3. Channels and conventions

- Order (21): FP1-F7, F7-T7, T7-P7, P7-O1, FP1-F3, F3-C3, C3-P3, P3-O1, FP2-F4, F4-C4, C4-P4, P4-O2, FP2-F8, F8-T8, T8-P8, P8-O2, FZ-CZ, CZ-PZ, T7-FT9, FT9-FT10, FT10-T8.
- Each channel is **first electrode minus second**, so it is reference-free. Orientation is fixed to the surface normal. Units are V per A*m.
- Dropped from the 23 EDF channels: `P7-T7` (sign-flipped copy of `T7-P7`) and `T8-P8-1` (duplicate), by `build_Hbip.py`.
- `channels` in the npz files are upper-case (`FP1-F7`); `electrodes` use `Fp1`.

## 4. Source space and patch

- fsaverage template (adult), `lateraloccipital-lh`, ico-5 source space, mindist 5 mm. The patch location is a design choice, not data-based.
- Flattening: geodesic distances on the surface, embedded in 2D, then cut to the 30 x 20 mm rectangle. Result: **73 sources** (220 sources lie inside the 25 mm geodesic disc).
- Saved patch: x from 1.5 to 29.9 mm, y from 0.0 to 19.7 mm. Nearest-neighbour source spacing: min 1.84, mean 2.45, max 3.04 mm.
- The pulse is 1.85 mm wide, about the same as the source spacing. This is why `A` averages over area. Treat sub-spacing detail as not recoverable.
- Flattening distortion (flat distance divided by on-surface distance): edges min 0.92, mean 1.08, max 1.27; all pairs, 5th to 95th percentile 0.92 to 1.12. So flat distances are good to roughly +/-10% typically and up to about +27% on single edges.
- No source lies within the first 1.5 mm of x. The sources there own those nodes through nearest-source assignment (Voronoi), so the first 1.5 mm is represented only by the nearest sources.

## 5. The default wave

- The 2D wave is **exactly uniform in y** at every saved time (max deviation across y = 0). The stimulus is a full-width strip: u = 1 for x from 0 to about 2.9 mm and all y at t = 0.
- So the default 2D dataset contains only x-profile information. Recovering a 1D u(x, t) is sufficient, and y-structure cannot be tested with this run.
- Amplitude is still settling in the first 60 s (max u 1.00 at 0 s, 0.83 at 60 s, 0.82 at 100 s). Cut or flag t < 60 s.
- u undershoots below rest behind the pulse (min -0.246 at 100 s, -0.216 from 200 s on; Member 3 reports about -0.2 to -0.3). **Do not clip u to [0, 1].**
- The model recovers faster than the data (decay 28 s vs 41 s), so the pulse tail at the sensors is somewhat short.
- Front speed is from literature, not verified by Member 2 yet.

## 6. Noise

- Source: 3 seizure-free files of one subject (chb01_01, chb01_02, chb01_08), pediatric epilepsy data used only as background.
- Anti-aliased decimation 256 Hz to 0.5 Hz (three stages of 8), 60 s trimmed from each end. This removes most normal EEG content, so the noise is band-limited.
- Noise is correlated across channels. Standard deviation is about 2.2 to 2.5 uV on the left channels but 10.8 uV (P4-O2) and 12.0 uV (P8-O2) on the right (values from `noise_cov.npz`; the per-file mean printed by `make_noise.py` is lower, 8.7 and 10.5 uV). The right-side channels carry very little usable information.
- `make_sensors` draws one random file and one random contiguous 271-sample segment per call. With only 3 files, repeated draws overlap.
- Slow drift, measured on the 3 noise files (18 windows of 540 s): the linear trend over 540 s has median 0.02 to 0.06 uV and 95th percentile 0.09 to 0.20 uV on the left channels, and 1.3 to 1.5 uV (95th percentile) on P8-O2 and P4-O2. A linear trend explains under 0.5% of the variance (median). This is small next to the signal at scales 1e-7 and 1e-8, but not at 1e-9, where the peak clean signal is only 0.13 to 0.58 uV. The EDF prefilter field in the `chb01_01` header is empty, so the files give no information about hardware filtering (empty is not the same as unfiltered). If the recordings were high-pass filtered, real slow drift is under-represented here.

## 7. Signal-to-noise (upper bounds)

Defined as peak |clean| for t >= 100 s divided by sqrt(diag C). SNR scales linearly with `scale`.

| Channel | Peak clean (uV, 1e-7) | SNR 1e-7 | SNR 1e-8 | SNR 1e-9 |
|---|---|---|---|---|
| P3-O1 | 42.6 | 19.6 | 1.96 | 0.20 |
| P7-O1 | 58.0 | 24.8 | 2.48 | 0.25 |
| T7-P7 | 37.6 | 16.4 | 1.64 | 0.16 |
| C3-P3 | 19.7 | 8.6 | 0.86 | 0.09 |
| P8-O2 | 15.2 | 1.3 | 0.13 | 0.01 |
| P4-O2 | 12.9 | 1.2 | 0.12 | 0.01 |
| Joint, whitened norm of the clean signal | | 38.9 | 3.9 | 0.39 |

For t >= 60 s the peaks are much larger (P3-O1 118.1, P7-O1 88.7, T7-P7 53.1, C3-P3 33.6, P8-O2 20.3, P4-O2 16.0 uV at 1e-7). The sensor signal is therefore strongest at 60 to 100 s and about 2 to 3 times lower later. Cause not tested.

## 8. How much is recoverable

- `H_bip` (21x73): at most 19 independent directions by construction (21 channels minus 2 electrode loops). The 19th is at float32 noise level, so numerical rank is 18 at 1e-6 and 10 at 1e-3 (relative to the largest singular value).
- x-collapsed operator (21x301): numerical rank 16 at 1e-6 and 7 at 1e-3.
- Noise-aware (whitened, scale 1e-7): top singular values 56, 23, 5.3, 5.1, 0.39, 0.08. **4 directions above 1, 5 above 0.3**, stable across eigenvalue floors 1e-6, 1e-4, 1e-2. Because the values scale linearly with `scale`, at 1e-8 only 2 are above 1 (4 above 0.3), and at 1e-9 none are above 1 (1 above 0.3).
- Position along x is hard to read from the channel pattern alone. For shifts of 1 to 5 mm the patterns are 97.5% similar or more (cosine similarity, strip operator). Timing and amplitude must carry localisation.
- The strip response is not trustworthy near x = 30 mm (patch edge; P7-O1 changes sign). The default front only reaches about 25.6 mm by 540 s.

## 9. Do not

- Do not clip `u`.
- Do not use the same discretisation for generating data and for inversion (inverse crime). Data are generated on the 0.1 mm grid; invert on a different, coarser representation.
- Do not assume `H_bip` is full rank, and do not invert it without regularisation.
- Do not quote SNR without the scale it was computed at, and do not treat 1e-7 as physiological.
- Do not use P8-O2 and P4-O2 as informative channels.
- Do not describe this as a CSD, EEG-CSD or clinical result.

## 10. Not verified

- FT9 and FT10 rows come from the same forward pipeline as the other electrodes but have not been independently checked (the loop check cancels them algebraically).
- The u-to-volts scale is illustrative. Whether 1e-7 A*m per source is plausible has not been checked against a source.
- Whether the sensor-signal decline after 60 s is due to the geometry or something else has not been tested.
- Adult template head versus pediatric EEG: location of sources should not be over-interpreted (D4). S3 template license is unconfirmed.
- The EDF prefilter field (read from the `chb01_01` header only) is empty, so it is unknown whether very slow drift was filtered out before the data were stored. The CHB-MIT documentation may say; this has not been checked.
- Drift has been measured in the noise (section 6) but no inversion stress test with added drift has been run. That needs Member 5's solver.
