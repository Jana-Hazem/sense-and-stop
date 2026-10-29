# Interface Note: PDE Model Output (Member 3 → Members 4, 5, 6)

**Version:** 0.2 (calibrated to S1 preliminary data) · **Module:** `src/model` · **Status:** draft for review
If you need a change to anything below, open an issue; do not edit arrays after the fact.

## 1. Calling the model

```python
from src.model import (simulate, load_config, with_resolution, save_result, load_result,
                       pulse_metrics, pulse_speed, calibrate_excitable)

cfg = load_config("configs/model_default.yaml")       # excitable, 1D by default
res = simulate(cfg)                                   # uncontrolled, primary code path
res = simulate(cfg, control=g)                        # Member 6 only: g(x, t) in 1D, g(X, Y, t) in 2D
save_result(res, "outputs/run_001.npz")               # arrays + JSON metadata in one file
```

For 2D set `cfg["model"]["dimension"] = 2`.

## 2. Returned dictionary

| Key | Shape | Unit | Meaning |
|---|---|---|---|
| `u` | **1D:** `[nt, nx]` · **2D:** `[nt, ny, nx]` | dimensionless | Hidden state. Rest = 0; pulse peak ≈ 0.81 (excitable) |
| `v` | same as `u`, or `None` | dimensionless | Recovery variable (excitable only) |
| `t` | `[nt]` | s | `0, 2, 4, …, 540` (exact multiples of `output_dt`) |
| `x`, `y` | `[nx]`, `[ny]` or `None` | mm | Node coordinates; nodes lie **on** the boundaries |
| `units`, `params`, `meta` | dict | | Units; parameters used; solver info and a full copy of the config |

**Axis order is time, then y, then x: `u[t, y, x]`.** `u[:, iy, ix]` is the time series at one node.

**Default sizes and cost** (dx = 0.1 mm, T = 540 s, output every 2 s):

| Run | `u` shape | Run time | Memory for `u` (and the same again for `v`) |
|---|---|---|---|
| 1D | 271 × 301 | ~1 s | < 1 MB |
| 2D (30 × 20 mm) | 271 × 201 × 301 | ~33 s | 131 MB |

Cost scales with the number of grid nodes: doubling both Lx and Ly gives about 4× the time and memory.

## 3. What the model represents (all from S1 preliminary data and the literature)

| Quantity | Value | Source |
|---|---|---|
| Temporal FWHM at one point | 42.8 s | S1-measured (rat, 3 animals) |
| Front speed | 2.6 mm/min = 0.0433 mm/s | Literature-derived, **not yet verified** |
| Pulse width (speed × FWHM) | **1.85 mm** | Derived |
| Sensitivity range for width | 1.4 to 6.4 mm | Literature broad speed range 2–9 mm/min |

## 4. Things downstream code MUST respect

1. **Do not clip or constrain `u` to [0, 1] for the excitable model.** `u` undershoots to about −0.2 to −0.3 behind the pulse. A [0,1] box constraint would make the true state unreachable. (Fisher–KPP stays within [0,1].)
2. **Units are mm and s.** Convert explicitly if H or the EEG pipeline uses metres, ms or Hz.
3. **Grid nodes include both boundary points.** Define whether the source map is node-centred (as here) or cell-centred.
4. **Parameters come from rat data plus unverified literature speed.** Never describe them as human or as final.
5. **Do not edit D and k by hand.** Re-calibrate with `calibrate_excitable(cfg, fwhm_s, speed_mm_per_min)`; the pulse is fragile in (a, ε, γ).
6. **Resolution:** the front needs dx ≤ ~0.1 mm. Do not coarsen `dx` in the PDE solve to match sensors; map the fine field to the sensors instead (Sec. 5). For the inversion grid, check `calibrate_excitable(...)["suggested_max_dx_mm"]`.

## 5. For Member 4: domain size, mapping and open questions

**The rat data cannot set the domain size**; it must come from your head-model patch. The model accepts any size.

**Proposal (to agree):**
- **2D flat patch 30 mm (x, direction of travel) × 20 mm (y)**; 1D line of 30 mm for debugging only.
- The pulse (1.85 mm wide) travels about 23 mm in the 540 s run, so a 30 mm patch holds the whole wave life.
- If you need a larger patch, tell me; a 60 × 40 mm patch would take about 2 min and ~1 GB per run.

**Questions for you:**
1. Which area of the template cortex does the patch represent (e.g., occipital), and what are its dimensions in mm?
2. What is the spacing of your source space? **The pulse is 1.85 mm wide, which is likely smaller than a typical source spacing (a few mm).** Proposal: H should act on a *source-averaged* field: for each source point, average the PDE nodes in its area (a sparse averaging matrix A), so the full map is H·A.
3. 1D or 2D for the first end-to-end test? (Recommendation: 2D for anything that touches H.)
4. Do you want the wide-pulse scenario (6.4 mm, 9 mm/min) as a second test case? It is one call to `calibrate_excitable`.

## 6. For Member 5

- Reconstruction target: `u[t, (y,) x]` on the axes above. For the excitable model, decide whether `v` is reconstructed (default: not observed, not reconstructed).
- Avoid the inverse crime: generate the truth on a fine grid with one scheme and invert on a coarser grid with an independent scheme:

```python
truth = simulate(with_resolution(cfg, dx=0.05, dt=0.05, method="rk4"))     # 601 nodes
inv   = simulate(with_resolution(cfg, dx=0.1,  method="solve_ivp"))        # 301 nodes
ref_on_inv_grid = truth["u"][:, ::2]     # fine nodes coincide with coarse nodes
```

- For this pair the noise-free numerical gap between the two solutions is **2.9 % relative RMS (max abs 0.057)**. This is the error floor your solver sees from discretisation alone, before noise or model mismatch.
- Use `pulse_metrics(result, x_probe)` on reconstructions to compare FWHM, rise and decay with the truth in the same definitions as the S1 data.

## 7. For Member 6 (control extension)

- `simulate(cfg, control=g)` subtracts `g·u` from u_t; `g` returns an array of the grid shape and should be bounded and localised. Tested: a strong localised g blocks the wave.
- Because the pulse is only 1.85 mm wide, a control region needs to be at least that wide; the search space for location/width must be on the mm scale.
- The control search, objective and constraints belong in `src/control/`.

## 8. Reproducibility

Configuration lives in `configs/model_default.yaml` (each value tagged with its provenance); a copy is stored in every result's `meta["config"]`. The solver is deterministic. Run the checks with `python -m pytest tests/ -v` (21 tests, about 20 s).
