# Interface Note: Member 6B — Hypothetical Wave Control

**Version:** 1.0 · **Module:** `src/control/` · **Status:** complete (depends on Members 3 & 5)

> **HYPOTHETICAL**: This module is a numerical thought experiment. The control term `g(x,t)` has no
> confirmed physical or clinical counterpart. All results are labelled accordingly and must not be
> presented as clinical guidance.

---

## 1. Control Formulation

### Controlled PDE

Member 6B extends Member 3's excitable model with a localised suppression term:

```
u_t  =  D u_xx  +  k·(u(1−u)(u−a) − v)  −  g(x,t)·u        [Member 6B]
v_t  =  k·ε·(u − γ·v)
```

Setting `g ≡ 0` recovers Member 3's model **identically** (same solver, same initial state, same grid).

### Control Term

```
g(x,t) = g0 · zone(x; xc, w, e) · gate(t; ts, dur, ramp)
```

| Component | Formula | Range |
|-----------|---------|-------|
| `zone(x; xc, w, e)` | ½·[tanh((x−(xc−w/2))/e) − tanh((x−(xc+w/2))/e)] | [0, 1] |
| `gate(t; ts, dur, ramp)` | smoothstep((t−ts)/ramp) · smoothstep((ts+dur−t)/ramp) | [0, 1] |
| `g0` | maximum removal rate [1/s] | [0, 2] |

**Hypothetical biophysics**: `g(x,t)·u` removes wave activity locally in space and time — analogous
to focal optogenetic inhibition, cooling, or drug application. The linear `g·u` term preserves
physiological state bounds (u cannot go below 0 due to this term alone).

### Objective Function

```
J(P) = D(P)/D₀  +  w · E(P)/E_ref        (minimise; J=1 means no improvement)
```

| Symbol | Definition | Notes |
|--------|-----------|-------|
| `D(P)` | ∫_prot ∫_t max(u(x,t;P), 0)^p dx dt | wave damage in protected region |
| `D₀` | D(P=0) on the **plant** | uncontrolled baseline |
| `E(P)` | ∫∫ g(x,t;P)² dx dt | control effort (energy) |
| `E_ref` | g0_max² · width · dur_max | maximum-effort normaliser |
| `w` | `objective.effort_weight` (default 0.1) | trades suppression vs effort |
| `p` | `objective.damage_exponent` (default 2.0) | quadratic → penalises peaks more |

### Search Box (Constraints)

| Parameter | Lower bound | Upper bound |
|-----------|------------|-------------|
| `xc` (location, mm) | p_front_obs + w/2 + margin_front | prot_lo − w/2 − margin_protect |
| `ts` (start, s) | t_obs | start_max_s |
| `dur` (duration, s) | duration_s[0] | duration_s[1] |
| `g0` (strength, 1/s) | 0.0 | strength_per_s[1] |

`p_front_obs` is the estimated wave-front position at `t_obs` derived from the reconstruction.

---

## 2. Input Contracts

### Member 3: PDE model config

**File**: `configs/model_default.yaml` (or `--model-config PATH`)

**Fields used by Member 6B** (via `src/control/model.py::load_model()`):

| Field | Used for |
|-------|---------|
| `model.type == "excitable"` | Required; raises error otherwise |
| `model.dimension == 1` | Required; control is 1D only |
| `model.excitable.{D, k, a, eps, gamma}` | PDE parameters |
| `domain.dx` | Stability check for RK4 sub-step |
| `time.dt` | Used for final `apply()` evaluation (accuracy run) |
| `initial_condition.*` | Fallback if `u_truth[0]` absent from reconstruction |

### Member 5: reconstruction file

**File**: `reconstruction.npz` — resolved from pattern list in `configs/control.yaml::inputs.reconstruction`

The resolver checks patterns in order and uses the most recently modified match. Priority:
1. `outputs/inverse/run_*/reconstruction*.npz` (full format, preferred)
2. `outputs/reconstruction*.npz`
3. `outputs/member5_recon.npz`

**Two accepted key layouts** (auto-detected):

| Key | Shape | Required | Comment |
|-----|-------|----------|---------|
| `t` | `(nt,)` float64 | ✅ | Time axis (s) |
| `x` | `(nx,)` float64 | ✅ | Space axis (mm) |
| `u_truth` **or** `truth` | `(nt, nx)` float64 | ✅ | Ground truth (initial state + solver check only) |
| `rec_<method>_<scale>` | `(nt, nx)` float64 | ≥1 | Full format: `outputs/inverse/run_*/reconstruction.npz` |
| `rec_<scale>` | `(nt, nx)` float64 | ≥1 | Compact format: `outputs/member5_recon.npz` |

**How `u_truth` is used**:
- `u_truth[0]`: initial state for the PDE (identical to what Member 3 simulated)
- `u_truth[all]`: solver verification only — relative L2 vs Member 6B's re-integration of the nominal model; logged to `control_metadata.json::solver_check_relL2`
- **Never used for planning decisions**

### Member 5: event estimate file

**File**: `event_estimate.json` — resolved from `configs/control.yaml::inputs.event_estimate`

**Three accepted layouts** (auto-detected by `src/control/events.py`):

```jsonc
// Layout A: per-seed arrivals (inverse/run_*/event_estimate.json)
{ "per_scale": { "1e-08": { "arrival_s_per_seed": [360.1, 361.5, 358.9] } } }

// Layout B: Member 5 actual format
{ "per_scale": { "1e-08": { "B": { "arrival_est_s": 362.0 } } } }

// Layout C: older format
{ "per_scale": { "1e-08": { "arrival_s_at_target": 361.5 } } }
```

Optional top-level keys:

| Key | Used for | Fallback |
|-----|---------|---------|
| `target_x_mm` | Front-arrival measurement point | `event.target_x_mm` in config |
| `truth_arrival_s` | Ground-truth arrival for error calculation | Plant solver arrival time |

---

## 3. Output Contract (for Member 6A)

All outputs written to `outputs/control/run_<YYYYMMDD_HHMMSS>/`.

### CSV files

| File | Description |
|------|-------------|
| `control_results.csv` | Per-scenario table (see column reference below) |
| `control_summary.csv` | Grouped stats for `rec` and `event` sources |
| `timing_robustness.csv` | Oracle plan applied to timing-shifted waves |

**`control_results.csv` column reference:**

| Column | Type | Description |
|--------|------|-------------|
| `scenario` | str | Unique key: `none`, `oracle`, `rec:<method>:<scale>`, `event:<scale>:<seed>` |
| `source` | str | `none` / `oracle` / `rec` / `event` |
| `method` | str | Reconstruction method (from key suffix) or `full_state` |
| `scale` | float | Noise scale (from Member 5) |
| `seed` | int | Seed index for event estimates |
| `status` | str | `ok` / `no_detection` / `infeasible` |
| `J_true` | float | Full-horizon cost on the **plant** (not planning model) |
| `damage_ratio` | float | D(P) / D₀ on plant |
| `effort_norm` | float | E(P) / E_ref on plant |
| `peak_ratio` | float | max(u_prot; controlled) / max(u_prot; uncontrolled) |
| `blocked` | bool | True if `peak_ratio < block_fraction` |
| `xc, ts, dur, g0` | float | Optimal control parameters (`NaN` if status ≠ ok) |
| `J_plan` | float | Cost on **planning** model at optimal P |
| `delta_s` | float | Clock offset (plant arrival − nominal arrival, seconds) |
| `p_front_obs` | float | Estimated front position (mm) at t_obs |
| `regret_vs_oracle` | float | J_true − J_true(oracle) |

### JSON + YAML metadata

| File | Description |
|------|-------------|
| `control_metadata.json` | Full provenance: CLI args, config dict, solver check relL2, input hashes, Python/library versions, git commit |
| `control_used.yaml` | Exact copy of `configs/control.yaml` used for this run |

### Figures (PNG, 160 dpi)

| File | Description |
|------|-------------|
| `fig_control_spacetime.png` | CSD-diverging-colourmap space-time plots for `none`, `oracle`, top-2 `rec`, top `event` |
| `fig_control_cost.png` | J-scatter by scenario with source-coloured points and median lines |
| `fig_search_history.png` | DE + NM-polish convergence (best J vs batch evaluation) for up to 8 scenarios |
| `fig_control_timing_robustness.png` | Oracle plan J vs timing shift; Member 5 error range shaded |

---

## 4. Scenarios

| Scenario | Planner model | Clock offset δ | Information source |
|----------|--------------|---------------|-------------------|
| `none` | — | — | Baseline: uncontrolled plant, J=1 by definition |
| `oracle` | plant | 0 | Upper bound: planner has perfect plant and timing |
| `rec:<m>:<s>` | nominal | Δ from front tracking in `rec_<m>_<s>` | Member 5 reconstruction |
| `event:<s>:<i>` | nominal | Δ from arrival estimate | Member 5 event estimates |

**Information isolation**: the planner never reads `u_truth` for planning. `rec_*` arrays are used
only via `estimate.delta_from_rec()` to extract a timing offset `delta`; then the planner runs
its own nominal model forward from that inferred state.

---

## 5. Speed Notes

### Warm-start optimisation (implemented in `search.py::optimise`)

The uncontrolled phase [0, t_obs − δ] is simulated **once** per scenario (K=1) via `model.warm_cache()`.
All DE population evaluations and NM polish steps start from the cached `(u_warm, v_warm)` state.

```
Time saved = (n_gen + n_polish) × (i_obs / nt) × dt_per_step
           ≈ (10 + 25) × (27%) × 0.5s ≈ 4.7s per scenario
```

### Early-exit on blocking

If the planning objective `J < 1e-4` after DE (wave fully blocked on planning model), Nelder-Mead
polish is **skipped**. This saves ~15s per scenario when the control is very effective.

### Typical runtimes (Intel CPU, single thread, grid 271×301)

| Mode | Approx. time |
|------|-------------|
| `--quick` (popsize=4, maxiter=5, polish=15) | **3–5 min** |
| Default (popsize=4, maxiter=10, polish=25) | **8–15 min** |
| `--no-event` (skip event scenarios) | saves ~40% |

### Tuning levers

```yaml
# configs/control.yaml
search:
  popsize: 4          # K = popsize × 4; more = wider search, slower
  maxiter: 10         # DE generations before polish
  polish_maxiter: 25  # NM steps (skipped when J < 1e-4)
  dt_target_s: 0.5    # search step (stability: actual dt ≈ 0.094s min)
```

---

## 6. Running

```bash
# Default run (auto-detects reconstruction + event estimate)
python src/control/run_control.py

# Fast smoke test
python src/control/run_control.py --quick --no-event

# Override specific files
python src/control/run_control.py \
  --reconstruction outputs/inverse/run_*/reconstruction.npz \
  --event-estimate outputs/inverse/run_*/event_estimate.json

# Specific reconstruction methods only
python src/control/run_control.py --methods B_model B_gauss

# As a module from any directory
python -m src.control.run_control --quick
```

---

## 7. Code Map

```
src/control/
├── model.py       PDE solver (simulate, warm_cache, load_model, initial_state)
├── objective.py   Design container, cost function, search bounds
├── search.py      Differential evolution + Nelder-Mead polish (optimise)
├── estimate.py    Front tracking and clock-offset estimation from reconstruction
├── events.py      Adapter for Member 5 event_estimate.json layouts
└── run_control.py Entry point, scenario orchestration, figures, provenance

configs/control.yaml   Control config (formulation, search params, I/O paths)
docs/member6B_interface_note.md   This document
```

---

*All results labelled HYPOTHETICAL. Model calibrated to preliminary rat CSD data (see
`configs/model_default.yaml`). No clinical claim is made or implied.*
