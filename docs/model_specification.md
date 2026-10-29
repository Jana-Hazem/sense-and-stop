# Model Specification: Reduced PDE Models of a CSD-Inspired Wave

**Owner:** Member 3 · **Status:** v0.2 (calibrated to S1 preliminary data) · **Code:** `src/model/` · **Config:** `configs/model_default.yaml`

## 0. Changes from v0.1

- Excitable model re-calibrated to Member 2's S1 measurements (pulse FWHM 42.8 s) and the literature speed (2.6 mm/min). v0.1 values were illustrative and gave a pulse about 4× too long and too wide.
- Excitable model is now the **default (primary)**; Fisher–KPP is kept as the front baseline.
- Grid refined from dx = 0.25 mm to **0.1 mm**, because the calibrated front is much sharper.
- New tool `calibrate_excitable()` computes D and k from target FWHM and speed, and `pulse_metrics()` measures the same features as Member 2.
- Tests increased from 16 to 21 and now include agreement with the S1 features.

## 1. Scope statement

The models are **reduced, phenomenological reaction–diffusion equations**, not ionic CSD simulators: no extracellular K⁺, Na⁺/Cl⁻/Ca²⁺ dynamics, pumps, glial buffering, cell swelling, oxygen, perfusion, or neurovascular coupling. Their purpose is to provide a hidden field u(x,t) with realistic *timing and scale* for testing an inverse-reconstruction method. The S1 data are **rat** recordings; nothing here is human data, and nothing here establishes anything about human migraine.

## 2. State variables (units: mm, s)

- **u(x,t)**: dimensionless excitation/depolarisation state. 0 = rest. Peak ≈ 0.81 (excitable) or → 1 (Fisher–KPP). Normalised, so the unknown amplifier gain of S1 does not matter.
- **v(x,t)**: dimensionless recovery variable (excitable model only); phenomenological, not a measured quantity.

## 3. Equations

**Excitable (FitzHugh–Nagumo type), primary:**

    u_t = D ∇²u + k [ u(1−u)(u−a) − v ]
    v_t = k ε ( u − γ v )

Role: traveling pulse with recovery, as recommended in the literature review (Secs. 2.2.2, 2.4).

**Fisher–KPP, baseline:**

    u_t = D ∇²u + r u(1−u)

Role: front-reconstruction baseline only. It has no recovery, no threshold and no refractoriness. Asymptotic speed c* = 2√(Dr), approached from below.

**Optional control hook (Member 6):** `simulate(cfg, control=g)` subtracts g(x,t)·u from u_t. Default `None`, so the primary path is unaffected.

## 4. Domain, boundary and initial conditions

| Item | Choice | Status |
|---|---|---|
| Domain | 1D: x ∈ [0, 30] mm. 2D: [0, 30] × [0, 20] mm rectangle | **Design choice, to be agreed with Member 4.** Rat data cannot set it |
| Boundary | No-flux (Neumann), second-order mirror ghost nodes | Design choice |
| Initial condition | `edge_step`: u = ½(1 − tanh((x − x₀)/w)), x₀ = 3 mm, w = 0.5 mm; v = 0 | Launches a wave in +x |
| Other IC types | `gaussian` (point stimulus); `kpp_exact_wave` (verification only) | |

## 5. Numerical method

- **Space:** second-order central finite differences, uniform grid, nodes on both boundaries. **dx = 0.1 mm** (must be ≤ 0.5·√(D/k) = 0.096 mm to resolve the front; convergence is tested).
- **Time:** RK4 (default), forward Euler, or `scipy.integrate.solve_ivp` (independent reference for avoiding the inverse crime). Internal dt = 0.1 s.
- **Stability guard:** with Λ = D·dt·(1/dx² + 1/dy²), require Λ ≤ 0.5 (Euler) or ≤ 0.69 (RK4), times a safety factor of 0.9; dt is reduced automatically. Default Λ = 0.19 (1D) and 0.37 (2D).
- **Output frames** are at exact multiples of `output_dt` = 2 s (about 21 frames per FWHM).

## 6. Parameters and provenance

| Parameter | Value | Unit | Provenance |
|---|---|---|---|
| Target temporal FWHM | 42.8 (animal means 36–56) | s | **S1-measured**, rat, 3 animals, 15 waves, PRELIMINARY |
| Target speed | 2.6 (range 2.3–2.9) | mm/min | **Literature-derived**, not yet verified by reading the papers |
| Implied pulse width | 1.85 (range 1.7–2.1) | mm | Derived: speed × FWHM |
| Excitable D | 0.018561 | mm²/s | Calibrated |
| Excitable k | 0.503758 | 1/s | Calibrated (k = T₀/FWHM, T₀ = 21.56) |
| Excitable a, ε, γ | 0.1, 0.01, 2.0 | – | Design choice (shape of the pulse) |
| Fisher–KPP D, r | 0.01296, 0.036 | mm²/s, 1/s | Hand-calibrated: c* = 0.0432 mm/s, front length √(D/r) = 0.6 mm |
| Lx, Ly | 30, 20 | mm | Design choice |

**Calibration method.** For the excitable model, time scales as 1/k and length as √(D/k). A dimensionless run (D = k = 1) gives the pulse's FWHM T₀ and speed c₀; then k = T₀/FWHM and D = (c/c₀)²/k. Call `calibrate_excitable(cfg, fwhm_s, speed_mm_per_min)`; do not edit D and k by hand.

## 7. Model versus S1 data

| Feature | S1 (animal-mean range) | Model | Note |
|---|---|---|---|
| FWHM at a point | 42.8 s (36–56) | 42.8 s | Matches **by construction** (it is the calibration target) |
| Speed | 2.3–2.9 mm/min (literature) | 2.60 mm/min | By construction |
| Pulse width | 1.7–2.1 mm (derived) | 1.85 mm | By construction |
| Rise 10 % → peak | 35.2 s (25–46) | 36.3 s | **Independent check, agrees** |
| Peak → 10 % | 40.7 s (28–47) | 28.2 s | **Independent check, at the low edge**: model recovers faster than the data |

Only the last two rows are independent checks. With three animals (R6 differs from R4) this is a plausibility check, not validation.

## 8. Sensitivity scenarios (literature broad range, FWHM kept at 42.8 s)

Generated with `calibrate_excitable`. For each, set `dx` ≤ the value shown and adjust `T` to the crossing time.

| Speed | Pulse width | D (mm²/s) | max dx | Time to cross 27 mm |
|---|---|---|---|---|
| 2.0 mm/min | 1.43 mm | 0.0110 | 0.074 mm | 810 s |
| **2.6 mm/min (baseline)** | **1.85 mm** | **0.0186** | **0.096 mm** | **623 s** |
| 9.0 mm/min | 6.42 mm | 0.2224 | 0.332 mm | 180 s |

## 9. Verification (`tests/test_model.py`, 21 tests, all passing)

- **Exact solution:** Fisher–KPP traveling wave (Ablowitz–Zeppetella) reproduced; error falls by exactly 4.00× when dx is halved (second order). Max error 4.9×10⁻⁴ at dx = 0.1 mm.
- **Convergence:** second order in space (RK4); first order in time (Euler).
- **Independent solver:** RK4 and `solve_ivp` agree (max difference < 5×10⁻³).
- **Boundary:** pure diffusion with no-flux walls conserves mass to 10⁻¹⁰.
- **Stability guard** works; 2D planar wave equals 1D to 10⁻¹⁰.
- **Excitable model:** pulse recovers behind itself; FWHM within 3 % of 42.8 s, speed in 2.3–2.9 mm/min, width 1.6–2.2 mm; halving dx changes FWHM and speed by < 2 %.
- **Measurement functions:** exact check on a Gaussian pulse.
- **Calibration:** a different target (50 s, 4 mm/min) is reproduced within 5 %; it raises an error if no pulse exists.
- **Control hook:** g = 0 changes nothing; a strong g blocks the wave.
- **Known mismatch recorded as a test:** model decay10 is shorter than the S1 mean.

## 10. Known limitations (for the later limitations phase)

1. **Data basis:** rat, 3 animals, preliminary; R6 differs from R4. Speed is literature-derived and unverified, since the data have one electrode.
2. **Decay too fast** (28 s vs 41 s). Pulse shape is only partly matched.
3. **Recovery undershoot:** u falls to about −0.2 to −0.3 behind the pulse and relaxes over a couple of minutes. S1 does not constrain this; it is a typical FitzHugh–Nagumo feature, not a physiological prediction.
4. **Sensitivity:** the pulse is fragile in (a, ε, γ): ε ≥ 0.015 or a ≥ 0.15 made the wave die in scans.
5. **Phenomenological:** no ionic, glial, metabolic or vascular mechanisms; flat rectangular domain with no folding or anisotropy.
6. **Domain size is a design choice**, not data.
7. Fisher–KPP has no recovery and is a front-only baseline.

## 11. Open decisions (for Member 1's decision log)

1. Primary model: excitable (recommended by the literature review, now default) versus Fisher–KPP. Both are available.
2. Domain size and 1D/2D, to be agreed with Member 4 (see the interface note, Sec. 5).
3. Final targets once Member 2 verifies the two speed papers and adds animals: re-run `calibrate_excitable` and the tests.
