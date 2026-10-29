"""
Sense & Stop - calibration of the excitable model to empirical pulse features (Member 3)

Why this exists
---------------
Member 2's data give the temporal width (FWHM) of the CSD pulse at one electrode and the
literature gives a front speed. The excitable model u_t = D lap(u) + k[u(1-u)(u-a) - v],
v_t = k eps (u - gamma v) has a simple scaling property:

    * time is measured in units of 1/k      ->  temporal FWHM  F = T0 / k
    * length is measured in sqrt(D/k)       ->  speed          c = c0 * sqrt(D k)

where T0 and c0 are the dimensionless FWHM and speed of the pulse for the shape parameters
(a, eps, gamma).  So: measure (T0, c0) once with D = k = 1, then solve for k and D.

Functions
---------
pulse_metrics(result, x_probe)   -> fwhm / rise / decay at one location (same definitions as Member 2)
pulse_speed(result, x1, x2)      -> speed from peak arrival times at two locations
calibrate_excitable(config, fwhm_s, speed_mm_per_min) -> (new_config, info)
"""
from __future__ import annotations

import copy

import numpy as np

from .solver import simulate

__all__ = ["pulse_metrics", "pulse_speed", "calibrate_excitable"]


def _cross_time(t, s, level, i_peak, direction):
    """Time at which s crosses `level` going away from the peak (direction -1 = before, +1 = after).
    Linear interpolation. Returns nan if no crossing inside the record."""
    i = i_peak
    while 0 <= i + direction < len(s):
        j = i + direction
        if s[j] <= level:
            frac = (s[i] - level) / (s[i] - s[j])
            return t[i] + frac * (t[j] - t[i])
        i = j
    return np.nan


def _peak_time(t, s):
    """Peak time with parabolic refinement (finer than the output spacing)."""
    i = int(np.argmax(s))
    if 0 < i < len(s) - 1:
        y0, y1, y2 = s[i - 1], s[i], s[i + 1]
        den = y0 - 2 * y1 + y2
        if den != 0:
            return t[i] + 0.5 * (y0 - y2) / den * (t[1] - t[0])
    return t[i]


def pulse_metrics(result: dict, x_probe: float, min_peak: float = 0.3):
    """Temporal features of the pulse at x_probe (mm), using Member 2's definitions:
       fwhm_s        full width at half of the peak height
       rise_10pct_s  time from the 10 % level (before the peak) to the peak
       decay_10pct_s time from the peak to the 10 % level (after the peak)
    The baseline is the value at t = 0 (rest). Works for 1D (u[t,x]) and 2D (mid-line in y).
    Returns None if no pulse of height >= min_peak passes the probe."""
    u, x, t = result["u"], result["x"], result["t"]
    if u.ndim == 3:
        u = u[:, u.shape[1] // 2, :]
    s = u[:, int(np.argmin(np.abs(x - x_probe)))]
    s = s - s[0]
    pk = s.max()
    if pk < min_peak:
        return None
    ip = int(np.argmax(s))
    tp = _peak_time(t, s)
    t_half_a = _cross_time(t, s, 0.5 * pk, ip, -1)
    t_half_b = _cross_time(t, s, 0.5 * pk, ip, +1)
    t_10_a = _cross_time(t, s, 0.1 * pk, ip, -1)
    t_10_b = _cross_time(t, s, 0.1 * pk, ip, +1)
    return {"peak": float(pk), "t_peak_s": float(tp),
            "fwhm_s": float(t_half_b - t_half_a),
            "rise_10pct_s": float(tp - t_10_a),
            "decay_10pct_s": float(t_10_b - tp)}


def pulse_speed(result: dict, x1: float, x2: float) -> float:
    """Speed (mm/s) from peak arrival times at x1 < x2 (use points away from the stimulus
    and the far boundary)."""
    m1, m2 = pulse_metrics(result, x1), pulse_metrics(result, x2)
    if m1 is None or m2 is None:
        raise ValueError("no pulse at one of the probe locations")
    return (x2 - x1) / (m2["t_peak_s"] - m1["t_peak_s"])


def calibrate_excitable(config: dict, fwhm_s: float, speed_mm_per_min: float):
    """Return (new_config, info) with D and k set so that, for the config's shape parameters
    (a, eps, gamma), the pulse has temporal FWHM `fwhm_s` at a point and front speed
    `speed_mm_per_min`.  Implied spatial width = speed * fwhm.

    Raises RuntimeError if no sustained pulse exists for those shape parameters."""
    cfg = copy.deepcopy(config)
    p = cfg["model"]["excitable"]

    # dimensionless reference run (D = k = 1); lengths in sqrt(D/k), time in 1/k
    ref = copy.deepcopy(cfg)
    ref["model"].update(type="excitable", dimension=1)
    ref["model"]["excitable"].update(D=1.0, k=1.0)
    ref["domain"].update(Lx=250.0, dx=0.5, boundary="neumann")
    ref["time"].update(T=420.0, output_dt=0.5, dt=0.05)
    ref["initial_condition"].update(type="edge_step", x0=25.0, width=2.0, amplitude=1.0)
    ref["solver"]["method"] = "rk4"
    r = simulate(ref)
    m1, m2 = pulse_metrics(r, 150.0), pulse_metrics(r, 200.0)
    if m1 is None or m2 is None:
        raise RuntimeError("no sustained pulse for these (a, eps, gamma); change the shape parameters")
    T0 = m2["fwhm_s"]                                   # dimensionless FWHM
    c0 = 50.0 / (m2["t_peak_s"] - m1["t_peak_s"])       # dimensionless speed

    c = speed_mm_per_min / 60.0                         # mm/s
    k = T0 / fwhm_s                                     # 1/s
    D = (c / c0) ** 2 / k                               # mm^2/s   (c = c0 sqrt(D k))
    p["D"], p["k"] = float(D), float(k)
    info = {"T0_dimensionless": T0, "c0_dimensionless": c0,
            "D_mm2_per_s": D, "k_per_s": k,
            "length_scale_mm": float(np.sqrt(D / k)),
            "implied_pulse_width_mm": c * fwhm_s,
            "suggested_max_dx_mm": float(0.5 * np.sqrt(D / k)),
            "pulse_shape_dimensionless": {kk: m2[kk] for kk in ("fwhm_s", "rise_10pct_s", "decay_10pct_s", "peak")}}
    return cfg, info
