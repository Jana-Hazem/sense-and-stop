"""
Sense & Stop - Control Strategy Module (Member 6)
Hypothetical localized intervention g(x, y, t) to suppress CSD wave propagation.
"""
from __future__ import annotations
import copy
import numpy as np
from typing import Dict, Tuple, Optional, Callable


class GaussianControl:
    """
    Spatiotemporal Gaussian intervention control term g(x, y, t).
    Formulation:
      g(x, y, t) = g0 * exp(-((x-x0)^2 + (y-y0)^2) / (2 * sigma_r^2)) * exp(-(t - t_mid)^2 / (2 * sigma_t^2))
    """
    def __init__(self, x0: float, y0: float, t0: float, duration: float, g0: float, sigma_r: float = 2.0):
        self.x0 = x0
        self.y0 = y0
        self.t0 = t0
        self.t_end = t0 + duration
        self.t_mid = t0 + duration / 2.0
        self.g0 = g0
        self.sigma_r = sigma_r
        # sigma_t selected so pulse is active primarily inside [t0, t_end]
        self.sigma_t = max(duration / 4.0, 1e-3)

    def __call__(self, X: np.ndarray, Y: Optional[np.ndarray], t: float) -> np.ndarray:
        # Temporal window gate
        if t < self.t0 or t > self.t_end:
            return np.zeros_like(X)

        t_factor = np.exp(-((t - self.t_mid) ** 2) / (2 * self.sigma_t ** 2))
        
        if Y is None:  # 1D case
            spatial = np.exp(-((X - self.x0) ** 2) / (2 * self.sigma_r ** 2))
        else:          # 2D case
            spatial = np.exp(-(((X - self.x0) ** 2 + (Y - self.y0) ** 2) / (2 * self.sigma_r ** 2)))

        return self.g0 * t_factor * spatial


def evaluate_control_efficacy(baseline_res: dict, controlled_res: dict) -> dict:
    """
    Quantifies the reduction in wave amplitude and spatial coverage after control.
    """
    u_base = baseline_res["u"]
    u_ctrl = controlled_res["u"]
    
    peak_base = float(np.max(u_base))
    peak_ctrl = float(np.max(u_ctrl))
    
    # Measure active area/time where depolarization state u > 0.3
    active_base = int(np.sum(u_base > 0.3))
    active_ctrl = int(np.sum(u_ctrl > 0.3))
    
    suppression_ratio = 1.0 - (active_ctrl / max(active_base, 1))
    amplitude_reduction = 1.0 - (peak_ctrl / max(peak_base, 1e-6))

    return {
        "baseline_peak": peak_base,
        "controlled_peak": peak_ctrl,
        "amplitude_reduction_pct": amplitude_reduction * 100.0,
        "suppression_ratio_pct": suppression_ratio * 100.0,
        "stopped_successfully": bool(suppression_ratio > 0.70)
    }


def bounded_control_search(
    simulate_fn: Callable,
    config: dict,
    target_x: float,
    target_y: float,
    trigger_time: float,
    g0_range: Tuple[float, float] = (0.5, 3.0),
    duration_range: Tuple[float, float] = (10.0, 40.0),
    n_steps: int = 5
) -> Tuple[GaussianControl, dict, dict]:
    """
    Bounded parameter search over intervention strength (g0) and duration (dt)
    to stop wave with minimal energy input.
    """
    baseline_result = simulate_fn(config, control=None)
    best_ctrl = None
    best_metrics = {"suppression_ratio_pct": -1.0}
    best_result = None

    g0_vals = np.linspace(g0_range[0], g0_range[1], n_steps)
    dur_vals = np.linspace(duration_range[0], duration_range[1], n_steps)

    for g0 in g0_vals:
        for dur in dur_vals:
            ctrl = GaussianControl(x0=target_x, y0=target_y, t0=trigger_time, duration=dur, g0=g0)
            res = simulate_fn(config, control=ctrl)
            metrics = evaluate_control_efficacy(baseline_result, res)

            if metrics["suppression_ratio_pct"] > best_metrics["suppression_ratio_pct"]:
                best_metrics = metrics
                best_ctrl = ctrl
                best_result = res

    return best_ctrl, best_metrics, best_result