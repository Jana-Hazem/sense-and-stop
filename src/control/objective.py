"""Control objective, design container and search bounds."""
from dataclasses import dataclass
import numpy as np
from .model import Ctrl


@dataclass
class Design:
    width: float; edge: float; ramp: float
    p: float; effort_weight: float; e_ref: float; block_fraction: float
    margin_front: float; margin_protect: float
    t_obs: float; start_max: float; dur: tuple; g: tuple
    prot_lo: float; prot_hi: float


def design_from_cfg(cfg, T):
    iv, b, ob = cfg["intervention"], cfg["bounds"], cfg["objective"]; pr = cfg["protected_region"]
    dur = tuple(float(v) for v in b["duration_s"]); g = tuple(float(v) for v in b["strength_per_s"])
    sm = float(b["start_max_s"]) if b.get("start_max_s") is not None else float(T - dur[0])
    if dur[0] < 2 * float(iv["ramp_s"]): raise ValueError(f"min duration {dur[0]} s must be >= 2*ramp_s = {2 * float(iv['ramp_s'])} s (the gate would never reach full strength)")
    e_ref = g[1] ** 2 * float(iv["width_mm"]) * dur[1]                 # effort of the strongest, longest, full-width control
    return Design(float(iv["width_mm"]), float(iv["edge_mm"]), float(iv["ramp_s"]), float(ob["damage_exponent"]), float(ob["effort_weight"]), e_ref,
                  float(ob["block_fraction"]), float(iv["margin_front_mm"]), float(iv["margin_protect_mm"]), float(b["t_obs_s"]), sm, dur, g,
                  float(pr["x_min"]), float(pr["x_max"]))


def make_ctrl(P, d):
    """P: (4, K) rows = xc, ts, dur, g0."""
    P = np.asarray(P, float); P = P[:, None] if P.ndim == 1 else P
    return Ctrl(P[0], P[1], P[2], P[3], d.width, d.edge, d.ramp)


def cost(out, D0, d):
    """J = damage ratio + effort_weight * normalised effort. Returns (J, damage_ratio, effort_norm)."""
    dr = out["damage"] / D0; en = out["effort"] / d.e_ref
    return dr + d.effort_weight * en, dr, en


def search_bounds(d, p_front_obs):
    """Box for (xc, ts, dur, g0). The zone must lie downstream of the front at t_obs and upstream of the protected region.
    Returns None if no admissible location exists."""
    lo = p_front_obs + d.width / 2 + d.margin_front; hi = d.prot_lo - d.width / 2 - d.margin_protect
    if hi <= lo: return None
    return [(lo, hi), (d.t_obs, max(d.start_max, d.t_obs + 1e-6)), d.dur, d.g]
