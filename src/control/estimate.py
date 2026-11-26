"""Turn a reconstruction (or an arrival estimate) into a clock offset delta for the planner.

The planner owns a NOMINAL uncontrolled trajectory (Member 3 model). A reconstruction says where the wave front is, which
fixes how far the real wave is shifted in time from the nominal one:  nominal_time = real_time - delta.
The controller never sees the truth: only the nominal model and the estimate."""
from dataclasses import dataclass
import numpy as np


def peak_pos(U, x):
    i = U.argmax(1); j = np.clip(i, 1, U.shape[1] - 2); r = np.arange(len(i)); a, b, c = U[r, j - 1], U[r, j], U[r, j + 1]
    den = a - 2 * b + c; off = np.clip(np.where(np.abs(den) > 1e-12, 0.5 * (a - c) / den, 0.0), -1, 1)
    return np.where(i == j, x[j] + off * (x[1] - x[0]), x[i])


def arrival_time(pos, t, level):
    k = np.where(pos >= level)[0]
    if not len(k) or k[0] == 0: return float("nan")
    i = k[0]; return float(t[i - 1] + (level - pos[i - 1]) / (pos[i] - pos[i - 1]) * (t[i] - t[i - 1]))


@dataclass
class Nominal:
    t: np.ndarray; x: np.ndarray; pos: np.ndarray; amp: float; tt: np.ndarray; pp: np.ndarray    # tt,pp: settled monotone front curve
    def at(self, t_nom): return float(np.interp(t_nom, self.t, self.pos))
    def time_of(self, p): return float(np.interp(p, self.pp, self.tt))


def make_nominal(u_nom, x, t, tol=0.2):
    pos = peak_pos(u_nom, x); n = len(t); sp = np.gradient(pos, t); ref = np.median(sp[2 * n // 3:])
    bad = np.where(np.abs(sp - ref) > tol * abs(ref))[0]; i0 = int(min(bad[-1] + 1, n - 2)) if len(bad) else 0
    tt, pp = t[i0:], np.maximum.accumulate(pos[i0:]); keep = np.r_[True, np.diff(pp) > 1e-9]
    return Nominal(t, x, pos, float(u_nom[n // 2:].max()), tt[keep], pp[keep])


def delta_from_rec(rec, x, t, nom, t_obs, detect_frac, window, min_points):
    """Fit the front position of the reconstruction over [t_obs-window, t_obs]; delta = t_obs - (nominal time at that position)."""
    sel = (t >= t_obs - window) & (t <= t_obs); det = sel & (rec.max(1) >= detect_frac * nom.amp)
    if det.sum() < min_points: return None
    pos = peak_pos(rec, x); c = np.polyfit(t[det], pos[det], 1); p_hat = float(np.polyval(c, t_obs))
    if not (nom.pp[0] <= p_hat <= nom.pp[-1]): return None
    return dict(delta=t_obs - nom.time_of(p_hat), p_hat=p_hat, n_det=int(det.sum()))


def delta_from_arrival(arr_est, arr_truth_json, a_nom):
    """Convert a retrospective arrival estimate into a nominal-clock offset.

    The estimate is referenced only to the arrival truth carried in the event
    file.  Do not use a plant arrival computed by the control study: that would
    leak unavailable plant information into the estimator.
    """
    return arr_est - a_nom
