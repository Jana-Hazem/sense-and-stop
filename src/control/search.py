"""Bounded, ROBUST search over intervention location, start time, duration and strength.

Stage 1: differential evolution (global, population simulated as ONE batch).  Stage 2: Nelder-Mead polish.

Why robust?  With a plain objective the optimiser pushes the strength g0 down until the wave is JUST blocked
(knife edge): a 2 s timing error then makes the plan fail.  Here every candidate is simulated for an ENSEMBLE of cases
(timing offsets around the estimate, optional model factors) and the cost uses the WORST damage ratio of the ensemble:
    J = max_over_cases(damage/D0) + w * mean_over_cases(effort/e_ref)
so a plan only looks good if it blocks the wave in all cases."""
import numpy as np
from scipy.optimize import differential_evolution, minimize
from .model import simulate
from .objective import make_ctrl


def optimise(cases, grid, u0, prot, d, bounds, s, seed_offset=0):
    """cases: list of dict(model=Model, D0=float, offset=float) = planner's belief set (offset = clock offset in s)."""
    groups = {}
    for c in cases: groups.setdefault(id(c["model"]), dict(model=c["model"], D0=c["D0"], off=[]))["off"].append(float(c["offset"]))
    hist, best = [], [np.inf]

    def J_batch(P):
        P = np.asarray(P, float); P = P[:, None] if P.ndim == 1 else P; K = P.shape[1]; drs, ens = [], []
        for g in groups.values():
            M = len(g["off"]); out = simulate(g["model"], grid, u0, prot, make_ctrl(np.tile(P, (1, M)), d), clock_offset=np.repeat(g["off"], K),
                                              p=d.p, dt_target=s["dt_target_s"])
            drs.append((out["damage"] / g["D0"]).reshape(M, K)); ens.append((out["effort"] / d.e_ref).reshape(M, K))
        dr = np.vstack(drs).max(0); en = np.vstack(ens).mean(0); J = dr + d.effort_weight * en
        best[0] = min(best[0], float(J.min())); hist.append(best[0]); return J

    r = differential_evolution(J_batch, bounds, seed=int(s["seed"]) + seed_offset, popsize=int(s["popsize"]), maxiter=int(s["maxiter"]), tol=float(s["tol"]),
                                polish=False, vectorized=True, updating="deferred", init="sobol")
    x, J = r.x, float(r.fun)
    if int(s.get("polish_maxiter", 0)) > 0:
        lo, hi = np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds])
        span = hi - lo
        # Polish scaled to each parameter range (especially g0 vs ts/dur) via [0, 1] normalization
        def obj_norm(w):
            z = lo + np.clip(w, 0.0, 1.0) * span
            return float(J_batch(z)[0])
        w0 = np.clip((x - lo) / np.maximum(span, 1e-12), 0.0, 1.0)
        nm = minimize(obj_norm, w0, method="Nelder-Mead", options=dict(maxiter=int(s["polish_maxiter"]), xatol=1e-3, fatol=1e-6))
        x_polished = lo + np.clip(nm.x, 0.0, 1.0) * span
        if nm.fun < J: x, J = x_polished, float(nm.fun)
    return dict(P=x, J_plan=J, history=hist, nfev_batches=len(hist))
