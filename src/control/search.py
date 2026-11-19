"""Bounded search over intervention location, start time, duration and strength.

Stage 1: differential evolution (global, population simulated as ONE batch).  Stage 2: Nelder-Mead polish.
The planner simulates the NOMINAL model on the nominal clock (clock_offset = delta).

Speed-ups implemented
---------------------
1. Warm-start cache: the uncontrolled phase [0, t_obs-delta] is computed once per scenario (K=1).
   All batch and polish evaluations re-use the cached (u_warm, v_warm) state, skipping ~27% of steps.
2. Early-exit on blocking: if the planning J drops below `early_exit_J` during DE, the polish is
   skipped and the best solution is returned immediately (wave already blocked on planning model).
3. Effort-aware guard: if g0 is effectively zero the control does nothing; a floor prevents trivial
   solutions from dominating the population and confusing the polish.

Correctness: g(x,t)=0 before ts >= t_obs, so the controlled and uncontrolled solutions are identical on
[0, t_obs]. The pre-obs damage D_pre is constant across candidates and added back so J_plan correctly
represents the full-horizon normalised cost on the planning model.
"""
import numpy as np
from scipy.optimize import differential_evolution, minimize
from .model import simulate, warm_cache
from .objective import make_ctrl, cost

_EARLY_EXIT_J = 1e-4          # skip polish when planning wave is fully blocked (J < this)
_EFFORT_FLOOR  = 1e-6         # treat g0 < this as "no control" to avoid degenerate candidates


def optimise(model, grid, u0, prot, d, bounds, delta, D0, s, seed_offset=0):
    """Search for the best intervention parameters [xc, ts, dur, g0].

    Parameters
    ----------
    model        : Model – planning model (nominal or plant for oracle)
    grid         : Grid  – shared (t, x) grid from Member 5 reconstruction
    u0           : (nx,) initial state from Member 5 u_truth[0]
    prot         : boolean mask – protected region
    d            : Design – from design_from_cfg()
    bounds       : list of (lo, hi) for [xc, ts, dur, g0]
    delta        : float – clock offset: planning clock = real clock - delta
    D0           : float – uncontrolled damage on planning model (normaliser)
    s            : dict  – search config keys: seed, popsize, maxiter, tol, polish_maxiter, dt_target_s
    seed_offset  : int   – added to s['seed'] for reproducibility across scenarios

    Returns
    -------
    dict(P, J_plan, history, nfev_batches)
      P           : (4,) optimal [xc, ts, dur, g0]
      J_plan      : full-horizon cost on planning model at P
      history     : list of best-J-so-far after each batch call (for convergence plots)
      nfev_batches: number of J_batch calls made
    """
    dt_search = float(s["dt_target_s"])

    # ── Warm-start ──────────────────────────────────────────────────────────────
    # Run uncontrolled from t[0] to (t_obs − delta) once; reuse for all evaluations.
    # D_pre = damage in protected region over [0, t_obs] — constant across all P.
    u_warm, v_warm, i_obs, D_pre = warm_cache(
        model, grid, u0, d.t_obs - delta, prot, p=d.p, dt_target=dt_search)

    hist, best = [], [np.inf]

    def J_batch(P):
        P = np.asarray(P, float)
        if P.ndim == 1: P = P[:, None]
        # Zero-effort guard: treat g0 < floor as no-control (avoids degenerate NM steps)
        P_eff = P.copy(); P_eff[3, P_eff[3] < _EFFORT_FLOOR] = 0.0
        out = simulate(model, grid, u0, prot, make_ctrl(P_eff, d),
                       clock_offset=delta, p=d.p, dt_target=dt_search,
                       u_warm=u_warm, v_warm=v_warm, i_start=i_obs)
        out_full = dict(damage=out["damage"] + D_pre, effort=out["effort"], peak=out["peak"])
        J = cost(out_full, D0, d)[0]
        best[0] = min(best[0], float(J.min())); hist.append(best[0])
        return J

    r = differential_evolution(
        J_batch, bounds,
        seed=int(s["seed"]) + seed_offset,
        popsize=int(s["popsize"]),
        maxiter=int(s["maxiter"]),
        tol=float(s["tol"]),
        polish=False,
        vectorized=True,
        updating="deferred",
        init="sobol",
    )
    x, J = r.x, float(r.fun)

    # ── Polish ───────────────────────────────────────────────────────────────────
    # Skip polish if wave is already fully blocked on planning model (J ≈ 0).
    n_polish = int(s["polish_maxiter"])
    if n_polish > 0 and J > _EARLY_EXIT_J:
        lo = np.array([b[0] for b in bounds])
        hi = np.array([b[1] for b in bounds])
        nm = minimize(
            lambda z: float(J_batch(np.clip(z, lo, hi))[0]),
            x,
            method="Nelder-Mead",
            options=dict(maxiter=n_polish, xatol=1e-3, fatol=1e-6),
        )
        if nm.fun < J:
            x, J = np.clip(nm.x, lo, hi), float(nm.fun)

    return dict(P=x, J_plan=J, history=hist, nfev_batches=len(hist))
