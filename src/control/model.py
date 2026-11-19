"""Controlled excitable (FitzHugh-Nagumo type) PDE, batch solver.

    u_t = D u_xx + k*( u(1-u)(u-a) - v ) - g(x,t) u        (g = 0: identical to Member 3's model)
    v_t = k*eps*( u - gamma v )
    g(x,t) = g0 * zone(x; xc, width) * gate(t; ts, dur)    bounded, localised, HYPOTHETICAL intervention

Control term specification
--------------------------
The hypothetical control g(x,t) is a smooth, bounded, spatially localised removal term:

    g(x,t) = g0 * zone(x; xc, w) * gate(t; ts, dur, ramp)

where
    zone(x; xc, w)   = tanh-smoothed top-hat of width w centred at xc          [dimensionless, in [0,1]]
    gate(t; ts, dur) = smoothstep ramp-on at ts, ramp-off at ts+dur             [dimensionless, in [0,1]]
    g0               = maximum removal rate  (units: 1/s, positive)
    xc, ts, dur      = intervention centre (mm), start time (s), duration (s)   [all bounded]

Biophysical interpretation (HYPOTHETICAL): g0 u represents a local, time-limited quench of wave activity,
e.g. via cooling, optogenetic inhibition or pharmaceutical suppression. The term is linear in u so it cannot
drive u below zero; for the excitable model this is consistent with physiological state bounds.

The PDE with control is identical to Member 3's model when g0 = 0 or outside the zone/gate support.

Parameters D,k,a,eps,gamma come from the Member 3 config; grid, time axis and initial state come from Member 5's
reconstruction.npz (t, x, u_truth[0]) so nothing is typed in.  K candidate controls are simulated at once (batch axis).

Speed-up: warm_simulate() runs the solver in two phases:
  Phase 1 (K=1):   t[0] → t_obs  (no control, shared across all candidates)
  Phase 2 (K=pop): t_obs → t[-1] (each candidate diverges here)
Call warm_cache() once per scenario to get (u_warm, v_warm), then pass them to warm_simulate().
"""
from dataclasses import dataclass
from pathlib import Path
import numpy as np, yaml


@dataclass
class Model:
    D: float; k: float; a: float; eps: float; gamma: float
    def scaled(self, d=1.0, k=1.0): return Model(self.D * d, self.k * k, self.a, self.eps, self.gamma)


@dataclass
class Grid:
    x: np.ndarray; t: np.ndarray
    @property
    def dx(self): return float(self.x[1] - self.x[0])
    @property
    def out_dt(self): return float(self.t[1] - self.t[0])


@dataclass
class Ctrl:
    xc: np.ndarray; ts: np.ndarray; dur: np.ndarray; g0: np.ndarray      # (K,) each
    width: float; edge: float; ramp: float


def load_model(path):
    cfg = yaml.safe_load(Path(path).read_text()); m = cfg["model"]
    if m.get("type") != "excitable": raise ValueError("control module expects model.type == excitable")
    if int(m.get("dimension", 1)) != 1: raise ValueError("control module is 1D")
    e = m["excitable"]
    return Model(**{q: float(e[q]) for q in ("D", "k", "a", "eps", "gamma")}), cfg


def initial_state(cfg, grid, u_truth0=None):
    """u0 = first row of Member 5's u_truth (exactly Member 3's initial state); else rebuilt from the config."""
    if u_truth0 is not None: return np.asarray(u_truth0, float), "reconstruction.npz:u_truth[0]"
    ic = cfg["initial_condition"]; x0, w = float(ic["x0"]), float(ic["width"]); amp = float(ic.get("amplitude", 1.0))
    return amp * 0.5 * (1 - np.tanh((grid.x - x0) / w)), "config edge_step (tanh) - verify against Member 3"


def smoothstep(z): z = np.clip(z, 0.0, 1.0); return z * z * (3 - 2 * z)
def zone(x, xc, w, edge): return 0.5 * (np.tanh((x[None, :] - (xc[:, None] - w / 2)) / edge) - np.tanh((x[None, :] - (xc[:, None] + w / 2)) / edge))
def gate(t, ts, dur, ramp): return smoothstep((t - ts) / ramp) * smoothstep((ts + dur - t) / ramp)


def _nsub(od, dt_target, dx, D):
    """Number of RK4 sub-steps per output step (stability-aware)."""
    ns = max(int(np.ceil(od / dt_target)), 1)
    while od / ns > 0.5 * 2.78 * dx * dx / (4 * D): ns += 1  # RK4 diffusion stability
    return ns


def _rk4_step(u, v, f, dt, tr, clock):
    """One RK4 step."""
    k1u, k1v = f(u, v, tr + clock)
    k2u, k2v = f(u + .5 * dt * k1u, v + .5 * dt * k1v, tr + .5 * dt + clock)
    k3u, k3v = f(u + .5 * dt * k2u, v + .5 * dt * k2v, tr + .5 * dt + clock)
    k4u, k4v = f(u + dt * k3u, v + dt * k3v, tr + dt + clock)
    return u + dt / 6 * (k1u + 2 * k2u + 2 * k3u + k4u), v + dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)


def simulate(model, grid, u0, prot, ctrl=None, clock_offset=0.0, p=2.0, store=False, dt_target=0.1,
             u_warm=None, v_warm=None, i_start=0):
    """Returns dict(damage, effort, peak, [u]) with batch axis K.

    prot: boolean mask over x (protected region). clock_offset: real_time = simulation_time + clock_offset.
    damage = integral over time and protected region of max(u,0)^p ; effort = integral of g^2 over x and t.

    Warm-start interface: pass u_warm (K, nx) and v_warm (K, nx) together with i_start (the index in grid.t
    where the warm state was captured). Damage and effort accumulation starts from i_start; steps 0..i_start-1
    are not re-computed. This avoids duplicating the shared uncontrolled phase across many candidates.
    """
    x, t, dx, od = grid.x, grid.t, grid.dx, grid.out_dt
    K = 1 if ctrl is None else len(ctrl.xc); nt, nx = len(t), len(x)
    nsub = _nsub(od, dt_target, dx, model.D)
    dt = od / nsub
    # Initialise state
    if u_warm is not None:
        u = np.asarray(u_warm, float); v = np.asarray(v_warm, float)
        if u.ndim == 1: u = np.tile(u, (K, 1))
        if v.ndim == 1: v = np.tile(v, (K, 1))
    else:
        u = np.tile(u0, (K, 1)); v = np.zeros_like(u)
        i_start = 0
    GS = None if ctrl is None else ctrl.g0[:, None] * zone(x, ctrl.xc, ctrl.width, ctrl.edge)
    D, k, a, eps, gam, dx2 = model.D, model.k, model.a, model.eps, model.gamma, dx * dx
    gt = (lambda tr: GS * gate(tr, ctrl.ts, ctrl.dur, ctrl.ramp)[:, None]) if GS is not None else (lambda tr: 0.0)

    def f(u, v, tr):
        lap = np.empty_like(u); lap[:, 1:-1] = u[:, 2:] - 2 * u[:, 1:-1] + u[:, :-2]
        lap[:, 0] = 2 * (u[:, 1] - u[:, 0]); lap[:, -1] = 2 * (u[:, -2] - u[:, -1]); lap /= dx2     # no-flux
        return D * lap + k * (u * (1 - u) * (u - a) - v) - gt(tr) * u, k * eps * (u - gam * v)

    dmg, eff, peak = np.zeros(K), np.zeros(K), np.full(K, -np.inf)
    U = np.empty((K, nt, nx)) if store else None
    for n in range(i_start, nt):
        wt = 0.5 if n in (i_start, nt - 1) else 1.0; tr0 = t[n] + clock_offset
        dmg += wt * od * dx * (np.maximum(u[:, prot], 0.0) ** p).sum(1); peak = np.maximum(peak, u[:, prot].max(1))
        if GS is not None: eff += wt * od * dx * (gt(tr0) ** 2).sum(1)
        if store: U[:, n] = u
        if n == nt - 1: break
        tn = t[n]
        for _ in range(nsub):
            u, v = _rk4_step(u, v, f, dt, tn, clock_offset)
            tn += dt
    out = dict(damage=dmg, effort=eff, peak=peak)
    if store:
        if i_start > 0: out["u"] = U[:, i_start:]   # only the post-warm portion is stored
        else: out["u"] = U
    return out


def warm_cache(model, grid, u0, t_obs, prot, p=2.0, dt_target=0.1):
    """Run the uncontrolled (g=0) solver from t[0] to t_obs.
    Returns (u_warm, v_warm, i_obs, D_pre) where:
      u_warm, v_warm : state arrays (nx,) at t[i_obs]
      i_obs          : index in grid.t just at or after t_obs
      D_pre          : damage integral over [0, t[i_obs]] in protected region
    This state can be reused across all candidates.
    """
    t = grid.t; x = grid.x; dx = grid.dx; od = grid.out_dt
    i_obs = int(np.searchsorted(t, t_obs, side="right")) - 1
    i_obs = max(0, min(i_obs, len(t) - 2))
    nsub = _nsub(od, dt_target, dx, model.D)
    dt = od / nsub
    D, k, a, eps, gam, dx2 = model.D, model.k, model.a, model.eps, model.gamma, dx * dx

    def f(u, v, _tr):
        lap = np.empty_like(u); lap[:, 1:-1] = u[:, 2:] - 2 * u[:, 1:-1] + u[:, :-2]
        lap[:, 0] = 2 * (u[:, 1] - u[:, 0]); lap[:, -1] = 2 * (u[:, -2] - u[:, -1]); lap /= dx2
        return D * lap + k * (u * (1 - u) * (u - a) - v), k * eps * (u - gam * v)

    u = u0[None, :].copy(); v = np.zeros_like(u)
    D_pre = 0.0
    for n in range(i_obs):
        wt = 0.5 if n in (0, i_obs - 1) else 1.0
        D_pre += wt * od * dx * (np.maximum(u[0, prot], 0.0) ** p).sum()
        tn = t[n]
        for _ in range(nsub):
            u, v = _rk4_step(u, v, f, dt, tn, 0.0)
            tn += dt
    return u[0], v[0], i_obs, D_pre          # (nx,), (nx,), int, float
