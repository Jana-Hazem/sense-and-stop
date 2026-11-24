"""Controlled excitable (FitzHugh-Nagumo type) PDE, batch solver.

    u_t = D u_xx + k*( u(1-u)(u-a) - v ) - g(x,t) u        (g = 0: identical to Member 3's model)
    v_t = k*eps*( u - gamma v )
    g(x,t) = g0 * zone(x; xc, width) * gate(t; ts, dur)    bounded, localised, HYPOTHETICAL intervention

Parameters D,k,a,eps,gamma come from the Member 3 config; grid, time axis and initial state come from Member 5's
reconstruction.npz (t, x, u_truth[0]) so nothing is typed in.  K candidate controls are simulated at once (batch axis).
"""
from dataclasses import dataclass
from pathlib import Path
import numpy as np, yaml


@dataclass
class Model:
    D: float; k: float; a: float; eps: float; gamma: float
    def scaled(self, D=None, k=None, a=None, eps=None, gamma=None, **kwargs):
        """Return a multiplicatively perturbed plant without changing self."""
        d_val = kwargs.get("d", kwargs.get("D", D if D is not None else 1.0))
        k_val = kwargs.get("k", kwargs.get("K", k if k is not None else 1.0))
        a_val = kwargs.get("a", kwargs.get("A", a if a is not None else 1.0))
        eps_val = kwargs.get("eps", kwargs.get("EPS", eps if eps is not None else 1.0))
        gam_val = kwargs.get("gamma", kwargs.get("gam", kwargs.get("GAM", gamma if gamma is not None else 1.0)))
        return Model(self.D * float(d_val), self.k * float(k_val), self.a * float(a_val), self.eps * float(eps_val), self.gamma * float(gam_val))


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


def simulate(model, grid, u0, prot, ctrl=None, clock_offset=0.0, p=2.0, store=False, dt_target=0.1):
    """Returns dict(damage, effort, peak, [u]) with batch axis K.
    prot: boolean mask over x (protected region). clock_offset: real_time = simulation_time + clock_offset; a scalar, or an array of length K (one offset per candidate).
    damage = integral over time and protected region of max(u,0)^p ; effort = integral of g^2 over x and t."""
    x, t, dx, od = grid.x, grid.t, grid.dx, grid.out_dt
    K = 1 if ctrl is None else len(ctrl.xc); nt, nx = len(t), len(x)
    nsub = max(int(np.ceil(od / dt_target)), 1)
    while od / nsub > 0.5 * 2.78 * dx * dx / (4 * model.D): nsub += 1            # RK4 diffusion stability
    dt = od / nsub
    u = np.tile(u0, (K, 1)); v = np.zeros_like(u)
    GS = None if ctrl is None else ctrl.g0[:, None] * zone(x, ctrl.xc, ctrl.width, ctrl.edge)
    D, k, a, eps, gam, dx2 = model.D, model.k, model.a, model.eps, model.gamma, dx * dx
    gt = (lambda tr: GS * gate(tr, ctrl.ts, ctrl.dur, ctrl.ramp)[:, None]) if GS is not None else (lambda tr: 0.0)

    def f(u, v, tr):
        lap = np.empty_like(u); lap[:, 1:-1] = u[:, 2:] - 2 * u[:, 1:-1] + u[:, :-2]
        lap[:, 0] = 2 * (u[:, 1] - u[:, 0]); lap[:, -1] = 2 * (u[:, -2] - u[:, -1]); lap /= dx2     # no-flux
        return D * lap + k * (u * (1 - u) * (u - a) - v) - gt(tr) * u, k * eps * (u - gam * v)

    dmg, eff, peak = np.zeros(K), np.zeros(K), np.full(K, -np.inf); U = np.empty((K, nt, nx)) if store else None
    for n in range(nt):
        wt = 0.5 if n in (0, nt - 1) else 1.0; tr0 = t[n] + clock_offset
        dmg += wt * od * dx * (np.maximum(u[:, prot], 0.0) ** p).sum(1); peak = np.maximum(peak, u[:, prot].max(1))
        if GS is not None: eff += wt * od * dx * (gt(tr0) ** 2).sum(1)
        if store: U[:, n] = u
        if n == nt - 1: break
        tn = t[n]
        for _ in range(nsub):
            tr = tn + clock_offset
            k1u, k1v = f(u, v, tr); k2u, k2v = f(u + .5 * dt * k1u, v + .5 * dt * k1v, tr + .5 * dt)
            k3u, k3v = f(u + .5 * dt * k2u, v + .5 * dt * k2v, tr + .5 * dt); k4u, k4v = f(u + dt * k3u, v + dt * k3v, tr + dt)
            u = u + dt / 6 * (k1u + 2 * k2u + 2 * k3u + k4u); v = v + dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v); tn += dt
    # This is deliberately the last-time-step maximum over the entire domain,
    # not the historical/protected-region peak used by the damage objective.
    out = dict(damage=dmg, effort=eff, peak=peak, final_max_u_domain=u.max(1))
    if store: out["u"] = U
    return out
