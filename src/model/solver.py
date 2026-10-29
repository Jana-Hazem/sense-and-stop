"""
Sense & Stop - PDE model (Member 3)

Reduced, phenomenological reaction-diffusion models of a CSD-inspired wave.
These are NOT ionic CSD simulators (no K+, glia, pumps, oxygen, perfusion).

Models
------
fisher_kpp : u_t = D lap(u) + r u (1 - u)
excitable  : u_t = D lap(u) + k ( u (1-u) (u-a) - v )
             v_t = k eps ( u - gamma v )
Optional control hook (Member 6): subtract g(x[,y],t) * u from u_t.

Public API
----------
simulate(config, control=None) -> dict
load_config(path)              -> dict
with_resolution(config, dx=None, dy=None, dt=None, method=None) -> dict
measure_front_speed(result, level=0.5, t_window=None) -> float
"""
from __future__ import annotations

import copy
from typing import Callable, Optional

import numpy as np
import yaml
from scipy.integrate import solve_ivp

__all__ = ["simulate", "load_config", "with_resolution", "measure_front_speed",
           "save_result", "load_result"]
__version__ = "0.1.0"


# --------------------------------------------------------------------------- config
def load_config(path: str) -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f)


def with_resolution(config: dict, dx=None, dy=None, dt=None, method=None) -> dict:
    """Return a deep copy of config with a different grid / step / integrator.
    Use this to make the 'fine truth grid' vs 'coarse inversion grid' pair."""
    cfg = copy.deepcopy(config)
    if dx is not None:
        cfg["domain"]["dx"] = dx
        if dy is None:
            cfg["domain"]["dy"] = dx
    if dy is not None:
        cfg["domain"]["dy"] = dy
    if dt is not None:
        cfg["time"]["dt"] = dt
    if method is not None:
        cfg["solver"]["method"] = method
    return cfg


# --------------------------------------------------------------------------- grid
def _build_grid(cfg: dict):
    dom = cfg["domain"]
    dim = int(cfg["model"]["dimension"])
    if dim not in (1, 2):
        raise ValueError("model.dimension must be 1 or 2")
    dx = float(dom["dx"])
    nx = int(round(dom["Lx"] / dx)) + 1
    x = np.linspace(0.0, dom["Lx"], nx)          # nodes include both boundaries
    dx = x[1] - x[0]
    if dim == 1:
        return dim, x, None, dx, None
    dy = float(dom["dy"])
    ny = int(round(dom["Ly"] / dy)) + 1
    y = np.linspace(0.0, dom["Ly"], ny)
    dy = y[1] - y[0]
    return dim, x, y, dx, dy


def _laplacian(f: np.ndarray, dx: float, dy: Optional[float], dim: int) -> np.ndarray:
    """Second-order central FD Laplacian with no-flux (Neumann) boundaries.
    Ghost nodes mirror the neighbour inside the domain (np.pad 'reflect'),
    which is the 2nd-order accurate discretisation of du/dn = 0."""
    if dim == 1:
        p = np.pad(f, 1, mode="reflect")
        return (p[2:] - 2.0 * p[1:-1] + p[:-2]) / dx ** 2
    p = np.pad(f, 1, mode="reflect")                 # f has shape (ny, nx)
    lap_x = (p[1:-1, 2:] - 2.0 * p[1:-1, 1:-1] + p[1:-1, :-2]) / dx ** 2
    lap_y = (p[2:, 1:-1] - 2.0 * p[1:-1, 1:-1] + p[:-2, 1:-1]) / dy ** 2
    return lap_x + lap_y


# --------------------------------------------------------------------------- initial condition
def _initial_condition(cfg, dim, x, y):
    ic = cfg["initial_condition"]
    amp = float(ic.get("amplitude", 1.0))
    w = float(ic.get("width", 1.0))
    if ic["type"] == "edge_step":
        prof = 0.5 * (1.0 - np.tanh((x - ic["x0"]) / w))
        u0 = prof if dim == 1 else np.tile(prof, (len(y), 1))
    elif ic["type"] == "kpp_exact_wave":
        # Ablowitz-Zeppetella exact travelling wave of u_t = D u_xx + r u(1-u):
        #   u(x,t) = [1 + exp(sqrt(r/(6D)) (x - x0 - c t))]^-2,  c = 5 sqrt(D r / 6).
        # Used ONLY for solver verification (1D, fisher_kpp).
        p = cfg["model"]["fisher_kpp"]
        lam = np.sqrt(p["r"] / (6.0 * p["D"]))
        u0 = (1.0 + np.exp(lam * (x - ic["x0"]))) ** -2
        return u0
    elif ic["type"] == "gaussian":
        cx, cy = ic["center"]
        if dim == 1:
            u0 = np.exp(-((x - cx) ** 2) / (2 * w ** 2))
        else:
            X, Y = np.meshgrid(x, y)
            u0 = np.exp(-((X - cx) ** 2 + (Y - cy) ** 2) / (2 * w ** 2))
    else:
        raise ValueError(f"unknown initial_condition.type: {ic['type']}")
    return amp * u0


# --------------------------------------------------------------------------- right-hand sides
def _make_rhs(cfg, dim, x, y, dx, dy, control: Optional[Callable]):
    mtype = cfg["model"]["type"]
    if dim == 2:
        X, Y = np.meshgrid(x, y)

    def g_term(u, t):
        if control is None:
            return 0.0
        g = control(X, Y, t) if dim == 2 else control(x, t)
        return g * u

    if mtype == "fisher_kpp":
        p = cfg["model"]["fisher_kpp"]
        D, r = float(p["D"]), float(p["r"])

        def rhs(t, state):
            u = state[0]
            du = D * _laplacian(u, dx, dy, dim) + r * u * (1.0 - u) - g_term(u, t)
            return np.stack([du])

        return rhs, D, 1

    if mtype == "excitable":
        p = cfg["model"]["excitable"]
        D, k, a, eps, gam = (float(p[s]) for s in ("D", "k", "a", "eps", "gamma"))

        def rhs(t, state):
            u, v = state
            du = D * _laplacian(u, dx, dy, dim) + k * (u * (1 - u) * (u - a) - v) - g_term(u, t)
            dv = k * eps * (u - gam * v)
            return np.stack([du, dv])

        return rhs, D, 2

    raise ValueError(f"unknown model.type: {mtype}")


# --------------------------------------------------------------------------- integrators
def _step_euler(rhs, t, s, dt):
    return s + dt * rhs(t, s)


def _step_rk4(rhs, t, s, dt):
    k1 = rhs(t, s)
    k2 = rhs(t + dt / 2, s + dt / 2 * k1)
    k3 = rhs(t + dt / 2, s + dt / 2 * k2)
    k4 = rhs(t + dt, s + dt * k3)
    return s + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


# --------------------------------------------------------------------------- main entry point
def simulate(config, control: Optional[Callable] = None) -> dict:
    """Run the forward model. `config` is a dict or a path to a YAML file.

    control : optional callable g(x, t) (1D) or g(X, Y, t) (2D) returning an array
              of the grid shape.  Default None = uncontrolled primary code path.

    Returns dict with keys
      u      : ndarray [t, x] (1D) or [t, y, x] (2D), dimensionless, 0 = rest, 1 = elevated
      v      : recovery variable, same shape (excitable only; else None)
      t, x, y: coordinate arrays (s, mm, mm); y is None in 1D
      units  : dict of units
      params : model parameters used
      meta   : solver metadata (dt actually used, stability number, versions, config copy)
    """
    cfg = load_config(config) if isinstance(config, str) else copy.deepcopy(config)
    if cfg["domain"].get("boundary", "neumann") != "neumann":
        raise NotImplementedError("only 'neumann' (no-flux) boundaries are implemented")

    dim, x, y, dx, dy = _build_grid(cfg)
    rhs, D, nvar = _make_rhs(cfg, dim, x, y, dx, dy, control)

    # ---- time grid: output frames are exact; internal dt divides output_dt ----------
    T = float(cfg["time"]["T"])
    out_dt = float(cfg["time"]["output_dt"])
    n_out = int(round(T / out_dt))
    t_out = np.arange(n_out + 1) * out_dt
    dt_target = float(cfg["time"]["dt"])
    method = cfg["solver"]["method"]

    inv_h2 = 1.0 / dx ** 2 + (1.0 / dy ** 2 if dim == 2 else 0.0)
    stab_limit = {"euler": 0.5, "rk4": 0.69}.get(method, None)
    if stab_limit is not None:
        dt_max = cfg["solver"].get("safety", 0.9) * stab_limit / (D * inv_h2)
        dt_target = min(dt_target, dt_max)
    n_sub = int(np.ceil(out_dt / dt_target - 1e-12))
    dt = out_dt / n_sub
    stab_number = D * dt * inv_h2                      # must be <= 0.5 (Euler) / 0.69 (RK4)

    # ---- initial state ---------------------------------------------------------------
    u0 = _initial_condition(cfg, dim, x, y)
    state0 = np.stack([u0] + ([np.zeros_like(u0)] if nvar == 2 else []))
    frames = np.empty((n_out + 1,) + state0.shape)
    frames[0] = state0

    if method in ("euler", "rk4"):
        step = _step_euler if method == "euler" else _step_rk4
        s, t = state0.copy(), 0.0
        for i in range(1, n_out + 1):
            for _ in range(n_sub):
                s = step(rhs, t, s, dt)
                t += dt
            frames[i] = s
    elif method == "solve_ivp":
        shape = state0.shape
        sol = solve_ivp(
            lambda t, z: rhs(t, z.reshape(shape)).ravel(),
            (0.0, t_out[-1]), state0.ravel(), t_eval=t_out,
            method=cfg["solver"].get("ivp_method", "RK45"),
            rtol=float(cfg["solver"].get("ivp_rtol", 1e-6)),
            atol=float(cfg["solver"].get("ivp_atol", 1e-8)),
        )
        if not sol.success:
            raise RuntimeError(f"solve_ivp failed: {sol.message}")
        frames = sol.y.T.reshape((n_out + 1,) + shape)
        dt, stab_number = float("nan"), float("nan")   # adaptive stepping
    else:
        raise ValueError(f"unknown solver.method: {method}")

    if not np.all(np.isfinite(frames)):
        raise RuntimeError("non-finite values: reduce dt or check parameters")

    mtype = cfg["model"]["type"]
    return {
        "u": frames[:, 0],
        "v": frames[:, 1] if nvar == 2 else None,
        "t": t_out,
        "x": x,
        "y": y,
        "units": {"x": "mm", "y": "mm", "t": "s", "u": "dimensionless (0=rest, 1=elevated)",
                  "array_order": "u[t, x] in 1D; u[t, y, x] in 2D"},
        "params": copy.deepcopy(cfg["model"][mtype]) | {"model_type": mtype, "dimension": dim},
        "meta": {"solver": method, "dt_internal_s": dt, "dx_mm": dx, "dy_mm": dy,
                 "stability_number": stab_number, "control_active": control is not None,
                 "solver_version": __version__, "config": cfg},
    }


# --------------------------------------------------------------------------- diagnostics
def measure_front_speed(result: dict, level: float = 0.5, t_window=None) -> float:
    """Leading-front speed (mm/s): linear fit of the position where u first crosses
    `level` (rightmost crossing, along the mid-line in 2D), over t_window=(t0, t1)."""
    u, x, t = result["u"], result["x"], result["t"]
    if u.ndim == 3:
        u = u[:, u.shape[1] // 2, :]
    pos, times = [], []
    for i, ti in enumerate(t):
        if t_window and not (t_window[0] <= ti <= t_window[1]):
            continue
        above = np.where(u[i] >= level)[0]
        if len(above) == 0 or above[-1] >= len(x) - 1:
            continue
        j = above[-1]
        # linear interpolation of the crossing between node j and j+1
        frac = (u[i, j] - level) / (u[i, j] - u[i, j + 1])
        pos.append(x[j] + frac * (x[j + 1] - x[j]))
        times.append(ti)
    if len(pos) < 3:
        raise ValueError("front not found in the requested window")
    return float(np.polyfit(times, pos, 1)[0])


# --------------------------------------------------------------------------- I/O
def save_result(result: dict, path: str) -> None:
    """Save a simulate() result to a single .npz (arrays) with a JSON metadata string."""
    import json
    meta = {k: result[k] for k in ("units", "params", "meta")}
    np.savez_compressed(
        path, u=result["u"], v=np.array([]) if result["v"] is None else result["v"],
        t=result["t"], x=result["x"], y=np.array([]) if result["y"] is None else result["y"],
        metadata=json.dumps(meta, default=float))


def load_result(path: str) -> dict:
    """Inverse of save_result; returns the same dict layout as simulate()."""
    import json
    z = np.load(path, allow_pickle=False)
    meta = json.loads(str(z["metadata"]))
    return {"u": z["u"], "v": z["v"] if z["v"].size else None, "t": z["t"], "x": z["x"],
            "y": z["y"] if z["y"].size else None, **meta}
