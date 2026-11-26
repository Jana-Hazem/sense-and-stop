"""Supplementary checks for the retrospective hypothetical control study.

This intentionally keeps slow confirmation checks out of run_control.py.  It
expects a completed control run's control_metadata.json, reuses the exact plan
and model configuration, and writes separate, audit-friendly CSV files.

Example:
  python scripts/control_checks.py --metadata results/.../control_metadata.json
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
import sys
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.control.model import Grid, load_model, simulate
from src.control.objective import design_from_cfg, make_ctrl, search_bounds
from src.control.estimate import make_nominal
from src.control.search import optimise


def args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--metadata", required=True, help="control_metadata.json from a completed run")
    p.add_argument("--config", default=None, help="control.yaml; defaults to control_used.yaml beside metadata")
    p.add_argument("--model-config", default=None, help="model YAML; defaults to the recorded input")
    p.add_argument("--out", default=None, help="directory for CSVs; defaults to the completed-run directory")
    p.add_argument("--convergence-iterations", type=int, default=None, help="DE iterations for convergence rerun; defaults to 2x config")
    return p.parse_args(argv)


def initial_state(model_cfg, x):
    ic = model_cfg["initial_condition"]
    return float(ic.get("amplitude", 1.0)) * .5 * (1 - np.tanh((x - float(ic["x0"])) / float(ic["width"])))


def blocked(out, baseline_peak, block_fraction, extinct_u):
    peak_ratio = float(out["peak"][0] / max(baseline_peak, 1e-12))
    extinct = bool(out["final_max_u_domain"][0] < extinct_u)
    return peak_ratio, extinct, bool(peak_ratio < block_fraction and extinct)


def main(argv=None):
    a = args(argv)
    meta_path = Path(a.metadata).resolve(); meta = json.loads(meta_path.read_text(encoding="utf-8"))
    out = Path(a.out).resolve() if a.out else meta_path.parent
    out.mkdir(parents=True, exist_ok=True)
    cfg_path = Path(a.config).resolve() if a.config else meta_path.parent / "control_used.yaml"
    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    model_path = Path(a.model_config).resolve() if a.model_config else Path(meta["inputs"]["model"])
    nominal, model_cfg = load_model(model_path)
    dom, tm = model_cfg["domain"], model_cfg["time"]
    dx, out_dt, T = map(float, (dom["dx"], tm["output_dt"], tm["T"]))
    x = np.arange(0.0, float(dom["Lx"]) + .5 * dx, dx); t = np.arange(0.0, T + .5 * out_dt, out_dt)
    grid = Grid(x, t); u0 = initial_state(model_cfg, x); d = design_from_cfg(cfg, float(t[-1]))
    prot = (x >= d.prot_lo) & (x <= d.prot_hi); dt = float(tm.get("dt", .1))
    P0 = np.array([meta["plan"][q] for q in ("xc", "ts", "dur", "g0")], float)
    onset = float(meta.get("plant_onset_shift_s", 0.0)); extinct_u = float(cfg["objective"].get("extinct_u", .1))
    plant_cfg = cfg["plant"]
    plant = nominal.scaled(**{q: float(plant_cfg.get(f"{q}_factor", 1.0)) for q in ("D", "k", "a", "eps", "gamma")})
    base = simulate(plant, grid, u0, prot, None, p=d.p, dt_target=dt)

    # 1. Strength margin at the two pre-declared factors.
    strength = []
    for factor in cfg.get("planning", {}).get("strength_factors", [.5, .75]):
        P = P0.copy(); P[3] *= float(factor)
        r = simulate(plant, grid, u0, prot, make_ctrl(P[:, None], d), clock_offset=onset, p=d.p, dt_target=dt)
        pr, extinct, ok = blocked(r, float(base["peak"][0]), d.block_fraction, extinct_u)
        strength.append(dict(factor=float(factor), g0=float(P[3]), duration_s=float(P[2]), dose_g0_s=float(P[2] * P[3]), peak_ratio="below 10%" if pr < .1 else f"{100 * pr:.1f}%", final_max_u_domain=float(r["final_max_u_domain"][0]), wave_extinct=extinct, blocked=ok))
    pd.DataFrame(strength).to_csv(out / "check_strength_margin.csv", index=False)

    # 2. All declared parameter perturbations, including a, eps, and gamma.
    cases = []
    for raw in cfg.get("evaluation", {}).get("plant_cases", []):
        fac = {q: 1.0 for q in ("D", "k", "a", "eps", "gamma")}; fac.update({q: float(v) for q, v in raw.items()})
        m = plant.scaled(**fac); b = simulate(m, grid, u0, prot, None, p=d.p, dt_target=dt)
        self_extinguishing = bool(b["final_max_u_domain"][0] < extinct_u)
        if self_extinguishing:
            cases.append(dict(**fac, self_extinguishing=True, blocked=False, peak_ratio="not applicable", final_max_u_domain=float(b["final_max_u_domain"][0]))); continue
        r = simulate(m, grid, u0, prot, make_ctrl(P0[:, None], d), clock_offset=onset, p=d.p, dt_target=dt)
        pr, extinct, ok = blocked(r, float(b["peak"][0]), d.block_fraction, extinct_u)
        cases.append(dict(**fac, self_extinguishing=False, blocked=ok, peak_ratio="below 10%" if pr < .1 else f"{100 * pr:.1f}%", final_max_u_domain=float(r["final_max_u_domain"][0]), wave_extinct=extinct))
    case_df = pd.DataFrame(cases); case_df.to_csv(out / "check_plant_cases.csv", index=False)
    pd.DataFrame([dict(n_plant_cases=len(case_df), n_self_extinguishing=int(case_df.self_extinguishing.sum()) if len(case_df) else 0)]).to_csv(out / "check_plant_case_summary.csv", index=False)

    # 3. Multi-seed plans and 4. a doubled-iteration convergence rerun.
    p_obs = make_nominal(simulate(nominal, grid, u0, prot, None, p=d.p, store=True, dt_target=dt)["u"][0], x, t).at(d.t_obs)
    bounds = search_bounds(d, p_obs)
    S = float(meta["plan"]["halfwidth_s"])
    bounds[1] = (min(max(bounds[1][0], d.t_obs + S), bounds[1][1] - 1e-6), bounds[1][1])
    factors = [tuple(map(float, q)) for q in cfg["planning"].get("model_factors", [[1, 1]])]
    planning_cases = []
    for fd, fk in factors:
        m = nominal.scaled(fd, fk); b = simulate(m, grid, u0, prot, None, p=d.p, dt_target=dt)
        planning_cases += [dict(model=m, D0=float(b["damage"][0]), offset=o) for o in (-S, 0.0, S)]
    s = dict(cfg["search"])
    seed_rows = []
    for seed in s.get("seeds", [s["seed"]]):
        r = optimise(planning_cases, grid, u0, prot, d, bounds, s, seed_offset=int(seed))
        seed_rows.append(dict(seed=int(seed), J_plan=float(r["J_plan"]), xc=float(r["P"][0]), ts=float(r["P"][1]), duration_s=float(r["P"][2]), g0=float(r["P"][3]), dose_g0_s=float(r["P"][2] * r["P"][3]), nfev_batches=int(r["nfev_batches"])))
    pd.DataFrame(seed_rows).to_csv(out / "check_multiseed_plans.csv", index=False)
    long_s = dict(s); long_s["maxiter"] = int(a.convergence_iterations or 2 * int(s["maxiter"]))
    long_r = optimise(planning_cases, grid, u0, prot, d, bounds, long_s, seed_offset=int(s.get("seeds", [s["seed"]])[0]))
    pd.DataFrame([dict(seed=int(s.get("seeds", [s["seed"]])[0]), maxiter=long_s["maxiter"], J_plan=float(long_r["J_plan"]), xc=float(long_r["P"][0]), ts=float(long_r["P"][1]), duration_s=float(long_r["P"][2]), g0=float(long_r["P"][3]), dose_g0_s=float(long_r["P"][2] * long_r["P"][3]), nfev_batches=int(long_r["nfev_batches"]))]).to_csv(out / "check_convergence_rerun.csv", index=False)


if __name__ == "__main__":
    main()
