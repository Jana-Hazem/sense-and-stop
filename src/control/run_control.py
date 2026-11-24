"""Member 6B v3 - hypothetical wave-control study (FAST). Called by run_pipeline.py or alone:

    python -m src.control.run_control --reconstruction <reconstruction.npz> --event-estimate <event_estimate.json>
    (flags: --quick  --no-event  --confirm-2d  --tradeoff  --no-mismatch)

Idea (why it is fast): the wave is uniform in y, so the control is a STRIPE g(x,t) = g0*zone(x)*gate(t) that covers the whole
wave front (1D model). The uncontrolled wave travels at constant speed in the control region, so a plan found for "timing error 0"
is shifted in time by the estimated timing offset:  plan(delta) = (xc, ts0 + delta, dur, g0).
We search ONCE for a robust plan (worst case over timing offsets +-S and model factors D,k), then apply it to every estimate
(Member 5 rec_* and event_*), and re-test all plans in single batches. Typical run: a few minutes, not an hour.

Scenarios: none | no_sensing (fixed P0) | oracle | rec:<method>:<scale> | event:<scale>:<i>
RETROSPECTIVE, exploratory, hypothetical numerical study; no clinical claim.
Reconstruction and event estimates are calculated from data through t_obs only.
"""
import argparse, glob, hashlib, json, math, platform, re, shutil, subprocess, sys, time
from pathlib import Path
import numpy as np, pandas as pd, yaml, scipy
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
from src.control.model import Grid, load_model, simulate
from src.control.objective import design_from_cfg, make_ctrl, cost, search_bounds
from src.control.estimate import make_nominal, peak_pos, arrival_time, delta_from_rec, delta_from_arrival
from src.control.events import load_events
from src.control.search import optimise
try:
    from src.model import load_config as m3_load, simulate as m3_simulate       # Member 3's own solver (reference + 2D check)
except Exception:                                                                 # control can also run without it
    m3_load = m3_simulate = None

METHOD_OF_LETTER = {"A": "A_smooth", "B": "B_model"}


def parse_args(argv=None):
    ap = argparse.ArgumentParser(); a = ap.add_argument
    a("--config", default="configs/control.yaml"); a("--model-config", default="configs/model_default.yaml"); a("--master-config", default="configs/config.yaml")
    a("--reconstruction", required=True); a("--event-estimate", default=None); a("--out", default=None)
    a("--no-event", action="store_true"); a("--quick", action="store_true"); a("--no-mismatch", action="store_true")
    a("--confirm-2d", action="store_true", help="re-run the oracle plan in Member 3's 2D solver (about 1 min)")
    a("--tradeoff", action="store_true", help="also plan for several timing half-widths S and write tradeoff.csv (slower)")
    a("--allow-solver-mismatch", action="store_true")
    return ap.parse_args(argv)


def rooted(p):
    if p is None: return None
    path = Path(p)
    path = path if path.is_absolute() else ROOT / path
    if not path.exists():
        parts = list(path.parts)
        if "latest" in parts:
            idx = parts.index("latest")
            base = Path(*parts[:idx])
            if base.exists():
                runs = sorted([d for d in base.glob("run_*") if d.is_dir()], key=lambda d: d.stat().st_mtime)
                if runs:
                    cand = runs[-1] / Path(*parts[idx + 1:])
                    if cand.exists(): return cand
        if "*" in str(p):
            m = sorted(glob.glob(str(ROOT / p if not p.is_absolute() else p)))
            if m: return Path(m[-1])
    return path


def load_rec(path):
    """Understands Member 5 files of both layouts. Returns t, x, truth|None, {name: (method, scale, array[t, x])}."""
    z = np.load(path, allow_pickle=True); f = z.files
    t = np.asarray(z["t" if "t" in f else "time_s"], float); x = np.asarray(z["x" if "x" in f else "x_mm"], float)
    truth = next((np.asarray(z[k], float) for k in ("u_truth", "truth") if k in f), None); recs = {}
    for k in f:
        m = re.fullmatch(r"rec_(.+)_([0-9.]+e[-+]?\d+)", k) or None
        if m: recs[f"{m.group(1)}:{m.group(2)}"] = (m.group(1), m.group(2), np.asarray(z[k], float)); continue
        m = re.fullmatch(r"([0-9.]+e[-+]?\d+)_([AB])", k)                                   # '1e-09_B' layout
        if m: recs[f"{METHOD_OF_LETTER[m.group(2)]}:{m.group(1)}"] = (METHOD_OF_LETTER[m.group(2)], m.group(1), np.asarray(z[k], float))
    if not recs: raise SystemExit(f"ERROR: no reconstruction arrays found in {path}; keys: {f}")
    return t, x, truth, recs


class Stripe:
    """Callable g(x,t) / g(X,Y,t) for Member 3's solver: the same stripe zone x gate used by the batch solver (uniform in y)."""
    def __init__(self, P, d): self.P, self.d = P, d
    def __call__(self, X, Y=None, t=None):
        if t is None: t, Y = Y, None
        from src.control.model import zone, gate
        xc, ts, dur, g0 = self.P; x1 = np.asarray(X)[0] if np.ndim(X) == 2 else np.asarray(X)
        z = zone(x1, np.array([xc]), self.d.width, self.d.edge)[0]; gt = float(gate(np.array(t), np.array(ts), np.array(dur), self.d.ramp))
        g = g0 * gt * z; return np.broadcast_to(g, np.shape(X)).copy()


def main(argv=None):
    a = parse_args(argv); cfg = yaml.safe_load(rooted(a.config).read_text(encoding="utf-8")); T0 = time.time()
    s = dict(cfg["search"]); est = cfg["estimation"]; pl = cfg.get("planning", {}); ev_cfg = cfg.get("evaluation", {}); notes = []
    if a.quick: s.update(popsize=4, maxiter=8, polish_maxiter=10)
    p_model = rooted(a.model_config); nominal, mcfg = load_model(p_model)
    plant_cfg = cfg["plant"]
    plant = nominal.scaled(**{q: float(plant_cfg.get(f"{q}_factor", 1.0)) for q in ("D", "k", "a", "eps", "gamma")})
    onset_rng = np.random.default_rng(int(plant_cfg.get("onset_seed", 0)))
    onset_shift = float(onset_rng.uniform(*map(float, plant_cfg.get("onset_shift_s", [0.0, 0.0]))))
    rt, rx, truth, recs = load_rec(rooted(a.reconstruction))
    dom, tm = mcfg["domain"], mcfg["time"]; dx, od, Tend = float(dom["dx"]), float(tm["output_dt"]), float(tm["T"])
    x = np.arange(0.0, float(dom["Lx"]) + 0.5 * dx, dx); t = np.arange(0.0, Tend + 0.5 * od, od); grid = Grid(x, t)     # model grid, NOT the coarse inverse grid
    d = design_from_cfg(cfg, float(t[-1])); prot = (x >= d.prot_lo) & (x <= d.prot_hi); dt_final = float(tm.get("dt", 0.1))
    if not prot.any(): raise SystemExit("ERROR: protected region is outside the domain")
    # ---- initial state and solver check against Member 3's own 1D solver
    if m3_simulate is not None:
        c1 = m3_load(str(p_model)); c1["model"]["dimension"] = 1; r1 = m3_simulate(c1); u_ref = np.asarray(r1["u"], float); u0 = u_ref[0]; u0src = "Member 3 solver (1D) u[0]"
        if len(r1["x"]) != len(x) or len(r1["t"]) != len(t): raise SystemExit(f"ERROR: Member 3 grid {len(r1['x'])}x{len(r1['t'])} differs from the control grid {len(x)}x{len(t)}")
    else: u_ref = None; ic = mcfg["initial_condition"]; u0 = float(ic.get("amplitude", 1)) * 0.5 * (1 - np.tanh((x - float(ic["x0"])) / float(ic["width"]))); u0src = "config edge_step"
    nm_run = simulate(nominal, grid, u0, prot, None, p=d.p, store=True, dt_target=dt_final); pl_run = simulate(plant, grid, u0, prot, None, p=d.p, store=True, dt_target=dt_final)
    Unm, Upl = nm_run["u"][0], pl_run["u"][0]; nom = make_nominal(Unm, x, t); pos_pl = peak_pos(Upl, x)
    verify = float(np.linalg.norm(Unm - u_ref) / np.linalg.norm(u_ref)) if u_ref is not None else float("nan")
    print(f"solver check vs Member 3 (1D, g=0): rel L2 = {verify:.2e} | initial state: {u0src}")
    if np.isfinite(verify) and verify > 0.05:
        if not a.allow_solver_mismatch: raise SystemExit(f"ERROR: solver check {verify:.2e} > 0.05; the control solver does not match Member 3's model")
        notes.append(f"WARNING solver check {verify:.2e} > 0.05")
    D0_pl, D0_nm = float(pl_run["damage"][0]), float(nm_run["damage"][0]); peak_pl = float(pl_run["peak"][0])
    if D0_pl < 1e-8: raise SystemExit("ERROR: the uncontrolled wave never reaches the protected region")
    # ---- events and timing estimates (delta: real = nominal + delta)
    ev = None; p_evt = None
    if not a.no_event and a.event_estimate:
        p_evt = rooted(a.event_estimate); ev = load_events(p_evt, str(cfg.get("event", {}).get("method_key", "B")), float(cfg.get("target", {}).get("x_mm", 18.0)))
        for w in ev["warnings"]: notes.append("event_estimate: " + w); print("[WARNING] event_estimate:", w)
    tx = ev["target_x_mm"] if ev else float(cfg.get("target", {}).get("x_mm", 18.0)); a_nom, a_pl = arrival_time(nom.pos, t, tx), arrival_time(pos_pl, t, tx)
    truth_arr = ev["truth_arrival_s"] if ev and ev["truth_arrival_s"] is not None else a_pl
    scen = []                                                                      # (name, source, method, scale, seed, delta, status, extra)
    for key, (method, sc, arr) in recs.items():
        name = f"rec:{method}:{sc}"; r = delta_from_rec(arr, rx, rt, nom, d.t_obs, float(est["detect_fraction"]), float(est["front_window_s"]), int(est["min_points"]))
        if r is None: scen.append((name, "rec", method, float(sc), 0, np.nan, "no_detection", {}))
        else: scen.append((name, "rec", method, float(sc), 0, r["delta"], "ok", dict(p_hat_est=r["p_hat"])))
    if ev:
        for sc, vals in ev["per_scale"].items():
            for i, v in enumerate(vals):
                name = f"event:{sc}:{i}"
                if not np.isfinite(v): scen.append((name, "event", "B_model", float(sc), i, np.nan, "no_detection", {})); continue
                scen.append((name, "event", "B_model", float(sc), i, delta_from_arrival(v, truth_arr, a_nom), "ok", dict(arrival_err_s=v - truth_arr)))
            for j in range(ev.get("n_failed", {}).get(sc, 0)):                      # failed fits count as "no estimate -> no control" (not hidden)
                scen.append((f"event:{sc}:failed{j}", "event", "B_model", float(sc), 1000 + j, np.nan, "no_detection", {}))
    dl = np.array([q[5] for q in scen if q[6] == "ok"]); dmax = float(np.max(np.abs(dl))) if len(dl) else 0.0
    # ---- robust plan, found ONCE
    hw = pl.get("timing_halfwidth_s")
    if hw in (None, "auto"): raise SystemExit("ERROR: planning.timing_halfwidth_s must be a pre-declared numeric value, not auto")
    S = float(hw)
    factors = [tuple(map(float, f)) for f in (pl.get("model_factors", [[1, 1]]) if not a.quick else [[1, 1]])]
    print(f"{len(scen)} retrospective scenarios (data through t_obs={d.t_obs:g} s) | largest |delta| = {dmax:.1f} s (not used to choose S) | pre-declared S = {S:g} s | model factors {factors} | random plant onset shift {onset_shift:+.2f} s")
    p_obs = float(nom.at(d.t_obs)); models = {}
    def cases_for(Sx):
        out = []
        for fd, fk in factors:
            if (fd, fk) not in models: m = nominal.scaled(fd, fk); models[(fd, fk)] = (m, float(simulate(m, grid, u0, prot, None, p=d.p, dt_target=dt_final)["damage"][0]))
            m, D0m = models[(fd, fk)]; out += [dict(model=m, D0=D0m, offset=float(o)) for o in sorted({-Sx, 0.0, Sx})]
        return out
    def plan(Sx):
        bnd = search_bounds(d, p_obs)
        if bnd is None: raise SystemExit("ERROR: no admissible control location between the front at t_obs and the protected region")
        bnd[1] = (min(max(bnd[1][0], d.t_obs + Sx), bnd[1][1] - 1e-6), bnd[1][1]); t1 = time.time()
        runs = [optimise(cases_for(Sx), grid, u0, prot, d, bnd, s, seed_offset=int(seed)) for seed in s.get("seeds", [s["seed"]])]
        # A plan is feasible only when every planning case both has a protected
        # peak below 10% and leaves no wave front alive at the end of the run.
        extinct_u = float(cfg["objective"].get("extinct_u", 0.1))
        for r in runs:
            checks = []
            for c in cases_for(Sx):
                o = simulate(c["model"], grid, u0, prot, make_ctrl(np.asarray(r["P"], float)[:, None], d), clock_offset=c["offset"], p=d.p, dt_target=dt_final)
                checks.append(bool(o["peak"][0] / max(float(simulate(c["model"], grid, u0, prot, None, p=d.p, dt_target=dt_final)["peak"][0]), 1e-12) < d.block_fraction and o["final_max_u_domain"][0] < extinct_u))
            r["feasible"] = all(checks)
        feasible = [r for r in runs if r["feasible"]]
        if not feasible:
            notes.append("WARNING: no multi-seed candidate satisfied every planning case; reporting the lowest-J candidate, not a feasible plan")
        chosen = min(feasible or runs, key=lambda r: r["J_plan"])
        chosen["g0_feasible_range"] = ([float(r["P"][3]) for r in feasible] or [float(chosen["P"][3])])
        chosen["seconds"] = time.time() - t1; chosen["n_seeds"] = len(runs); return chosen
    R = plan(S); P0 = np.array(R["P"], float); print(f"robust plan: xc {P0[0]:.2f} mm, ts {P0[1]:.1f} s, duration {P0[2]:.1f} s, dose {P0[2] * P0[3]:.4f}, feasible g0 range [{min(R['g0_feasible_range']):.4f}, {max(R['g0_feasible_range']):.4f}] /s | timing tolerance ±{S:g} s | {R['n_seeds']} seeds | {R['seconds']:.0f} s")
    tradeoff = []
    if a.tradeoff:
        for Sx in pl.get("tradeoff_halfwidths_s", [0, 5, 10, 20, 30]):
            Rx = plan(float(Sx)); Px = Rx["P"]; ox = simulate(plant, grid, u0, prot, make_ctrl(np.array(Px, float), d), p=d.p, dt_target=dt_final)
            tradeoff.append(dict(S_s=float(Sx), xc=Px[0], ts=Px[1], dur=Px[2], g0=Px[3], effort_norm=float(ox["effort"][0] / d.e_ref), J_plan=Rx["J_plan"]))
    # ---- one plan per scenario: same plan, shifted in time by the estimated offset
    P_oracle = P0.copy(); P_oracle[1] += onset_shift  # oracle knows the true random onset; no-sensing does not
    plans = {"no_sensing": P0.copy(), "oracle": P_oracle}; rows = [dict(scenario="none", source="none", method=None, scale=np.nan, seed=-1, status="ok")]
    rows.append(dict(scenario="no_sensing", source="no_sensing", method="fixed_P0", scale=np.nan, seed=-1, status="ok", delta_s=0.0, xc=P0[0], ts=P0[1], dur=P0[2], g0=P0[3], dose_g0_s=P0[2] * P0[3], J_plan=R["J_plan"]))
    rows.append(dict(scenario="oracle", source="oracle", method="full_state", scale=np.nan, seed=-1, status="ok", delta_s=onset_shift, xc=P_oracle[0], ts=P_oracle[1], dur=P_oracle[2], g0=P_oracle[3], dose_g0_s=P_oracle[2] * P_oracle[3], J_plan=R["J_plan"]))
    for name, src, method, sc, seed, dlt, status, extra in scen:
        base = dict(scenario=name, source=src, method=method, scale=sc, seed=seed, delta_s=dlt, status=status, **extra)
        if status != "ok": rows.append(base); continue
        P = P0.copy(); P[1] = P0[1] + dlt
        if P[1] < d.t_obs:
            rows.append({**base, "status": "ts_before_t_obs"})
            continue
        plans[name] = P; rows.append({**base, "xc": P[0], "ts": P[1], "dur": P[2], "g0": P[3], "dose_g0_s": P[2] * P[3]})
    names = list(plans); Pm = np.array([plans[n] for n in names]).T
    def run_case(model, shift=0.0, store=False):
        r0 = simulate(model, grid, u0, prot, None, p=d.p, dt_target=dt_final); o = simulate(model, grid, u0, prot, make_ctrl(Pm, d), clock_offset=shift, p=d.p, dt_target=dt_final)
        J, dr, en = cost(o, float(r0["damage"][0]), d); return J, dr, en, o["peak"] / float(r0["peak"][0]), o["final_max_u_domain"]
    J, dr, en, pk, final_u = run_case(plant, onset_shift)
    extinct_u = float(cfg["objective"].get("extinct_u", 0.1)); blocked = (pk < d.block_fraction) & (final_u < extinct_u)
    res = {n: dict(J_true=J[i], damage_ratio=dr[i], effort_norm=en[i], peak_ratio="below 10%" if pk[i] < d.block_fraction else f"{100 * pk[i]:.1f}%", final_max_u_domain=final_u[i], wave_extinct=bool(final_u[i] < extinct_u), blocked=bool(blocked[i])) for i, n in enumerate(names)}
    df = pd.DataFrame(rows)
    for c in ("J_true", "damage_ratio", "effort_norm", "peak_ratio", "final_max_u_domain", "wave_extinct", "blocked"):
        df[c] = [res[n][c] if n in res else ({"J_true": 1.0, "damage_ratio": 1.0, "effort_norm": 0.0, "peak_ratio": "100.0%", "final_max_u_domain": np.nan, "wave_extinct": False, "blocked": False}[c]) for n in df.scenario]
    # ---- re-tests (one batch per case)
    rb, mm = [], []
    if not a.no_mismatch:
        for sh in [float(v) for v in cfg.get("robustness", {}).get("timing_shifts_s", list(range(-40, 41, 5)))]:
            J_, dr_, en_, pk_, fu_ = run_case(plant, onset_shift + sh); rb += [dict(scenario=n, shift_s=sh, damage_ratio=float(dr_[i]), peak_ratio_report="below 10%" if pk_[i] < d.block_fraction else f"{100 * pk_[i]:.1f}%", final_max_u_domain=float(fu_[i]), wave_extinct=bool(fu_[i] < extinct_u), blocked=bool(pk_[i] < d.block_fraction and fu_[i] < extinct_u)) for i, n in enumerate(names)]
        self_extinguished = 0
        for pc in ev_cfg.get("plant_cases", []):
            fac = {q: 1.0 for q in ("D", "k", "a", "eps", "gamma")}; fac.update({q: float(v) for q, v in pc.items()})
            m = plant.scaled(**fac); base = simulate(m, grid, u0, prot, None, p=d.p, dt_target=dt_final)
            if float(base["final_max_u_domain"][0]) < extinct_u:
                self_extinguished += 1; notes.append(f"plant case {fac} self-extinguishes; excluded from controlled comparisons"); continue
            if float(base["damage"][0]) < 1e-8: notes.append(f"plant case {fac}: wave never reaches protected region; skipped"); continue
            J_, dr_, en_, pk_, fu_ = run_case(m, onset_shift); mm += [dict(scenario=n, **{f"{q}_factor": fac[q] for q in fac}, damage_ratio=float(dr_[i]), peak_ratio_report="below 10%" if pk_[i] < d.block_fraction else f"{100 * pk_[i]:.1f}%", final_max_u_domain=float(fu_[i]), wave_extinct=bool(fu_[i] < extinct_u), blocked=bool(pk_[i] < d.block_fraction and fu_[i] < extinct_u)) for i, n in enumerate(names)]
    rb, mm = pd.DataFrame(rb), pd.DataFrame(mm); tol_cfg = cfg.get("robustness", {}).get("tolerance_s", "auto"); tol = S if tol_cfg == "auto" else float(tol_cfg)
    if len(rb):                                    # plan = P0 shifted by delta, so the total timing error of the plan is (shift - delta): count shifts with |shift - delta| <= tol
        dmap = dict(zip(df.scenario, df.delta_s.fillna(0.0))); rb["eff_err_s"] = onset_shift + rb.shift_s - rb.scenario.map(dmap)
        df = df.merge(rb[rb.eff_err_s.abs() <= tol].groupby("scenario").blocked.mean().rename("blocked_frac_timing").reset_index(), on="scenario", how="left")
    if len(mm): df = df.merge(mm.groupby("scenario").blocked.mean().rename("blocked_frac_mismatch").reset_index(), on="scenario", how="left")
    # A deliberately weaker-plan check: how much of the selected dose margin
    # remains at 50% and 75% of the chosen strength?
    strength_rows = []
    baseline = simulate(plant, grid, u0, prot, None, p=d.p, dt_target=dt_final)
    for factor in pl.get("strength_factors", [0.5, 0.75]):
        Pweak = P0.copy(); Pweak[3] *= float(factor)
        oo = simulate(plant, grid, u0, prot, make_ctrl(Pweak[:, None], d), clock_offset=onset_shift, p=d.p, dt_target=dt_final)
        pr = float(oo["peak"][0] / baseline["peak"][0]); extinct = bool(oo["final_max_u_domain"][0] < extinct_u)
        strength_rows.append(dict(factor=float(factor), g0=float(Pweak[3]), duration_s=float(Pweak[2]), dose_g0_s=float(Pweak[2] * Pweak[3]), peak_ratio_report="below 10%" if pr < d.block_fraction else f"{100 * pr:.1f}%", final_max_u_domain=float(oo["final_max_u_domain"][0]), wave_extinct=extinct, blocked=bool(pr < d.block_fraction and extinct)))
    # ---- outputs
    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()); out_base = rooted(a.out or (cfg["output"].get("quick_dir") if a.quick else cfg["output"]["dir"]))
    out = out_base / f"run_{stamp}"; out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "control_results.csv", index=False); sub = df[df.source.isin(["no_sensing", "oracle", "rec", "event"])]
    agg = dict(n=("J_true", "size"), blocked_frac=("blocked", "mean"), damage_ratio_med=("damage_ratio", "median"), damage_ratio_max=("damage_ratio", "max"), effort_med=("effort_norm", "median"),
               no_plan=("status", lambda v: float((v != "ok").mean())), delta_med_s=("delta_s", "median"), delta_absmax_s=("delta_s", lambda v: float(np.nanmax(np.abs(v))) if np.isfinite(v).any() else np.nan),
               event_bias_s=("arrival_err_s", "mean"), event_spread_sd_s=("arrival_err_s", "std"))
    for c in ("blocked_frac_timing", "blocked_frac_mismatch"):
        if c in sub: agg[c] = (c, "mean")
    g = sub.groupby(["source", "method", "scale"], dropna=False).agg(**agg).reset_index(); g.to_csv(out / "control_summary.csv", index=False); print("\n" + g.assign(scale=g.scale.map(lambda v: f"{v:.0e}" if pd.notna(v) and np.isfinite(v) else "baseline")).round(4).to_string())
    win = None
    if len(rb):
        ro = rb[rb.scenario == "oracle"].sort_values("shift_s"); ok0 = ro[ro.shift_s == 0]
        if len(ok0) and bool(ok0.blocked.iloc[0]):
            lo = hi = 0.0; sh = list(ro.shift_s); bl = dict(zip(ro.shift_s, ro.blocked))
            for q in sorted([v for v in sh if v < 0], reverse=True):
                if bl[q]: lo = q
                else: break
            for q in sorted([v for v in sh if v > 0]):
                if bl[q]: hi = q
                else: break
            win = [lo, hi]; print(f"plan blocks the wave when the real wave is shifted by {lo:g} .. +{hi:g} s (planned +-{S:g} s)")
    if len(rb): rb.to_csv(out / "timing_robustness_all.csv", index=False); rb[rb.scenario == "oracle"].drop(columns="scenario").to_csv(out / "timing_robustness.csv", index=False)
    if len(mm): mm.to_csv(out / "control_mismatch.csv", index=False)
    if not a.no_mismatch:
        pd.DataFrame([dict(n_plant_cases=len(ev_cfg.get("plant_cases", [])), n_self_extinguishing=self_extinguished)]).to_csv(out / "plant_case_summary.csv", index=False)
    pd.DataFrame(strength_rows).to_csv(out / "strength_margin.csv", index=False)
    if tradeoff: pd.DataFrame(tradeoff).to_csv(out / "tradeoff.csv", index=False)
    # ---- figures
    # Keep the A-method outcome visible even when it produces no usable
    # estimate: it is explicitly a no-detection, not an omitted scenario.
    has_a = any(q[1] == "rec" and q[2].startswith("A_") and q[6] != "ok" for q in scen)
    sel = ["none", "no_sensing", "oracle"] + (["no_detection_A"] if has_a else []) + [n for n in names if n.startswith("rec:B_model")][:1] + [n for n in names if n.startswith("event")][:1]; fig, axs = plt.subplots(1, len(sel), figsize=(3.6 * len(sel), 3.8), squeeze=False)
    for ax, n in zip(axs[0], sel):
        if n in ("none", "no_detection_A"): U, P = Upl, None
        else: P = plans[n]; U = simulate(plant, grid, u0, prot, make_ctrl(P, d), clock_offset=onset_shift, p=d.p, store=True, dt_target=dt_final)["u"][0]
        im = ax.imshow(U.T, origin="lower", aspect="auto", extent=[t[0], t[-1], x[0], x[-1]], vmin=-0.3, vmax=1); ax.axhspan(d.prot_lo, d.prot_hi, fc="none", ec="w", ls="--", lw=.8)
        if P is not None: ax.add_patch(Rectangle((P[1] - onset_shift, P[0] - d.width / 2), P[2], d.width, fc="none", ec="r", lw=1.5))
        label = "A method: no detection" if n == "no_detection_A" else n
        ax.set_title(f"{label}\n(n=1)", fontsize=8); ax.set_xlabel("time (s)")
    axs[0][0].set_ylabel("x (mm)"); fig.colorbar(im, ax=axs[0].tolist(), label="u"); fig.savefig(out / "fig_control_spacetime.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    if len(rb):
        fig, ax = plt.subplots(figsize=(5.5, 3.8)); rb["grp"] = rb.scenario.map(lambda n: n if n in ("oracle", "no_sensing") else ":".join(n.split(":")[:2]) if n.startswith("rec") else "event")
        for gname, gg in rb.groupby("grp"): fr = gg.groupby("shift_s").blocked.mean(); ax.plot(fr.index, fr.values, "o-", label=f"{gname}, pooled across scales (n={gg.scenario.nunique()})")
        ax.axvspan(-S, S, alpha=.12, label=f"planned tolerance +-{S:g} s"); ax.set_xlabel("extra timing shift of the real wave (s)"); ax.set_ylabel("fraction of plans that block the wave"); ax.set_ylim(-.05, 1.05); ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(1.01, 1))
        ax.set_title("Timing robustness (lines pooled across scales)")
        fig.savefig(out / "fig_control_timing_robustness.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    cats = [(f"{r.source}:{r.method}:{r.scale:g} (n={r.n})" if np.isfinite(r.scale) else f"{r.source} (n={r.n})", r) for r in g.itertuples()]
    if cats:
        fig, ax = plt.subplots(figsize=(1.0 * len(cats) + 2, 3.6)); ax.bar(range(len(cats)), [r.blocked_frac for _, r in cats]); ax.set_xticks(range(len(cats))); ax.set_xticklabels([c for c, _ in cats], rotation=40, ha="right", fontsize=7)
        ax.set_ylim(0, 1.05); ax.set_ylabel("fraction of plans that block the wave"); fig.savefig(out / "fig_control_blocked.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    # ---- optional 2D confirmation with Member 3's solver (the stripe is uniform in y)
    conf = None
    if a.confirm_2d and m3_simulate is not None:
        c2 = m3_load(str(p_model)); c2["model"]["dimension"] = 2; b2 = m3_simulate(c2); k2 = m3_simulate(c2, control=Stripe(P0, d)); ub, uc = np.asarray(b2["u"]), np.asarray(k2["u"]); px = np.asarray(b2["x"])
        m2 = (px >= d.prot_lo) & (px <= d.prot_hi); conf = dict(peak_ratio_protected_2d=float(uc[:, :, m2].max() / ub[:, :, m2].max()), min_over_y_ratio=float(np.min([uc[:, j, m2].max() for j in range(uc.shape[1])]) / ub[:, :, m2].max()),
              max_over_y_ratio=float(np.max([uc[:, j, m2].max() for j in range(uc.shape[1])]) / ub[:, :, m2].max()), blocked_2d=bool(uc[:, :, m2].max() / ub[:, :, m2].max() < d.block_fraction))
        print("2D confirmation (Member 3 solver, stripe control):", conf)
    # ---- provenance
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16] if p and Path(p).is_file() else "n/a"
    def git(*c):
        try: r = subprocess.run(["git", "-C", str(ROOT), *c], capture_output=True, text=True); return r.stdout.strip() if r.returncode == 0 else "n/a"
        except OSError: return "n/a"
    status = git("status", "--porcelain", "--", ".", ":(exclude)outputs", ":(exclude)results")
    if status not in ("", "n/a"): notes.append("git working tree is NOT clean: results are not reproducible from the commit alone")
    notes += [f"RETROSPECTIVE sensing: reconstruction and event estimates use data through t_obs={d.t_obs:g} s only; event scenarios use the arrival at target x", f"one open-loop plan tolerates about ±{S:g} s of timing error and ±10% speed error; plant cases {ev_cfg.get('plant_cases')}",
              "control is a STRIPE uniform in y (the simulated wave is uniform in y); a localised 2D Gaussian cannot stop a planar wave"]
    shutil.copy(rooted(a.config), out / "control_used.yaml")
    json.dump(dict(args=vars(a), notes=notes, solver_check_relL2=verify, initial_state=u0src, timing_window_s=win, n_failed_fits=ev.get("n_failed") if ev else None, plan=dict(xc=P0[0], ts=P0[1], dur=P0[2], g0=P0[3], dose_g0_s=P0[2] * P0[3], g0_feasible_range=R["g0_feasible_range"], halfwidth_s=S, J_plan=R["J_plan"], n_seeds=R["n_seeds"]), plant_onset_shift_s=onset_shift, strength_margin=strength_rows, self_extinguishing_plant_cases=self_extinguished if not a.no_mismatch else None, confirm_2d=conf, seconds_total=time.time() - T0,
                   inputs=dict(model=str(p_model), reconstruction=str(rooted(a.reconstruction)), event=str(p_evt)), hashes=dict(model=sha(p_model), reconstruction=sha(rooted(a.reconstruction)), event=sha(p_evt)),
                   versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, pandas=pd.__version__, platform=platform.platform()), git=dict(commit=git("rev-parse", "HEAD"), status=status),
                   hypothetical="control is a hypothetical numerical study; no clinical claim"), open(out / "control_metadata.json", "w", encoding="utf-8"), indent=2, default=str)
    print(f"ALL DONE in {time.time() - T0:.0f} s -> {out}"); return out


if __name__ == "__main__":
    main()
