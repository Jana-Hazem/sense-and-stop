"""Member 6B - hypothetical wave-control study.   Run from anywhere:

    python -m src.control.run_control                   (all paths from configs/control.yaml)
    python src/control/run_control.py --quick           (small search, 1 rec-seed: smoke test)
    python src/control/run_control.py --no-event        (skip event_estimate.json scenarios)

Scenarios (all scored on the PLANT with J = damage_ratio + w*effort_norm):
  none                   no control, J = 1 by construction
  oracle                 idealised benchmark: planner knows the plant and exact timing
  rec:<method>:<scale>   planner timing from Member 5's reconstruction (reconstruction.npz)
  event:<scale>:<i>      planner timing from Member 5's arrival estimate (event_estimate.json)

PRIMARY REPORTING METRICS: blocked (peak_ratio < block_fraction) and damage_ratio.
J itself should not be used to rank plans when the wave is fully blocked -- tiny effort
differences dominate and "negative regret" (a rec plan beats oracle) is a numerical artefact,
not a real gain.

KNOWN LIMITATIONS (fully noted in control_metadata.json::notes)
  1. Causality: reconstructions use the whole 540 s record; a real controller at t_obs would
     see only 0-t_obs data.  Results are optimistic.
  2. Inverse crime: plant == nominal by default (D/k_factor = 1).  Addressed partially by
     the plant_variants mismatch grid in robustness.plant_variants.
  3. Protected region (20-28 mm): wave enters at ~393 s and only partly crosses before T=540 s.
"""
import argparse, glob, hashlib, json, platform, shutil, subprocess, sys, time, zlib
from pathlib import Path
import numpy as np, pandas as pd, yaml, scipy
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
from src.control.model import Grid, load_model, initial_state, simulate
from src.control.objective import design_from_cfg, make_ctrl, cost, search_bounds
from src.control.estimate import make_nominal, peak_pos, arrival_time, delta_from_rec, delta_from_arrival
from src.control.events import load_events
from src.control.search import optimise


def _is_scale(s):
    """Return True if s looks like a floating-point scale (e.g. '1e-09', '0.001')."""
    try: float(s); return True
    except ValueError: return False


def parse_args(argv=None):
    ap = argparse.ArgumentParser(); a = ap.add_argument
    a("--config", default=str(ROOT / "configs" / "control.yaml"))
    a("--root", default=str(ROOT))
    a("--model-config", default=None)
    a("--reconstruction", default=None)
    a("--event-estimate", default=None)
    a("--out", default=None)
    a("--methods", nargs="*", default=None, help="reconstruction methods to include (default: all rec_*)")
    a("--event-seeds", type=int, default=None, help="max arrival estimates per scale")
    a("--no-event", action="store_true")
    a("--quick", action="store_true", help="small search + 1 rec-seed: smoke test only")
    return ap.parse_args(argv)


def resolve(root, p):
    """Resolve one path or a list of glob patterns relative to root; most-recently-modified wins."""
    pats = [p] if isinstance(p, (str, Path)) else list(p)
    for pat in pats:
        q = Path(pat) if Path(pat).is_absolute() else Path(root) / pat
        hits = [Path(h) for h in glob.glob(str(q), recursive=True) if Path(h).is_file()]
        if hits: return max(hits, key=lambda h: h.stat().st_mtime)
    near = sorted({str(h) for h in glob.glob(str(Path(root) / "outputs" / "**" / "*.np[zy]"), recursive=True)})[:12]
    raise SystemExit("ERROR: file not found. Tried: " + "; ".join(str(Path(root) / q) for q in pats)
                     + ("\n  .npz/.npy files under outputs/: " + ", ".join(near) if near else
                        "\n  nothing under outputs/ yet: run Member 5's script first")
                     + "\n  or pass explicitly: --reconstruction path/to/file.npz")


def seed_offset(name): return zlib.crc32(name.encode("utf-8")) % 10000


def main(argv=None):  # noqa: C901
    a = parse_args(argv)
    root = Path(a.root).resolve()
    cfg  = yaml.safe_load(Path(a.config).read_text(encoding="utf-8"))
    s    = dict(cfg["search"])
    est  = cfg["estimation"]
    evc  = cfg.get("event", {})
    rb_cfg = cfg.get("robustness", {})

    # -- Quick-mode overrides ----------------------------------------------------
    if a.quick:
        s.update(popsize=4, maxiter=8, polish_maxiter=20, n_rec_seeds=1)
        print("[quick mode] popsize=4 maxiter=8 polish=20 n_rec_seeds=1")

    n_rec_seeds = int(s.get("n_rec_seeds", 1))

    # -- Resolve input files -----------------------------------------------------
    p_model = resolve(root, a.model_config or cfg["inputs"]["model_config"])
    p_rec   = resolve(root, a.reconstruction or cfg["inputs"]["reconstruction"])
    notes   = []
    p_evt   = None
    ev      = None
    if not a.no_event:
        try:
            p_evt = resolve(root, a.event_estimate or cfg["inputs"]["event_estimate"])
            ev    = load_events(p_evt, str(evc.get("method_key", "B")), float(evc.get("target_x_mm", 15.0)))
            for w in ev["warnings"]: notes.append("event_estimate: " + w); print("[WARNING] event_estimate:", w)
            if not ev["per_scale"]:
                notes.append("event_estimate.json found but no usable arrival estimates")
                print("[WARNING] no usable arrival estimates in", p_evt)
        except SystemExit as e:
            print(f"[note] event scenarios skipped ({e})")
            notes.append(f"event scenarios skipped: {e}")

    # -- Models and grid ---------------------------------------------------------
    nominal_model, mcfg = load_model(p_model)
    plant_model = nominal_model.scaled(float(cfg["plant"]["D_factor"]), float(cfg["plant"]["k_factor"]))
    z = np.load(p_rec, allow_pickle=True)
    _truth_key = "u_truth" if "u_truth" in z.files else ("truth" if "truth" in z.files else None)
    if _truth_key is None:
        raise SystemExit(f"ERROR: reconstruction.npz has neither 'u_truth' nor 'truth'. Keys: {list(z.files)}")
    t, x, truth = np.asarray(z["t"], float), np.asarray(z["x"], float), np.asarray(z[_truth_key], float)
    grid  = Grid(x, t)
    u0, u0src = initial_state(mcfg, grid, truth[0])
    d     = design_from_cfg(cfg, float(t[-1]))
    prot  = (x >= d.prot_lo) & (x <= d.prot_hi)
    if not prot.any(): raise SystemExit("ERROR: protected region is outside the domain")
    dt_final = float(mcfg.get("time", {}).get("dt", 0.1))
    out_dir  = (Path(a.out) if a.out else root / cfg["output"]["dir"]) / time.strftime("run_%Y%m%d_%H%M%S")
    out_dir.mkdir(parents=True, exist_ok=True)

    # -- Reference runs and solver verification ----------------------------------
    _ref = lambda m: simulate(m, grid, u0, prot, None, p=d.p, store=True, dt_target=dt_final)
    nm_run, pl_run = _ref(nominal_model), _ref(plant_model)
    Unm, Upl = nm_run["u"][0], pl_run["u"][0]
    nom    = make_nominal(Unm, x, t)
    pos_pl = peak_pos(Upl, x)
    m_     = (t >= nom.tt[0])[:, None] & np.ones_like(truth, bool)
    verify = float(np.linalg.norm((Unm - truth)[m_]) / np.linalg.norm(truth[m_]))
    print(f"solver check: relL2(nominal vs u_truth) = {verify:.2e} | initial state: {u0src}")
    if verify > 0.05:
        notes.append(f"WARNING solver relL2 = {verify:.2e} > 0.05: solver or initial state differs from Member 3")

    D0_nom, D0_pl = float(nm_run["damage"][0]), float(pl_run["damage"][0])
    if D0_pl < 1e-8: raise SystemExit("ERROR: uncontrolled wave never reaches protected region")
    peak_pl = float(pl_run["peak"][0])

    tx       = ev["target_x_mm"] if ev else float(evc.get("target_x_mm", 15.0))
    a_nom    = arrival_time(nom.pos, t, tx)
    a_pl     = arrival_time(pos_pl, t, tx)
    truth_arr = ev["truth_arrival_s"] if ev and ev["truth_arrival_s"] is not None else a_pl
    if ev and not (np.isfinite(a_nom) and np.isfinite(a_pl)):
        msg = f"WARNING front never reaches target x={tx:g} mm in nominal/plant run: event scenarios infeasible"
        notes.append(msg); print("[WARNING]", msg)
    if ev:
        print(f"event layout: {ev['layout']} | target x = {tx:g} mm | truth arrival {truth_arr:.1f} s")
    print(f"protected region {d.prot_lo}-{d.prot_hi} mm | D0_plant = {D0_pl:.3f} | "
          f"arrival @{tx:g} mm: nominal {a_nom:.1f} s, plant {a_pl:.1f} s | effort_ref {d.e_ref:.3g}")

    # Known-limitation notes
    notes.append("CAUSALITY: reconstructions are fitted to the full 540 s record; a real controller "
                 f"at t_obs={d.t_obs:.0f} s would have only 0-{d.t_obs:.0f} s data. Results are optimistic.")
    notes.append(f"INVERSE CRIME: plant D_factor={cfg['plant']['D_factor']}, k_factor={cfg['plant']['k_factor']}. "
                 "Plant == nominal model; mismatch grid (robustness.plant_variants) partially addresses this.")
    arr_prot_lo = arrival_time(pos_pl, t, d.prot_lo)
    notes.append(f"PROTECTED REGION: uncontrolled wave enters at ~{arr_prot_lo:.0f} s "
                 f"and reaches {d.prot_hi} mm only partially before T={t[-1]:.0f} s. "
                 "D0 reflects a partial crossing; controlled runs prevent entry.")

    # -- Scenario helpers --------------------------------------------------------
    rows, trajs, hists = [], {}, {}

    def apply_plan(P, model=None, clock_offset=0.0):
        """Evaluate plan P (shape (4,)) on `model` at timing offset `clock_offset`."""
        m = model if model is not None else plant_model
        P1 = np.asarray(P, float).ravel()          # ensure (4,) -- make_ctrl reshapes to (4,1)
        o = simulate(m, grid, u0, prot, make_ctrl(P1, d),
                     p=d.p, store=True, dt_target=dt_final, clock_offset=clock_offset)
        J, dr, en = cost(o, D0_pl, d)
        return dict(J=float(J[0]), dr=float(dr[0]), en=float(en[0]),
                    pk=float(o["peak"][0] / peak_pl), u=o["u"][0])


    def none_row(name, source, method, scale, seed, status, extra=None):
        rows.append(dict(scenario=name, source=source, method=method, scale=scale, seed=seed,
                         status=status, J_true=1.0, damage_ratio=1.0, effort_norm=0.0,
                         peak_ratio=1.0, blocked=False, **(extra or {})))

    def run_scenario(name, source, method, scale, seed, planner_model, delta, D0_plan, p_front_obs, extra=None):
        base = dict(scenario=name, source=source, method=method, scale=scale, seed=seed,
                    delta_s=delta, p_front_obs=p_front_obs, **(extra or {}))
        bnd = search_bounds(d, p_front_obs) if np.isfinite(delta) and np.isfinite(p_front_obs) else None
        if bnd is None:
            rows.append({**base, "status": "infeasible", "J_true": 1.0, "damage_ratio": 1.0,
                         "effort_norm": 0.0, "peak_ratio": 1.0, "blocked": False})
            return
        t0 = time.time()
        r  = optimise(planner_model, grid, u0, prot, d, bnd, delta, D0_plan, s, seed_offset=seed_offset(name))
        P  = r["P"]
        ap = apply_plan(P)
        elapsed = time.time() - t0
        blocked = bool(ap["pk"] < d.block_fraction)
        rows.append({**base, "status": "ok", "xc": P[0], "ts": P[1], "dur": P[2], "g0": P[3],
                     "J_plan": r["J_plan"], "J_true": ap["J"], "damage_ratio": ap["dr"],
                     "effort_norm": ap["en"], "peak_ratio": ap["pk"], "blocked": blocked,
                     "nfev_batches": r["nfev_batches"]})
        trajs[name] = (ap["u"], P)
        hists[name] = r["history"]
        status_flag = "[BLOCKED]" if blocked else "[not blocked]"
        print(f"  {name:<36} delta={delta:+.1f}s | xc={P[0]:.2f} ts={P[1]:.1f} dur={P[2]:.1f} "
              f"g0={P[3]:.4f} | dr={ap['dr']:.4f} {status_flag} [{elapsed:.0f}s]", flush=True)

    # -- Scenario loop -----------------------------------------------------------
    none_row("none", "none", None, np.nan, -1, "ok")
    trajs["none"] = (Upl, None)

    # Oracle: plant model, perfect timing
    run_scenario("oracle", "oracle", "full_state", np.nan, -1,
                 plant_model, 0.0, D0_pl, float(np.interp(d.t_obs, t, pos_pl)))

    # Reconstruction-based scenarios: n_rec_seeds optimizer seeds per (method, scale)
    for key in [k for k in z.files if k.startswith("rec_")]:
        suffix = key[4:]
        parts  = suffix.rsplit("_", 1)
        if len(parts) == 2 and _is_scale(parts[1]):
            method, scale_str = parts
        else:
            method, scale_str = "unknown", suffix
        if a.methods and method not in a.methods: continue
        try: scale_f = float(scale_str)
        except ValueError: scale_f = float("nan")
        r_est = delta_from_rec(np.asarray(z[key], float), x, t, nom, d.t_obs,
                               float(est["detect_fraction"]), float(est["front_window_s"]),
                               int(est["min_points"]))
        if r_est is None:
            none_row(f"rec:{method}:{scale_str}", "rec", method, scale_f, 0, "no_detection")
            print(f"  rec:{method}:{scale_str:<24} no front detected -> no plan (J = 1)")
            continue
        for si in range(n_rec_seeds):
            name_si = f"rec:{method}:{scale_str}" if n_rec_seeds == 1 else f"rec:{method}:{scale_str}:s{si}"
            run_scenario(name_si, "rec", method, scale_f, si,
                         nominal_model, r_est["delta"], D0_nom,
                         nom.at(d.t_obs - r_est["delta"]),
                         dict(p_hat_est=r_est["p_hat"]))

    # Event-based scenarios
    if ev:
        for sc, vals in ev["per_scale"].items():
            try: scf = float(sc)
            except ValueError: scf = np.nan
            max_ev = a.event_seeds or (2 if a.quick else None)
            for i, arr in enumerate(vals[:max_ev]):
                name = f"event:{sc}:{i}"
                if not np.isfinite(arr):
                    none_row(name, "event", "B_model", scf, i, "no_detection")
                    print(f"  {name:<36} arrival NaN -> no plan"); continue
                dl = delta_from_arrival(arr, truth_arr, a_pl, a_nom)
                run_scenario(name, "event", "B_model", scf, i,
                             nominal_model, dl, D0_nom,
                             nom.at(d.t_obs - dl),
                             dict(arrival_err_s=float(arr - truth_arr)))

    # -- Colour palette (shared by all subsequent sections) ----------------------
    SRC_COLOR = {"none": "#888", "oracle": "#1a9641", "rec": "#2171b5", "event": "#d94801"}

    # -- Timing + mismatch robustness: every successful plan ---------------------
    shifts        = [float(sh) for sh in rb_cfg.get("timing_shifts_s", [0])]
    pv_cfgs       = rb_cfg.get("plant_variants", [{"D_factor": 1.0, "k_factor": 1.0}])
    pv_models     = [(pv, nominal_model.scaled(float(pv["D_factor"]), float(pv["k_factor"])))
                     for pv in pv_cfgs]
    pv_D0         = {}   # D0 for each variant (computed lazily)

    def _pv_key(pv): return "D%s_k%s" % (pv['D_factor'], pv['k_factor'])

    # Compute uncontrolled D0 for every plant variant at search speed (run once per variant)
    dt_rb = float(s["dt_target_s"])          # robustness grid uses search-speed RK4 (4x faster than dt_final)
    for pv, pv_model in pv_models:
        pk = _pv_key(pv)
        ref_pv = simulate(pv_model, grid, u0, prot, None, p=d.p, store=False, dt_target=dt_rb)
        pv_D0[pk] = float(ref_pv["damage"][0])
        if pv_D0[pk] < 1e-8: pv_D0[pk] = D0_pl  # fallback if variant wave misses prot region

    rb_rows, mm_rows = [], []
    ok_rows = [r for r in rows if r.get("status") == "ok" and r.get("scenario") != "none"]
    n_rb_evals = len(ok_rows) * (len(shifts) + len(pv_models))
    print("\n-- Robustness grid (%d plans x %d timing shifts + %d plant variants = ~%d evals at dt=%.1fs) --"
          % (len(ok_rows), len(shifts), len(pv_models), n_rb_evals, dt_rb))
    for row in ok_rows:
        P = np.array([row["xc"], row["ts"], row["dur"], row["g0"]])
        sc_name = row["scenario"]

        # Timing robustness on nominal plant (search-speed RK4)
        for sh in shifts:
            o = simulate(plant_model, grid, u0, prot, make_ctrl(P, d),
                         clock_offset=float(sh), p=d.p, store=False, dt_target=dt_rb)
            J, dr, en = cost(o, D0_pl, d)
            rb_rows.append(dict(scenario=sc_name, source=row["source"], method=row.get("method"),
                                scale=row.get("scale"), seed=row.get("seed"),
                                shift_s=sh, J_plant=float(J[0]), damage_ratio=float(dr[0]),
                                blocked=bool(o["peak"][0] / peak_pl < d.block_fraction)))

        # Model-mismatch robustness (search-speed RK4)
        for pv, pv_model in pv_models:
            pk = _pv_key(pv)
            o = simulate(pv_model, grid, u0, prot, make_ctrl(P, d),
                         p=d.p, store=False, dt_target=dt_rb)
            J, dr, en = cost(o, pv_D0[pk], d)
            mm_rows.append(dict(scenario=sc_name, source=row["source"], method=row.get("method"),
                                scale=row.get("scale"), seed=row.get("seed"),
                                D_factor=pv["D_factor"], k_factor=pv["k_factor"],
                                J_plant=float(J[0]), damage_ratio=float(dr[0]),
                                blocked=bool(o["peak"][0] / peak_pl < d.block_fraction)))

    # -- Build DataFrames and CSVs -----------------------------------------------
    df = pd.DataFrame(rows)
    jo = df.loc[df.scenario == "oracle", "J_true"]
    if len(jo): df["regret_vs_oracle"] = df.J_true - float(jo.iloc[0])
    df.to_csv(out_dir / "control_results.csv", index=False)

    sub = df[df.source.isin(["rec", "event"])].copy()
    if len(sub):
        g = (sub.groupby(["source", "method", "scale"])
             .agg(n=("J_true", "size"),
                  blocked_frac=("blocked", "mean"),
                  damage_ratio_med=("damage_ratio", "median"),
                  damage_ratio_max=("damage_ratio", "max"),
                  J_med=("J_true", "median"),
                  J_q25=("J_true", lambda v: v.quantile(.25)),
                  J_q75=("J_true", lambda v: v.quantile(.75)),
                  no_plan=("status", lambda v: float((v != "ok").mean())))
             .reset_index())
        # Add n_reconstructions column: rec rows are optimiser restarts of ONE reconstruction.
        # Only event rows use genuinely different M5 estimates per seed.
        g["n_reconstructions"] = g["source"].apply(lambda src: 1 if src == "rec" else None)
        g["note"] = g["source"].apply(
            lambda src: "n = optimiser restarts of 1 reconstruction (not independent noise draws)"
            if src == "rec" else "n = independent M5 arrival estimates (true replicates)")
        g.to_csv(out_dir / "control_summary.csv", index=False)
        # Format scale in scientific notation for display (round(4) would show 0.0 for 1e-09)
        g_disp = g[["source","method","scale","n","n_reconstructions",
                    "blocked_frac","damage_ratio_med","damage_ratio_max","no_plan"]].copy()
        g_disp["scale"] = g_disp["scale"].apply(lambda s_: "%g" % s_ if pd.notna(s_) else "nan")
        print("\n-- Summary by (source, method, scale) --")
        print("  NOTE: rec n = optimiser restarts of 1 reconstruction; event n = independent M5 estimates")
        print(g_disp.to_string(index=False))

    if rb_rows:
        rb_df = pd.DataFrame(rb_rows)
        rb_df.to_csv(out_dir / "plan_robustness.csv", index=False)
        # timing_robustness.csv: oracle-only subset (backward-compatible interface contract)
        orb = rb_df[rb_df.scenario == "oracle"]
        if len(orb):
            orb[["shift_s","damage_ratio","blocked"]].to_csv(
                out_dir / "timing_robustness.csv", index=False)
            print("\n-- Oracle plan timing robustness --")
            print(orb[["shift_s","damage_ratio","blocked"]].round(4).to_string(index=False))

    if mm_rows:
        mm_df = pd.DataFrame(mm_rows)
        mm_df.to_csv(out_dir / "mismatch_robustness.csv", index=False)
        # mismatch_summary: mean blocked over variants (not worst-case — "worst" was mislabelled)
        # blocked_all_variants = True only if plan blocks under every D/k variant tested
        wc = (mm_df.groupby("scenario")
              .agg(blocked_frac_over_variants=("blocked", "mean"),
                   blocked_all_variants=("blocked", lambda v: bool(v.all())),
                   damage_ratio_worst=("damage_ratio", "max"))
              .reset_index())
        wc.to_csv(out_dir / "mismatch_summary.csv", index=False)
        print("\n-- Mismatch robustness: blocked_frac_over_variants (D/k variants)")
        print("  blocked_all_variants = True only when every D/k variant is blocked")
        print(wc.round(4).to_string(index=False))

    # ==========================================================================
    # FIGURES -- publication-quality
    # ==========================================================================
    CMAP = "RdBu_r"

    # -- Fig 1: space-time maps -------------------------------------------------
    names_st = ([n for n in ("none", "oracle") if n in trajs]
                + [n for n in trajs if n.startswith("rec:")][:2]
                + [n for n in trajs if n.startswith("event:")][:1])
    ncols = max(len(names_st), 1)
    fig_st, axs_st = plt.subplots(1, ncols, figsize=(3.5 * ncols, 4.2), squeeze=False)
    fig_st.suptitle("Hypothetical CSD wave control -- space-time (u field)\n"
                    "[HYPOTHETICAL: no clinical claim]", fontsize=8, y=1.02)
    for ax, n in zip(axs_st[0], names_st):
        U, P = trajs[n]
        im = ax.imshow(U.T, origin="lower", aspect="auto",
                       extent=[t[0], t[-1], x[0], x[-1]], cmap=CMAP, vmin=-0.3, vmax=1.0)
        ax.axhline(d.prot_lo, c="white", ls="--", lw=0.8)
        ax.axhline(d.prot_hi, c="white", ls="--", lw=0.8)
        # Uncontrolled front trajectory (white dotted)
        ax.plot(t, peak_pos(trajs["none"][0], x), "w:", lw=0.9, alpha=0.5)
        if P is not None:
            ax.add_patch(Rectangle((P[1], P[0] - d.width / 2), P[2], d.width,
                                   fc="none", ec="#ff4444", lw=1.8, zorder=5))
            ax.text(P[1] + P[2] / 2, P[0] + d.width / 2 + 0.5, "ctrl",
                    ha="center", va="bottom", fontsize=5, color="#ff4444")
        src = n.split(":")[0]
        ax.set_title(n, fontsize=7, color=SRC_COLOR.get(src, "k"), pad=3)
        ax.set_xlabel("time (s)", fontsize=8); ax.tick_params(labelsize=7)
    axs_st[0][0].set_ylabel("x (mm)", fontsize=8)
    cb = fig_st.colorbar(im, ax=axs_st[0].tolist(), shrink=0.85, pad=0.02)
    cb.set_label("u (dimensionless)", fontsize=8); cb.ax.tick_params(labelsize=7)
    fig_st.savefig(out_dir / "fig_control_spacetime.png", dpi=160, bbox_inches="tight")
    plt.close(fig_st)

    # -- Fig 2: damage_ratio scatter (PRIMARY metric) ---------------------------
    rr = df[df.source == "rec"].copy()
    ee_ = df[df.source == "event"].copy()
    # Group: for rec with multiple seeds, group by (method, scale)
    _cats_rec = sorted({(m, sc) for m, sc in zip(rr.method, rr.scale) if pd.notna(sc)},
                       key=lambda ms: (ms[0], ms[1]))
    _cats_ev  = sorted({sc for sc in ee_.scale.dropna()})
    cats = ["none", "oracle"] + [f"rec:{m}:{sc:g}" for m, sc in _cats_rec] + [f"event:{sc:g}" for sc in _cats_ev]

    fig_c, ax_c = plt.subplots(figsize=(max(7, 1.5 * len(cats)), 4.8))
    rng = np.random.default_rng(42)
    for i, c in enumerate(cats):
        src = c.split(":")[0]
        col = SRC_COLOR.get(src, "#555")
        if c in ("none", "oracle"):
            v = df[df.scenario == c].damage_ratio.values
            bk = df[df.scenario == c].blocked.values
        elif c.startswith("rec:"):
            m_p, sc_p = c.split(":")[1], c.split(":")[2]
            mask = (rr.method == m_p) & rr.scale.apply(
                lambda s_: np.isclose(s_, float(sc_p)) if pd.notna(s_) else False)
            v  = rr[mask].damage_ratio.values
            bk = rr[mask].blocked.values
        else:
            sc_p = float(c.split(":")[1])
            m2 = ee_.scale.apply(lambda s_: np.isclose(s_, sc_p) if pd.notna(s_) else False)
            v  = ee_[m2].damage_ratio.values
            bk = ee_[m2].blocked.values
        jit = rng.uniform(-0.2, 0.2, len(v))
        # Distinguish blocked (filled) from not-blocked (open)
        ax_c.scatter(i + jit[bk], v[bk], s=28, color=col, alpha=0.9, zorder=3, marker="o")
        ax_c.scatter(i + jit[~bk], v[~bk], s=28, color=col, alpha=0.9, zorder=3,
                     marker="o", facecolors="none", linewidths=1.5)
        if len(v):
            med = float(np.median(v))
            ax_c.hlines(med, i - 0.3, i + 0.3, colors=col, lw=2.5, zorder=4)
            ax_c.text(i, med + 0.015, f"{med:.3f}", ha="center", va="bottom", fontsize=5.5, color=col)

    ax_c.axhline(1, c="#aaa", ls="--", lw=1.0)
    ax_c.axhline(0, c="#ccc", ls=":", lw=0.7)
    ax_c.set_xticks(range(len(cats)))
    ax_c.set_xticklabels(cats, rotation=40, ha="right", fontsize=6)
    ax_c.set_ylabel("damage_ratio  D(P)/D0  on plant\n(0 = fully blocked, 1 = no improvement)", fontsize=8)
    ax_c.set_title("Hypothetical control comparison -- damage_ratio (PRIMARY metric)\n"
                   "filled=blocked  open=not blocked  bar=median  [HYPOTHETICAL]", fontsize=8)
    from matplotlib.lines import Line2D
    leg = ([Line2D([0],[0], marker="o", color=v, label=k, ls="", ms=7)
            for k, v in SRC_COLOR.items() if k in {"none","oracle","rec","event"}]
           + [Line2D([0],[0], marker="o", color="k", ls="", ms=7, label="blocked (filled)"),
              Line2D([0],[0], marker="o", color="k", ls="", ms=7, label="not blocked (open)",
                     mfc="none", mew=1.5)])
    ax_c.legend(handles=leg, fontsize=6, loc="upper right", ncol=2)
    ax_c.tick_params(labelsize=7)
    fig_c.savefig(out_dir / "fig_control_damage_ratio.png", dpi=160, bbox_inches="tight")
    plt.close(fig_c)

    # -- Fig 3: optimiser convergence -------------------------------------------
    if hists:
        fig_h, ax_h = plt.subplots(figsize=(6.5, 3.8))
        for ni, (n, h) in enumerate(list(hists.items())[:10]):
            src = n.split(":")[0]; col = SRC_COLOR.get(src, "#555")
            ax_h.plot(range(1, len(h)+1), h, color=col, alpha=0.75,
                      lw=1.8 if ni < 2 else 1.0, label=n)
        ax_h.axhline(1, c="#aaa", ls="--", lw=0.8, label="no control")
        ax_h.set_xlabel("Batch evaluation (each = vectorised population x 4 params)", fontsize=8)
        ax_h.set_ylabel("Best planning J so far", fontsize=8)
        ax_h.set_title("Optimiser convergence (DE + Nelder-Mead polish)\n"
                        f"popsize={s['popsize']} maxiter={s['maxiter']} polish={s['polish_maxiter']}", fontsize=8)
        ax_h.legend(fontsize=5.5, ncol=2); ax_h.tick_params(labelsize=7)
        fig_h.savefig(out_dir / "fig_search_history.png", dpi=160, bbox_inches="tight")
        plt.close(fig_h)

    # -- Fig 4: timing robustness -- ALL plans -----------------------------------
    if rb_rows:
        rb_df = pd.DataFrame(rb_rows)
        # Actual delta range from this run (not just event-estimate errors)
        all_deltas = [row.get("delta_s", 0.0) for row in rows
                      if row.get("status") == "ok" and row.get("scenario") != "none"
                      and np.isfinite(row.get("delta_s", np.nan))]
        delta_lo = min(all_deltas) if all_deltas else 0.0
        delta_hi = max(all_deltas) if all_deltas else 0.0

        fig_rb, ax_rb = plt.subplots(figsize=(7, 4.0))
        # Plot each plan's robustness curve, colour by source
        for sc_name in rb_df.scenario.unique():
            sub_rb = rb_df[rb_df.scenario == sc_name].sort_values("shift_s")
            src = sc_name.split(":")[0]
            col = SRC_COLOR.get(src, "#555")
            lw  = 2.0 if sc_name == "oracle" else 0.8
            ax_rb.plot(sub_rb.shift_s, sub_rb.damage_ratio, color=col, lw=lw,
                       alpha=0.9 if sc_name == "oracle" else 0.45,
                       label=sc_name if sc_name == "oracle" else "_nolegend_")
            # Star where blocked
            blk = sub_rb[sub_rb.blocked]
            if len(blk):
                ax_rb.scatter(blk.shift_s, blk.damage_ratio, marker="*", s=50,
                              color=col, zorder=5, alpha=0.6)

        ax_rb.axhline(1, c="#aaa", ls="--", lw=1.0, label="no control")
        ax_rb.axhline(0, c="#ccc", ls=":", lw=0.7)
        # Actual delta range band
        if delta_lo < delta_hi or abs(delta_lo) + abs(delta_hi) > 0.5:
            ax_rb.axvspan(delta_lo, delta_hi, alpha=0.12, color="#2171b5",
                          label=f"actual delta range [{delta_lo:.0f},{delta_hi:.0f}] s")
        ax_rb.axvspan(-15, 15, alpha=0.07, color="#888", label="+/-15 s tolerance")
        ax_rb.set_xlabel("Timing shift applied to real wave (s)  [+ve = wave later than planned]", fontsize=8)
        ax_rb.set_ylabel("damage_ratio on plant  (0=blocked, 1=no improvement)", fontsize=8)
        ax_rb.set_title("Plan timing robustness -- all plans tested\n"
                        "thick=oracle, thin=rec/event, star=blocked  [HYPOTHETICAL]", fontsize=8)
        from matplotlib.lines import Line2D as _L2D
        leg2 = [_L2D([0],[0], color=v, label=k, lw=2) for k, v in SRC_COLOR.items()]
        ax_rb.legend(handles=leg2, fontsize=7, loc="upper left")
        ax_rb.tick_params(labelsize=7)
        fig_rb.savefig(out_dir / "fig_plan_timing_robustness.png", dpi=160, bbox_inches="tight")
        plt.close(fig_rb)

    # -- Fig 5: model-mismatch heatmap -----------------------------------------
    if mm_rows:
        mm_df = pd.DataFrame(mm_rows)
        # Pivot: scenario (rows) x variant (cols) -> damage_ratio
        pv_labels = [f"Dx{pv['D_factor']} kx{pv['k_factor']}" for pv in pv_cfgs]
        ok_scenarios = [r["scenario"] for r in ok_rows]
        heat = np.full((len(ok_scenarios), len(pv_cfgs)), np.nan)
        for ri, sc in enumerate(ok_scenarios):
            for ci, pv in enumerate(pv_cfgs):
                sub_mm = mm_df[(mm_df.scenario == sc) &
                                mm_df.D_factor.apply(lambda d_: np.isclose(d_, pv["D_factor"])) &
                                mm_df.k_factor.apply(lambda k_: np.isclose(k_, pv["k_factor"]))]
                if len(sub_mm): heat[ri, ci] = float(sub_mm.damage_ratio.iloc[0])

        # Only plot if more than 1 variant and more than 1 scenario
        if heat.shape[1] > 1 and heat.shape[0] > 0:
            # Clip row labels for readability
            sc_labels = [s_[:28] for s_ in ok_scenarios]
            fig_mm, ax_mm = plt.subplots(figsize=(max(5, len(pv_cfgs) * 1.6),
                                                   max(3, len(ok_scenarios) * 0.35 + 1.5)))
            im_mm = ax_mm.imshow(heat, aspect="auto", cmap="RdYlGn_r", vmin=0, vmax=1)
            ax_mm.set_xticks(range(len(pv_labels))); ax_mm.set_xticklabels(pv_labels, fontsize=7, rotation=25, ha="right")
            ax_mm.set_yticks(range(len(sc_labels))); ax_mm.set_yticklabels(sc_labels, fontsize=5.5)
            ax_mm.set_title("Model-mismatch robustness: damage_ratio\n"
                             "green=blocked, red=wave survived  [HYPOTHETICAL]", fontsize=8)
            plt.colorbar(im_mm, ax=ax_mm, label="damage_ratio", shrink=0.85)
            fig_mm.savefig(out_dir / "fig_mismatch_robustness.png", dpi=160, bbox_inches="tight")
            plt.close(fig_mm)

    # -- Provenance -------------------------------------------------------------
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16] if p and Path(p).is_file() else "n/a"
    def git(*c):
        try: r = subprocess.run(["git", "-C", str(root), *c], capture_output=True, text=True)
        except OSError: return "n/a"
        return r.stdout.strip() if r.returncode == 0 else "n/a"
    status = git("status", "--porcelain", "--", ".", ":(exclude)outputs")
    if status not in ("", "n/a"):
        notes.append("REPRODUCIBILITY: git working tree is not clean -- commit before the final run")

    shutil.copy(a.config, out_dir / "control_used.yaml")
    shutil.copy(p_model,  out_dir / ("model_" + Path(p_model).name))
    vers = dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                pandas=pd.__version__, matplotlib=matplotlib.__version__, pyyaml=yaml.__version__,
                platform=platform.platform())
    meta = dict(
        args=vars(a), config=cfg, quick=a.quick, notes=notes,
        solver_check_relL2=verify, initial_state=u0src,
        nominal_model=vars(nominal_model), plant_model=vars(plant_model),
        event=(dict(layout=ev["layout"], target_x_mm=tx, truth_arrival_s=truth_arr) if ev else None),
        inputs=dict(model_config=str(p_model), reconstruction=str(p_rec), event_estimate=str(p_evt)),
        hashes=dict(model=sha(p_model), reconstruction=sha(p_rec), event=sha(p_evt)),
        search_params=s,
        versions=vers,
        git=dict(commit=git("rev-parse", "HEAD"), status=status),
        hypothetical="All results are from a hypothetical numerical study. No clinical claim is made.",
    )
    with open(out_dir / "control_metadata.json", "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2, default=str)

    print("ALL DONE ->", out_dir)
    return out_dir


if __name__ == "__main__":
    main()
