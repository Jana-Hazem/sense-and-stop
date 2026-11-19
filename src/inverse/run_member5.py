"""Member 5 - inverse reconstruction of the hidden wave u(x,t) from simulated scalp EEG (1D, Fisher/excitable pulse).

IMPORTABLE MODULE + COMMAND LINE.        python run_member5.py --data <data/derived> [--out <outputs/inverse>] [--quick]
Pipeline use:  inp = load_inputs(args); Y = simulate_observations(inp, scale, noise_window, idx);
               out = reconstruct(inp, Y, idx, fold, scale, method="B_model");  m = evaluate(inp, out["rec"])

Everything is READ from Member 3 / Member 4 files, nothing is typed in:
  u_2d_default.npy (M3 truth [t,y,x]) | flat_patch.npz (dx_mm) | sensors_default.npz (dt_s, scale_Am, clean_V check)
  A_matrix.npz + H_bip.npz or H_matrix.npz (M4 operator) | noise_2s.npz (M4 noise files and channel names)
  model config + solver import (M3) for the pulse template.  noise_cov.npz is NOT read: the covariance is re-estimated
  from noise_2s.npz on the files other than the held-out one, so the whitening always matches the noise files in use.
If Member 4 / 3 change H, A, the noise files or the grid, rerun with no code change.  Missing or ambiguous inputs STOP the run
(pass --allow-fallbacks to continue; every fallback is then written to metadata.json -> notes).

Model: y(t) = scale * G u(t) + eps(t),  G = A^T H^T,  1D operator = G summed over y.  Whitening with the fold covariance.
Methods (inversion never sees the truth, except where a method name says oracle):
  A_smooth / A_transport_fit / A_transport_lit / A_transport_truth(oracle)   Tikhonov, lambda by the discrepancy principle
  B_model  travelling pulse, template from the PDE solver at a DIFFERENT resolution (and optional different config)
  B_gauss  travelling pulse, Gaussian shape
  B has per-channel offset+slope nuisance terms (projected out), amplitude >= 0 is checked, Jacobian-based sd is reported.
Scores: relL2 = main metric on the settled part (t >= late_t0, x <= x_max); relL2_fit on t >= fit_t0.
"""
import argparse, json, subprocess, sys, time, hashlib, shutil, platform
from pathlib import Path
from types import SimpleNamespace as NS
import numpy as np, pandas as pd, scipy, scipy.sparse as sp
from scipy.sparse.linalg import splu
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import least_squares
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2] if len(Path(__file__).resolve().parents) > 2 else Path.cwd()

def parse_args(argv=None):
    ap = argparse.ArgumentParser()
    a = ap.add_argument
    a("--data", default=str(ROOT / "data" / "derived")); a("--out", default=str(ROOT / "outputs" / "inverse"))
    a("--u2d", default=None, help="truth file (default <data>/u_2d_default.npy)")
    a("--config", default=str(ROOT / "configs" / "model_default.yaml"), help="M3 model config")
    a("--template-config", default=None, help="different M3 config for the template (default: --config)")
    a("--template-dx-factor", type=float, default=2.0, help="template solved at dx * factor (avoids the inverse crime)")
    a("--template-method", default="solve_ivp"); a("--template-at", type=float, default=0.5, help="template snapshot: peak at this fraction of its domain")
    a("--gen-data", default=None, help="folder with perturbed A/H(/noise) used ONLY to generate the data")
    a("--allow-fallbacks", action="store_true"); a("--allow-overlap", action="store_true")
    a("--seeds", type=int, default=None, help="default: all disjoint noise windows"); a("--quick", action="store_true")
    a("--dx", type=float, default=None); a("--dt", type=float, default=None)
    a("--scales", type=float, nargs="+", default=None); a("--scale-factors", type=float, nargs="+", default=[0.01, 0.1, 1.0])
    a("--scale-assumed-factor", type=float, default=1.0, help="inversion assumes scale x this (scale-error test)")
    a("--drop", nargs="*", default=None, help="manual channel drop; default: per window from the baseline noise; 'none' keeps all")
    a("--bad-ratio", type=float, default=3.0, help="drop a channel whose baseline sd exceeds this x the median sd")
    a("--single", default=None)
    a("--exclude-files", nargs="*", default=[], help="noise files to leave out (names from the 'files' key, e.g. chb01_23, or 0-based indices)")
    a("--no-mask", action="store_true", help="ignore 'artifact_mask' in noise_2s.npz even if present")
    a("--max-flag-frac", type=float, default=0.05, help="skip a noise window if more than this share of its samples is flagged in any used channel")
    a("--fit-t0", type=float, default=None, help="s; default: start-up end of the TEMPLATE run")
    a("--late-t0", type=float, default=None, help="s; scoring start (default: truth pulse speed settled)")
    a("--settle-tol", type=float, default=0.03); a("--x-max", type=float, default=None, help="mm; default: furthest truth front")
    a("--dxc-per-fwhm", type=float, default=4.0); a("--dxc", type=float, default=None)
    a("--modes", nargs="+", default=["smooth", "transport_fit", "transport_lit"], choices=["smooth", "transport_fit", "transport_lit", "transport_truth"])
    a("--c-lit", type=float, default=None, help="mm/min literature speed for transport_lit (default: from config if present)")
    a("--lam-min", type=float, default=1e-8); a("--lam-max", type=float, default=1e4); a("--nlam", type=int, default=9)
    a("--ridge", type=float, default=1e-3); a("--tau", type=float, default=1.05)
    a("--c-min", type=float, default=0.5); a("--c-max", type=float, default=6.0)
    a("--x0-min", type=float, default=None); a("--x0-max", type=float, default=None)
    a("--s-min", type=float, default=0.5); a("--s-max", type=float, default=2.0)
    a("--s-min-gauss", type=float, default=0.1); a("--s-max-gauss", type=float, default=3.0)
    a("--nuisance", default="offset+slope", choices=["none", "offset", "offset+slope"]); a("--ar1", action="store_true")
    a("--h-to-volt", default="auto"); a("--detect-frac", type=float, default=0.3)
    a("--target-x", type=float, default=None); a("--mismatch-d", type=float, default=1.15); a("--drift-uv", type=float, default=5.0)
    a("--cut-sens", type=float, nargs="*", default=None)
    args = ap.parse_args(argv)
    if args.quick and args.seeds is None: args.seeds = 3
    return args

# ================================================================== loading (strict: errors instead of silent fallbacks)
def _fail(inp_notes, msg, strict):
    if strict: raise SystemExit("ERROR: " + msg + "   (use --allow-fallbacks to continue and record it)")
    inp_notes.append("FALLBACK: " + msg); print("[fallback]", msg)

def _find(data, name, alts, notes):
    p = Path(data) / name
    if p.exists(): return p
    for alt in alts:
        if (Path(data) / alt).exists(): notes.append(f"{name} not found, used {alt}"); return Path(data) / alt
    c = [q for q in Path(data).glob(Path(name).stem + "*" + Path(name).suffix) if "backup" not in q.name.lower()]
    if len(c) == 1: notes.append(f"{name} not found, used {c[0].name}"); return c[0]
    raise SystemExit(f"ERROR: cannot resolve {name} in {data}; candidates: {[q.name for q in c]}")

def _npz(p):
    z = np.load(p, allow_pickle=True); return {k: z[k] for k in z.files}
scalar = lambda d, k: float(np.asarray(d[k]).ravel()[0]) if k in d and np.asarray(d[k]).size else None

def load_A(p):
    d = _npz(p)
    if "indptr" in d:
        fmt = str(np.asarray(d.get("format", "csr")).item()).lower().replace("b'", "").replace("'", "")
        return (sp.csc_matrix if "csc" in fmt else sp.csr_matrix)((d["data"], d["indices"], d["indptr"]), shape=tuple(int(s) for s in d["shape"]))
    return sp.csr_matrix(d[next(iter(d))])

def load_hrows(p, nsrc, notes):
    """Returns (function name -> bipolar H row over sources, format string)."""
    d = _npz(p)
    if "H_bip" in d:
        Hm, nm, fmt = np.asarray(d["H_bip"], float), [str(c) for c in d["channels"]], "bipolar"
    elif "H" in d:
        Hm, nm, fmt = np.asarray(d["H"], float), [str(c) for c in d["channels"]], "monopolar"
    else: raise SystemExit(f"ERROR: {p.name} has neither 'H_bip' nor 'H' (keys {list(d)})")
    if Hm.shape[1] != nsrc and Hm.shape[0] == nsrc: Hm = Hm.T; notes.append(f"{p.name}: H transposed to (channels, sources)")
    if Hm.shape[1] != nsrc: raise SystemExit(f"ERROR: H {Hm.shape} does not connect to A with {nsrc} sources")
    alias = {"T3": "T7", "T4": "T8", "T5": "P7", "T6": "P8"}; nrm = lambda e: alias.get(e.strip().upper(), e.strip().upper())
    idx = {nrm(n): i for i, n in enumerate(nm)}
    def row(name):
        if fmt == "bipolar" and name in nm: return Hm[nm.index(name)]
        parts = [q for q in name.replace(" ", "").split("-") if q][:2]
        if fmt == "bipolar": raise KeyError(name)
        return Hm[idx[nrm(parts[0])]] - Hm[idx[nrm(parts[1])]] if len(parts) == 2 else Hm[idx[nrm(parts[0])]]
    return row, fmt

def operator(data, names, u, A, notes, strict, label):
    """G0 (grid nodes x channels) = A^T H^T for the requested channel names; the folder's own A and H."""
    row, fmt = load_hrows(_find(data, "H_bip.npz", ["H_matrix.npz"], notes), A.shape[0], notes)
    rows, got, skipped = [], [], []
    for n in names:
        try: rows.append(row(n)); got.append(n)
        except KeyError: skipped.append(n)
    if skipped: notes.append(f"{label}: channels without an H row (skipped): {skipped}")
    return np.asarray(A.T @ np.array(rows).T), got, fmt

def load_inputs(args):
    notes, strict = [], not args.allow_fallbacks
    data = Path(args.data); inp = NS(args=args, notes=notes, strict=strict, data=data)
    pA = _find(data, "A_matrix.npz", [], notes); inp.A = A = load_A(pA)
    pn = _find(data, "noise_2s.npz", [], notes); dn = _npz(pn)
    nv = np.asarray(dn["noise_V"], float)
    if "channels" not in dn: raise SystemExit("ERROR: noise_2s.npz has no 'channels' list")
    # --- optional artifact mask (Member 4): same layout as noise_V (files, channels, samples); True = flagged sample
    amask = None
    if "artifact_mask" in dn and not args.no_mask:
        amask = np.asarray(dn["artifact_mask"]).astype(bool)
        if amask.shape != nv.shape and amask.ndim == 3 and amask.transpose(0, 2, 1).shape == nv.shape: amask = amask.transpose(0, 2, 1)
        if amask.shape != nv.shape: raise SystemExit(f"ERROR: artifact_mask shape {amask.shape} does not match noise_V {nv.shape}")
    elif "artifact_mask" not in dn: notes.append("noise_2s.npz has no artifact_mask: no artifact screening of noise windows")
    fnames = [str(q) for q in dn["files"]] if "files" in dn else [f"file{i}" for i in range(nv.shape[0])]
    drop_f = set()
    for e in args.exclude_files:
        if e.isdigit() and int(e) < len(fnames): drop_f.add(int(e))
        else:
            hit = [i for i, q in enumerate(fnames) if e.lower() in q.lower()]
            if len(hit) != 1: raise SystemExit(f"ERROR: --exclude-files '{e}' matches {len(hit)} files in {fnames}")
            drop_f.add(hit[0])
    keep_f = [i for i in range(nv.shape[0]) if i not in drop_f]
    nv = nv[keep_f]; fnames = [fnames[i] for i in keep_f]; amask = amask[keep_f] if amask is not None else None
    inp.fnames, inp.excluded = fnames, [str(q) for q in dn["files"][sorted(drop_f)]] if "files" in dn and drop_f else sorted(drop_f)
    if nv.shape[0] < 2: raise SystemExit("ERROR: need >= 2 noise files after exclusion; got %d" % nv.shape[0])
    inp.nfile, nch, inp.nsamp = nv.shape
    noise_names = [str(c) for c in dn["channels"]]
    pu = Path(args.u2d) if args.u2d else _find(data, "u_2d_default.npy", [], notes)
    u = np.load(pu).astype(float); u = u[:, None, :] if u.ndim == 2 else u
    inp.NT, inp.NY, inp.NX = u.shape; inp.u = u
    if A.shape[1] != inp.NY * inp.NX: raise SystemExit(f"ERROR: A has {A.shape[1]} nodes, u is {inp.NY}x{inp.NX}")
    dfp = _npz(_find(data, "flat_patch.npz", [], notes)); dsd = _npz(_find(data, "sensors_default.npz", [], notes))
    if "nx" in dfp and (int(dfp["nx"]), int(dfp["ny"])) != (inp.NX, inp.NY): raise SystemExit("ERROR: u grid differs from flat_patch.npz")
    dts = [v for v in (scalar(dsd, "dt_s"), scalar(dn, "dt_s")) if v]
    if len(dts) == 2 and abs(dts[0] - dts[1]) > 1e-9 * dts[0]: _fail(notes, f"dt differs between sensors ({dts[0]}) and noise ({dts[1]})", strict)
    inp.DT = args.dt or (dts[0] if dts else None); inp.DX = args.dx or scalar(dfp, "dx_mm")
    if inp.DT is None or inp.DX is None: raise SystemExit("ERROR: dt or dx not in the files; pass --dt/--dx")
    inp.t = np.arange(inp.NT) * inp.DT; inp.x = np.arange(inp.NX) * inp.DX
    sc_file = scalar(dsd, "scale_Am")
    inp.SCALES = args.scales or (list(np.asarray(dsd["scales_Am"], float)) if "scales_Am" in dsd else ([sc_file * f for f in args.scale_factors] if sc_file else None))
    if inp.SCALES is None: raise SystemExit("ERROR: no dipole scale (pass --scales or provide scale_Am in sensors_default.npz)")
    if args.quick: inp.SCALES = [inp.SCALES[len(inp.SCALES) // 2]]
    # channels and operators
    G0, names, inp.hfmt = operator(data, noise_names, u, A, notes, strict, "main")
    cols = np.array([noise_names.index(n) for n in names]); inp.names = names
    inp.G0 = G0; inp.G1 = G0.reshape(inp.NY, inp.NX, -1).sum(axis=0)
    YU0 = u.reshape(inp.NT, -1) @ G0; inp.hinfo = dict(h=1.0, mismatch=None, source="documented V per A*m")
    if str(args.h_to_volt).lower() == "auto" and "clean_V" in dsd and sc_file:
        cs = [str(c) for c in dsd["channels"]]; cv = np.asarray(dsd["clean_V"], float)
        cv = cv.T if (cv.shape[0] == len(cs) and cv.shape[0] != inp.NT) else cv
        sel = [(i, cs.index(n)) for i, n in enumerate(names) if n in cs]
        if len(sel) >= 2 and cv.shape[0] == inp.NT:
            i, j = map(list, zip(*sel)); p_, c_ = sc_file * YU0[:, i], cv[:, j]
            h = float((p_ * c_).sum() / (p_ * p_).sum()); err = float(np.linalg.norm(c_ - h * p_) / np.linalg.norm(c_))
            inp.hinfo = dict(h=h, mismatch=err, source=f"fitted to clean_V ({len(i)} channels)")
            if err > 0.05: _fail(notes, f"clean_V does not match A^T H^T u (mismatch {err:.2e}); h=1 used", strict); inp.hinfo["h"] = 1.0
        else: _fail(notes, "clean_V not comparable; h=1 used", strict)
    elif str(args.h_to_volt).lower() != "auto": inp.hinfo = dict(h=float(args.h_to_volt), mismatch=None, source="--h-to-volt")
    else: _fail(notes, "no clean_V/scale_Am to calibrate h; h=1 used", strict)
    inp.HV = inp.hinfo["h"]; inp.G1 = inp.G1 * inp.HV; inp.YU = YU0 * inp.HV
    inp.YU_gen = inp.YU; inp.gen_info = None                       # generation operator (perturbed if --gen-data)
    if args.gen_data:
        gA = load_A(_find(Path(args.gen_data), "A_matrix.npz", [], notes)); gG0, gnames, gfmt = operator(Path(args.gen_data), names, u, gA, notes, strict, "gen")
        if gnames != names: raise SystemExit("ERROR: --gen-data operator lacks some inversion channels")
        if gA.shape[1] != A.shape[1]: raise SystemExit("ERROR: gen A grid differs")
        inp.YU_gen = u.reshape(inp.NT, -1) @ gG0; inp.gen_info = dict(folder=args.gen_data, format=gfmt)
    inp.NZ = [nv[f][cols].T * 1e6 for f in range(inp.nfile)]       # (samples, channels) in uV
    inp.COV = {}
    for f in range(inp.nfile):
        o = [g for g in range(inp.nfile) if g != f]
        Xs = []
        for g in o:                                                # covariance only from samples with no flag in any used channel
            ok = ~amask[g][cols].any(0) if amask is not None else np.ones(inp.nsamp, bool)
            Xs.append(nv[g][cols][:, ok] - nv[g][cols][:, ok].mean(1, keepdims=True))
        X = np.concatenate(Xs, axis=1) * 1e6
        inp.COV[f] = X @ X.T / (X.shape[1] - 1)
    inp.kmax = inp.nsamp // (2 * inp.NT)                           # disjoint (data, baseline) window pairs per file
    # valid window pairs: both windows have <= max_flag_frac flagged samples in any used channel (Member 4 mask)
    inp.wins, inp.skipped_win = [], 0
    for k in range(inp.kmax):
        for f in range(inp.nfile):
            s1 = k * 2 * inp.NT
            fr = 0.0 if amask is None else float(amask[f][cols][:, s1:s1 + 2 * inp.NT].any(0).mean())
            if fr <= args.max_flag_frac: inp.wins.append((f, s1))
            else: inp.skipped_win += 1
    inp.maskinfo = dict(used=amask is not None, flagged_share_per_channel=None if amask is None else dict(zip(names, np.round(amask[:, cols].mean((0, 2)), 5).tolist())),
                        windows_skipped=inp.skipped_win, max_flag_frac=args.max_flag_frac)
    if len(inp.wins) < 2: raise SystemExit("ERROR: fewer than 2 clean noise windows left")
    inp.nseeds_max = len(inp.wins)
    inp.nseeds = args.seeds or inp.nseeds_max
    if inp.nseeds > inp.nseeds_max and not args.allow_overlap:
        raise SystemExit(f"ERROR: {inp.nseeds} seeds requested but only {inp.nseeds_max} disjoint windows exist (use --allow-overlap)")
    inp.sig_ch = np.linalg.norm(inp.G1, axis=0)
    # scoring geometry from the TRUTH (scoring only; the inversion never uses these)
    x, t = inp.x, inp.t; truth = u.mean(axis=1); inp.truth = truth
    pos = peak_pos(truth, x); inp.pos_t = pos; inp.XMAX = args.x_max or float(pos.max())
    inp.FWHM = fwhm_of(truth[inp.NT // 2], inp.DX)
    inp.LATE_T0 = args.late_t0 if args.late_t0 is not None else settle_time(pos, t, inp.FWHM, args.settle_tol)
    inp.tl = t >= inp.LATE_T0; inp.xi = x <= inp.XMAX
    inp.mk = inp.tl & (pos <= inp.XMAX - inp.FWHM)
    inp.TRUE_SPEED = float(np.polyfit(t[inp.mk], pos[inp.mk], 1)[0] * 60)
    inp.TARGET_X = args.target_x or 0.6 * x[-1]
    inp.TRUE_ARR = crossing(pos, t, inp.TARGET_X)
    inp.X0MIN = args.x0_min if args.x0_min is not None else -0.4 * x[-1]; inp.X0MAX = args.x0_max if args.x0_max is not None else 0.4 * x[-1]
    # template (Member 3 solver, independent resolution) -> shape, width, fit start, inversion grid
    inp.tmpl, inp.FW_T, FIT_T0_T, inp.tmpl_src = build_template(inp)
    inp.FIT_T0 = args.fit_t0 if args.fit_t0 is not None else FIT_T0_T
    inp.tfit = t >= inp.FIT_T0
    inp.DXC = args.dxc or max(inp.FW_T / args.dxc_per_fwhm, inp.DX)
    inp.gauss = lambda zz: np.exp(-0.5 * (zz * 2.3548 / inp.FW_T) ** 2)
    xc = np.arange(0, x[-1] + 1e-9, inp.DXC); inp.xc = xc
    inp.Pm = np.column_stack([np.interp(x, xc, e) for e in np.eye(len(xc))])
    inp.c_lit = args.c_lit or find_in_config(args.config, ("speed_mm_per_min", "speed_mm_min", "front_speed_mm_per_min"))
    inp._linv = {}
    return inp

# ================================================================== geometry helpers
def peak_pos(U, x):
    """Sub-grid position of the maximum of each row (parabolic interpolation)."""
    i = U.argmax(1); j = np.clip(i, 1, U.shape[1] - 2); r = np.arange(len(i))
    a, b, c = U[r, j - 1], U[r, j], U[r, j + 1]; den = a - 2 * b + c
    off = np.where(np.abs(den) > 1e-12, 0.5 * (a - c) / den, 0.0); off = np.clip(off, -1, 1)
    return np.where((i == j), x[j] + off * (x[1] - x[0]), x[i])

def fwhm_of(prof, dx):
    ab = np.where(prof >= 0.5 * prof.max())[0]; return max((ab[-1] - ab[0] + 1) * dx, 2 * dx)

def crossing(pos, t, level):
    k = np.where(pos >= level)[0]
    if not len(k) or k[0] == 0: return float("nan")
    i = k[0]; return float(t[i - 1] + (level - pos[i - 1]) / (pos[i] - pos[i - 1]) * (t[i] - t[i - 1]))

def settle_time(pos, t, fwhm, tol):
    """First time after which the local pulse speed stays within tol of its final value."""
    n = len(t); h = n // 2; cf = np.polyfit(t[h:], pos[h:], 1)[0]
    nw = max(int(round(2 * fwhm / max(abs(cf), 1e-9) / (t[1] - t[0]))), 4)
    if n - nw < 3: return float(t[0])
    sp_ = np.array([np.polyfit(t[k:k + nw], pos[k:k + nw], 1)[0] for k in range(n - nw)])
    bad = np.where(np.abs(sp_ - cf) > tol * abs(cf))[0]
    return float(t[min(bad[-1] + 1, n - 1)]) if len(bad) else float(t[0])

def find_in_config(path, keys):
    p = Path(path)
    if not p.exists(): return None
    try:
        import yaml; cfg = yaml.safe_load(p.read_text())
    except Exception: return None
    def walk(o):
        if isinstance(o, dict):
            for k, v in o.items():
                if k in keys and isinstance(v, (int, float)): return float(v)
                r = walk(v)
                if r is not None: return r
    return walk(cfg)

def member3_run(cfg_path, dx, method, dmult=None):
    sys.path.insert(0, str(ROOT))
    from src.model import simulate, load_config, with_resolution
    c = with_resolution(load_config(cfg_path), dx=dx, method=method)
    if dmult: c["model"]["dimension"] = 1; c["model"]["excitable"]["D"] = float(c["model"]["excitable"]["D"]) * dmult
    r = simulate(c); U = np.asarray(r["u"]); return (U[:, 0, :] if U.ndim == 3 else U), np.asarray(r["x"])

def build_template(inp):
    a = inp.args
    try:
        U, xs = member3_run(a.template_config or a.config, inp.DX * a.template_dx_factor, a.template_method); src = "member3_solver"
        src += f" dx={inp.DX * a.template_dx_factor:g} cfg={Path(a.template_config or a.config).name}"
        tt = np.arange(U.shape[0]) * inp.DT
    except Exception as e:
        _fail(inp.notes, f"Member 3 solver not usable ({repr(e)[:90]}); template taken from the TRUTH (optimistic)", inp.strict)
        U, xs, tt, src = inp.truth, inp.x, inp.t, "truth_snapshot_OPTIMISTIC"
    pos = peak_pos(U, xs); k = int(np.argmin(np.abs(pos - a.template_at * xs[-1])))
    prof = U[k]; xp = xs - pos[k]; fw = fwhm_of(prof, float(np.median(np.diff(xs))))
    fit_t0 = settle_time(pos, tt, fw, a.settle_tol)
    return (lambda zz: np.interp(zz, xp, prof, left=0.0, right=0.0)), fw, fit_t0, src

# ================================================================== whitening, windows, observations
def channels(inp, case, f, nzb):
    a = inp.args; n = len(inp.names)
    if case == "all": return np.arange(n)
    if case == "single":
        nsd = np.concatenate([inp.NZ[g] for g in range(inp.nfile) if g != f]).std(0); return np.array([int(np.argmax(inp.sig_ch / nsd))])
    if a.drop is not None:
        if [d.lower() for d in a.drop] == ["none"]: return np.arange(n)
        return np.array([i for i, nm in enumerate(inp.names) if nm.upper() not in {d.upper() for d in a.drop}])
    sd = nzb.std(0); ok = (sd <= a.bad_ratio * np.median(sd)) & (inp.sig_ch > 0)       # decided on the BASELINE window only
    return np.where(ok)[0]

def linv(inp, f, idx):
    key = (f, tuple(idx))
    if key not in inp._linv:
        C = inp.COV[f][np.ix_(idx, idx)]; L = cholesky(C + 1e-3 * np.mean(np.diag(C)) * np.eye(len(idx)), lower=True)
        inp._linv[key] = solve_triangular(L, np.eye(len(idx)), lower=True)
    return inp._linv[key]

def window(inp, seed):
    if seed < len(inp.wins): f, s1 = inp.wins[seed]                 # clean disjoint (data, baseline) pairs
    else: f = seed % inp.nfile; s1 = int(np.random.default_rng(seed).integers(0, inp.nsamp - 2 * inp.NT + 1))   # only with --allow-overlap
    return f, inp.NZ[f][s1:s1 + inp.NT], inp.NZ[f][s1 + inp.NT:s1 + 2 * inp.NT]

def simulate_observations(inp, scale, nzd, idx):
    """Sensor data in uV: scale * (generation operator) u + measured noise window (generation operator = --gen-data if given)."""
    return scale * 1e6 * inp.YU_gen[:, idx] + nzd[:, idx]

def lag1(inp, f, idx):
    w = np.concatenate([inp.NZ[g] for g in range(inp.nfile) if g != f])[:, idx] @ linv(inp, f, idx).T
    return float(np.mean(w[:-1] * w[1:]) / np.mean(w ** 2))

# ================================================================== Method B: travelling pulse (nuisance projected out)
def make_B(inp, Y, idx, f, scale, shape, tm):
    a = inp.args; tt = inp.t[tm]; Li = linv(inp, f, idx); m = len(idx); nt = tt.size
    Gw = inp.G1[:, idx] * scale * a.scale_assumed_factor * 1e6
    rho = lag1(inp, f, idx) if a.ar1 else 0.0
    T = np.eye(nt) - rho * np.eye(nt, k=-1); T[0, 0] = np.sqrt(1 - rho ** 2)
    tf = lambda M: T @ M
    yw = tf(Y[tm] @ Li.T).ravel(); ttc = tt - tt.mean(); cols = []
    nuis = a.nuisance if m >= 4 else "none"          # offset/slope are not identifiable from < 4 channels (they absorb the pulse)
    if nuis != "none":
        for i in range(m):
            cols.append(tf(np.outer(np.ones(nt), Li[:, i])).ravel())
            if nuis == "offset+slope": cols.append(tf(np.outer(ttc, Li[:, i])).ravel())
    if cols:
        Q, _ = np.linalg.qr(np.column_stack(cols)); proj = lambda v: v - Q @ (Q.T @ v)
    else: proj = lambda v: v
    yp = proj(yw)
    def feat(th, shp): return proj(tf(((shp((inp.x[None, :] - th[0] - th[1] * tt[:, None]) / th[2])) @ Gw) @ Li.T).ravel())
    def resid(th, shp):
        f_ = feat(th, shp); return yp - (f_ @ yp / max(f_ @ f_, 1e-30)) * f_
    return resid, feat, yp, len(cols)

def fit_B(inp, Y, idx, f, scale, shape, lb_s, ub_s, tm=None, fixed_s=None):
    a = inp.args; tm = inp.tfit if tm is None else tm
    resid, feat, yp, ncols = make_B(inp, Y, idx, f, scale, shape, tm)
    LB = np.array([inp.X0MIN, a.c_min / 60, lb_s]); UB = np.array([inp.X0MAX, a.c_max / 60, ub_s])
    c0s = np.geomspace(a.c_min / 60 * 1.5, a.c_max / 60 / 1.5, 3); best, nfail = None, 0
    for x0 in (inp.X0MIN / 3, inp.X0MAX / 3):
        for c in c0s:
            for s in (lb_s + 0.15 * (ub_s - lb_s), lb_s + 0.6 * (ub_s - lb_s)):
                try: r = least_squares(lambda th: resid(th, shape), [x0, c, s], bounds=(LB, UB), x_scale=[1, .02, .5])
                except Exception: nfail += 1; continue
                if best is None or r.cost < best.cost: best = r
    if best is None: raise RuntimeError("all optimiser starts failed in fit_B")
    th = best.x; f_ = feat(th, shape); amp = float(f_ @ yp / (f_ @ f_)); tol = 1e-3 * (UB - LB)
    bound = bool(np.any(np.isclose(th, LB, rtol=0, atol=tol) | np.isclose(th, UB, rtol=0, atol=tol)))
    n_eff = max(yp.size - 4 - ncols, 1); s2 = float(np.sum(best.fun ** 2) / n_eff)
    try: sd = np.sqrt(np.diag(s2 * np.linalg.pinv(best.jac.T @ best.jac)))
    except Exception: sd = np.full(3, np.nan)
    rec = amp * shape((inp.x[None, :] - th[0] - th[1] * inp.t[:, None]) / th[2])
    if amp < 0: rec = np.zeros_like(rec)
    return dict(rec=rec, th=th, amp=amp, rms=float(np.sqrt(np.mean(best.fun ** 2))), at_bound=bound, neg_amp=bool(amp < 0), nfail=nfail,
                sd_x0=float(sd[0]), sd_c_mm_min=float(sd[1] * 60), sd_s=float(sd[2]), resid=resid, LB=LB, UB=UB)

# ================================================================== Method A: Tikhonov (normal equations, sparse LU)
def prep_A(inp):
    nc = len(inp.xc); nt2 = int(inp.tfit.sum())
    d1 = lambda n: sp.diags([-np.ones(n - 1), np.ones(n - 1)], [0, 1], shape=(n - 1, n)).tocsr()
    I_t, I_c = sp.identity(nt2).tocsr(), sp.identity(nc).tocsr()
    R = sp.vstack([sp.kron(I_t, d1(nc)), sp.kron(d1(nt2), I_c)]).tocsr()
    inp.A_ = NS(nc=nc, nt2=nt2, d1=d1, I_t=I_t, I_c=I_c, RtR=(R.T @ R).tocsr(), Idn=sp.identity(nt2 * nc).tocsr(),
                Dc1=sp.csr_matrix(np.gradient(np.eye(nc), axis=0)), LAMS=np.geomspace(inp.args.lam_min, inp.args.lam_max, inp.args.nlam))

def fit_A(inp, Y, idx, f, scale, mode, c_mm_s, nzb):
    if not hasattr(inp, "A_"): prep_A(inp)
    q = inp.A_; a = inp.args; Li = linv(inp, f, idx)
    Gw = inp.G1[:, idx] * scale * a.scale_assumed_factor * 1e6
    M = (inp.Pm.T @ Gw) @ Li.T; Yw = Y[inp.tfit] @ Li.T; s0 = np.linalg.norm(M, 2) ** 2
    AtA = sp.kron(q.I_t, sp.csr_matrix(M @ M.T)).tocsr(); rhs = (Yw @ M.T).ravel()
    if mode == "smooth": RtR = q.RtR
    else:
        cr = c_mm_s * inp.DT / inp.DXC; avg = sp.diags([np.ones(q.nt2 - 1), np.ones(q.nt2 - 1)], [0, 1], shape=(q.nt2 - 1, q.nt2)).tocsr() * 0.5
        Rt = (sp.kron(q.d1(q.nt2), q.I_c) + cr * sp.kron(avg, q.Dc1)).tocsr(); RtR = (Rt.T @ Rt + 0.1 * q.RtR).tocsr()
    RtR = RtR + a.ridge * q.Idn; out = []
    for lam in q.LAMS:
        xs = splu((AtA + lam * s0 * RtR).tocsc(), permc_spec="MMD_AT_PLUS_A").solve(rhs)
        X = xs.reshape(q.nt2, q.nc); rec = np.zeros((inp.NT, inp.NX)); rec[inp.tfit] = X @ inp.Pm.T
        out.append((lam, float(np.sqrt(np.mean((X @ M - Yw) ** 2))), evaluate(inp, rec)["relL2"], rec))
    sig = float(np.sqrt(np.mean((nzb[:, idx] @ Li.T) ** 2)))                       # sigma from the adjacent baseline window
    ok = [o for o in out if o[1] <= a.tau * sig]
    pick = max(ok, key=lambda o: o[0]) if ok else out[-1]                          # none fits -> most regularised (no wave found)
    orc = min(out, key=lambda o: o[2])
    return dict(rec=pick[3], lam=pick[0], rms=pick[1], dp_fallback=not ok, relL2_oracle=orc[2], lam_oracle=orc[0], sweep=[(o[0], o[1], o[2]) for o in out])

# ================================================================== public interface for the pipeline
def evaluate(inp, rec):
    a = inp.args; d = lambda m, xs_: np.linalg.norm((rec - inp.truth)[m][:, xs_]) / np.linalg.norm(inp.truth[m][:, xs_])
    det = rec[inp.tfit].max() >= a.detect_frac * inp.truth.max()
    pr = inp.x[rec.argmax(1)]; e = np.abs(np.where(rec.max(1) >= a.detect_frac * inp.truth.max(), pr, np.nan) - inp.pos_t)[inp.mk]
    return dict(relL2=float(d(inp.tl, inp.xi)), relL2_fit=float(d(inp.tfit, inp.xi)), no_detect=bool(not det),
                front_err_mm=float(np.nanmean(e)) if np.isfinite(e).any() else np.nan)

def reconstruct(inp, Y, idx, f, scale, method="B_model", c_mm_min=None, nzb=None):
    a = inp.args
    if method == "B_model": return fit_B(inp, Y, idx, f, scale, inp.tmpl, a.s_min, a.s_max)
    if method == "B_gauss": return fit_B(inp, Y, idx, f, scale, inp.gauss, a.s_min_gauss, a.s_max_gauss)
    mode = method[2:]; c = (c_mm_min or inp.TRUE_SPEED) / 60
    return fit_A(inp, Y, idx, f, scale, "smooth" if mode == "smooth" else "transport", c, nzb)

# ================================================================== experiment driver
def run(inp):
    a = inp.args; OUT = Path(a.out) / time.strftime("run_%Y%m%d_%H%M%S"); OUT.mkdir(parents=True, exist_ok=True)
    t, x, NT = inp.t, inp.x, inp.NT; sc_mid = inp.SCALES[len(inp.SCALES) // 2]
    print(f"A {inp.A.shape} | H {inp.hfmt} | channels {len(inp.names)} | noise files {inp.nfile} | disjoint windows {inp.nseeds_max} (using {inp.nseeds})")
    print(f"dx={inp.DX} mm dt={inp.DT} s | scales {inp.SCALES} | h={inp.HV:.6g} (mismatch {inp.hinfo['mismatch']})")
    print(f"TRUTH (scoring only): speed {inp.TRUE_SPEED:.3f} mm/min | FWHM {inp.FWHM:.2f} mm | late_t0 {inp.LATE_T0:.0f} s | x_max {inp.XMAX:.1f} mm | arrival at {inp.TARGET_X:.1f} mm = {inp.TRUE_ARR:.0f} s")
    print(f"INVERSION (no truth): template {inp.tmpl_src} | template FWHM {inp.FW_T:.2f} mm | fit_t0 {inp.FIT_T0:.0f} s | inversion dx {inp.DXC:.2f} mm")
    rows, store, sweeps = [], {}, {}
    nm_of = lambda mode: "A_" + mode
    for cn in ("main", "all", "single"):
        for scale in inp.SCALES:
            for seed in range(inp.nseeds):
                f, nzd, nzb = window(inp, seed); idx = channels(inp, cn, f, nzb); Y = simulate_observations(inp, scale, nzd, idx)
                base = dict(case=cn, fold=f, scale=scale, seed=seed, n_ch=len(idx)); thB = None
                for mn in ("B_model", "B_gauss"):
                    o = reconstruct(inp, Y, idx, f, scale, mn); th = o["th"]
                    if mn == "B_model": thB = th
                    rows.append({**base, "method": mn, **evaluate(inp, o["rec"]), "speed_mm_min": th[1] * 60, "speed_err": th[1] * 60 - inp.TRUE_SPEED,
                                 "arrival_s": (inp.TARGET_X - th[0]) / th[1], "x0": th[0], "s": th[2], "amp": o["amp"], "amp_x_s": o["amp"] * th[2], "rms": o["rms"],
                                 "at_bound": o["at_bound"], "neg_amp": o["neg_amp"], "sd_speed_jac": o["sd_c_mm_min"], "sd_x0_jac": o["sd_x0"], "nfail": o["nfail"]})
                    if seed == 0 and cn == "main": store[(mn, scale)] = o["rec"]
                for mode in a.modes:
                    if mode == "transport_lit" and not inp.c_lit: continue
                    c = {"smooth": None, "transport_fit": thB[1] * 60, "transport_lit": inp.c_lit, "transport_truth": inp.TRUE_SPEED}[mode]
                    o = reconstruct(inp, Y, idx, f, scale, nm_of(mode), c, nzb)
                    rows.append({**base, "method": nm_of(mode), **evaluate(inp, o["rec"]), "lam": o["lam"], "rms": o["rms"], "dp_fallback": o["dp_fallback"],
                                 "relL2_oracle": o["relL2_oracle"]})
                    if seed == 0 and cn == "main": store[(nm_of(mode), scale)] = o["rec"]; sweeps[(nm_of(mode), scale)] = o["sweep"]
            print(f"done {cn} scale {scale:.0e}", flush=True)
    df = pd.DataFrame(rows); df.to_csv(OUT / "results.csv", index=False)
    # cut sensitivity of the pulse fit
    cuts = a.cut_sens or [round(inp.FIT_T0 * m_) for m_ in (0.25, 0.5, 1.0, 1.5, 2.0)]; cs = []
    for cut in cuts:
        for seed in range(inp.nseeds):
            f, nzd, nzb = window(inp, seed); idx = channels(inp, "main", f, nzb); Y = simulate_observations(inp, sc_mid, nzd, idx)
            o = fit_B(inp, Y, idx, f, sc_mid, inp.tmpl, a.s_min, a.s_max, tm=t >= cut)
            cs.append(dict(fit_start_s=cut, seed=seed, speed_err=o["th"][1] * 60 - inp.TRUE_SPEED, **evaluate(inp, o["rec"])))
    cs = pd.DataFrame(cs); cs.to_csv(OUT / "cut_sensitivity.csv", index=False)
    print("\ncut sensitivity (B_model, scale %.0e, median)" % sc_mid); print(cs.groupby("fit_start_s")[["relL2", "relL2_fit", "speed_err"]].median().round(3).to_string())
    # extra tests: model mismatch (D x factor), drift (nuisance on/off)
    extra = []
    def score(Y, ut, scale, idx, f, test, seed, sp_true, nuis):
        a.nuisance, keep = nuis, a.nuisance; o = fit_B(inp, Y, idx, f, scale, inp.tmpl, a.s_min, a.s_max); a.nuisance = keep
        d = lambda m: float(np.linalg.norm((o["rec"] - ut)[m][:, inp.xi]) / np.linalg.norm(ut[m][:, inp.xi]))
        extra.append(dict(test=test, nuisance=nuis, scale=scale, seed=seed, speed_err=o["th"][1] * 60 - sp_true, at_bound=o["at_bound"], relL2=d(inp.tl), relL2_fit=d(inp.tfit)))
    try:
        um, _ = member3_run(a.config, inp.DX, "solve_ivp", dmult=a.mismatch_d); assert um.shape == (NT, inp.NX), um.shape
        pm = peak_pos(um, x); mm = inp.tl & (pm <= inp.XMAX - inp.FWHM); spm = float(np.polyfit(t[mm], pm[mm], 1)[0] * 60)
        print("mismatch truth D x %.2f: speed %.3f mm/min" % (a.mismatch_d, spm))
        for scale in inp.SCALES:
            for seed in range(inp.nseeds):
                f, nzd, nzb = window(inp, seed); idx = channels(inp, "main", f, nzb)
                score(scale * 1e6 * (um @ inp.G1)[:, idx] + nzd[:, idx], um, scale, idx, f, "mismatch_D", seed, spm, a.nuisance)
    except Exception as e: _fail(inp.notes, f"mismatch test skipped: {repr(e)[:100]}", inp.strict)
    for scale in inp.SCALES:
        for seed in range(inp.nseeds):
            f, nzd, nzb = window(inp, seed); idx = channels(inp, "main", f, nzb); rg = np.random.default_rng(1000 + seed)
            drift = (t / t[-1])[:, None] * rg.normal(0, a.drift_uv, len(idx))[None, :]
            Y = scale * 1e6 * (inp.truth @ inp.G1)[:, idx] + nzd[:, idx] + drift
            for nuis in ("none", a.nuisance if a.nuisance != "none" else "offset+slope"): score(Y, inp.truth, scale, idx, f, "drift_%guV" % a.drift_uv, seed, inp.TRUE_SPEED, nuis)
    ex = pd.DataFrame(extra); ex.to_csv(OUT / "extra_tests.csv", index=False)
    if len(ex): print("\nextra tests (median)"); print(ex.groupby(["test", "nuisance", "scale"])[["relL2", "relL2_fit", "speed_err", "at_bound"]].median().round(3).to_string())
    # noise-only: all disjoint windows (data and baseline of every pair), every scale
    fa = {}
    for scale in inp.SCALES:
        n_f, n_t = 0, 0
        for seed in range(inp.nseeds_max):
            f, nzd, nzb = window(inp, seed)
            for w1, w2 in ((nzd, nzb), (nzb, nzd)):
                idx = channels(inp, "main", f, w2); o = fit_B(inp, w1[:, idx], idx, f, scale, inp.tmpl, a.s_min, a.s_max); n_t += 1
                n_f += int((not o["at_bound"]) and (not o["neg_amp"]) and o["rec"][inp.tfit].max() >= a.detect_frac * inp.truth.max())
        fa["%.0e" % scale] = f"{n_f}/{n_t}"
    print("noise-only false alarms:", fa)
    # profile along the width s (seed 0, middle scale)
    f, nzd, nzb = window(inp, 0); idx = channels(inp, "main", f, nzb); Y = simulate_observations(inp, sc_mid, nzd, idx)
    o0 = fit_B(inp, Y, idx, f, sc_mid, inp.tmpl, a.s_min, a.s_max); resid = o0["resid"]; prof = []
    for s in np.linspace(a.s_min, a.s_max, 11):
        r = least_squares(lambda p: resid([p[0], p[1], s], inp.tmpl), o0["th"][:2], bounds=(o0["LB"][:2], o0["UB"][:2]), x_scale=[1, .02])
        prof.append(dict(s=s, cost=float(2 * r.cost), speed_mm_min=float(r.x[1] * 60), x0=float(r.x[0])))
    pd.DataFrame(prof).to_csv(OUT / "profile_s.csv", index=False)
    # summary
    q = lambda p: (lambda v: v.quantile(p))
    agg = dict(n=("relL2", "size"), relL2_med=("relL2", "median"), relL2_q25=("relL2", q(.25)), relL2_q75=("relL2", q(.75)), relL2_fit_med=("relL2_fit", "median"),
               front_mm_med=("front_err_mm", "median"), no_detect_frac=("no_detect", "mean"))
    for c_, v_ in (("speed_err", "speed_err"),):
        agg.update(speed_bias=(v_, "median"), speed_err_q25=(v_, q(.25)), speed_err_q75=(v_, q(.75)), speed_sd=("speed_mm_min", "std"))
    for c_ in ("at_bound", "neg_amp", "dp_fallback", "relL2_oracle", "sd_speed_jac", "amp_x_s"):
        if c_ in df: agg[c_ + ("_med" if c_ in ("relL2_oracle", "sd_speed_jac", "amp_x_s") else "_frac")] = (c_, "median" if c_ in ("relL2_oracle", "sd_speed_jac", "amp_x_s") else "mean")
    sm = df.groupby(["case", "scale", "method"]).agg(**agg).reset_index(); sm.to_csv(OUT / "summary.csv", index=False); print(sm.round(3).to_string())
    # outputs for the pipeline / Member 6
    est = dict(true_speed_mm_min=inp.TRUE_SPEED, target_x_mm=inp.TARGET_X, truth_arrival_s=inp.TRUE_ARR, arrival_definition="pulse peak reaches target_x; B_model: (target_x - x0)/c",
               noise_only_false_alarms=fa, per_scale={})
    for sc in inp.SCALES:
        d = df[(df.case == "main") & (df.scale == sc) & (df.method == "B_model")]
        est["per_scale"][f"{sc:.0e}"] = dict(n=int(len(d)), speed_mm_min_mean=float(d.speed_mm_min.mean()), speed_mm_min_sd=float(d.speed_mm_min.std()), at_bound_frac=float(d.at_bound.mean()),
            arrival_s_median=float(d.arrival_s.median()), arrival_err_s_median=float(d.arrival_s.median() - inp.TRUE_ARR),
            arrival_s_per_seed=[float(v) for v in d.arrival_s], arrival_err_s_per_seed=[float(v - inp.TRUE_ARR) for v in d.arrival_s])
    json.dump(est, open(OUT / "event_estimate.json", "w"), indent=2)
    np.savez_compressed(OUT / "reconstruction.npz", t=t, x=x, u_truth=inp.truth, **{f"rec_{k[0]}_{k[1]:.0e}": v for k, v in store.items()})
    shutil.copy(__file__, OUT / "script_used.py")
    for p in {Path(a.config), Path(a.template_config) if a.template_config else Path(a.config)}:
        if p.exists(): shutil.copy(p, OUT / ("config_" + p.name))
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16] if Path(p).is_file() else "n/a"
    git = lambda *c: (lambda r: r.stdout.strip() if r.returncode == 0 else "n/a")(subprocess.run(["git", "-C", str(ROOT), *c], capture_output=True, text=True))
    meta = dict(args=vars(a), notes=inp.notes, versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, pandas=pd.__version__),
                git=dict(commit=git("rev-parse", "HEAD"), status=git("status", "--porcelain")), h=inp.hinfo, h_format=inp.hfmt, gen_data=inp.gen_info,
                derived=dict(dx=inp.DX, dt=inp.DT, scales=inp.SCALES, channels=inp.names, truth_speed=inp.TRUE_SPEED, late_t0=inp.LATE_T0, x_max=inp.XMAX, fwhm_truth=inp.FWHM,
                             fit_t0=inp.FIT_T0, template=inp.tmpl_src, template_fwhm=inp.FW_T, dxc=inp.DXC, target_x=inp.TARGET_X, truth_arrival=inp.TRUE_ARR,
                             nseeds=inp.nseeds, nseeds_max=inp.nseeds_max, c_lit=inp.c_lit,
                             noise_files=inp.fnames, noise_files_excluded=inp.excluded, artifact_mask=inp.maskinfo,
                             lag1_whitened=[lag1(inp, f_, np.arange(len(inp.names))) for f_ in range(inp.nfile)]),
                noise_protocol="covariance from the other files; data and baseline windows are disjoint windows of the held-out file; bad channels from the baseline window",
                hashes={n: sha(_find(inp.data, n, ["H_matrix.npz"], [])) for n in ("A_matrix.npz", "H_bip.npz", "noise_2s.npz", "sensors_default.npz", "flat_patch.npz")},
                hash_u2d=sha(Path(a.u2d) if a.u2d else inp.data / "u_2d_default.npy"))
    json.dump(meta, open(OUT / "metadata.json", "w"), indent=2, default=str)
    # figures
    meths = [m_ for m_ in ("A_smooth", "A_transport_fit", "A_transport_lit", "A_transport_truth", "B_model", "B_gauss") if (df.method == m_).any()]
    fig, axs = plt.subplots(1, len(inp.SCALES), figsize=(4.8 * len(inp.SCALES), 4.2), sharey=True, squeeze=False)
    for ax, sc in zip(axs[0], inp.SCALES):
        d = df[(df.case == "main") & (df.scale == sc)]
        for i, m_ in enumerate(meths):
            v = d[d.method == m_].relL2.dropna().values; ax.scatter(i + np.random.default_rng(0).uniform(-.12, .12, len(v)), v, s=18)
            if len(v): ax.hlines(np.median(v), i - .25, i + .25, colors="k", lw=2)
        ax.axhline(1, c="grey", ls="--"); ax.set_yscale("log"); ax.set_xticks(range(len(meths))); ax.set_xticklabels(meths, rotation=30, fontsize=7)
        ax.set_title(f"{sc:.0e} A·m (n={inp.nseeds})"); ax.set_ylabel("relative L2 (t >= %.0f s), dashed = zero" % inp.LATE_T0)
    fig.savefig(OUT / "fig_error_by_method.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    if sweeps:
        fig, axs = plt.subplots(1, len(sweeps), figsize=(4 * len(sweeps), 3.8), squeeze=False)
        for ax, ((mn, sc), sw) in zip(axs[0], sweeps.items()):
            sw = np.array(sw); ax.loglog(sw[:, 0], sw[:, 2], "o-", label="rel L2"); ax.loglog(sw[:, 0], sw[:, 1], "s-", label="whitened misfit"); ax.set_title(f"{mn}, {sc:.0e} (seed 0)", fontsize=8); ax.set_xlabel("lambda (relative)")
        axs[0][0].legend(fontsize=7); fig.savefig(OUT / "fig_lambda_sweep.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    pan = [("truth", inp.truth)] + [(k, store[(k, sc_mid)]) for k in meths if (k, sc_mid) in store]
    fig, axs = plt.subplots(1, len(pan), figsize=(3.6 * len(pan), 3.8), squeeze=False)
    for ax, (ttl, r) in zip(axs[0], pan):
        im = ax.imshow(r[:, inp.xi].T, origin="lower", aspect="auto", extent=[t[0], t[-1], 0, x[inp.xi][-1]], vmin=-0.3, vmax=1); ax.set_xlabel("time (s)"); ax.set_title(ttl, fontsize=8)
        ax.axvline(inp.FIT_T0, c="w", lw=.6, ls="--")
    axs[0][0].set_ylabel("x (mm)"); fig.colorbar(im, ax=axs[0].tolist(), label="u"); fig.savefig(OUT / "fig_spacetime.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    cm = cs.groupby("fit_start_s")[["relL2", "speed_err"]].median(); fig, axs = plt.subplots(1, 2, figsize=(9, 3.6))
    axs[0].plot(cm.index, cm.relL2, "o-"); axs[0].set_ylabel("relative L2"); axs[1].plot(cm.index, cm.speed_err, "o-"); axs[1].axhline(0, c="grey", ls="--"); axs[1].set_ylabel("speed error (mm/min)")
    for ax in axs: ax.set_xlabel("fit start time (s)")
    fig.savefig(OUT / "fig_cut_sensitivity.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    print("ALL DONE ->", OUT); return OUT

def main(argv=None):
    run(load_inputs(parse_args(argv)))

if __name__ == "__main__":
    main()