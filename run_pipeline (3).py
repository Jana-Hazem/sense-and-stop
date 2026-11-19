"""Sense & Stop pipeline (repo root). Every path and setting comes from configs/config.yaml; none is typed here.

    python run_pipeline.py                         # all stages
    python run_pipeline.py --stages forward        # one stage (reuses the files earlier stages saved)
    python run_pipeline.py --quick                 # small runs in inverse and control (smoke test)

Stages, in order, and who owns them:
    model    Member 3  PDE -> u[t, y, x]          saved as <data_dir>/<forward.files.u_default>
    forward  Member 4  check_grid, bipolar H aligned with the noise channels, sensors per scale, saved in <data_dir>
    inverse  Member 5  python -m src.inverse.run_member5 (NOT modified by the pipeline) -> <outputs_dir>/<inverse.subdir>/run_<time>/
                       then "latest" copies are exported to <outputs_dir>/<inverse.subdir>/latest/ for the control stage
    control  Member 6  src.control.run_control (HYPOTHETICAL study) -> output.dir of configs/control.yaml

Each stage checks that the files it needs exist and stops with the exact path if one is missing.
The script moves to the repo root first, so it can be started from any folder.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
os.chdir(ROOT)                                    # every relative path in the configs is relative to the repo root

import numpy as np

from src.model import simulate
from src.forward.sensors import load_cfg, data_dir, model_cfg, check_grid, make_sensors
from src.inverse.run_member5 import load_A, load_hrows, _find   # Member 5's own loaders: forward and inverse build H the same way

ALL_STAGES = ("model", "forward", "inverse", "control")
# file names run_member5.py looks for in the data folder (it does not read them from the config)
M5_NAMES = {"A": "A_matrix.npz", "noise": "noise_2s.npz", "patch": "flat_patch.npz", "H_bip": "H_bip.npz", "u_default": "u_2d_default.npy"}
M5_EXTRA = ("sensors_default.npz",)               # also required by run_member5.py (dt_s, scale_Am, clean_V)


def need(path, who):
    p = Path(path)
    if not p.is_file():
        raise SystemExit(f"[{who}] required file not found: {p.resolve()}")
    return p


def latest_dir(cfg):
    return Path(cfg["paths"]["outputs_dir"]) / cfg["inverse"].get("subdir", "inverse") / "latest"


def check_names(cfg):
    files = cfg["forward"]["files"]
    for k, want in M5_NAMES.items():
        if files.get(k) != want:
            print(f"[warning] config forward.files.{k} = {files.get(k)!r}, but run_member5.py looks for {want!r}")


def bipolar_for_noise(cfg, D):
    """H rows for the noise channels. Returns A, H (channels x sources), kept indices, noise npz, names, skipped, format."""
    files = cfg["forward"]["files"]
    A = load_A(need(D / files["A"], "forward"))
    nz = np.load(need(D / files["noise"], "forward"), allow_pickle=True)
    notes = []
    row, fmt = load_hrows(_find(D, M5_NAMES["H_bip"], ["H_matrix.npz"], notes), A.shape[0], notes)
    names = [str(c) for c in nz["channels"]]
    rows, kept, skipped = [], [], []
    for j, nm in enumerate(names):                # one H row per NOISE channel; channels without a row are dropped from H and noise together
        try:
            rows.append(row(nm)); kept.append(j)
        except KeyError:
            skipped.append(nm)
    if not rows:
        raise SystemExit("[forward] no noise channel has an H row. Noise channels: %s" % names)
    if fmt == "bipolar" and skipped:              # run_member5.py matches bipolar names exactly (case-sensitive)
        print("[warning] H file is bipolar and these noise channels have no row: %s. Names must match noise_2s.npz exactly (same case)." % skipped)
    for n_ in notes:
        print("note:", n_)
    return A, np.array(rows), kept, nz, names, skipped, fmt


def stage_model(cfg, D):
    files = cfg["forward"]["files"]
    t0 = time.time()
    res = simulate(model_cfg(cfg))
    u = np.asarray(res["u"])
    np.save(D / files["u_default"], u.astype(np.float32))
    print("model: u %s, min %.3f, max %.3f, %.0f s -> %s" % (u.shape, u.min(), u.max(), time.time() - t0, D / files["u_default"]))
    return res


def stage_forward(cfg, D, res):
    files, fwd = cfg["forward"]["files"], cfg["forward"]
    patch = np.load(need(D / files["patch"], "forward"), allow_pickle=True)
    A, H_bip, kept, noise_npz, names, skipped, fmt = bipolar_for_noise(cfg, D)
    noise_dt, noise_V = float(noise_npz["dt_s"]), noise_npz["noise_V"][:, kept, :]
    if res is None:                               # model stage skipped: use the saved u (only shape and time step are checkable)
        u = np.load(need(D / files["u_default"], "forward")).astype(np.float64)
        nx, ny, dx = int(patch["nx"]), int(patch["ny"]), float(patch["dx_mm"])
        res = {"u": u, "x": np.arange(nx) * dx, "y": np.arange(ny) * dx, "t": np.arange(u.shape[0]) * noise_dt}
    check_grid(res, patch, noise_dt)
    print("forward: grid check PASS")
    print("forward: H format %s | %d of %d noise channels used | no H row: %s" % (fmt, len(kept), len(names), skipped))
    u = np.asarray(res["u"])
    rng = np.random.default_rng(fwd["seed"])
    for scale in fwd["scales_am"]:
        y, clean = make_sensors(u, H_bip, A, scale, noise_V, rng, fwd["noise_gain"])
        out = D / ("sensors_%.0e.npz" % scale)
        np.savez(out, sensors=y, clean=clean, scale_am=scale, channels=np.array([names[j] for j in kept]))
        print("forward: scale %.0e A*m -> sensors %s -> %s" % (scale, y.shape, out))


def export_latest(cfg, run_dir):
    """Copy the files Member 6 reads to <...>/latest/. event_estimate.json keeps only VALID arrival estimates
    (fit not on a bound, amplitude > 0), taken from results.csv; the untouched original is saved as event_estimate_raw.json."""
    import pandas as pd
    L = latest_dir(cfg); L.mkdir(parents=True, exist_ok=True)
    names = cfg["inverse"]["outputs"]
    for n in ("results.csv", "summary.csv", "metadata.json", names["reconstruction"]):
        if (run_dir / n).is_file():
            shutil.copy(run_dir / n, L / n)
    raw = run_dir / names["event_estimate"]
    shutil.copy(raw, L / "event_estimate_raw.json")
    ev = json.load(open(raw))
    df = pd.read_csv(run_dir / "results.csv")
    d = df[(df.case == "main") & (df.method == "B_model")]
    for sc, e in ev["per_scale"].items():
        m = d[np.isclose(d.scale.astype(float), float(sc), rtol=1e-6)]
        arr = e.get("arrival_s_per_seed", [])
        if len(m) != len(arr):
            print(f"[warning] scale {sc}: results.csv has {len(m)} rows but event file has {len(arr)} arrivals; left unfiltered")
            continue
        ok = ~(m.at_bound.astype(bool).values | m.neg_amp.astype(bool).values)
        truth_arr = float(ev["truth_arrival_s"])
        e["arrival_s_per_seed"] = [float(v) for v, k in zip(arr, ok) if k]
        e["arrival_err_s_per_seed"] = [float(v) - truth_arr for v, k in zip(arr, ok) if k]
        e["n_failed"] = int((~ok).sum()); e["n_valid"] = int(ok.sum())
    ev["note"] = "arrival lists keep only valid fits (not at a bound, amplitude > 0); n_failed counts the others; raw file: event_estimate_raw.json. Fits use the whole record (retrospective, not causal)."
    json.dump(ev, open(L / names["event_estimate"], "w"), indent=2)
    (L / "latest_run.txt").write_text(str(run_dir.resolve()) + "\n")
    print("inverse: latest copies in", L.resolve())
    return L


def stage_inverse(cfg, D, quick):
    inv = cfg["inverse"]
    files = cfg["forward"]["files"]
    for n in (files["A"], files["noise"], files["patch"], files["u_default"], *M5_EXTRA):
        need(D / n, "inverse")
    _find(D, M5_NAMES["H_bip"], ["H_matrix.npz"], [])          # stops with the candidates if neither H file exists
    check_names(cfg)
    base = Path(cfg["paths"]["outputs_dir"]) / inv.get("subdir", "inverse")
    base.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-m", "src.inverse.run_member5", "--data", str(D), "--out", str(base), "--config", cfg["paths"]["model_config"]]
    if inv.get("target_x_mm"):
        cmd += ["--target-x", str(inv["target_x_mm"])]
    if inv.get("seeds"):
        cmd += ["--seeds", str(inv["seeds"])]
    cmd += [str(a) for a in (inv.get("extra_args") or [])]
    if quick:
        cmd.append("--quick")
    print("inverse:", " ".join(cmd))
    t0 = time.time()
    if subprocess.run(cmd, cwd=ROOT).returncode != 0:
        raise SystemExit("[inverse] src/inverse/run_member5.py failed (see its message above)")
    runs = sorted(p for p in base.glob("run_*") if p.is_dir() and p.stat().st_mtime >= t0 - 2)
    if not runs:
        raise SystemExit(f"[inverse] no new run_<time> folder found in {base.resolve()}")
    run_dir = runs[-1]
    for name in inv["outputs"].values():
        need(run_dir / name, "inverse")
    print("inverse: run folder", run_dir.resolve())
    export_latest(cfg, run_dir)


def stage_control(cfg, quick):
    ctl = cfg["control"]
    if not ctl.get("enabled", True):
        print("control: disabled in config.yaml (control.enabled: false)")
        return
    L = latest_dir(cfg)
    names = cfg["inverse"]["outputs"]
    rec = need(L / names["reconstruction"], "control (run the inverse stage first)").resolve()
    argv = ["--config", ctl["config"], "--model-config", cfg["paths"]["model_config"], "--reconstruction", str(rec)]
    if ctl["use_event_estimate"]:
        argv += ["--event-estimate", str(need(L / names["event_estimate"], "control").resolve())]
    else:
        argv.append("--no-event")
    argv += [str(a) for a in (ctl.get("extra_args") or [])]
    if quick:
        argv.append("--quick")
    print("control (HYPOTHETICAL study):", " ".join(argv))
    try:
        from src.control.run_control import main      # imported here so the other stages do not need the control code
    except ImportError as e:
        raise SystemExit(f"[control] cannot import src.control.run_control: {e}")
    main(argv)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stages", nargs="+", default=list(ALL_STAGES), choices=ALL_STAGES)
    ap.add_argument("--config", default="configs/config.yaml")
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()

    cfg = load_cfg(a.config)
    D = data_dir(cfg)
    res = None
    for s in ALL_STAGES:                          # always in pipeline order
        if s not in a.stages:
            continue
        if s == "model":
            res = stage_model(cfg, D)
        elif s == "forward":
            stage_forward(cfg, D, res)
        elif s == "inverse":
            stage_inverse(cfg, D, a.quick)
        elif s == "control":
            stage_control(cfg, a.quick)
    print("pipeline done:", ", ".join(s for s in ALL_STAGES if s in a.stages))


if __name__ == "__main__":
    main()
