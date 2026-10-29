"""
csd_wave_features_v4.py - Steps 6-7 of the Member 2 task (S1 rat CSD recordings).

Place at: src/data/csd_wave_features_v4.py   (next to quicklook_s1.py, which it imports)
Usage:    python csd_wave_features_v4.py "<path to .smr>"
Needs:    pip install neo matplotlib scipy

Reads ONE recording (raw file is never modified). Everything is saved in a folder called
"csd_results" NEXT TO THIS SCRIPT (so you always know where it is), zipped into one file,
and that folder is opened automatically on Windows:
  <name>_v4_waves.csv        one row per detected CSD wave
  <name>_v4_exclusions.csv   time windows excluded, with reasons
  <name>_v4_waves_QC.png     picture of every detected wave - CHECK THESE BY EYE
  <name>_v4_summary.txt      the text printed on screen
  <name>_v4_results.zip      all four files above in one zip

Method (all thresholds are in PARAMS below):
  - DC and Doppler are averaged to 1 Hz.
  - Excluded: saturated samples, time before the amplifier note + settle time, time after "cull".
  - DC baseline = rolling median (10 min); waves = peaks of the baseline-corrected, smoothed DC.
  - rise_s / decay_s follow Pomierny et al. 2025 (min within 100 s before the peak -> peak;
    peak -> min within 350 s after). These are noise-sensitive, so also reported:
    fwhm_s (width at half amplitude) and rise_10pct_s / decay_10pct_s (time from/to the
    10 % amplitude crossing). Prefer fwhm_s and the 10 % times for model timescales.
  - Narrow spikes (KCl-application artifacts) are removed by a 15 s running median first.
  - flag = near_start_settling: wave is within 10 min of the analysis start, so the baseline
    may not have settled; treat with caution / exclude from summary statistics.
  - Markers: the Keyboard channel holds ASCII key codes; typed text is rebuilt from them.
    Text starting 'zol' = injection, 'cul' = end of experiment, 'mov' = amplifier adjusted.
"""
import csv
import os
import sys
import zipfile
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter, uniform_filter1d
from scipy.signal import find_peaks

OUT_TAG = "v4"   # added to every output file name, so new results never overwrite old ones

PARAMS = dict(
    baseline_win_s=601,   # rolling-median window for DC baseline (s, odd)
    despike_s=15,         # running-median window that removes narrow spikes, e.g. KCl artifacts (s)
    smooth_s=5,           # smoothing window for DC (s)
    prom_k=6.0,           # wave prominence >= prom_k * noise sigma ...
    prom_min_v=0.03,      # ... and at least this many ADC volts
    min_gap_s=180,        # minimum spacing between waves (s)
    settle_s=60,          # extra seconds excluded after the amplifier note
    sat_v=4.99,           # |V| at or above this = ADC saturation
    sat_pad_s=10,         # padding around saturated samples (s)
    gap_s=5.0,            # keystrokes closer than this belong to one typed note
)


# ---------- helpers -------------------------------------------------------
def to_1hz(t, v, factor=10):
    n = (len(v) // factor) * factor
    return t[:n].reshape(-1, factor).mean(1), v[:n].reshape(-1, factor).mean(1)


def decode_markers(markers, gap=5.0):
    """markers: [(time_s, '109'), ...] -> [(t_start, t_end, 'text'), ...]"""
    clusters, cur, last = [], [], None
    for t, lab in markers:
        if last is not None and t - last > gap:
            clusters.append(cur)
            cur = []
        cur.append((t, lab))
        last = t
    if cur:
        clusters.append(cur)
    out = []
    for c in clusters:
        txt = "".join(chr(int(l)) if l.strip().isdigit() and 31 < int(l) < 127 else "?" for _, l in c)
        out.append((c[0][0], c[-1][0], txt))
    return out


def runs(mask, t):
    """contiguous True runs -> [(t_start, t_end)]"""
    out, start = [], None
    for i, m in enumerate(mask):
        if m and start is None:
            start = i
        if not m and start is not None:
            out.append((t[start], t[i - 1]))
            start = None
    if start is not None:
        out.append((t[start], t[-1]))
    return out


def wave_features(t, sm, k):
    amp = sm[k]
    lo = max(0, k - 100)
    i1 = lo + int(np.argmin(sm[lo:k + 1]))
    hi = min(len(sm), k + 351)
    i2 = k + int(np.argmin(sm[k:hi]))
    half = amp / 2.0
    left = k
    while left > 0 and sm[left] > half:
        left -= 1
    right = k
    while right < len(sm) - 1 and sm[right] > half:
        right += 1
    tenth = 0.1 * amp
    l10 = k
    while l10 > 0 and sm[l10] > tenth:
        l10 -= 1
    r10 = k
    while r10 < len(sm) - 1 and sm[r10] > tenth:
        r10 += 1
    return dict(amp=float(amp), rise=float(t[k] - t[i1]), decay=float(t[i2] - t[k]),
                fwhm=float(t[right] - t[left]),
                rise10=float(t[k] - t[l10]), decay10=float(t[r10] - t[k]))


# ---------- core analysis (no file access; testable) ------------------------
def analyze(dc10, dop10, t10, markers, p=PARAMS):
    t, dc = to_1hz(t10, dc10)
    _, dop = to_1hz(t10, dop10)
    notes = decode_markers(markers, p["gap_s"])
    note = lambda pre: next((n for n in notes if n[2].lower().startswith(pre)), None)
    n_mov, n_inj, n_end = note("mov"), note("zol"), note("cul")

    excl = []
    valid = np.ones_like(t, dtype=bool)
    sat = (np.abs(dc) >= p["sat_v"]) | (np.abs(dop) >= p["sat_v"])
    sat = uniform_filter1d(sat.astype(float), size=2 * p["sat_pad_s"] + 1) > 0
    valid &= ~sat
    for a, b in runs(sat, t):
        excl.append((a, b, "ADC saturation (DC or Doppler at +-5 V)"))
    start = t[0]
    if n_mov:
        start = n_mov[1] + p["settle_s"]
        excl.append((t[0], start, "before/around amplifier adjustment (note: %r) + settle time" % n_mov[2]))
    else:
        start = t[0] + 120
        excl.append((t[0], start, "initial 2 min (no amplifier note found)"))
    valid &= t >= start
    end = n_end[0] if n_end else t[-1]
    if n_end:
        excl.append((end, t[-1], "after end-of-experiment note (%r)" % n_end[2]))
    valid &= t <= end
    inj_t = n_inj[0] if n_inj else None

    dcf = np.interp(t, t[valid], dc[valid])
    base = median_filter(dcf, size=p["baseline_win_s"], mode="nearest")
    sm = uniform_filter1d(median_filter(dcf - base, size=p["despike_s"], mode="nearest"),
                          size=p["smooth_s"], mode="nearest")
    med = np.median(sm[valid])
    sigma = 1.4826 * np.median(np.abs(sm[valid] - med))
    prom = max(p["prom_k"] * sigma, p["prom_min_v"])
    masked = np.where(valid, sm, med)
    peaks, _ = find_peaks(masked, prominence=prom, distance=p["min_gap_s"])

    dopf = np.interp(t, t[valid], dop[valid])
    dsm = uniform_filter1d(dopf, size=11, mode="nearest")
    rows = []
    for n, k in enumerate(peaks, 1):
        f = wave_features(t, sm, k)
        a, b = max(0, k - 30), min(len(t), k + 300)
        kd = a + int(np.argmax(dsm[a:b]))
        c0, c1 = max(0, k - 180), max(1, k - 30)
        dbase = float(np.median(dsm[c0:c1]))
        rows.append(dict(
            wave_id=n,
            phase=("pre_injection" if (inj_t is None or t[k] < inj_t) else "post_injection"),
            dc_peak_time_s=round(float(t[k]), 1),
            dc_amp_V=round(f["amp"], 4), rise_s=f["rise"], decay_s=f["decay"], fwhm_s=f["fwhm"],
            rise_10pct_s=f["rise10"], decay_10pct_s=f["decay10"],
            doppler_peak_time_s=round(float(t[kd]), 1),
            dc_to_doppler_delay_s=round(float(t[kd] - t[k]), 1),
            doppler_amp_V=round(float(dsm[kd] - dbase), 3),
            noise_sigma_V=round(float(sigma), 5), prominence_threshold_V=round(float(prom), 4),
            flag=("near_start_settling" if t[k] - start < p["baseline_win_s"] else "")))
    info = dict(notes=notes, injection_s=inj_t, valid_start_s=start, valid_end_s=end)
    return rows, excl, info, (t, sm, dsm, peaks)


# ---------- script entry ------------------------------------------------------
def main(path_str):
    from quicklook_s1 import open_reader, read_markers, read_channels_downsampled
    path = Path(path_str)
    stem = path.stem
    base_dir = Path(__file__).resolve().parent / "csd_results"
    out_dir = base_dir / f"{stem}_{OUT_TAG}"
    out_dir.mkdir(parents=True, exist_ok=True)
    lines = []

    def log(msg=""):
        print(msg)
        lines.append(str(msg))

    reader = open_reader(path)
    reader.parse_header()
    markers = read_markers(reader)
    ch = read_channels_downsampled(reader, target_hz=10.0)
    rows, excl, info, (t, sm, dsm, peaks) = analyze(ch["CSD-DC"]["v"], ch["doppler"]["v"],
                                                    ch["CSD-DC"]["t"], markers)
    log("Typed notes found:")
    for a, b, txt in info["notes"]:
        log(f"  {a:8.1f}-{b:8.1f} s  {txt!r}")
    log(f"\nInjection at {info['injection_s']} s | analysis window "
        f"{info['valid_start_s']:.0f}-{info['valid_end_s']:.0f} s")
    log(f"\nDetected {len(rows)} waves:")
    for r in rows:
        log("  #{wave_id:<2} {phase:<14} t={dc_peak_time_s:>7.1f}s amp={dc_amp_V:.3f}V rise={rise_s:.0f}s "
            "decay={decay_s:.0f}s fwhm={fwhm_s:.0f}s on10={rise_10pct_s:.0f}s off10={decay_10pct_s:.0f}s "
            "delay_to_doppler={dc_to_doppler_delay_s:.0f}s {flag}".format(**r))

    f_waves = out_dir / f"{stem}_{OUT_TAG}_waves.csv"
    f_excl = out_dir / f"{stem}_{OUT_TAG}_exclusions.csv"
    f_png = out_dir / f"{stem}_{OUT_TAG}_waves_QC.png"
    f_txt = out_dir / f"{stem}_{OUT_TAG}_summary.txt"
    if rows:
        with open(f_waves, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    with open(f_excl, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "start_s", "end_s", "reason"])
        for a, b, why in excl:
            w.writerow([path.name, round(a, 1), round(b, 1), why])

    ncol = 4
    nrow = max(1, int(np.ceil(len(peaks) / ncol)))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 2.8 * nrow), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    for ax, k, r in zip(axes.ravel(), peaks, rows):
        ax.axis("on")
        a, b = max(0, k - 400), min(len(t), k + 400)
        ax.plot(t[a:b], sm[a:b], lw=1)
        ax.set_title(f"#{r['wave_id']} {r['phase'][:4]} {t[k]/60:.1f} min {r['flag']}", fontsize=9)
        ax2 = ax.twinx()
        ax2.plot(t[a:b], dsm[a:b], color="tab:orange", lw=0.8, alpha=0.7)
        ax.axvline(t[k], color="r", alpha=0.3)
    fig.suptitle(f"{path.name}: DC (blue, baseline-corrected V) and Doppler (orange)")
    fig.tight_layout()
    fig.savefig(f_png, dpi=110)
    plt.close(fig)

    zip_path = base_dir / f"{stem}_{OUT_TAG}_results.zip"
    log(f"\nAll results are in: {out_dir}")
    log(f"ZIP file (send me this one): {zip_path}")
    with open(f_txt, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
        for f in (f_waves, f_excl, f_png, f_txt):
            if f.exists():
                z.write(f, arcname=f.name)
    try:
        if hasattr(os, "startfile"):      # Windows: open the folder for you
            os.startfile(str(base_dir))
    except Exception:
        pass


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
