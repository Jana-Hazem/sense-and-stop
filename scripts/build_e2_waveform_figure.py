#!/usr/bin/env python
"""
E2 waveform figure for rat S1 recordings (Sense & Stop).

RUN:
    python build_e2_waveform_figure.py

PROJECT:
    D:\\sense-and-stop-data

RAW DATA:
    D:\\sense-and-stop-data\\raw\\S1-repod-csd\\

The script:
- Reads raw .smr files only (read-only).
- Processes recordings in chunks.
- Block-averages signals to approximately 1 Hz.
- Generates E2_waveforms.png and E2_waveforms.pdf.
- Writes the figure caption.
- Updates figure_index.csv.
- Reports mismatches without changing wave selections.
"""

import argparse
import csv
import re
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


# ============================================================
# 1. PROJECT CONFIGURATION
# ============================================================

ROOT = Path(r"D:\sense-and-stop-data").resolve()

RAW_DIR = ROOT / "raw" / "S1-repod-csd"
DATA_DIR = ROOT / "results" / "E2_rat_comparison" / "data"
FIG_DIR = ROOT / "results" / "E2_rat_comparison" / "figures"
DOCS_DIR = ROOT / "docs"
REPORTS_DIR = ROOT / "reports"


# ============================================================
# 2. RAT CONFIGURATION
# ============================================================

RATS = {
    "R4": ("Z_R4", "26032024_Z_R4.smr"),
    "R5": ("Z_R5", "26032024_Z_R5.S2R.smr"),
    "R6": ("Z_R6", "27032024_Z_R6.smr"),
}

CHANNELS = {
    "dc": "CSD-DC",
    "doppler": "doppler",
}

OUT_HZ = 1.0
WINDOW_S = 200.0
BASELINE_S = (-200.0, -60.0)

TASK_EXPECTED_TOTAL = {
    "R4": 10,
    "R5": 13,
}

EXPECTED_BASELINE = {
    "R4": [2, 3, 4, 5, 6],
    "R5": [6, 7, 8, 9],
    "R6": [6, 7, 8, 9, 10, 11],
}


# ============================================================
# 3. FIGURE STYLE
# ============================================================

STYLE = {
    "font.family": "sans-serif",
    "font.size": 9,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "figure.dpi": 100,
    "pdf.fonttype": 42,
}

C_DC = "#1f4e79"
C_DOP = "#d9822b"
C_BASE = "#111111"
C_EXCL = "#8a8a8a"
C_POST = "#c0392b"
C_SHADE = "#bdbdbd"


# ============================================================
# 4. VALIDATE PROJECT PATHS
# ============================================================

def validate_paths():
    required = [ROOT, RAW_DIR, DATA_DIR, DOCS_DIR]

    for folder in required:
        if not folder.exists():
            raise FileNotFoundError(
                f"\nRequired folder not found:\n{folder}\n"
                "Please check your project structure."
            )

    selection_file = DOCS_DIR / "s1_wave_selection.csv"

    if not selection_file.is_file():
        raise FileNotFoundError(
            f"Selection file not found:\n{selection_file}"
        )

    print("\nPROJECT PATHS")
    print("-" * 55)
    print("Project root:", ROOT)
    print("Raw data:", RAW_DIR)
    print("Processed data:", DATA_DIR)
    print("Selection file:", selection_file)
    print("\nProject paths verified successfully.")


# ============================================================
# 5. FILE HELPERS
# ============================================================

def find_one(folder, pattern):
    hits = sorted(folder.glob(pattern))

    if len(hits) != 1:
        raise FileNotFoundError(
            f"Expected exactly one file matching {pattern!r}\n"
            f"Folder: {folder}\n"
            f"Found: {hits}"
        )

    return hits[0]


def load_waves(path):
    rows = []

    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            rows.append({
                "wave_id": int(r["wave_id"]),
                "phase": r["phase"].strip(),
                "t_dc": float(r["dc_peak_time_s"]),
                "t_dop": float(r["doppler_peak_time_s"]),
                "delay": float(r["dc_to_doppler_delay_s"]),
                "fwhm": float(r["fwhm_s"]),
                "flag": (r.get("flag") or "").strip(),
            })

    return rows


def load_selection(path):
    selection = {}

    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            key = (r["animal"].strip(), int(r["wave_id"]))
            selection[key] = r["decision"].strip().lower()

    return selection


def load_exclusions(path):
    exclusions = []

    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            exclusions.append((
                float(r["start_s"]),
                float(r["end_s"]),
                r.get("reason", "")
            ))

    return exclusions


def load_markers(path):
    """Read typed notes from *_summary.txt."""
    pattern = re.compile(
        r"^\s*([0-9]+(?:\.[0-9]+)?)\s*-\s*"
        r"([0-9]+(?:\.[0-9]+)?)\s*s\s+(\S.*?)\s*$"
    )

    markers = []
    unknown = []

    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)

        if not match:
            continue

        note = match.group(3)

        if not note or note[0] not in "'\"":
            continue

        text = note[1:-1]
        low = text.lower()

        if low.startswith("kc"):
            kind = "kcl"
        elif low.startswith("zol"):
            kind = "inj"
        elif low.startswith("cul"):
            kind = "end"
        elif low.startswith("mov"):
            kind = "amp"
        else:
            kind = "other"
            unknown.append(text)

        markers.append((float(match.group(1)), kind, text))

    return markers, unknown


def check_gitignore():
    path = ROOT / ".gitignore"

    lines = (
        [
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
        ]
        if path.exists()
        else []
    )

    checks = {
        "data/raw/": any(
            line.rstrip("/").lstrip("/") == "data/raw"
            for line in lines
        ),
        "*.smr": "*.smr" in lines,
    }

    for name, passed in checks.items():
        print(
            f"{name:<15}: "
            f"{'OK' if passed else 'MISSING - add before committing'}"
        )

    return all(checks.values())


# ============================================================
# 6. RAW SPIKE2 FILE ACCESS
# ============================================================

def open_reader(path):
    from neo.rawio import Spike2RawIO

    if not path.is_file():
        raise FileNotFoundError(
            f"Raw recording not found:\n{path}\n"
            "Check the filename and raw data directory."
        )

    reader = Spike2RawIO(filename=str(path))
    reader.parse_header()
    return reader


def channel_info(reader, name):
    signals = reader.header["signal_channels"]
    streams = reader.header["signal_streams"]

    names = [str(channel["name"]) for channel in signals]

    if name not in names:
        raise KeyError(
            f"Channel {name!r} not found.\n"
            f"Available channels: {names}"
        )

    channel = signals[names.index(name)]

    stream_index = [
        str(stream["id"]) for stream in streams
    ].index(str(channel["stream_id"]))

    same_stream = [
        str(item["name"])
        for item in signals
        if str(item["stream_id"]) == str(channel["stream_id"])
    ]

    fs = float(
        reader.get_signal_sampling_rate(stream_index=stream_index)
    )

    return (
        stream_index,
        same_stream.index(name),
        fs,
        float(channel["gain"]),
        float(channel["offset"]),
    )


def read_decimated(
    reader,
    name,
    t0=None,
    t1=None,
    out_hz=OUT_HZ,
    chunk=4_000_000,
):
    """
    Read one channel in chunks and block-average it.
    Raw files are never modified.
    """
    stream_index, channel_index, fs, gain, offset = channel_info(
        reader, name
    )

    total = reader.get_signal_size(
        block_index=0,
        seg_index=0,
        stream_index=stream_index
    )

    i0 = 0 if t0 is None else max(0, int(round(t0 * fs)))
    i1 = total if t1 is None else min(total, int(round(t1 * fs)))

    if i1 <= i0:
        raise ValueError(f"Invalid time range: {t0} to {t1}")

    factor = max(1, int(round(fs / out_hz)))
    step = max(factor, (chunk // factor) * factor)

    parts = []

    for start in range(i0, i1, step):
        stop = min(start + step, i1)

        raw = reader.get_analogsignal_chunk(
            block_index=0,
            seg_index=0,
            i_start=start,
            i_stop=stop,
            stream_index=stream_index,
        )

        signal = (
            raw[:, channel_index].astype("float64") * gain + offset
        )

        usable = (len(signal) // factor) * factor

        if usable:
            averaged = signal[:usable].reshape(
                -1, factor
            ).mean(axis=1)
            parts.append(averaged)

    if not parts:
        raise ValueError(f"No data read from channel {name}")

    y = np.concatenate(parts)

    t = (
        i0 / fs
        + (np.arange(len(y)) + 0.5) * factor / fs
    )

    return t, y, fs


# ============================================================
# 7. PLOT HELPERS
# ============================================================

def valid_ylim(t, y, exclusions, lo_p=0.5, hi_p=99.5, pad=0.12):
    mask = np.ones_like(t, dtype=bool)

    for a, b, _ in exclusions:
        mask &= ~((t >= a) & (t <= b))

    values = y[mask] if mask.sum() > 10 else y
    lo, hi = np.percentile(values, [lo_p, hi_p])
    span = max(hi - lo, 1e-6)

    return lo - pad * span, hi + pad * span


def draw_marker_lines(ax, markers):
    for time, kind, _ in markers:
        if kind == "kcl":
            ax.axvline(
                time / 60,
                color="#555555",
                lw=0.6,
                ls=(0, (3, 3)),
                zorder=1
            )
        elif kind == "inj":
            ax.axvline(
                time / 60,
                color=C_POST,
                lw=1.3,
                zorder=2
            )
        elif kind == "end":
            ax.axvline(
                time / 60,
                color="black",
                lw=1.3,
                zorder=2
            )


# ============================================================
# 8. BUILD FIGURE
# ============================================================

def build_figure(D, reps, out_png, out_pdf):
    plt.rcParams.update(STYLE)

    n = len(RATS)

    fig = plt.figure(figsize=(11, 14))

    outer = fig.add_gridspec(
        2,
        1,
        height_ratios=[3.3, 1.5],
        left=0.08,
        right=0.97,
        top=0.925,
        bottom=0.055,
        hspace=0.20
    )

    # PANEL A: FULL RECORDINGS
    ga = GridSpecFromSubplotSpec(
        n, 1, subplot_spec=outer[0], hspace=0.32
    )

    first = None
    marker_counts = {}

    for i, rat in enumerate(RATS):
        d = D[rat]

        sub = GridSpecFromSubplotSpec(
            2, 1, subplot_spec=ga[i], hspace=0.07
        )

        ax_dc = fig.add_subplot(sub[0], sharex=first)

        if first is None:
            first = ax_dc

        ax_dop = fig.add_subplot(sub[1], sharex=first)

        ax_dc.plot(
            d["t_dc"] / 60,
            d["dc"],
            color=C_DC,
            lw=0.5,
            zorder=3
        )

        ax_dop.plot(
            d["t_dop"] / 60,
            d["dop"],
            color=C_DOP,
            lw=0.5,
            zorder=3
        )

        for ax in (ax_dc, ax_dop):
            for a, b, _ in d["excl"]:
                ax.axvspan(
                    a / 60,
                    b / 60,
                    color=C_SHADE,
                    alpha=0.40,
                    lw=0,
                    zorder=0
                )

            draw_marker_lines(ax, d["markers"])

        lo, hi = valid_ylim(d["t_dc"], d["dc"], d["excl"])
        ax_dc.set_ylim(lo, hi)

        ax_dop.set_ylim(
            *valid_ylim(d["t_dop"], d["dop"], d["excl"])
        )

        span = hi - lo

        counts = {
            "base": 0,
            "preexcl": 0,
            "post": 0,
            "offscale": 0
        }

        for w in d["waves"]:
            y = (
                float(np.interp(w["t_dc"], d["t_dc"], d["dc"]))
                + 0.07 * span
            )

            yc = min(
                max(y, lo + 0.04 * span),
                hi - 0.04 * span
            )

            counts["offscale"] += int(yc != y)

            kw = {
                "ls": "none",
                "ms": 5.5,
                "zorder": 6,
                "clip_on": False
            }

            if w["cls"] == "base":
                ax_dc.plot(
                    w["t_dc"] / 60,
                    yc,
                    marker="o",
                    mfc=C_BASE,
                    mec=C_BASE,
                    **kw
                )
            elif w["cls"] == "preexcl":
                ax_dc.plot(
                    w["t_dc"] / 60,
                    yc,
                    marker="o",
                    mfc="white",
                    mec=C_EXCL,
                    mew=1.0,
                    **kw
                )
            else:
                ax_dc.plot(
                    w["t_dc"] / 60,
                    yc,
                    marker="^",
                    mfc=C_POST,
                    mec=C_POST,
                    **kw
                )

            counts[w["cls"]] += 1

        marker_counts[rat] = counts

        ax_dc.set_ylabel("DC (V)")
        ax_dop.set_ylabel("Doppler (V)")

        ax_dc.set_title(rat, loc="left", fontweight="bold", pad=3)

        plt.setp(ax_dc.get_xticklabels(), visible=False)

        if i < n - 1:
            plt.setp(ax_dop.get_xticklabels(), visible=False)
        else:
            ax_dop.set_xlabel(
                "Time since start of recording (min)"
            )

    tmax = max(D[r]["t_dc"][-1] for r in RATS) / 60
    first.set_xlim(0, tmax * 1.01)

    # PANEL B: REPRESENTATIVE WAVES
    gb = GridSpecFromSubplotSpec(
        1, n, subplot_spec=outer[1], wspace=0.38
    )

    for i, rat in enumerate(RATS):
        r = reps[rat]

        sub = GridSpecFromSubplotSpec(
            2, 1, subplot_spec=gb[i], hspace=0.10
        )

        ax1 = fig.add_subplot(sub[0])
        ax2 = fig.add_subplot(sub[1], sharex=ax1)

        for ax, key, color in (
            (ax1, "dc", C_DC),
            (ax2, "dop", C_DOP)
        ):
            ax.plot(r["t"], r[key], color=color, lw=1.1)
            ax.axvline(0, color=C_DC, lw=0.9, ls="--")
            ax.axvline(r["delay"], color=C_DOP, lw=0.9, ls="--")
            ax.axhline(0, color="#bbbbbb", lw=0.5, zorder=0)

        ax1.set_xlim(-WINDOW_S, WINDOW_S)

        ax1.set_title(
            f"{rat} wave {r['wave_id']} "
            f"({r['t_peak'] / 60:.1f} min)",
            pad=3
        )

        ax1.set_ylabel("DC (V)")
        ax2.set_ylabel("Doppler (V)")
        ax2.set_xlabel("Time from DC peak (s)")

        plt.setp(ax1.get_xticklabels(), visible=False)

        lo2, hi2 = ax2.get_ylim()
        ya = hi2 - 0.12 * (hi2 - lo2)

        ax2.set_ylim(
            lo2,
            hi2 + 0.1 * (hi2 - lo2)
        )

        ax2.annotate(
            "",
            xy=(r["delay"], ya),
            xytext=(0, ya),
            arrowprops={
                "arrowstyle": "<->",
                "color": "black",
                "lw": 0.8
            }
        )

        ax2.text(
            r["delay"] / 2,
            ya + 0.03 * (hi2 - lo2),
            f"delay {r['delay']:.0f} s",
            ha="center",
            va="bottom",
            fontsize=8
        )

    # LEGEND
    handles = [
        Line2D(
            [], [], marker="o", ls="none",
            mfc=C_BASE, mec=C_BASE,
            label="baseline wave, kept"
        ),
        Line2D(
            [], [], marker="o", ls="none",
            mfc="white", mec=C_EXCL,
            label="pre-injection wave, excluded"
        ),
        Line2D(
            [], [], marker="^", ls="none",
            mfc=C_POST, mec=C_POST,
            label="post-injection wave"
        ),
        Line2D(
            [], [], color="#555555", lw=0.8,
            ls=(0, (3, 3)),
            label="KCl application (typed)"
        ),
        Line2D(
            [], [], color=C_POST, lw=1.3,
            label="zolpidem injection (typed)"
        ),
        Line2D(
            [], [], color="black", lw=1.3,
            label="end note (typed)"
        ),
        Patch(
            facecolor=C_SHADE,
            alpha=0.4,
            label="excluded interval"
        )
    ]

    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.985),
        ncol=4,
        frameon=False,
        fontsize=8
    )

    fig.text(0.012, 0.955, "(a)", fontsize=12, fontweight="bold")
    fig.text(0.012, 0.31, "(b)", fontsize=12, fontweight="bold")

    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)

    return marker_counts


# ============================================================
# 9. WRITE FIGURE CAPTION
# ============================================================

def write_caption(path, D, reps, counts):
    reps_txt = "; ".join(
        f"{rat}: wave {reps[rat]['wave_id']} "
        f"at {reps[rat]['t_peak'] / 60:.1f} min "
        f"(FWHM {reps[rat]['fwhm']:.0f} s, "
        f"median of kept baseline waves "
        f"{reps[rat]['median']:.1f} s, "
        f"DC to Doppler delay "
        f"{reps[rat]['delay']:.0f} s)"
        for rat in RATS
    )

    n_base = sum(counts[rat]["base"] for rat in RATS)

    caption = f"""**Figure E2_waveforms. Waveforms of the three rat recordings (R4, R5, R6).**

All three anesthetized rats received zolpidem during the recording; cortical spreading depolarization (CSD) was evoked with KCl. Each recording has one DC electrode channel and one Doppler (blood flow) channel.

**(a)** Full recordings, one rat per row, DC (top) and Doppler (bottom), against time since the start of the recording in minutes (shared time axis). Traces are block averages at about 1 Hz of the raw signal, in ADC volts before amplifier gain; the y-axes show the 0.5 to 99.5 percentile range of the valid data, so large artifacts are off scale. Vertical dashed lines mark typed KCl applications, the red line the typed zolpidem injection, the black line the typed end note (times from the hand-typed notes, so approximate). Grey shading marks intervals excluded by the extraction (ADC saturation, settling before the amplifier adjustment, time after the end note). Each detected wave is marked at its DC peak (marker drawn slightly above the trace; off-scale peaks are drawn at the axis edge): filled circle = baseline (pre-injection) wave kept (n = {n_base}: R4 {counts['R4']['base']}, R5 {counts['R5']['base']}, R6 {counts['R6']['base']}); open circle = pre-injection wave excluded (selection in docs/s1_wave_selection.csv, not a new rule); triangle = post-injection wave, a separate drug-condition group (all post-injection waves are shown, because the selection file covers pre-injection waves only).

**(b)** One representative baseline wave per rat, DC and Doppler in a window of +/- {WINDOW_S:.0f} s around the DC peak. Both traces are shifted so that their median over {BASELINE_S[0]:.0f} to {BASELINE_S[1]:.0f} s before the DC peak is zero. Dashed lines mark the DC peak and the Doppler peak (times from the waves table); the arrow shows the delay between them. **These waves are representative and were not selected for appearance.** Rule: among the kept baseline waves of a rat, the one whose FWHM is closest to the median FWHM of that rat's kept baseline waves (ties: the earlier wave). Selected: {reps_txt}.
"""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(caption, encoding="utf-8")


# ============================================================
# 10. UPDATE FIGURE INDEX
# ============================================================

def update_figure_index(path, rel_script, rel_data):
    default = [
        "figure_id",
        "description",
        "script",
        "data_file",
        "seed",
        "commit",
        "owner",
        "status"
    ]

    header = default
    rows = []

    if path.exists():
        with open(path, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            header = reader.fieldnames or default
            rows = list(reader)

    old = next(
        (
            row for row in rows
            if row.get("figure_id") == "E2_waveforms"
        ),
        {}
    )

    new = {
        "figure_id": "E2_waveforms",
        "description": (
            "Full recordings (DC and Doppler) of rats "
            "R4, R5, R6 with markers, excluded intervals "
            "and wave classes; one representative "
            "baseline wave per rat"
        ),
        "script": rel_script,
        "data_file": rel_data,
        "seed": "n/a (no randomness)",
        "commit": old.get("commit", ""),
        "owner": old.get("owner") or "Member 2",
        "status": old.get("status") or "draft"
    }

    rows = [
        row for row in rows
        if row.get("figure_id") != "E2_waveforms"
    ]
    rows.append(new)

    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=header,
            extrasaction="ignore"
        )
        writer.writeheader()

        for row in rows:
            writer.writerow({
                key: row.get(key, "")
                for key in header
            })


# ============================================================
# 11. MAIN FUNCTION
# ============================================================

def main():
    # IMPORTANT: global declarations must come before ROOT is used
    # anywhere inside this function.
    global ROOT, RAW_DIR, DATA_DIR
    global FIG_DIR, DOCS_DIR, REPORTS_DIR

    parser = argparse.ArgumentParser(
        description="Build E2 rat waveform figure"
    )

    parser.add_argument(
        "--root",
        default=str(ROOT),
        help="Repository root directory"
    )

    # Normal VS Code / PowerShell execution.
    args = parser.parse_args()

    ROOT = Path(args.root).resolve()

    RAW_DIR = ROOT / "raw" / "S1-repod-csd"
    DATA_DIR = ROOT / "results" / "E2_rat_comparison" / "data"
    FIG_DIR = ROOT / "results" / "E2_rat_comparison" / "figures"
    DOCS_DIR = ROOT / "docs"
    REPORTS_DIR = ROOT / "reports"

    # Validate folders
    print("\n" + "=" * 60)
    print("E2 WAVEFORM FIGURE GENERATION")
    print("=" * 60)

    validate_paths()

    print("\nCHECKING GITIGNORE")
    check_gitignore()

    selection = load_selection(DOCS_DIR / "s1_wave_selection.csv")

    D = {}
    problems = []

    # Load processed files
    print("\nLOADING WAVE TABLES AND METADATA")

    for rat, (animal, raw_name) in RATS.items():
        print(f"\nLoading {rat}...")

        waves = load_waves(
            find_one(DATA_DIR, f"*{rat}*_waves.csv")
        )

        markers, unknown = load_markers(
            find_one(DATA_DIR, f"*{rat}*_summary.txt")
        )

        exclusions = load_exclusions(
            find_one(DATA_DIR, f"*{rat}*_exclusions.csv")
        )

        if unknown:
            problems.append(
                f"{rat}: unrecognised typed notes: {unknown}"
            )

        # Assign wave classes.
        for w in waves:
            if w["phase"] == "post_injection":
                w["cls"] = "post"
            else:
                decision = selection.get((animal, w["wave_id"]))

                if decision is None:
                    problems.append(
                        f"{rat} wave {w['wave_id']}: "
                        "missing from selection file"
                    )
                    decision = "excluded"

                w["cls"] = (
                    "base" if decision == "kept" else "preexcl"
                )

        pre_ids = {
            w["wave_id"]
            for w in waves
            if w["phase"] == "pre_injection"
        }

        selection_ids = {
            wid
            for (an, wid) in selection
            if an == animal
        }

        if pre_ids != selection_ids:
            problems.append(
                f"{rat}: wave IDs differ between "
                f"CSV {sorted(pre_ids)} and "
                f"selection {sorted(selection_ids)}"
            )

        D[rat] = {
            "waves": waves,
            "markers": markers,
            "excl": exclusions,
            "raw": RAW_DIR / raw_name
        }

        print(
            f"  Waves: {len(waves)} | "
            f"Markers: {len(markers)} | "
            f"Excluded intervals: {len(exclusions)}"
        )

    # Open raw recordings
    print("\nOPENING RAW RECORDINGS")

    readers = {}
    sampling_rates = {}

    for rat in RATS:
        print(f"Opening {rat}...")

        readers[rat] = open_reader(D[rat]["raw"])

        sampling_rates[rat] = {
            key: channel_info(readers[rat], channel)[2]
            for key, channel in CHANNELS.items()
        }

        print(
            f"  DC: {sampling_rates[rat]['dc']:.0f} Hz"
            f" | Doppler: {sampling_rates[rat]['doppler']:.0f} Hz"
        )

    if len({
        tuple(sorted(rates.items()))
        for rates in sampling_rates.values()
    }) != 1:
        problems.append(
            f"Sampling rates differ between rats: {sampling_rates}"
        )

    # Read and decimate complete recordings
    print("\nREADING AND DECIMATING RECORDINGS")

    for rat in RATS:
        print(f"\nProcessing {rat}...")

        D[rat]["t_dc"], D[rat]["dc"], _ = read_decimated(
            readers[rat],
            CHANNELS["dc"]
        )

        D[rat]["t_dop"], D[rat]["dop"], _ = read_decimated(
            readers[rat],
            CHANNELS["doppler"]
        )

        print(f"  DC samples: {len(D[rat]['dc']):,}")
        print(f"  Doppler samples: {len(D[rat]['dop']):,}")

    # Select representative baseline waves
    print("\nSELECTING REPRESENTATIVE WAVES")

    reps = {}

    for rat in RATS:
        baseline = [
            w for w in D[rat]["waves"]
            if w["cls"] == "base"
        ]

        if not baseline:
            raise ValueError(f"{rat} has no kept baseline waves.")

        median_fwhm = float(np.median([
            w["fwhm"] for w in baseline
        ]))

        representative = min(
            baseline,
            key=lambda w: (
                abs(w["fwhm"] - median_fwhm),
                w["wave_id"]
            )
        )

        t0 = representative["t_dc"]

        # Read only the representative wave window.
        td, yd, _ = read_decimated(
            readers[rat],
            CHANNELS["dc"],
            t0 - WINDOW_S - 5,
            t0 + WINDOW_S + 5
        )

        tp, yp, _ = read_decimated(
            readers[rat],
            CHANNELS["doppler"],
            t0 - WINDOW_S - 5,
            t0 + WINDOW_S + 5
        )

        tt = np.arange(-WINDOW_S, WINDOW_S + 1)

        dc_signal = np.interp(tt + t0, td, yd)
        doppler_signal = np.interp(tt + t0, tp, yp)

        baseline_mask = (
            (tt >= BASELINE_S[0])
            & (tt <= BASELINE_S[1])
        )

        if not baseline_mask.any():
            raise ValueError(f"Empty baseline window for {rat}.")

        reps[rat] = {
            "wave_id": representative["wave_id"],
            "t_peak": t0,
            "delay": representative["delay"],
            "fwhm": representative["fwhm"],
            "median": median_fwhm,
            "t": tt,
            "dc": (
                dc_signal
                - np.median(dc_signal[baseline_mask])
            ),
            "dop": (
                doppler_signal
                - np.median(doppler_signal[baseline_mask])
            )
        }

        print(
            f"  {rat}: wave {representative['wave_id']} "
            f"| FWHM = {representative['fwhm']:.1f} s "
            f"| median = {median_fwhm:.1f} s"
        )

    # Generate figure
    print("\nBUILDING FIGURE")

    FIG_DIR.mkdir(parents=True, exist_ok=True)

    out_png = FIG_DIR / "E2_waveforms.png"
    out_pdf = FIG_DIR / "E2_waveforms.pdf"

    counts = build_figure(D, reps, out_png, out_pdf)

    # Write caption
    caption_path = REPORTS_DIR / "E2_waveforms_caption.md"

    write_caption(caption_path, D, reps, counts)

    # Update figure index
    index_path = ROOT / "results" / "figure_index.csv"

    update_figure_index(
        index_path,
        "scripts/build_e2_waveform_figure.py",
        "raw/S1-repod-csd/*.smr; "
        "results/E2_rat_comparison/data/*_waves.csv; "
        "docs/s1_wave_selection.csv"
    )

    # CROSS-CHECK
    print("\n" + "=" * 60)
    print("CROSS-CHECK")
    print("=" * 60)

    print(
        f"{'Rat':<6}"
        f"{'CSV':>8}"
        f"{'Marked':>9}"
        f"{'Pre':>7}"
        f"{'Post':>7}"
        f"{'Kept':>7}"
        f"{'Excluded':>10}"
        f"{'Off-scale':>11}"
    )

    total_kept = 0

    for rat in RATS:
        c = counts[rat]
        waves = D[rat]["waves"]

        marked = c["base"] + c["preexcl"] + c["post"]
        total_kept += c["base"]

        npre = sum(
            1 for w in waves
            if w["phase"] == "pre_injection"
        )

        print(
            f"{rat:<6}"
            f"{len(waves):>8}"
            f"{marked:>9}"
            f"{npre:>7}"
            f"{len(waves) - npre:>7}"
            f"{c['base']:>7}"
            f"{c['preexcl']:>10}"
            f"{c['offscale']:>11}"
        )

        if marked != len(waves):
            problems.append(
                f"{rat}: marked {marked} waves "
                f"but CSV contains {len(waves)}"
            )

        kept_ids = sorted(
            w["wave_id"]
            for w in waves
            if w["cls"] == "base"
        )

        if kept_ids != EXPECTED_BASELINE[rat]:
            problems.append(
                f"{rat}: kept baseline IDs {kept_ids} "
                f"differ from expected {EXPECTED_BASELINE[rat]}"
            )

        expected_total = TASK_EXPECTED_TOTAL.get(rat)

        if (
            expected_total is not None
            and expected_total != len(waves)
        ):
            problems.append(
                f"{rat}: expected {expected_total} waves, "
                f"CSV contains {len(waves)}"
            )

    print(
        f"\nTotal kept baseline waves: {total_kept} "
        f"(expected: 15)"
    )

    if total_kept != 15:
        problems.append(
            f"Total kept baseline waves = {total_kept}, expected 15"
        )

    print(
        "\nRepresentative waves:",
        {rat: reps[rat]["wave_id"] for rat in RATS}
    )

    # REPORT PROBLEMS
    print("\n" + "=" * 60)
    print("FINAL REPORT")
    print("=" * 60)

    if problems:
        print("\nMISMATCHES / PROBLEMS:")
        for problem in problems:
            print(" -", problem)
    else:
        print("\nNo mismatches found.")

    # OUTPUTS
    print("\n" + "=" * 60)
    print("FIGURE GENERATION COMPLETED")
    print("=" * 60)

    print("\nGenerated files:")
    print("\nPNG:")
    print(out_png)
    print("\nPDF:")
    print(out_pdf)
    print("\nCaption:")
    print(caption_path)
    print("\nFigure index:")
    print(index_path)

    print(
        "\nNext step: confirm the visual style against "
        "E2_fwhm_dotplot.png and add the commit hash "
        "after committing."
    )


# ============================================================
# 12. RUN
# ============================================================

if __name__ == "__main__":
    main()
