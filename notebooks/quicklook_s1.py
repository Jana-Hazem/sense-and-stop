"""
quicklook_s1.py - Step 3b of the Member 2 task: look at one S1 recording.

Place at: src/data/quicklook_s1.py
Usage:    python quicklook_s1.py "<path to .smr/.smrx>"
Needs:    pip install neo matplotlib

Does three things (read-only; the raw file is never modified):
 1. prints ALL Keyboard markers with their labels (what was typed at each time)
 2. prints min / max / mean of each channel (in ADC volts, before amplifier gain)
 3. saves <filename>_quicklook.png in the CURRENT folder: DC trace, doppler trace,
    and the marker times as vertical lines (signals averaged down to 10 Hz for plotting)
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def open_reader(path: Path):
    if path.suffix.lower() == ".smrx":
        from neo.rawio import CedRawIO
        return CedRawIO(filename=str(path))
    from neo.rawio import Spike2RawIO
    return Spike2RawIO(filename=str(path))


def read_markers(reader):
    out = []
    for i, ev in enumerate(reader.header["event_channels"]):
        n = reader.event_count(block_index=0, seg_index=0, event_channel_index=i)
        if not n:
            continue
        ts, _, labels = reader.get_event_timestamps(block_index=0, seg_index=0,
                                                    event_channel_index=i)
        secs = reader.rescale_event_timestamp(ts, dtype="float64")
        out.extend((float(s), str(l)) for s, l in zip(secs, labels))
    return sorted(out)


def read_channels_downsampled(reader, target_hz=10.0, chunk=5_000_000):
    """Return {channel_name: (time_s, volts)} averaged to ~target_hz, read in chunks."""
    sig = reader.header["signal_channels"]
    result = {}
    for s_idx, stream in enumerate(reader.header["signal_streams"]):
        fs = float(reader.get_signal_sampling_rate(stream_index=s_idx))
        names = [c["name"] for c in sig if c["stream_id"] == stream["id"]]
        gains = np.array([c["gain"] for c in sig if c["stream_id"] == stream["id"]], dtype="float64")
        offs = np.array([c["offset"] for c in sig if c["stream_id"] == stream["id"]], dtype="float64")
        factor = max(1, int(round(fs / target_hz)))
        n_total = reader.get_signal_size(block_index=0, seg_index=0, stream_index=s_idx)
        pieces, stats = [], None
        step = (chunk // factor) * factor
        for i0 in range(0, n_total, step):
            raw = reader.get_analogsignal_chunk(block_index=0, seg_index=0, i_start=i0,
                                                i_stop=min(i0 + step, n_total), stream_index=s_idx)
            volts = raw.astype("float32") * gains.astype("float32") + offs.astype("float32")
            mn, mx, sm = volts.min(0), volts.max(0), volts.sum(0, dtype="float64")
            stats = (mn, mx, sm, volts.shape[0]) if stats is None else (
                np.minimum(stats[0], mn), np.maximum(stats[1], mx), stats[2] + sm, stats[3] + volts.shape[0])
            k = (volts.shape[0] // factor) * factor
            if k:
                pieces.append(volts[:k].reshape(-1, factor, volts.shape[1]).mean(1))
        data = np.concatenate(pieces, axis=0)
        t = (np.arange(data.shape[0]) + 0.5) * factor / fs
        for j, name in enumerate(names):
            result[name] = dict(t=t, v=data[:, j], fs=fs, min=float(stats[0][j]),
                                max=float(stats[1][j]), mean=float(stats[2][j] / stats[3]))
    return result


def main(path_str):
    path = Path(path_str)
    reader = open_reader(path)
    reader.parse_header()

    markers = read_markers(reader)
    print(f"FILE: {path.name}\n\nKEYBOARD MARKERS ({len(markers)}): time_s  label")
    for t, lab in markers:
        print(f"  {t:9.1f}  {lab!r}")

    chans = read_channels_downsampled(reader)
    print("\nCHANNEL RANGES (ADC volts, before amplifier gain)")
    for name, d in chans.items():
        print(f"  {name:<10} fs={d['fs']:.0f} Hz  min={d['min']:.4f}  max={d['max']:.4f}  mean={d['mean']:.4f}")

    names = [n for n in ("CSD-DC", "doppler") if n in chans] or list(chans)[:2]
    fig, axes = plt.subplots(len(names), 1, figsize=(14, 3.2 * len(names)), sharex=True)
    axes = np.atleast_1d(axes)
    for ax, n in zip(axes, names):
        ax.plot(chans[n]["t"] / 60, chans[n]["v"], lw=0.6)
        ax.set_ylabel(f"{n} (V, 10 Hz avg)")
        for t, _ in markers:
            ax.axvline(t / 60, color="r", alpha=0.25, lw=0.8)
    axes[-1].set_xlabel("time (min)  -  red lines = keyboard markers")
    fig.suptitle(path.name)
    fig.tight_layout()
    out = Path.cwd() / f"{path.stem}_quicklook.png"
    fig.savefig(out, dpi=120)
    print(f"\nSaved plot: {out}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
