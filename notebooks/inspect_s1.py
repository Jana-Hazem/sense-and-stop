"""
inspect_s1.py - Step 3 of the Member 2 task: read-only header inspection of S1 files.

Place at: src/data/inspect_s1.py
Usage:    python src/data/inspect_s1.py data/raw/S1-repod-csd/<file>.smr[x]
Needs:    pip install neo        (.smrx files also need CED's `sonpy`)

Prints, for each file: duration, every analog channel (name, units, sampling rate,
gain/offset), and every event/marker channel (count + first timestamps).
It NEVER writes to the raw file. Copy the output into the data dictionary (S1 rows)
and cite "S1 file header" as the source.
"""
import sys
from pathlib import Path


def open_reader(path: Path):
    ext = path.suffix.lower()
    if ext == ".smrx":
        from neo.rawio import CedRawIO          # needs sonpy
        return CedRawIO(filename=str(path))
    from neo.rawio import Spike2RawIO            # .smr (32-bit CED format)
    return Spike2RawIO(filename=str(path))


def main(path_str: str) -> None:
    path = Path(path_str)
    reader = open_reader(path)
    reader.parse_header()

    t0 = reader.segment_t_start(block_index=0, seg_index=0)
    t1 = reader.segment_t_stop(block_index=0, seg_index=0)
    print(f"FILE: {path.name}")
    print(f"Segments: {reader.segment_count(0)} | t_start={t0:.3f} s | t_stop={t1:.3f} s "
          f"| duration={(t1 - t0) / 60:.1f} min\n")

    print("ANALOG CHANNELS")
    print(f"{'idx':>3}  {'name':<20} {'units':<8} {'fs (Hz)':>9} {'dtype':<8} {'gain':>10} {'offset':>10}")
    for i, ch in enumerate(reader.header["signal_channels"]):
        print(f"{i:>3}  {str(ch['name']):<20} {str(ch['units']):<8} "
              f"{float(ch['sampling_rate']):>9.2f} {str(ch['dtype']):<8} "
              f"{float(ch['gain']):>10.4g} {float(ch['offset']):>10.4g}")

    print("\nEVENT / MARKER CHANNELS")
    n_ev = len(reader.header["event_channels"])
    if n_ev == 0:
        print("  none found -> KCl / injection times are NOT in the file markers; check the paper/notes")
    for i, ev in enumerate(reader.header["event_channels"]):
        n = reader.event_count(block_index=0, seg_index=0, event_channel_index=i)
        line = f"  [{i}] {ev['name']}: {n} events"
        if n:
            ts, _, _ = reader.get_event_timestamps(block_index=0, seg_index=0,
                                                   event_channel_index=i, t_start=None, t_stop=None)
            secs = reader.rescale_event_timestamp(ts, dtype="float64")
            line += f" | first (s): {[round(float(s), 1) for s in secs[:6]]}"
        print(line)

    print("\nNEXT: record units, fs, channel roles and marker meaning in docs/data_dictionary.md (S1 rows).")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
