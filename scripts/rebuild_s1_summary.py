"""Rebuild the S1 per-animal summary from per-wave CSVs (replaces the lost merge script).
Rule: baseline = pre_injection, no near_start_settling flag, not in the exclusion list.
Run from the repo root:  python scripts\\rebuild_s1_summary.py
"""
import csv
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
WAVES = ROOT / "outputs" / "CSD_results"
EXCL = ROOT / "docs" / "s1_manual_exclusions_reconstructed.csv"
REF = WAVES / "csd_per_animal_summary_clean.csv"
OUT_T = ROOT / "results" / "E2_rat_comparison" / "tables"
OUT_D = ROOT / "results" / "E2_rat_comparison" / "data"
OUT_T.mkdir(parents=True, exist_ok=True)
OUT_D.mkdir(parents=True, exist_ok=True)

# short id -> (file stem, name used in the old summary)
ANIMALS = {
    "Z_R4": ("26032024_Z_R4", "26032024_Z_R4"),
    "Z_R5": ("26032024_Z_R5.S2R", "26032024_Z_R5_S2R"),
    "Z_R6": ("27032024_Z_R6", "27032024_Z_R6"),
}
GROUP = {"pre_injection": "baseline(pre_injection)",
         "post_injection": "post_injection(separate)"}
FEATS = ["fwhm_s", "rise_10pct_s", "decay_10pct_s",
         "dc_to_doppler_delay_s", "dc_amp_V", "doppler_amp_V"]
TOL = {"dc_amp_V": 0.0011, "doppler_amp_V": 0.0011}   # old summary is rounded
DEFAULT_TOL = 0.006


def sd(x):
    return st.stdev(x) if len(x) > 1 else float("nan")


with open(EXCL, newline="") as f:
    excluded = {(r["animal"], int(r["wave_id"])) for r in csv.DictReader(f)}
print(f"exclusion list: {len(excluded)} waves")

kept = []   # (animal, phase, row)
for aid, (stem, _) in ANIMALS.items():
    path = WAVES / f"{stem}_v4" / f"{stem}_v4_waves.csv"
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    n_drop = 0
    for r in rows:
        if r["flag"] == "near_start_settling" or (aid, int(r["wave_id"])) in excluded:
            n_drop += 1
            continue
        kept.append((aid, r["phase"], r))
    print(f"{aid}: {len(rows)} detected, {n_drop} dropped, {len(rows) - n_drop} kept")

# ---- per animal and group summary
summary = []
for aid in ANIMALS:
    for phase in ("pre_injection", "post_injection"):
        sel = [r for a, p, r in kept if a == aid and p == phase]
        row = {"animal": aid, "group": GROUP[phase], "n_waves": len(sel)}
        for ft in FEATS:
            vals = [float(r[ft]) for r in sel]
            row[ft + "_mean"] = st.mean(vals)
            row[ft + "_sd"] = sd(vals)
        summary.append(row)

with open(OUT_T / "s1_per_animal_rebuilt.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0]))
    w.writeheader()
    for r in summary:
        w.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})

with open(OUT_D / "s1_baseline_waves.csv", "w", newline="") as f:
    cols = ["animal", "wave_id", "dc_peak_time_s"] + FEATS
    w = csv.writer(f)
    w.writerow(cols)
    for a, p, r in kept:
        if p == "pre_injection":
            w.writerow([a, r["wave_id"], r["dc_peak_time_s"]] + [r[ft] for ft in FEATS])

# ---- compare with Member 2's summary
with open(REF, newline="") as f:
    ref = {(r["animal"], r["group"]): r for r in csv.DictReader(f)}
name_of = {aid: nm for aid, (_, nm) in ANIMALS.items()}
bad, total = 0, 0
for r in summary:
    old = ref[(name_of[r["animal"]], r["group"])]
    total += 1
    if int(old["n_waves"]) != r["n_waves"]:
        bad += 1
        print(f"MISMATCH {r['animal']} {r['group']} n: rebuilt {r['n_waves']} vs old {old['n_waves']}")
    for ft in FEATS:
        for stat in ("mean", "sd"):
            a, b = r[f"{ft}_{stat}"], float(old[f"{ft}_{stat}"])
            if abs(a - b) > TOL.get(ft, DEFAULT_TOL):
                bad += 1
                print(f"MISMATCH {r['animal']} {r['group']} {ft}_{stat}: rebuilt {a:.4f} vs old {b:.4f}")
print(f"\ncompared {total} animal-groups x {len(FEATS)} features x (mean, sd): "
      + ("ALL MATCH" if bad == 0 else f"{bad} mismatches"))

# ---- overall baseline numbers (how the headline values are formed)
base = [(a, r) for a, p, r in kept if p == "pre_injection"]
print(f"\nBaseline overall, n = {len(base)} waves from {len(ANIMALS)} animals:")
for ft in ("fwhm_s", "rise_10pct_s", "decay_10pct_s"):
    allv = [float(r[ft]) for _, r in base]
    am = [st.mean([float(r[ft]) for a, r in base if a == aid]) for aid in ANIMALS]
    print(f"  {ft:14s} mean of animal means {st.mean(am):6.2f} | "
          f"all waves pooled {st.mean(allv):6.2f} | wave SD {sd(allv):5.2f}")
print(f"\nsaved: {OUT_T / 's1_per_animal_rebuilt.csv'}\nsaved: {OUT_D / 's1_baseline_waves.csv'}")