"""E2 (rat S1 side): per-animal table, FWHM dot plot, selection-sensitivity table.

Inputs : data/s1_waves/*_waves.csv (39 detected waves),
         data/s1_waves/csd_per_animal_summary_clean.csv (cross-check only)
Outputs: results/E2_rat_comparison/tables/E2_s1_per_animal.csv
         results/E2_rat_comparison/tables/E2_selection_sensitivity.csv
         results/E2_rat_comparison/figures/E2_fwhm_dotplot.png / .pdf
Run from the repo root:  python scripts/build_e2_tables_and_dotplot.py

Baseline set = the 15 pre-injection waves listed in BASELINE_IDS (identified by
matching the per-animal summary file you verified; asserted below).
SEED only affects the horizontal jitter of the dots; no statistic depends on it.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]          # script lives in <repo>/scripts/
OUT = REPO / "results" / "E2_rat_comparison"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
SEED = 20261119
rng = np.random.default_rng(SEED)

FILES = {"R4": "26032024_Z_R4_v4_waves.csv",
         "R5": "26032024_Z_R5_S2R_v4_waves.csv",
         "R6": "27032024_Z_R6_v4_waves.csv"}
BASELINE_IDS = {"R4": [2, 3, 4, 5, 6], "R5": [6, 7, 8, 9], "R6": [6, 7, 8, 9, 10, 11]}
# input CSVs: first folder that contains ALL the expected files
DATA = next((d for d in [OUT / "data", REPO / "data" / "s1_waves"]
             if all((d / f).exists() for f in FILES.values())), None)
if DATA is None:
    raise SystemExit("Input CSVs not found. Put these files in results/E2_rat_comparison/data/:\n  "
                     + "\n  ".join(list(FILES.values()) + ["csd_per_animal_summary_clean.csv"]))
print("Reading inputs from:", DATA)
M = {"fwhm_s": "FWHM (s)", "rise_10pct_s": "Rise 10% (s)", "decay_10pct_s": "Decay 10% (s)"}

frames = []
for a, f in FILES.items():
    d = pd.read_csv(DATA / f)
    d.insert(0, "animal", a)
    frames.append(d)
w = pd.concat(frames, ignore_index=True)
w["baseline"] = [i in BASELINE_IDS[a] for a, i in zip(w.animal, w.wave_id)]
pre = w[w.phase == "pre_injection"]
base = w[w.baseline]
assert len(w) == 39 and len(base) == 15 and base.phase.eq("pre_injection").all()

# ---- cross-check against the verified summary file
ref = pd.read_csv(DATA / "csd_per_animal_summary_clean.csv")
ref = ref[ref.group.str.startswith("baseline")].copy()
ref["animal"] = ref.animal.str.extract(r"(R\d)")[0]
for a in FILES:
    g = base[base.animal == a]
    r = ref[ref.animal == a].iloc[0]
    assert len(g) == r.n_waves
    for c in M:
        assert abs(g[c].mean() - r[c + "_mean"]) < 0.01, (a, c)
        assert abs(g[c].std(ddof=1) - r[c + "_sd"]) < 0.01, (a, c)

# ---- Table: per animal + across-animal rows
am = base.groupby("animal")[list(M)].mean()
rows = []
for a in FILES:
    g = base[base.animal == a]
    r = {"row": a, "n_waves": len(g)}
    for c, lab in M.items():
        r[lab + " mean"] = round(g[c].mean(), 2)
        r[lab + " SD (within animal)"] = round(g[c].std(ddof=1), 2)
        r[lab + " min-max"] = f"{g[c].min():.0f}-{g[c].max():.0f}"
    rows.append(r)
summ = [("Mean of animal means (n=3 animals)", lambda c: round(am[c].mean(), 2)),
        ("SD across animals (n=3)", lambda c: round(am[c].std(ddof=1), 2)),
        ("Range across animals", lambda c: f"{am[c].min():.2f}-{am[c].max():.2f}")]
for name, fn in summ:
    r = {"row": name, "n_waves": ""}
    for c, lab in M.items():
        r[lab + " mean"] = fn(c)
    rows.append(r)
r = {"row": "Pooled over waves (n=15) - reference only", "n_waves": 15}
for c, lab in M.items():
    r[lab + " mean"] = round(base[c].mean(), 2)
    r[lab + " SD (within animal)"] = round(base[c].std(ddof=1), 2)
rows.append(r)
tab = pd.DataFrame(rows).fillna("")
tab.to_csv(OUT / "tables/E2_s1_per_animal.csv", index=False)


# ---- Sensitivity of the three targets to the selection rule
def summarize(sub, label):
    a_ = sub.groupby("animal")[list(M)].mean()
    out = {"rule": label, "n_waves": len(sub)}
    for c in M:
        out[M[c] + " (mean of animal means)"] = round(a_[c].mean(), 2)
    out["FWHM SD across animals"] = round(a_["fwhm_s"].std(ddof=1), 2)
    return out


sens = pd.DataFrame([
    summarize(base, "A. Baseline used for calibration (by-eye selection)"),
    summarize(pre[pre.flag.isna()], "B. All pre-injection, excluding only the near_start_settling flag"),
    summarize(pre, "C. All pre-injection waves, no exclusion"),
])
sens.to_csv(OUT / "tables/E2_selection_sensitivity.csv", index=False)

# ---- Dot plot
col = {"R4": "#0072B2", "R5": "#E69F00", "R6": "#009E73"}
order = list(FILES)
grand = am.fwhm_s.mean()
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1, 1.15]})


def panel(axx, data, show_excl):
    for i, a in enumerate(order):
        g = data[data.animal == a]
        kept = g[g.baseline]
        j = rng.uniform(-0.13, 0.13, len(kept))
        axx.scatter(i + j, kept.fwhm_s, s=46, color=col[a], edgecolor="k", lw=.5, zorder=3,
                    label="Retained baseline wave" if i == 0 else None)
        if show_excl:
            ex = g[~g.baseline]
            j2 = rng.uniform(-0.13, 0.13, len(ex))
            axx.scatter(i + j2, ex.fwhm_s, s=46, facecolors="none", edgecolors="0.35", lw=1.1, zorder=2,
                        label="Excluded pre-injection wave" if i == 0 else None)
        axx.hlines(kept.fwhm_s.mean(), i - .3, i + .3, color="k", lw=2.2, zorder=4,
                   label="Animal mean (retained)" if i == 0 else None)
    axx.axhline(grand, color="0.45", ls="--", lw=1, zorder=1)
    axx.set_xticks(range(3))
    axx.set_xticklabels([f"{a}\n(n={int((base.animal == a).sum())} retained)" for a in order])
    axx.set_xlim(-.6, 2.6)
    axx.set_ylabel("FWHM of DC wave (s)")
    axx.grid(axis="y", alpha=.25)


panel(ax[0], base, False)
ax[0].set_ylim(0, 85)
ax[0].text(2.58, grand + 1.2, f"mean of animal means = {grand:.2f} s", ha="right", va="bottom",
           fontsize=8, color="0.3")
ax[0].set_title("(a) Retained baseline waves (n = 15)", fontsize=10)
panel(ax[1], pre, True)
ax[1].set_ylim(0, 215)
ax[1].set_title("(b) All pre-injection waves (n = 26); excluded waves open", fontsize=10)
ax[1].legend(fontsize=8, loc="upper left", frameon=False)
fig.suptitle("S1 rat DC-wave FWHM per wave, grouped by animal (3 rats, 1 electrode; preliminary)", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, .95))
fig.savefig(OUT / "figures/E2_fwhm_dotplot.png", dpi=200)
fig.savefig(OUT / "figures/E2_fwhm_dotplot.pdf")

pd.set_option("display.width", 250, "display.max_columns", 30)
print(tab.to_string())
print()
print(sens.to_string())
print("\npost-injection waves:", int((w.phase == "post_injection").sum()),
      "| excluded pre-injection:", len(pre) - len(base))