from pathlib import Path
exps = ["E1_solver", "E2_rat_comparison", "E3_forward_operator",
        "E4_inversion", "E5_control", "E6_reproducibility"]
for e in exps:
    for sub in ("figures", "tables", "data"):
        p = Path("results", e, sub)
        p.mkdir(parents=True, exist_ok=True)
        (p / ".gitkeep").touch()
Path("results/tables").mkdir(parents=True, exist_ok=True)
idx = Path("results/figure_index.csv")
if not idx.exists():
    idx.write_text("figure_id,experiment,description,script,data_file,seed,commit,owner,status\n")
print("done")