"""Adapter for Member 5's event_estimate.json -> one normalised structure for the controller.

Two layouts are understood:
  Member 5 (actual):  per_scale[s]["B"]["arrival_est_s"]  = one number (or a list)      key "B" is configurable
  per-seed (planned): per_scale[s]["arrival_s_per_seed"]  = list of numbers
  older Member 5:     per_scale[s]["arrival_s_at_target"] = one number
Optional top-level keys: target_x_mm, truth_arrival_s. If absent the caller falls back to the config target and to
the plant's own arrival time. Nothing is skipped silently: problems are returned in `warnings`.
"""
import json
from pathlib import Path


def _floats(v):
    if v is None: return []
    if isinstance(v, (int, float)): v = [v]
    out = []
    for q in v:
        try: out.append(float(q))
        except (TypeError, ValueError): out.append(float("nan"))
    return out


def load_events(path, method_key="B", default_target_x=15.0):
    ev = json.loads(Path(path).read_text(encoding="utf-8"))
    per = ev.get("per_scale"); warnings = []
    if not isinstance(per, dict) or not per:
        return dict(target_x_mm=default_target_x, target_source="config", truth_arrival_s=None, per_scale={}, layout="none",
                    warnings=[f"{path}: no 'per_scale' dictionary; top-level keys = {sorted(ev)}"])
    scales, layouts, failed = {}, set(), {}
    for sc, e in per.items():
        e = e if isinstance(e, dict) else {}
        if "arrival_s_per_seed" in e: vals = _floats(e["arrival_s_per_seed"]); layouts.add("per_seed")
        elif isinstance(e.get(method_key), dict) and "arrival_est_s" in e[method_key]: vals = _floats(e[method_key]["arrival_est_s"]); layouts.add(f"member5[{method_key}]")
        elif "arrival_s_at_target" in e: vals = _floats(e["arrival_s_at_target"]); layouts.add("member5_old")
        else:
            warnings.append(f"scale {sc}: no 'arrival_s_per_seed' and no '{method_key}.arrival_est_s'; keys found = {sorted(e)}"); continue
        scales[str(sc)] = vals
        k = e.get("n_failed", (int(e["n"]) - len(vals)) if e.get("n") is not None else 0)       # fits Member 5 dropped (not a valid arrival)
        failed[str(sc)] = max(int(k), 0)
    tx = ev.get("target_x_mm"); truth = ev.get("truth_arrival_s")
    return dict(target_x_mm=float(tx) if tx is not None else float(default_target_x), target_source="json" if tx is not None else "config",
                truth_arrival_s=float(truth) if truth is not None else None, per_scale=scales, n_failed=failed, layout="+".join(sorted(layouts)) or "none", warnings=warnings)
