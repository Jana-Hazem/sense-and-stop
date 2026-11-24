"""Tests for the event_estimate.json adapter (no model config needed).   Run:  python tests/test_events.py"""
import json, sys, tempfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from src.control.events import load_events

def _w(obj):
    p = Path(tempfile.mkdtemp()) / "event_estimate.json"; p.write_text(json.dumps(obj), encoding="utf-8"); return p

def test_member5_layout_single_value():
    r = load_events(_w({"per_scale": {"1": {"B": {"arrival_est_s": 301.5}}, "0.5": {"B": {"arrival_est_s": 305.0}}}}), "B", 15.0)
    assert r["per_scale"] == {"1": [301.5], "0.5": [305.0]} and r["target_x_mm"] == 15.0 and r["target_source"] == "config" and r["truth_arrival_s"] is None and not r["warnings"]

def test_per_seed_layout_and_json_overrides():
    r = load_events(_w({"target_x_mm": 18, "truth_arrival_s": 300, "per_scale": {"1": {"arrival_s_per_seed": [299, 301, None]}}}), "B", 15.0)
    assert r["target_x_mm"] == 18.0 and r["target_source"] == "json" and r["truth_arrival_s"] == 300.0 and len(r["per_scale"]["1"]) == 3 and r["per_scale"]["1"][2] != r["per_scale"]["1"][2]

def test_older_member5_key():
    r = load_events(_w({"per_scale": {"1e-07": {"arrival_s_at_target": 310.0}}}), "B", 15.0)
    assert r["per_scale"] == {"1e-07": [310.0]}

def test_unknown_layout_is_reported_not_silent():
    r = load_events(_w({"per_scale": {"1": {"A": {"x": 1}}}}), "B", 15.0)
    assert r["per_scale"] == {} and r["warnings"]

def test_missing_per_scale_is_reported():
    assert load_events(_w({"foo": 1}), "B", 15.0)["warnings"]

if __name__ == "__main__":
    fails = 0
    for n, f in sorted(globals().items()):
        if n.startswith("test_"):
            try: f(); print("PASS", n)
            except Exception as e: fails += 1; print("FAIL", n, repr(e)[:200])
    sys.exit(1 if fails else 0)
