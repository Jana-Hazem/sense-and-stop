"""Verification tests for src/model (Member 3, v0.2).  Run:  python -m pytest tests/ -v

Groups:  A) numerical correctness (exact solution, convergence, solver agreement, BCs)
         B) agreement with Member 2's S1 data features (excitable model)
         C) calibration tool and measurement functions
         D) interface contract and control hook
"""
import os
import numpy as np
import pytest

from src.model import (simulate, load_config, with_resolution, measure_front_speed,
                       save_result, load_result, pulse_metrics, pulse_speed, calibrate_excitable)

CFG_PATH = os.path.join(os.path.dirname(__file__), "..", "configs", "model_default.yaml")


@pytest.fixture
def cfg():
    """Default config = excitable model, 1D (the primary test system)."""
    return load_config(CFG_PATH)


@pytest.fixture
def kpp(cfg):
    cfg["model"]["type"] = "fisher_kpp"
    return cfg


# =============================================================== A. numerical correctness
def test_kpp_stays_in_unit_interval(kpp):
    u = simulate(kpp)["u"]
    assert u.min() > -1e-6 and u.max() < 1 + 1e-6


def test_kpp_exact_travelling_wave(kpp):
    """Ablowitz-Zeppetella: u = [1+exp(sqrt(r/6D)(x-x0-ct))]^-2, c = 5 sqrt(D r/6) is an EXACT
    solution of Fisher-KPP. The numerical error must be small AND fall at 2nd order when dx is halved."""
    D, r_ = kpp["model"]["fisher_kpp"]["D"], kpp["model"]["fisher_kpp"]["r"]
    kpp["initial_condition"].update(type="kpp_exact_wave", x0=30.0)
    kpp["domain"]["Lx"] = 80.0
    kpp["time"]["T"] = 300.0
    c = 5.0 * np.sqrt(D * r_ / 6.0)
    lam = np.sqrt(r_ / (6.0 * D))
    errs = []
    for dx in (0.2, 0.1):
        r = simulate(with_resolution(kpp, dx=dx, dt=0.05))
        exact = (1.0 + np.exp(lam * (r["x"][None, :] - 30.0 - c * r["t"][:, None]))) ** -2
        errs.append(np.abs(r["u"] - exact).max())
    assert errs[1] < 1e-3                               # measured ~4.9e-4 at dx = 0.1 mm
    assert errs[0] / errs[1] == pytest.approx(4.0, rel=0.1)   # second order: error / 4 when dx / 2
    assert measure_front_speed(r, t_window=(50, 300)) == pytest.approx(c, rel=0.002)


def test_kpp_speed_below_cstar_and_increasing(kpp):
    """A pulled front from a steep initial condition approaches c* = 2 sqrt(D r) from BELOW."""
    D, r_ = kpp["model"]["fisher_kpp"]["D"], kpp["model"]["fisher_kpp"]["r"]
    cstar = 2 * np.sqrt(D * r_)
    r = simulate(kpp)
    c1 = measure_front_speed(r, t_window=(100, 300))
    c2 = measure_front_speed(r, t_window=(300, 500))
    assert c1 < c2 < cstar


def test_kpp_speed_approaches_cstar_at_late_time(kpp):
    kpp["domain"].update(Lx=120.0, dx=0.2, dy=0.2)
    kpp["time"].update(T=1800.0)
    r = simulate(kpp)
    D, r_ = kpp["model"]["fisher_kpp"]["D"], kpp["model"]["fisher_kpp"]["r"]
    cstar = 2 * np.sqrt(D * r_)
    early = measure_front_speed(r, t_window=(100, 400))
    late = measure_front_speed(r, t_window=(1000, 1600))
    assert late > early
    assert late == pytest.approx(cstar, rel=0.05)


def test_spatial_convergence_second_order(kpp):
    kpp["time"].update(T=200.0)
    ref = simulate(with_resolution(kpp, dx=0.05, dt=0.05))["u"]
    errs = []
    for dx, f in [(0.4, 8), (0.2, 4), (0.1, 2)]:
        u = simulate(with_resolution(kpp, dx=dx, dt=0.05))["u"]
        errs.append(np.abs(u - ref[:, ::f]).max())
    assert errs[0] > errs[1] > errs[2]
    assert np.log2(errs[1] / errs[2]) > 1.7            # ~2nd order in space


def test_temporal_convergence_euler_first_order(kpp):
    kpp["time"].update(T=100.0, output_dt=10.0)
    kpp["domain"]["dx"] = 0.2
    ref = simulate(with_resolution(kpp, dt=0.005, method="rk4"))["u"]
    errs = [np.abs(simulate(with_resolution(kpp, dt=dt, method="euler"))["u"] - ref).max()
            for dt in (0.8, 0.4, 0.2)]
    assert errs[0] > errs[1] > errs[2]
    assert np.log2(errs[1] / errs[2]) > 0.8            # ~1st order in time


def test_independent_solver_agrees(cfg):
    """Explicit RK4 vs scipy solve_ivp (adaptive): needed to avoid the inverse crime."""
    cfg["time"].update(T=300.0)
    a = simulate(with_resolution(cfg, method="rk4"))["u"]
    b = simulate(with_resolution(cfg, method="solve_ivp"))["u"]
    assert np.abs(a - b).max() < 5e-3


def test_stability_guard_reduces_dt(cfg):
    cfg["time"].update(T=50.0, dt=5.0, output_dt=5.0)  # absurdly large target step
    cfg["solver"]["method"] = "euler"
    r = simulate(cfg)
    assert r["meta"]["stability_number"] <= 0.5 + 1e-12
    assert np.isfinite(r["u"]).all()


def test_neumann_conserves_mass_for_pure_diffusion(kpp):
    kpp["model"]["fisher_kpp"]["r"] = 0.0
    kpp["initial_condition"].update(type="gaussian", center=[10.0, 0.0], width=2.0)
    r = simulate(kpp)
    mass = np.trapezoid(r["u"], r["x"], axis=1)
    assert np.abs(mass / mass[0] - 1).max() < 1e-10


def test_2d_planar_wave_equals_1d(cfg):
    cfg["time"].update(T=100.0)
    cfg["domain"].update(Lx=20.0, Ly=2.0, dx=0.5, dy=0.5)
    r1 = simulate(cfg)
    cfg2 = with_resolution(cfg)
    cfg2["model"]["dimension"] = 2
    r2 = simulate(cfg2)
    assert np.abs(r2["u"] - r1["u"][:, None, :]).max() < 1e-10


# =============================================================== B. agreement with S1 data
def _excitable_run(cfg):
    r = simulate(cfg)
    return r, pulse_metrics(r, 15.0), pulse_speed(r, 10.0, 20.0)


def test_excitable_pulse_matches_s1_timing_and_speed(cfg):
    """S1 (rat, PRELIMINARY): FWHM animal means 36-56 s (mean 42.8); rise10 25-46 s.
    Literature speed (narrow range) 2.3-2.9 mm/min; implied width 1.7-2.1 mm."""
    r, m, c = _excitable_run(cfg)
    assert 36.0 <= m["fwhm_s"] <= 56.0
    assert m["fwhm_s"] == pytest.approx(42.8, rel=0.03)
    assert 25.0 <= m["rise_10pct_s"] <= 46.0
    assert 2.3 <= c * 60 <= 2.9
    assert 1.6 <= c * m["fwhm_s"] <= 2.2


def test_excitable_decay_is_known_to_be_shorter_than_s1(cfg):
    """KNOWN MISMATCH, recorded on purpose: model decay10 ~28 s is at the LOW edge of the S1
    animal-mean range (28-47 s; mean 40.7 s). The model recovers faster than the data pulse."""
    _, m, _ = _excitable_run(cfg)
    assert 25.0 <= m["decay_10pct_s"] <= 47.0
    assert m["decay_10pct_s"] < 40.7


def test_excitable_pulse_recovers_and_is_bounded(cfg):
    r = simulate(cfg)
    u, x, t = r["u"], r["x"], r["t"]
    i = np.searchsorted(t, 400.0)
    peak_x = x[u[i].argmax()]
    behind = u[i][x < peak_x - 6.0]
    assert behind.max() < 0.1                     # no sustained elevated state: a pulse, not a front
    assert u.max() < 1.05 and -0.4 < u.min() < -0.1   # refractory undershoot (not confined to [0,1])


def test_excitable_resolution_converged_at_default_dx(cfg):
    """dx = 0.1 mm must resolve the front: halving dx changes FWHM and speed by < 2 %."""
    a = simulate(cfg)
    b = simulate(with_resolution(cfg, dx=0.05))
    ma, mb = pulse_metrics(a, 15.0), pulse_metrics(b, 15.0)
    assert ma["fwhm_s"] == pytest.approx(mb["fwhm_s"], rel=0.02)
    assert pulse_speed(a, 10, 20) == pytest.approx(pulse_speed(b, 10, 20), rel=0.02)


# =============================================================== C. calibration + measurement
def test_pulse_metrics_on_gaussian_pulse():
    """Exact check of the measurement definitions: Gaussian with sigma = 20 s."""
    t = np.arange(0.0, 400.0, 0.5)
    s = np.exp(-((t - 200.0) ** 2) / (2 * 20.0 ** 2))
    res = {"u": s[:, None], "x": np.array([5.0]), "t": t}
    m = pulse_metrics(res, 5.0)
    assert m["fwhm_s"] == pytest.approx(2.3548 * 20.0, rel=2e-3)
    assert m["rise_10pct_s"] == pytest.approx(np.sqrt(2 * np.log(10)) * 20.0, rel=2e-3)
    assert m["decay_10pct_s"] == pytest.approx(m["rise_10pct_s"], rel=1e-3)
    assert m["t_peak_s"] == pytest.approx(200.0, abs=0.1)


def test_calibration_hits_requested_targets(cfg):
    new, info = calibrate_excitable(cfg, fwhm_s=50.0, speed_mm_per_min=4.0)
    new = with_resolution(new, dx=min(0.1, info["suggested_max_dx_mm"]))
    new["time"].update(T=360.0)
    r, m, c = _excitable_run(new)
    assert m["fwhm_s"] == pytest.approx(50.0, rel=0.05)
    assert c * 60 == pytest.approx(4.0, rel=0.05)
    assert info["implied_pulse_width_mm"] == pytest.approx(4.0 / 60 * 50.0)


def test_calibration_raises_when_no_pulse_exists(cfg):
    cfg["model"]["excitable"].update(a=0.3, eps=0.2)      # no sustained pulse from the edge stimulus
    with pytest.raises(RuntimeError):
        calibrate_excitable(cfg, fwhm_s=43.0, speed_mm_per_min=2.6)


# =============================================================== D. contract + control
def test_output_contract_1d(kpp):
    r = simulate(kpp)
    nt, nx = len(r["t"]), len(r["x"])
    assert r["u"].shape == (nt, nx)
    assert r["v"] is None and r["y"] is None
    assert r["t"][0] == 0 and r["t"][-1] == kpp["time"]["T"]
    assert np.allclose(np.diff(r["t"]), kpp["time"]["output_dt"])
    assert r["x"][0] == 0 and r["x"][-1] == pytest.approx(kpp["domain"]["Lx"])
    for key in ("units", "params", "meta"):
        assert key in r
    assert r["meta"]["stability_number"] <= 0.69


def test_output_contract_2d_excitable(cfg):
    cfg["model"].update(dimension=2)
    cfg["domain"].update(Lx=10.0, Ly=5.0, dx=0.5, dy=0.5)
    cfg["time"].update(T=20.0)
    r = simulate(cfg)
    assert r["u"].shape == (len(r["t"]), len(r["y"]), len(r["x"]))     # u[t, y, x]
    assert r["v"].shape == r["u"].shape


def test_control_hook_default_off_and_effective(cfg):
    cfg["time"].update(T=400.0)
    base = simulate(cfg)
    zero = simulate(cfg, control=lambda x, t: 0.0 * x)
    assert np.abs(base["u"] - zero["u"]).max() < 1e-14
    block = simulate(cfg, control=lambda x, t: 1.0 * ((x > 10) & (x < 14)))
    assert block["meta"]["control_active"] is True
    assert block["u"][-1][block["x"] > 16].max() < 0.5 * base["u"][-1][base["x"] > 16].max()  # wave blocked


def test_save_load_roundtrip(cfg, tmp_path):
    cfg["time"]["T"] = 50.0
    r = simulate(cfg)
    f = str(tmp_path / "run.npz")
    save_result(r, f)
    r2 = load_result(f)
    assert np.array_equal(r["u"], r2["u"]) and np.array_equal(r["v"], r2["v"])
    assert r2["params"]["model_type"] == "excitable" and r2["units"]["t"] == "s"
