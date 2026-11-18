import os
import numpy as np
import yaml
from pathlib import Path
from scipy.sparse import load_npz

NT, NY, NX = 271, 201, 301
DT_S, DX_MM = 2.0, 0.1
SPEED_MM_PER_MIN = 2.6
FWHM_MM = 1.85
START_X_MM = 2.0
FWHM_TO_SIGMA = 2.3548
SCALES_AM = (1e-9, 1e-8, 1e-7)       # dipole moment per source at u = 1 (1, 10, 100 nAm)
V_TO_UV = 1e6
SEED = 0


def load_cfg(path="configs/config.yaml"):
    """Read the master config. Relative to the repo root: run scripts from there."""
    with open(path) as f:
        return yaml.safe_load(f)


def data_dir(cfg):
    """Folder with the derived files (H_bip.npz, A_matrix.npz, ...), kept outside Git.
    Override per machine with the environment variable SENSE_STOP_DATA."""
    d = Path(os.environ.get("SENSE_STOP_DATA", cfg["paths"]["data_dir"]))
    if not d.is_dir():
        raise FileNotFoundError(
            f"data folder not found: {d}. Set SENSE_STOP_DATA or paths.data_dir in configs/config.yaml")
    return d


def model_cfg(cfg):
    """M3's config with the forward side's override from config.yaml applied.
    model_default.yaml stays 1D for M3's tests; the dimension comes from config.yaml."""
    from src.model import load_config
    m = load_config(cfg["paths"]["model_config"])
    m["model"]["dimension"] = cfg["model"]["dimension"]
    return m


def check_grid(res, patch, noise_dt):
    """Fail loudly if M3's output does not match the grid that A was built for."""
    u, x, y, t = res["u"], res["x"], res["y"], res["t"]
    nx, ny, dx = int(patch["nx"]), int(patch["ny"]), float(patch["dx_mm"])
    assert u.ndim == 3, "forward side needs 2D u[t, y, x]; set model.dimension = 2"
    assert u.shape[1:] == (ny, nx), f"grid {u.shape[1:]} != A's grid {(ny, nx)}"
    assert np.allclose(x, np.arange(nx) * dx) and np.allclose(y, np.arange(ny) * dx), \
        "x/y nodes differ from A's grid"
    assert np.allclose(np.diff(t), noise_dt), "output_dt differs from the noise time step"


def make_sensors(u, H_bip, A, scale_am, noise_V, rng, noise_gain=1.0):
    """u: [t, y, x] (mm, s; not clipped). H_bip: [ch, src] V per A*m. A: [src, ny*nx].
    noise_V: [file, ch, time] from noise_2s.npz. Returns sensors [t, ch] in volts and the clean part."""
    nt = u.shape[0]
    clean = scale_am * (H_bip @ (A @ u.reshape(nt, -1).T)).T          # [t, ch]
    f = rng.integers(noise_V.shape[0])                                 # one file, one contiguous segment
    s = rng.integers(noise_V.shape[2] - nt + 1)
    noise = noise_V[f, :, s:s + nt].T                                  # [t, ch]
    return clean + noise_gain * noise, clean


if __name__ == '__main__':
    cfg = load_cfg()
    D = data_dir(cfg)
    files = cfg["forward"]["files"]
    H_bip = np.load(D / files["H_bip"], allow_pickle=True)['H_bip']
    A = load_npz(D / files["A"])
    noise_V = np.load(D / files["noise"], allow_pickle=True)['noise_V']
    rng = np.random.default_rng(SEED)

    t = np.arange(NT) * DT_S
    x = np.arange(NX) * DX_MM
    sigma = FWHM_MM / FWHM_TO_SIGMA
    front = START_X_MM + SPEED_MM_PER_MIN / 60.0 * t
    u1 = np.exp(-0.5 * ((x[None, :] - front[:, None]) / sigma) ** 2)  # [t, x]
    u = np.repeat(u1[:, None, :], NY, axis=1)                          # [t, y, x]
    print('u shape', u.shape, ' front reaches x =', round(front[-1], 1), 'mm (patch is 30 mm)')

    noise_std = noise_V.std(axis=2).mean(axis=0)
    for scale in SCALES_AM:
        y, clean = make_sensors(u, H_bip, A, scale, noise_V, rng)
        peak = np.abs(clean).max(axis=0)
        snr = peak / noise_std
        print('scale %.0e A*m: peak clean uV median %.3f max %.3f | peak/noise std: median %.3f max %.3f (%s)' %
              (scale, np.median(peak) * V_TO_UV, peak.max() * V_TO_UV,
               np.median(snr), snr.max(), 'channel ' + str(int(np.argmax(snr)))))
    print('sensors shape:', y.shape, '(time, channel)')