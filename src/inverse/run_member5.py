"""Member 5 v2: Tikhonov (Method A) + travelling-pulse (Method B), first 60 s removed,
noisy channels optionally dropped, coarse 0.3 mm solver grid, drift stress test, and plots.
Computational study, not clinical."""
import sys, json
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.sparse import load_npz
from scipy.sparse.linalg import lsmr
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import least_squares
import matplotlib

# ------------------ EDIT THIS ------------------
DATA = Path(r"D:\CUFE\PDE\Sense And Stop\src\data")
DROP = ("P8-O2", "P4-O2")          # set to () to keep all 21 channels
SHOW = True                        # True: also open the figures on screen; False: only save PNGs
# -----------------------------------------------
if not SHOW:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("outputs"); OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(DATA))
from sensors import make_sensors

NT, NY, NX = 271, 201, 301
DT, DX, DXC = 2.0, 0.1, 0.3
SCALES = (1e-9, 1e-8, 1e-7)
T_CUT, X_MAX = 60.0, 25.0
SEED, RIDGE = 123, 1e-4
LAMS = np.geomspace(1e-1, 1e4, 11)
DRIFTS_UV = (0.0, 0.5, 2.0, 5.0)    # total drift over 540 s, per channel, random sign
LB = np.array([-10.0, 0.02, 0.4]); UB = np.array([10.0, 0.20, 4.0])
TRUE_SPEED = 2.6                    # mm/min (Member 3 default)

# ---------- load ----------
H = np.load(DATA / "H_bip.npz", allow_pickle=True)["H_bip"]
A = load_npz(DATA / "A_matrix.npz").tocsr()
noise_V = np.load(DATA / "noise_2s.npz", allow_pickle=True)["noise_V"]
cov = np.load(DATA / "noise_cov.npz", allow_pickle=True)
names = [str(c) for c in cov["channels"]]
keep = np.array([n not in DROP for n in names])
kept_names = [n for n in names if n not in DROP]
print("channels kept:", keep.sum(), "| dropped:", [n for n in names if n in DROP])
C = (cov["C"] * 1e12)[np.ix_(keep, keep)]                              # uV^2
u = np.load(DATA / "u_2d_default.npy").astype(np.float64)              # [t, y, x]
t = np.arange(NT) * DT; x = np.arange(NX) * DX
truth = u[:, 0, :]                                                     # wave is uniform in y
G1_unit = A.T.dot(H.T).reshape(NY, NX, -1).sum(axis=0)[:, keep]       # [NX, m]
m = int(keep.sum())
L = cholesky(C + 1e-3 * np.mean(np.diag(C)) * np.eye(m), lower=True)
Linv = solve_triangular(L, np.eye(m), lower=True)

ti = t >= T_CUT; tt = t[ti]; nt2 = len(tt)
xc = np.arange(0, 30.0 + 1e-9, DXC); nc = len(xc)
P = np.column_stack([np.interp(x, xc, e) for e in np.eye(nc)])        # coarse -> fine interpolation [NX, nc]

def front_pos(rec):
    return x[rec.argmax(1)]

def front_err(rec, lo, hi):
    pk_t = front_pos(truth); pk_r = front_pos(rec)
    mk = (t >= lo) & (t < hi) & (pk_t <= 24.0)
    return float(np.mean(np.abs(pk_r - pk_t)[mk])) if mk.any() else float("nan")

def scores(rec):
    xi = x <= X_MAX
    rel = np.linalg.norm((rec - truth)[ti][:, xi]) / np.linalg.norm(truth[ti][:, xi])
    return float(rel), front_err(rec, 60, 300), front_err(rec, 300, 541)

# ---------- Method A: Tikhonov on coarse grid, t >= 60 s only ----------
def d1(n):
    return sp.diags([-np.ones(n - 1), np.ones(n - 1)], [0, 1], shape=(n - 1, n)).tocsr()
Dx = sp.kron(sp.identity(nt2), d1(nc)).tocsr()
Dt = sp.kron(d1(nt2), sp.identity(nc)).tocsr()
Id = sp.identity(nt2 * nc).tocsr()
nreg = Dx.shape[0] + Dt.shape[0] + Id.shape[0]

def fit_A(Yuv, Gw):
    M = solve_triangular(L, (P.T @ Gw).T, lower=True)                  # [m, nc]
    Ab = sp.kron(sp.identity(nt2), sp.csr_matrix(M)).tocsr()
    b = np.concatenate([solve_triangular(L, Yuv[ti].T, lower=True).T.ravel(), np.zeros(nreg)])
    best = None
    for lam in LAMS:                                                   # ascending: keeps largest lam with rms <= 1
        Ast = sp.vstack([Ab, np.sqrt(lam) * Dx, np.sqrt(lam) * Dt, np.sqrt(RIDGE) * Id]).tocsr()
        xs = lsmr(Ast, b, atol=1e-8, btol=1e-8, maxiter=3000)[0]
        rms = np.sqrt(np.mean((Ab @ xs - b[:Ab.shape[0]]) ** 2))
        if best is None or rms <= 1.0:
            best = (lam, xs, rms)
    lam, xs, rms = best
    rec = np.zeros((NT, NX)); rec[ti] = (P @ xs.reshape(nt2, nc).T).T   # no clipping
    return rec, lam, rms

# ---------- Method B: travelling pulse on coarse grid, optional per-channel drift nuisance ----------
def fit_B(Yuv, Gw, nuis):
    yw = solve_triangular(L, Yuv[ti].T, lower=True)                    # [m, nt2]
    ttc = tt - tt.mean()
    D = [np.outer(Linv[:, i], ttc).ravel() for i in range(m)] if nuis else []
    def design(th):
        sc = np.exp(-0.5 * ((xc[None, :] - th[0] - th[1] * tt[:, None]) / th[2]) ** 2)
        f = solve_triangular(L, ((sc @ P.T) @ Gw).T, lower=True)
        X = np.column_stack([f.ravel()] + D)
        coef = np.linalg.lstsq(X, yw.ravel(), rcond=None)[0]
        return X, coef
    def resid(th):
        X, coef = design(th)
        return yw.ravel() - X @ coef
    best = None
    for x0 in (0.0, 3.0):
        for c in (0.03, 0.0433, 0.08, 0.12):
            for s in (0.8, 1.5):
                try:
                    r = least_squares(resid, [x0, c, s], bounds=(LB, UB), x_scale=[1, .02, .5])
                except Exception:
                    continue
                if best is None or r.cost < best.cost:
                    best = r
    th = best.x; a = design(th)[1][0]
    sc = np.exp(-0.5 * ((xc[None, :] - th[0] - th[1] * t[:, None]) / th[2]) ** 2)
    rec = a * (sc @ P.T)
    return rec, th, a, float(np.sqrt(2 * best.cost / yw.size))

# ---------- run ----------
rng = np.random.default_rng(SEED)
res, recA, recB, recB_drift = {}, {}, {}, {}
sens = {}                                                              # sensor data kept for plotting
print("true: x0=2.0, c=2.6 mm/min")
for scale in SCALES:
    Y, clean = make_sensors(u, H, A, scale, noise_V, rng)
    Y = (Y * 1e6)[:, keep]; clean_uV = (clean * 1e6)[:, keep]
    Gw = G1_unit * scale * 1e6
    sens[scale] = (Y.copy(), clean_uV.copy())
    rA, lam, rmsA = fit_A(Y, Gw)
    recA[scale] = rA
    sA = scores(rA)
    print("\nscale %.0e | A (Tikhonov): lam %.3g rms %.2f relL2 %.3f front err early %.2f late %.2f mm"
          % ((scale, lam, rmsA) + sA))
    res[f"{scale:.0e}"] = {"A": dict(lam=float(lam), rms=float(rmsA), relL2=sA[0],
                                     front_early=sA[1], front_late=sA[2]), "B": []}
    drift_dir = rng.uniform(-1, 1, size=m)
    for dv in DRIFTS_UV:
        Yd = Y + dv * drift_dir[None, :] * (t[:, None] / t[-1])        # slow linear drift
        for nuis in (False, True):
            rec, th, a, rms = fit_B(Yd, Gw, nuis)
            s = scores(rec)
            if dv == 0.0 and not nuis:
                recB[scale] = rec
            recB_drift[(scale, dv, nuis)] = rec
            print("  B drift %.1f uV nuisance=%-5s | x0 %.2f c %.2f mm/min sigma %.2f a %.2f | rms %.2f | "
                  "relL2 %.3f | front err early %.2f late %.2f mm"
                  % (dv, nuis, th[0], th[1] * 60, th[2], a, rms, s[0], s[1], s[2]))
            res[f"{scale:.0e}"]["B"].append(dict(
                drift_uV=dv, nuisance=nuis, x0=float(th[0]), c_mm_min=float(th[1] * 60),
                sigma=float(th[2]), amp=float(a), rms=rms, relL2=s[0], front_early=s[1], front_late=s[2]))

json.dump(dict(dropped=list(DROP), dx_inv=DXC, t_cut=T_CUT, results=res),
          open(OUT / "member5_v2_results.json", "w"), indent=1)
print("\nsaved outputs/member5_v2_results.json")

# =====================================================================
# VISUALISATION
# =====================================================================
ext = [x[0], x[-1], t[0], t[-1]]
kw = dict(aspect="auto", origin="lower", cmap="viridis", vmin=-0.3, vmax=1, extent=ext)

def mark(a):
    a.axhline(T_CUT, c="w", lw=.7, ls="--"); a.axvline(X_MAX, c="w", lw=.7, ls="--")

# Figure 1: space-time maps (truth, Method A, Method B)
fig, ax = plt.subplots(3, 4, figsize=(17, 11), sharex=True, sharey=True)
for r, (lab, rec_d) in enumerate([("Truth", None), ("Method A (Tikhonov)", recA), ("Method B (pulse)", recB)]):
    for c, scale in enumerate((None,) + SCALES):
        a = ax[r, c]
        if c == 0:
            if r == 0:
                im = a.imshow(truth, **kw); a.set_title("TRUTH u(x,t)")
            else:
                a.axis("off")
            continue
        if rec_d is None:
            a.axis("off"); continue
        im = a.imshow(rec_d[scale], **kw)
        a.set_title("%s, scale %.0e" % (lab, scale), fontsize=9)
    for a in ax[r]:
        if a.axison: mark(a)
for a in ax[-1]: a.set_xlabel("x (mm)")
for a in ax[:, 0]: a.set_ylabel("t (s)")
fig.colorbar(im, ax=ax, shrink=.6, label="u (dimensionless)")
fig.suptitle("Wave reconstruction (dashed: 60 s cut and 25 mm edge). Dark = undershoot below 0, not clipped.")
fig.savefig(OUT / "fig1_spacetime_maps.png", dpi=130, bbox_inches="tight")

# Figure 2: where is the front? (position of the maximum versus time)
fig, ax = plt.subplots(1, 3, figsize=(16, 4.5), sharey=True)
for a, scale in zip(ax, SCALES):
    a.plot(t, front_pos(truth), "k", lw=2.5, label="truth")
    a.plot(t[ti], front_pos(recA[scale])[ti], "tab:red", lw=1.2, label="Method A")
    a.plot(t[ti], front_pos(recB[scale])[ti], "tab:blue", lw=2, label="Method B")
    a.axhline(X_MAX, c="grey", ls="--", lw=.8); a.axvline(T_CUT, c="grey", ls="--", lw=.8)
    a.set_title("front position, scale %.0e" % scale); a.set_xlabel("t (s)"); a.legend()
ax[0].set_ylabel("x of the maximum (mm)")
fig.savefig(OUT / "fig2_front_tracking.png", dpi=130, bbox_inches="tight")

# Figure 3: snapshots (the wave as a profile along x) at three times
times = (120, 300, 480)
fig, ax = plt.subplots(3, 3, figsize=(16, 10), sharex=True, sharey=True)
for r, scale in enumerate(SCALES):
    for c, tv in enumerate(times):
        k = int(round(tv / DT)); a = ax[r, c]
        a.plot(x, truth[k], "k", lw=2.5, label="truth")
        a.plot(x, recA[scale][k], "tab:red", lw=1.2, label="A")
        a.plot(x, recB[scale][k], "tab:blue", lw=2, label="B")
        a.axhline(0, c="grey", lw=.5); a.axvline(X_MAX, c="grey", ls="--", lw=.8)
        a.set_title("scale %.0e, t = %d s" % (scale, tv), fontsize=9)
        if r == 0 and c == 0: a.legend()
for a in ax[-1]: a.set_xlabel("x (mm)")
for a in ax[:, 0]: a.set_ylabel("u")
fig.savefig(OUT / "fig3_snapshots.png", dpi=130, bbox_inches="tight")

# Figure 4: drift stress test (speed and front error versus drift size)
fig, ax = plt.subplots(2, 3, figsize=(16, 8), sharex=True)
for c, scale in enumerate(SCALES):
    rows = res[f"{scale:.0e}"]["B"]
    for nuis, col, lab in [(False, "tab:red", "no drift term"), (True, "tab:green", "with drift term")]:
        sel = [r_ for r_ in rows if r_["nuisance"] == nuis]
        d = [r_["drift_uV"] for r_ in sel]
        ax[0, c].plot(d, [r_["c_mm_min"] for r_ in sel], "o-", c=col, label=lab)
        ax[1, c].plot(d, [r_["front_late"] for r_ in sel], "o-", c=col, label=lab + " (late)")
        ax[1, c].plot(d, [r_["front_early"] for r_ in sel], "s--", c=col, alpha=.5, label=lab + " (early)")
    ax[0, c].axhline(TRUE_SPEED, c="k", ls="--", label="true speed")
    ax[0, c].set_title("scale %.0e" % scale); ax[1, c].set_xlabel("drift over 540 s (uV)")
    if c == 0:
        ax[0, c].set_ylabel("recovered speed (mm/min)"); ax[1, c].set_ylabel("mean front error (mm)")
        ax[0, c].legend(fontsize=8); ax[1, c].legend(fontsize=7)
fig.suptitle("Drift stress test, Method B")
fig.savefig(OUT / "fig4_drift_test.png", dpi=130, bbox_inches="tight")

# Figure 5: what the sensors see (strongest scale): clean, noisy, and the noise level per channel
scale = SCALES[-1]
Yn, Yc = sens[scale]
fig, ax = plt.subplots(1, 3, figsize=(17, 5.5))
vmax = np.abs(Yc).max()
for a, arr, ttl in [(ax[0], Yc, "clean sensors (uV)"), (ax[1], Yn, "noisy sensors (uV)")]:
    im = a.imshow(arr.T, aspect="auto", origin="lower", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                  extent=[t[0], t[-1], -.5, m - .5])
    a.axvline(T_CUT, c="k", ls="--", lw=.8); a.set_yticks(range(m)); a.set_yticklabels(kept_names, fontsize=6)
    a.set_title("%s, scale %.0e" % (ttl, scale)); a.set_xlabel("t (s)"); plt.colorbar(im, ax=a, shrink=.8)
ax[2].barh(range(m), np.sqrt(np.diag(C)), color="tab:grey", label="noise std")
ax[2].barh(range(m), np.abs(Yc).max(0), height=.4, color="tab:orange", label="peak clean signal")
ax[2].set_yticks(range(m)); ax[2].set_yticklabels(kept_names, fontsize=6); ax[2].set_xlabel("uV"); ax[2].legend()
ax[2].set_title("signal versus noise per channel")
fig.savefig(OUT / "fig5_sensors.png", dpi=130, bbox_inches="tight")

print("saved figures: fig1_spacetime_maps, fig2_front_tracking, fig3_snapshots, fig4_drift_test, fig5_sensors")
if SHOW:
    plt.show()