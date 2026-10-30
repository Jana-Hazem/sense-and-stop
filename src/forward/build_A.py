# PROVISIONAL: paths below are hard-coded to the author's machine; set them before running.
import mne, numpy as np
from pathlib import Path
from scipy.sparse import coo_matrix, csr_matrix, save_npz
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree

SUBJECTS_DIR = Path(r'D:\mne_data')
OUT_DIR = Path(r'D:\sense_stop_scratch')
LABEL_NAME = 'lateraloccipital-lh'
DISC_RADIUS_MM = 25.0
PATCH_X_MM = 30.0
PATCH_Y_MM = 20.0
DX_MM = 0.1                       # PDE grid spacing (Member 3 default)
M_TO_MM = 1000.0
CONST_TOL = 1e-12
LINEAR_TEST_TOL_MM = 0.5          # A @ (x field) vs source x, loose sanity bound

OUT_DIR.mkdir(exist_ok=True)
NX = int(round(PATCH_X_MM / DX_MM)) + 1    # 301
NY = int(round(PATCH_Y_MM / DX_MM)) + 1    # 201

# ---- 1. Same selection and flattening as the clean run ----
src = mne.read_source_spaces(SUBJECTS_DIR/'fsaverage'/'bem'/'fsaverage-ico-5-src.fif', verbose=False)
labels = mne.read_labels_from_annot('fsaverage', 'aparc', subjects_dir=SUBJECTS_DIR, verbose=False)
lab = [l for l in labels if l.name == LABEL_NAME][0]
h = src[0]
vn = np.intersect1d(h['vertno'], lab.vertices)
P = h['rr'][vn] * M_TO_MM
idx = {v: i for i, v in enumerate(vn)}
edges = set()
for t in h['use_tris']:
    for a, b in ((0, 1), (1, 2), (2, 0)):
        if t[a] in idx and t[b] in idx:
            edges.add(tuple(sorted((idx[t[a]], idx[t[b]]))))
e = np.array(sorted(edges))
d = np.linalg.norm(P[e[:, 0]] - P[e[:, 1]], axis=1)
G = coo_matrix((np.r_[d, d], (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])),
               shape=(len(vn),) * 2).tocsr()
c = np.argmin(np.linalg.norm(P - P.mean(0), axis=1))
sel = np.where(dijkstra(G, indices=c) <= DISC_RADIUS_MM)[0]
Gs = G[sel][:, sel]
ci = np.where(sel == c)[0][0]
sel = sel[np.isfinite(dijkstra(Gs, indices=ci))]
Gs = G[sel][:, sel]
D = dijkstra(Gs, directed=False)
n = len(sel)
J = np.eye(n) - 1.0 / n
ev, evec = np.linalg.eigh(-0.5 * J @ (D ** 2) @ J)
o = np.argsort(ev)[::-1][:2]
X = evec[:, o] * np.sqrt(ev[o])

# ---- 2. Fix axis signs: flat x ~ anterior, flat y ~ superior ----
P_sel = P[sel]
for axis, scanner_dim in ((0, 1), (1, 2)):
    r = np.corrcoef(X[:, axis], P_sel[:, scanner_dim])[0, 1]
    print(f'flat axis {axis} vs scanner axis {scanner_dim}: corr {r:.2f}')
    if r < 0:
        X[:, axis] *= -1

# ---- 3. Keep sources inside the 30 x 20 rectangle ----
mask = (np.abs(X[:, 0]) <= PATCH_X_MM / 2) & (np.abs(X[:, 1]) <= PATCH_Y_MM / 2)
S_xy = X[mask] + np.array([PATCH_X_MM / 2, PATCH_Y_MM / 2])   # patch coords: x 0..30, y 0..20
S_vert = vn[sel][mask]                                         # vertex numbers in hemisphere
S_pos = P_sel[mask]                                            # 3D position, mm
S_nn = h['nn'][S_vert]                                         # unit surface normals
ns = len(S_vert)
print('sources in patch:', ns)

# ---- 4. Area-averaging matrix A (Voronoi cells on the PDE grid) ----
xs = np.arange(NX) * DX_MM
ys = np.arange(NY) * DX_MM
gx, gy = np.meshgrid(xs, ys)                 # shape [NY, NX]
nodes = np.c_[gx.ravel(), gy.ravel()]        # index = iy * NX + ix (x fastest)
dist, owner = cKDTree(S_xy).query(nodes)
counts = np.bincount(owner, minlength=ns)
assert (counts > 0).all(), 'a source owns no PDE node'
A = csr_matrix((1.0 / counts[owner], (owner, np.arange(len(nodes)))), shape=(ns, len(nodes)))

# ---- 5. Checks ----
row_sums = np.asarray(A.sum(axis=1)).ravel()
print('A shape:', A.shape, ' nnz:', A.nnz)
print('max |row sum - 1|:', np.abs(row_sums - 1).max())
print('A @ ones - ones, max abs:', np.abs(A @ np.ones(len(nodes)) - 1).max(), '(tol', CONST_TOL, ')')
cell_area = counts * DX_MM ** 2
print('cell area mm^2: min %.1f mean %.1f max %.1f' % (cell_area.min(), cell_area.mean(), cell_area.max()))
print('max node-to-owner distance mm: %.2f' % dist.max())
print('linear test, max |A@x - source x| mm: %.2f (bound %.1f)' %
      (np.abs(A @ nodes[:, 0] - S_xy[:, 0]).max(), LINEAR_TEST_TOL_MM))
print('normals unit length, max error:', np.abs(np.linalg.norm(S_nn, axis=1) - 1).max())

# ---- 6. Save outside the repo ----
np.savez(OUT_DIR / 'flat_patch.npz', vertno=S_vert, xy_mm=S_xy, pos_mm=S_pos, normals=S_nn,
         hemi=0, label=LABEL_NAME, dx_mm=DX_MM, nx=NX, ny=NY)
save_npz(OUT_DIR / 'A_matrix.npz', A)
print('saved to', OUT_DIR)