"""Speed benchmark for the optimise() function with warm-start."""
import sys, time; sys.path.insert(0, '.')
import numpy as np, yaml
from pathlib import Path
from src.control.model import Grid, load_model, initial_state, simulate
from src.control.objective import design_from_cfg, search_bounds
from src.control.search import optimise
from src.control.estimate import make_nominal, peak_pos

cfg = yaml.safe_load(Path('configs/control.yaml').read_text())
s_params = dict(cfg['search'])
s_params.update(popsize=4, maxiter=6, polish_maxiter=20)
print('Search params:', s_params)

model, mcfg = load_model(Path('configs/model_default.yaml'))
z = np.load('outputs/inverse/run_20261119_122348/reconstruction.npz', allow_pickle=True)
t, x, truth = np.asarray(z['t'], float), np.asarray(z['x'], float), np.asarray(z['u_truth'], float)
grid = Grid(x, t)
u0, _ = initial_state(mcfg, grid, truth[0])
d = design_from_cfg(cfg, float(t[-1]))
prot = (x >= d.prot_lo) & (x <= d.prot_hi)

# Reference run
ref = simulate(model, grid, u0, prot, None, p=d.p, store=True, dt_target=0.5)
Unm = ref['u'][0]
D0_nom = float(ref['damage'][0])
nom = make_nominal(Unm, x, t)
p_front = nom.at(d.t_obs - 0.0)
bnd = search_bounds(d, p_front)
print('D0_nom=%.3f, p_front=%.2fmm' % (D0_nom, p_front))

# Time optimise
t0 = time.time()
r = optimise(model, grid, u0, prot, d, bnd, 0.0, D0_nom, s_params, seed_offset=0)
t1 = time.time()
print('optimise time: %.1fs, J_plan=%.4f, nfev_batches=%d' % (t1-t0, r['J_plan'], r['nfev_batches']))
print('P=[xc=%.2f, ts=%.1f, dur=%.1f, g0=%.4f]' % tuple(r['P']))
