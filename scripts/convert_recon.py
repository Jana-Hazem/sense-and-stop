import numpy as np
z = np.load("outputs/member5_recon.npz", allow_pickle=True)
out = {"t": z["t"], "x": z["x"], "u_truth": z["truth"]}
for k in z.files:
    if k.startswith("rec_"):                 # rec_1e-09 -> rec_smooth_1e-09
        out["rec_smooth_" + k[4:]] = z[k]
np.savez_compressed("outputs/reconstruction.npz", **out)
print(sorted(out))