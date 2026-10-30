# PROVISIONAL: paths below are hard-coded to the author's machine; set them before running.
import mne, numpy as np
from pathlib import Path

SUBJECTS_DIR = Path(r'D:\mne_data')
SCRATCH = Path(r'D:\sense_stop_scratch')
EDF = Path(r'D:\chb_mit\chb01_01.edf')
SUBJECT = 'fsaverage'
MONTAGE = 'standard_1020'
SFREQ = 256.0
MINDIST_MM = 5.0
CHECK_SECONDS = 60
DROP_NAMES = ('P7-T7', 'T8-P8-1')     # sign-flipped copy; duplicate (checked below)
CASE_FIX = {'FP1': 'Fp1', 'FP2': 'Fp2', 'FZ': 'Fz', 'CZ': 'Cz', 'PZ': 'Pz'}
ELEC_NAMES_PER_BYTE = 16
UNIT_BYTES = 8
HEADER_BYTES = 256

# ---- 1. Units straight from the EDF header ----
with open(EDF, 'rb') as f:
    main = f.read(HEADER_BYTES)
    ns = int(main[252:256])
    sig = f.read(ns * HEADER_BYTES)
labels = [sig[i*ELEC_NAMES_PER_BYTE:(i+1)*ELEC_NAMES_PER_BYTE].decode().strip() for i in range(ns)]
unit_off = ns * (ELEC_NAMES_PER_BYTE + 80)
units = [sig[unit_off + i*UNIT_BYTES: unit_off + (i+1)*UNIT_BYTES].decode().strip() for i in range(ns)]
print('EDF physical units:', sorted(set(units)))

# ---- 2. Are the two T8-P8 channels identical? ----
raw = mne.io.read_raw_edf(EDF, preload=True, verbose=False)
names = list(raw.ch_names)
a, b = names.index('T8-P8-0'), names.index('T8-P8-1')
seg = raw.get_data(start=0, stop=int(CHECK_SECONDS * SFREQ))
print('T8-P8 duplicates: max abs diff %.3e V, std of channel %.3e V' %
      (np.abs(seg[a] - seg[b]).max(), seg[a].std()))
p7t7, t7p7 = names.index('P7-T7'), names.index('T7-P7')
print('P7-T7 + T7-P7 max abs (should be ~0): %.3e' % np.abs(seg[p7t7] + seg[t7p7]).max())

# ---- 3. 21-electrode forward ----
pairs = [n.replace('-0', '') for n in names if n not in DROP_NAMES]
pairs = [p for p in pairs]
elec = []
for p in pairs:
    for e in p.split('-'):
        e = CASE_FIX.get(e, e)
        if e not in elec:
            elec.append(e)
print('electrodes needed:', len(elec), elec)

bem_dir = SUBJECTS_DIR / SUBJECT / 'bem'
src = mne.read_source_spaces(bem_dir / 'fsaverage-ico-5-src.fif', verbose=False)
bem_file = sorted(bem_dir.glob('*-bem-sol.fif'))[0]
info = mne.create_info(elec, SFREQ, 'eeg')
info.set_montage(MONTAGE, verbose=False)
fwd = mne.make_forward_solution(info, bem_dir / 'fsaverage-trans.fif', src, bem_file,
                                meg=False, eeg=True, mindist=MINDIST_MM, verbose=False)
fwd = mne.convert_forward_solution(fwd, surf_ori=True, force_fixed=True, verbose=False)
patch = np.load(SCRATCH / 'flat_patch.npz', allow_pickle=True)
S_vert = patch['vertno']
lh = fwd['src'][0]['vertno']
rows = np.searchsorted(lh, S_vert)
assert (lh[rows] == S_vert).all()
H_el = fwd['sol']['data'][:, rows].astype(np.float64)       # [n_elec, 73]

# ---- 4. Derivation matrix M and bipolar H ----
M = np.zeros((len(pairs), len(elec)))
for r, p in enumerate(pairs):
    pos, neg = [CASE_FIX.get(e, e) for e in p.split('-')]
    M[r, elec.index(pos)] = 1.0
    M[r, elec.index(neg)] = -1.0
H_bip = M @ H_el
print('H_bip shape:', H_bip.shape, ' M shape:', M.shape)
s = np.linalg.svd(H_bip, compute_uv=False)
print('rel singular values:', np.array2string(s / s[0], precision=2, max_line_width=100))
for tol in (1e-2, 1e-3, 1e-4):
    print('effective rank at', tol, ':', int((s / s[0] > tol).sum()))
np.savez(SCRATCH / 'H_bip.npz', H_bip=H_bip, M=M, electrodes=np.array(elec),
         channels=np.array(pairs), vertno=S_vert,
         units='V per A*m, bipolar (reference-free), fixed normal orientation')
print('saved', SCRATCH / 'H_bip.npz', '| channels:', pairs)