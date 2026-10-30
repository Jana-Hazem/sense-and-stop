# PROVISIONAL: paths below are hard-coded to the author's machine; set them before running.
import mne, numpy as np
from pathlib import Path
from scipy.signal import decimate

CHB = Path(r'D:\chb_mit')
SCRATCH = Path(r'D:\sense_stop_scratch')
FILES = ['chb01_01.edf', 'chb01_02.edf', 'chb01_08.edf']
SFREQ = 256
PDE_DT_S = 2.0
TARGET_HZ = 1.0 / PDE_DT_S                         # 0.5 Hz
DECIM_STAGES = (8, 8, 8)                           # 256 -> 32 -> 4 -> 0.5 Hz (product 512)
DUP, SIGN_FLIP = 'T8-P8-1', 'P7-T7'
EDGE_TRIM_S = 60.0                                 # drop filter edge effects
HEADER_BYTES = 256
FIELD_BYTES = (16, 80, 8, 8, 8, 8, 8, 80, 8, 32)   # EDF per-signal header layout
PREFILTER_INDEX = 7

assert int(np.prod(DECIM_STAGES)) * TARGET_HZ == SFREQ

# ---- seizure check from the summary ----
lines = (CHB / 'chb01-summary.txt').read_text().splitlines()
for i, ln in enumerate(lines):
    if ln.startswith('File Name:') and ln.split(':')[1].strip() in FILES:
        nxt = [l for l in lines[i:i + 6] if 'Number of Seizures' in l]
        print(ln.strip(), '|', nxt[0].strip() if nxt else 'seizure line not found')

# ---- prefilter string from the header ----
with open(CHB / FILES[0], 'rb') as f:
    main = f.read(HEADER_BYTES)
    ns = int(main[252:256])
    sig = f.read(ns * HEADER_BYTES)
off = ns * sum(FIELD_BYTES[:PREFILTER_INDEX])
pre = {sig[off + i * 80: off + (i + 1) * 80].decode().strip() for i in range(ns)}
print('EDF prefilter field:', pre if pre != {''} else 'empty')

# ---- decimate ----
ref_names = list(np.load(SCRATCH / 'H_bip.npz', allow_pickle=True)['channels'])
series, raw_std = [], []
for fn in FILES:
    raw = mne.io.read_raw_edf(CHB / fn, preload=True, verbose=False)
    names = list(raw.ch_names)
    data = raw.get_data()                           # volts
    a = names.index('T8-P8-0')
    print(fn, 'T8-P8 duplicate max diff:', np.abs(data[a] - data[names.index(DUP)]).max())
    keep = [i for i, n in enumerate(names) if n not in (DUP, SIGN_FLIP)]
    got = [names[i].replace('-0', '') for i in keep]
    assert got == ref_names, (got, ref_names)
    x = data[keep]
    raw_std.append(x.std(axis=1))
    for q in DECIM_STAGES:
        x = decimate(x, q, axis=1, zero_phase=True)
    trim = int(EDGE_TRIM_S * TARGET_HZ)
    x = x[:, trim:-trim]
    series.append(x)
    print(fn, 'decimated shape', x.shape, ' duration s', x.shape[1] * PDE_DT_S)

noise = np.stack(series)                            # [file, channel, time]
raw_std = np.mean(raw_std, axis=0)
dec_std = noise.std(axis=2).mean(axis=0)
print('channel | raw std uV | decimated std uV | fraction')
for n, r, d in zip(ref_names, raw_std, dec_std):
    print('%-9s %9.2f %12.2f %10.3f' % (n, r * 1e6, d * 1e6, d / r))
print('overall decimated std uV: median %.2f max %.2f' % (np.median(dec_std) * 1e6, dec_std.max() * 1e6))
np.savez(SCRATCH / 'noise_2s.npz', noise_V=noise, channels=np.array(ref_names), dt_s=PDE_DT_S,
         files=np.array(FILES), note='anti-aliased decimation to 2 s; seizure-free chb01 files; edge-trimmed')
print('saved', SCRATCH / 'noise_2s.npz')