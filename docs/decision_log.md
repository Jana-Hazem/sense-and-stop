# Decision Log

## D1: Data source roles (date: 2026-10-28)

Three complementary sources are used in a computational study. They are not a synchronized dataset and do not constitute direct human CSD recordings.

## D2: Experimental CSD source: candidate located, verification pending

- Candidate: RepOD "Cortical spreading depolarisation electrophysiology data", DOI `10.18150/HIOPYS`, version 1.0, CC BY 4.0 (rat DC potential and laser Doppler recordings, KCl-induced CSD).
- Owner: Member 2. Gate: approximately two working days. Open a sample file; confirm usable CSD traces, units, timing, and electrode arrangement.
- Pass: use recordings to inform PDE domain scale and parameter bounds.
- Fail: record the blocker and use published speed ranges as literature-derived bounds, labeled as such, never as measured data.
- Prefer saline-group files for baseline wave properties; zolpidem is used to suppress CSD.

## D3: CHB-MIT EEG (PhysioNet version 1.0.0, DOI `10.13026/C2K01R`, ODC-By 1.0)

- Used only as measured background/noise. These are pediatric epilepsy-monitoring data.
- Never label them as CSD or migraine data. Exclude files listed in `RECORDS-WITH-SEIZURES`.
- Download only a few baseline files plus each case's `chbNN-summary.txt` montage file.

## D4: Head model: MNE fsaverage template

- Modeling scaffold, not subject-specific anatomy.
- Assumption: adult template versus pediatric EEG subjects (ages 1.5-22). Source locations must not be over-interpreted.
- Open item: confirm the license of the downloaded template data files. MNE-Python's BSD-3-Clause license covers the software, not necessarily the data.

## D5: No clinical claims

## D6 update: forward operator filled in (2026-10-30)

- Owner: Member 4. Interface touched: `u[t, y, x]` (Member 3), `H_bip` and `sensors` (Member 5). Details: `docs/interface_note_forward.md`.
- Open items from D6 now filled in:
  - Patch: `lateraloccipital-lh` of fsaverage (ico-5, mindist 5 mm). Flattening: geodesic distances, embedded in 2D, cut to 30 x 20 mm. Result: 73 sources (220 in the 25 mm geodesic disc), spacing 1.84 to 3.04 mm (mean 2.45).
  - Orientation: fixed to the surface normal.
  - `A`: sparse [73, 60501], rows sum to 1, node index `iy*301 + ix`.
  - `H_bip`: [21, 73], V per A*m, each channel is first electrode minus second.
  - Observation: `sensors[271, 21] = scale * H_bip @ (A @ u) + noise_gain * eps`, units V, `scale` illustrative (1e-9, 1e-8, 1e-7 A*m per source).
- Units and grid: mm and s, dx = 0.1 mm, output every 2 s, 0 to 540 s.
- Flattening distortion: flat/surface distance ratio on edges 0.92 to 1.27.
- Known limits: the default 2D wave is uniform in y; the first 60 s are still settling; the right-side channels (P4-O2, P8-O2) carry little information; `u` is not clipped; the same discretisation must not be used for generating data and for inversion.
- Not verified: FT9/FT10 rows, the u-to-volts scale, S3 data license.
- Remaining open (Member 3): grid spacing wording (0.1 mm default vs 0.096 mm in the work summary).
- Status: proposed; files for Member 5 are shared outside the repository.