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

No document, code, or commit message claims diagnosis, detection, or treatment of migraine.
## D6: PDE grid to cortical-source mapping (proposed, 2026-10-30)

- Owners: Member 3, Member 4. Touches: grid, units, H/y formats.
- Decision proposed:
  - PDE domain is a flat 2D patch, 30 mm (x) by 20 mm (y), stored as `u[t, y, x]` in mm and s. 1D is for debugging only.
  - The patch is placed on the occipital cortex of the MNE fsaverage template. This is a design choice, not data-based.
  - Sources inside the patch receive an area-averaged field, `s = A u_flat`, where `A` is a sparse area-averaging matrix (PDE nodes to sources). Full forward map: `y = H A u + epsilon`.
  - The PDE grid stays at the solver's `dx` (currently 0.1 mm) and is not coarsened to match sensors.
  - The first end-to-end test is 2D.
- Reason: the pulse is about 1.85 mm wide, narrower than the expected source spacing, so point sampling would alias it.
- Assumptions: a folded cortical patch is approximated as flat; fsaverage is an adult template while the EEG is pediatric (see D4); `u` is dimensionless, so any scale to sensor amplitude is illustrative.
- Not changed: `u` is not clipped (it undershoots to about -0.3 behind the pulse); no clinical claims.
- Open (Member 4, to fill after loading fsaverage): patch coordinates and flattening method; source spacing and number of sources; orientation convention (fixed normal vs free); shapes and units of `H` and `A`.
- Open (Member 3): confirm the grid spacing needed to resolve the front (0.1 mm default vs 0.096 mm stated in the work summary).
- Status: open until the items above are filled in.

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