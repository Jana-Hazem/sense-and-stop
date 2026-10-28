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
