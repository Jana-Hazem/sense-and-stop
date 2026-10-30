# Data Manifest

Raw data are never committed. Local copies live under `data/raw/<source-id>/` (git-ignored). This file records provenance only.

Status values: `listed` = seen on repository page; `downloaded` = local copy exists; `verified` = checksum matched. Page-level details below were reviewed on 2026-10-28. No data had been downloaded at the time this manifest was drafted.

## S1: Experimental CSD recordings (rat)

| Field | Value |
|---|---|
| Source ID | `S1-repod-csd` |
| Title | Cortical spreading depolarisation electrophysiology data |
| Author | Kucharczyk, Mateusz |
| Repository | RepOD |
| Identifier | [DOI 10.18150/HIOPYS](https://doi.org/10.18150/HIOPYS) |
| Version / date | 1.0, published 2026-09-21 |
| License | CC BY 4.0 |
| Content | Rat cortical DC potential and laser Doppler recordings during KCl-induced CSD; saline and zolpidem groups |
| Related paper | Pomierny et al., *Biomedicine & Pharmacotherapy* (2025), DOI [10.1016/j.biopha.2025.118320](https://doi.org/10.1016/j.biopha.2025.118320) |
| Format | Spike2 (`.smr`, `.smrx`); listed files are approximately 200-820 MB each; 12 files total |
| Local path | `data/raw/S1-repod-csd/` (not in Git) |
| Status | Listed; suitability gate pending (owner: Member 2) |
| Role | Informs model scale and parameter bounds only; not human data |
| Open items | Electrode count/layout unknown; prefer saline-group files; test a compatible reader for `.smrx` |

Files visible in the repository listing (MD5 values as listed by RepOD; verify after download):

| File | Group | Size | MD5 | Status |
|---|---|---:|---|---|
| `24072023_S_R1.smrx` | Saline R1 | 756.1 MB | `e0ccc4798dd9aa6fb49ff751fce85f55` | listed |
| `26062023_S_R3.smrx` | Saline R3 | 820.2 MB | `14d2f7a2c5dd674da27c55d835f1ead4` | listed |
| `26072023_S_R4.smrx` | Saline R4 | 710.2 MB | `d8650574be0769d031b740594b85223b` | listed |
| `28062023_S_R5.smrx` | Saline R5 | 706.9 MB | `9fda2962d2a1a4f007b0cd87c4bea833` | listed |
| `05052023_Z_R1.smrx` | Zolpidem R1 | 721.8 MB | `2184aac06827d3499f4e419b535c918c` | listed |
| `12052023_Z_R2.smrx` | Zolpidem R2 | 698.6 MB | `5081e95429405f91b4428e6565802db0` | listed |
| `16052023_Z_R3.smrx` | Zolpidem R3 | 745.6 MB | `2756b693dacd782da27c55d835f1ead4` | listed |
| `26032024_Z_R4.smr` | Zolpidem R4 | 202.7 MB | `2f0b02c9acfa7ec3d120ebc0701078f0` | listed |
| `26032024_Z_R5.S2R.smr` | Zolpidem R5 | 219.8 MB | `327c913cf50f5d7058bf006d357ab41f` | listed |
| `27032024_Z_R6.smr` | Zolpidem R6 | 223.3 MB | `51fc67e606fbf8833273bfed012fbbe7` | listed |
| Two additional files | Not yet reviewed | TBD | TBD | listed status/details pending |

Citation: Kucharczyk, M. (2026). *Cortical spreading depolarisation electrophysiology data* (Version 1.0). RepOD. [https://doi.org/10.18150/HIOPYS](https://doi.org/10.18150/HIOPYS)

## S2: CHB-MIT Scalp EEG Database

| Field | Value |
|---|---|
| Source ID | `S2-chbmit` |
| Repository | PhysioNet |
| Identifier | DOI [10.13026/C2K01R](https://doi.org/10.13026/C2K01R) |
| Version / date | 1.0.0, published 2010-06-09 |
| License | Open Data Commons Attribution License v1.0 (ODC-By) |
| Access | Open, subject to license terms |
| Content | Pediatric scalp EEG, 22 subjects / 23 cases, 10-20 montage, 256 Hz, 16-bit, mostly 23 channels |
| Total size | 42.6 GB (do not download in full) |
| Integrity / selection files | `SHA256SUMS.txt`, `RECORDS`, `RECORDS-WITH-SEIZURES`, `SUBJECT-INFO` |
| Local path | Local copy outside the repository (not in Git) |
| Status | Downloaded (selected files only); checksum vs `SHA256SUMS.txt`: git diff --stat (owner: Member 4) |
| Role | Measured background/noise only. Epilepsy-monitoring data, never CSD or migraine data |
| Selection rule | Download a few seizure-free baseline files plus each case's `chbNN-summary.txt`; exclude every file in `RECORDS-WITH-SEIZURES` |
| Selected files | `chb01_01.edf`, `chb01_02.edf`, `chb01_08.edf` plus `chb01-summary.txt` (case chb01; seizure-free per summary and absent from `RECORDS-WITH-SEIZURES`) |
| Processing | Anti-aliased decimation 256 Hz to 0.5 Hz (three stages of 8), 60 s trimmed from each end; stored outside the repository as `noise_2s.npz` |
| Prefilter | The EDF prefilter field in the `chb01_01` header is empty. This does not show whether hardware filtering was applied, so slow drift may be under-represented |

Citation: Guttag, J. (2010). *CHB-MIT Scalp EEG Database* (Version 1.0.0). PhysioNet. [https://doi.org/10.13026/C2K01R](https://doi.org/10.13026/C2K01R). Also cite the original database publication and Goldberger et al. (2000), as specified by PhysioNet.

## S3: Template head model / lead field

| Field | Value |
|---|---|
| Source ID | `S3-mne-fsaverage` |
| Provider | MNE-Python `fsaverage` template subject, 3-layer BEM, compatible 10-20 sensor positions |
| Version | MNE-Python 1.11.0; template fetched git status |
| License | MNE-Python software: BSD-3-Clause. License for the downloaded template data files is unconfirmed; check before use. |
| Local path | MNE data directory (outside repository) |
| Status | Listed; license check pending (owner: Member 4) |
| Role | Observation operator `H` scaffold. Adult template, not subject-specific anatomy |
| Open items | Adult template versus pediatric EEG subjects (see decision log D4); avoid over-interpreting source locations |

## Checksum procedure

- S2: compare downloaded files against `SHA256SUMS.txt` from the PhysioNet record.
- S1: compare downloaded files against the MD5 values above using `Get-FileHash -Algorithm MD5` (PowerShell), `certutil -hashfile <file> MD5` (Windows), or `md5sum <file>` (Git Bash/WSL).
- Record each local filename, download date, checksum result, and status here. Do not mark a file `verified` until its checksum matches.