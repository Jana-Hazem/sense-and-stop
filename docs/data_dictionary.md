# Data Dictionary

This table records each data field, its units, dimensions, and provenance. `TBD` means not yet verified and names the owner responsible for checking it. Replace TBD only using the file header, paper, or repository page, and record the source.

## 1. Conventions (project-wide)

| Item | Convention | Status |
|---|---|---|
| Time unit | Seconds (`s`) in derived arrays | Proposed; Members 3, 4, and 5 confirm |
| Space unit | Millimetres (`mm`) for the PDE domain | Proposed; Member 3 confirms |
| State `u(x,t)` | Dimensionless normalized excitation state | Per plan |
| Array layout | `(time, space)` for fields; `(time, channel)` for sensor data | Proposed; confirm at interface freeze |
| Missing values | `NaN`; never silently dropped; exclusions logged with reason | Per plan |
| Provenance | Every saved array carries source ID, version, and processing configuration | Per plan |

## 2. Source S1: RepOD rat CSD recordings (Owner: Member 2)

| Field | Meaning | Unit | Dimensions | Sampling rate | Provenance | Status |
|---|---|---|---|---|---|---|
| DC potential channels | Cortical DC potential during CSD | ADC volts ($\pm 5\text{ V}$, $0.0001526\text{ V/bit}$ before gain; normalize per animal) | 1 site (M1, right side) | $10{,}000\text{ Hz}$ (`CSD-DC`) | S1 file header and Pomierny et al. (2025) | Verified for `26032024_Z_R4.smr` |
| Laser Doppler channel | Capillary blood-flow measurement | ADC volts (flux in arbitrary units per paper) | 1 probe ($250\text{--}300\ \mu\text{m}$ rostral to DC) | $1{,}000\text{ Hz}$ (`doppler`) | S1 file header and Pomierny et al. (2025) | Verified for `26032024_Z_R4.smr` |
| Electrode layout | Recording positions on cortex | AP coordinates in mm (visual to M1 $\sim 6.5\text{ mm}$) | 1D spatial points | n/a | Paper methods | Partially verified (AP distance clear, DV/ML garbled in paper) |
| Group label | Saline or zolpidem | Categorical | Per file ($n=6$ per group, 12 files total) | n/a | File name / RepOD description | Known from page (2 unreviewed files mapped to saline R2 & R6) |
| KCl application schedule / event times | CSD induction schedule and recorded event times | Schedule reported as every 15 min; actual times typed at $1174, 2089, 3010, 3913, 4780, 5698, 6561, 7362\text{ s}$ | Event timestamps | n/a | Repository description / paper; confirm actual file timing (S1 Keyboard channel) | Verified for `26032024_Z_R4.smr` (typed ASCII key codes; markers don't align precisely with wave arrivals) |
| Missing data / artifacts | Gaps, drift, clipping | n/a | Time windows | n/a | Member 2 exclusion log | Verified (Excluded: $0\text{--}2\text{ min}$ insertion saturation, $0\text{--}18\text{ min}$ drift/amp offset shift, $79.7\text{ min}$ KCl spike, $>145.6\text{ min}$ post-mortem saturation) |

Known from the repository description: Spike2 format (`.smr` / `.smrx`), rat recordings, CSD induced with 1 M KCl (5 µL). Do not treat the stated application schedule as verified event timestamps in a particular recording.

## 3. Source S2: CHB-MIT EEG (Owner: Member 4)

| Field | Meaning | Unit | Dimensions | Sampling rate | Provenance | Status |
|---|---|---|---|---|---|---|
| EEG channels | Scalp potential, 10-20 montage | TBD; read from each EDF header | `(samples, approximately 23 channels)` | 256 Hz | PhysioNet record and EDF header | Rate / typical count known; units per file TBD |
| Resolution | ADC bit depth | bits | n/a | n/a | PhysioNet record | 16-bit |
| Dummy channels | Placeholder channels named `-` | n/a | Up to 5 per file | n/a | EDF header / PhysioNet description | Drop and record each drop |
| Channel names / montage | Electrode names and pairs | n/a | Per case | n/a | `chbNN-summary.txt` | TBD per selected case |
| Segment length | Recording duration | Seconds (`s`) | Per file | n/a | EDF header / PhysioNet record | Most files 1 h; some 2 h or 4 h |
| Seizure annotations | Start/end of seizures | Seconds from file start (`s`) | Per file | n/a | `RECORDS-WITH-SEIZURES` and summary files | Use only to exclude seizure files/intervals |
| Subject age / sex | Case metadata | Years / categorical | Per case | n/a | `SUBJECT-INFO` | Listed in source; record only for selected cases |

## 4. Source S3: MNE `fsaverage` template (Owner: Member 4)

| Field | Meaning | Unit | Dimensions | Coordinate frame | Provenance | Status |
|---|---|---|---|---|---|---|
| Source space | Cortical source locations | Metres (`m`); confirm on installed files | 20,484 vertices (ico-5, both hemispheres) in the cited tutorial output | MRI surface RAS, transformed to head frame for forward computation | MNE tutorial output | Count seen in tutorial; confirm on installation |
| BEM | Three-layer head conductor model | n/a | 5,120 triangles per layer in cited tutorial files | MRI / head; record transform | MNE `fsaverage` files | Confirm file and model used |
| Sensor positions | 10-20 / `standard_1005` montage | Metres (`m`); confirm on installed montage | Per channel | TBD; record MNE montage frame and transform to the forward model frame | MNE montage and project configuration | Confirm channel names and frame match CHB-MIT after any documented renaming |
| Lead field `H` | Maps source components to sensor space | TBD | TBD (typically sensor channels by source components) | Head | Computed by Member 4 | TBD |

The MNE tutorial warns against over-interpreting source locations when an individual MRI is unavailable. Do not assume a coordinate frame from the montage name; record the frame and transform used by the actual installed objects.

## 5. Derived fields (interfaces between members)

| Field | Meaning | Unit | Dimensions | Producer | Consumer | Status |
|---|---|---|---|---|---|---|
| `u(x,t)` | Hidden wave state | Dimensionless | `(n_t, n_x)` | Member 3 | Members 4 and 5 | TBD at grid freeze |
| Grid `x`, `t` | Named spatial and time grids | `mm`, `s` (proposed) | `n_x`, `n_t` | Member 3 | All | TBD |
| `H` | Forward operator in `y = H u + epsilon` | Per-entry units TBD | TBD; agree source-component convention | Member 4 | Member 5 | TBD |
| `y(t)` | Constructed sensor observations | TBD; read EEG units from EDF header and check consistency with `H u` | `(n_t, n_channels)` | Member 4 | Member 5 | TBD |
| `epsilon(t)` | Background/noise from CHB-MIT segments | Same as `y` | `(n_t, n_channels)` | Member 4 | Member 5 | TBD |
| Reconstruction | Estimated hidden state `u` | Dimensionless | `(n_t, n_x)` | Member 5 | Downstream analysis | TBD |

## Open interface question: 1D PDE grid to cortical sources

The PDE is one-dimensional in space, while the head model's sources lie on a cortical surface. Members 3 and 4 must agree and document how the 1D grid maps to source positions before the interface is frozen. Record the mapping rule, units, coordinate transforms, source orientation convention, and resulting array shapes. Until then, do not assume that `H` directly accepts the PDE grid.
