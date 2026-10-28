# Sense & Stop

PDE-based reconstruction of a hidden cortical spreading depolarization (CSD)-inspired
wave from sparse, noisy observations. Team: Neuroflux.

## Scope
- **Primary:** reconstruct the hidden spatiotemporal field u(x,t) from y(t) = H u(t) + ε.
- **Secondary (later):** hypothetical localized control term g(x,t), kept separate from the
  reconstruction code path.
- **Out of scope here:** evaluation, results, limitations (reserved for a later phase).

> This is a computational framework. It does not claim that scalp EEG records human CSD,
> and it makes no clinical claim. Data sources play distinct roles (see docs/).

## Data roles
| Source | Role | Boundary |
|---|---|---|
| Experimental CSD recordings (CRCNS / Dryad / Figshare) | Model scale and parameter bounds | Not human migraine data |
| PhysioNet CHB-MIT EEG | Empirical background/noise only | Epilepsy data; never CSD labels |
| Template head model / lead field | Observation operator H | Scaffold, not participant anatomy |

Raw data is never committed. See `data/MANIFEST.md`.

## Repository layout
`docs/` · `data/` (manifests only) · `configs/` · `src/{data,model,forward,inverse,control}/` ·
`notebooks/` · `outputs/`

## Team responsibilities
| Member | Name | Role |
|---|---|---|
| 1 | Jana H | Project and data lead |
| 2 | Suty | Experimental CSD data steward |
| 3 | Ganna | PDE and numerical model lead |
| 4 | Jana M | Forward sensing and EEG lead |
| 5 | Jana W | Inverse problem lead |
| 6 | Habiba | Integration, docs, control extension |

## Setup
Member 6 owns the environment file. Placeholder until it exists.

## Workflow
See `CONTRIBUTING.md`.
