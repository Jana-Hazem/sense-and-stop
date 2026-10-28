# Sense & Stop

**Neuroflux** is developing a computational framework to reconstruct a hidden, CSD-inspired cortical wave from sparse, noisy observations using partial differential equations (PDEs) and inverse methods.

## The problem

Cortical spreading depolarization (CSD) is a slowly propagating wave of neuronal and glial depolarization associated with migraine aura. Direct electrophysiological recordings of CSD during human migraine attacks have not been established; non-invasive human measurements are indirect. This makes it important to distinguish physiological motivation from what a computational model or EEG signal can actually show.

The technical problem studied here is an inverse one: given incomplete, noisy observations, can we recover the hidden spatiotemporal model state \(u(x,t)\)? In simplified form, observations are represented as

\[
y(t) = H[u(t)] + \epsilon(t),
\]

where \(H\) is the observation/forward operator and \(\epsilon\) represents measurement noise. Sensor coverage, sensor placement, noise, and regularization can all affect whether the reconstruction is identifiable and stable.

## Proposed computational approach

1. **Define a wave model.** Use a tractable reaction–diffusion PDE to generate a known spatiotemporal wave. The literature review cautions that scalar Fisher–KPP can serve as a traveling-front baseline, but it has no recovery state and does not represent key CSD physiology. A two-variable activator–recovery model is a candidate for the primary CSD-inspired test system. The team has not frozen this choice; see [`docs/decision_log.md`](docs/decision_log.md).
2. **Create partial observations.** Sample the simulated field at selected sensor locations and add controlled noise. The forward operator and its mapping from the project's 1D PDE grid to source/sensor space remain explicit design decisions.
3. **Reconstruct the hidden state.** First target trajectory-state reconstruction with model coefficients fixed. Treat initial-condition recovery as a separate task; defer parameter estimation until identifiability has been assessed.
4. **Keep later evaluation separate.** The current phase establishes the data, model, units, grid, and inverse-solver contracts. Evaluation, results, and limitations are reserved for a later phase.

The literature review proposes a bounded research gap: among the sources it examined, it did not identify a CSD-focused study jointly assessing how point-sensor number and placement, observation noise, and regularization affect state reconstruction while testing model or discretization mismatch. This is a gap in the reviewed sources, not a claim that no related work exists. The project contribution is framed as applying established inverse-problem methods in a CSD-motivated setting—not as inventing sparse-wave reconstruction or a clinical diagnostic method. Read the [literature review synthesis](reports/literature_review/literature_review.md) for the full rationale and references.

## Data selection and what each source can support

These sources have separate roles. They are not a synchronized dataset and must not be treated as if they were collected from the same people, animals, or experiment.

| Source | Why it was selected | Use in this project | What it cannot establish |
|---|---|---|---|
| [RepOD experimental CSD recordings](https://doi.org/10.18150/HIOPYS) | Rat cortical recordings during KCl-induced CSD, with saline and zolpidem groups | Inform temporal event shape and model scale/parameter bounds. The inspected sample has a named cortical DC channel plus Doppler and ECG channels; use saline recordings for baseline wave properties where verified. | Not human migraine data or spontaneous CSD. The inspected file has only one brain channel, so it cannot establish a spatial propagation profile. Electrode placement and acquisition metadata still need confirmation. |
| [PhysioNet CHB-MIT EEG](https://doi.org/10.13026/C2K01R) | Open pediatric scalp EEG with channel-level recordings and montage metadata | Empirical background/noise only, using selected seizure-free baseline files. Follow `RECORDS` and `RECORDS-WITH-SEIZURES`; do not download the full 42.6 GB database. | Never CSD or migraine labels. Seizure intervals are excluded and must not be used as CSD labels. |
| [MNE-Python fsaverage template](https://mne.tools/stable/auto_tutorials/forward/35_eeg_no_mri.html) | A template head model and lead-field scaffold for the observation operator \(H\) | Explore how model sources could map to sensor measurements. | Adult template, not participant-specific anatomy. It does not resolve the 1D PDE-to-cortical-source mapping; source locations must not be over-interpreted. Template-data licensing also remains to be confirmed. |

An initial inspection of one RepOD `.smr` file found that its header contains channel labels, sampling information, and voltage scaling, but not the electrode location, amplifier gain, filter settings, or the meaning of its keyboard markers. Channel names alone do not prove channel identity, and slow DC deflections are not confirmed CSD events without protocol and lab verification. See [`data/MANIFEST.md`](data/MANIFEST.md), [`docs/data_dictionary.md`](docs/data_dictionary.md), and the open gate in [`docs/decision_log.md`](docs/decision_log.md).

**Raw recordings are never committed.** Manifests and data dictionaries record provenance and verified metadata only.

## Scope and safeguards

- **Primary objective:** reconstruct the hidden model state from sparse, noisy observations.
- **Later extension:** explore a hypothetical localized control term \(g(x,t)\), kept separate from the reconstruction code path.
- This is a computational framework. It does not claim that scalp EEG records human CSD and makes no diagnosis, treatment, or clinical-efficacy claim.
- Experimental animal recordings inform model assumptions; CHB-MIT EEG supplies background/noise only; the template head model is a modeling scaffold.

## Repository layout

`docs/` · `data/` (manifests only) · `configs/` · `src/{data,model,forward,inverse,control}/` · `notebooks/` · `reports/` · `outputs/`

## Team

| Member | Name | Role |
|---|---|---|
| 1 | Jana H | Project and data lead |
| 2 | Suty | Experimental CSD data steward |
| 3 | Ganna | PDE and numerical model lead |
| 4 | Jana M | Forward sensing and EEG lead |
| 5 | Jana W | Inverse problem lead |
| 6 | Habiba | Integration, docs, control extension |

## Setup and contribution

The environment file is maintained by Member 6 and is not yet in place. Follow [`CONTRIBUTING.md`](CONTRIBUTING.md) for branching, pull requests, data-handling, and documentation rules.
