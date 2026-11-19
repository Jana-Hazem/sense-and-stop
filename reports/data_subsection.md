# Data

The study uses three public sources with different roles. They were recorded in
different species, subjects and settings and are not a synchronized dataset.
None of them is a cortical spreading depolarization (CSD) or migraine recording
from the same subject, and no result here is a clinical claim.

| | S1: rat CSD recordings | S2: CHB-MIT scalp EEG | S3: fsaverage head model |
|---|---|---|---|
| Role | Sets the scale (timing) of the simulated wave | Measured background noise | Forward model (head geometry, lead field) |
| Source | RepOD, DOI 10.18150/HIOPYS, version 1.0 | PhysioNet, version 1.0.0, DOI 10.13026/C2K01R | MNE-Python fsaverage template |
| License | CC BY 4.0 | ODC-By 1.0 | [FILL: unconfirmed, see Limitations] |
| What it is | Anesthetized rat, KCl-induced CSD, one DC electrode and one laser Doppler probe | Pediatric epilepsy-monitoring EEG, 10-20 montage | Adult template head, template BEM, ico-5 source space |
| Used here | 3 files: 26032024_Z_R4, 26032024_Z_R5.S2R, 27032024_Z_R6 (of 12 in the dataset) | 3 seizure-free files of one subject: chb01_01, chb01_02, chb01_08 | Patch in lateraloccipital-lh |
| Size | About 1 GB per file [CONFIRM: exact sizes] | 3 files; 1740 samples each after processing | 73 sources in the patch (220 in a 25 mm geodesic disc) |

## S1: rat recordings (scale of the wave)

S1 provides one DC potential channel (10 kHz) and one laser Doppler channel
(1 kHz) per animal, recorded during KCl-induced CSD in anesthetized rats. We used
three files, one per animal (R4, R5, R6). All three belong to the zolpidem group
and were injected during the recording. Only waves before the injection form the
baseline (15 waves: R4 5, R5 4, R6 6). The 10 retained waves after the injection are
kept as a separate group and are not part of the calibration targets.

Waves were detected from the DC channel, then described by FWHM, 10% rise time and
10% decay time. Event criteria, artifact rejection and exclusions are in the
Experiments text and in docs/s1_data_card.md. The result is preliminary: three
animals, one electrode, and no front speed measured in our data. The speed
(2.6 mm/min, range 2.3 to 2.9) is taken from the literature. The calibration
target is the FWHM (42.76 s, mean of three animal means). Rise (35.19 s) and
decay (40.68 s) serve as independent checks.

## S2: CHB-MIT EEG (background noise)

CHB-MIT is pediatric epilepsy-monitoring EEG (256 Hz, 16 bit, most files with 23
channels and one hour of data). It is used only as measured background. Files
listed as containing seizures are excluded. We used three seizure-free files from
one subject (chb01). Each was decimated from 256 Hz to a 2 s step (0.5 Hz) with
anti-aliasing, and about 60 s were trimmed from each end. This leaves 1740 samples
per file on the 2 s grid of the PDE output (3480 s). A decimation factor of 512
leaves only content below 0.25 Hz, so the noise is band-limited and most normal
EEG content is removed. The noise is correlated across channels, and the
right-side channels P4-O2 and P8-O2 are much noisier (10.8 and 12.0 uV) than the
left (about 2.2 to 2.5 uV). The processed noise array has 21 bipolar channels.
Subject age and sex: [FILL: from SUBJECT-INFO]. The filter type and exact trim are
[CONFIRM: M4's script].

## S3: fsaverage head model (forward operator)

The forward operator uses the MNE fsaverage template (adult). A cortical patch in
lateraloccipital-lh (ico-5 source space, minimum distance 5 mm) is flattened by
geodesic distances to a 30 x 20 mm rectangle. This gives 73 sources with a
nearest-neighbour spacing of 1.84 to 3.04 mm (mean 2.45 mm). Source orientation is
fixed to the surface normal, and the reference is the average reference. The patch
location is a design choice, not based on data. The patch is a flat approximation
of a folded cortex (distance distortion on edges 0.92 to 1.27). The BEM,
conductivities and montage file used are [CONFIRM: M4's head-model provenance
note].

The observation model is

    sensors[t, ch] = scale * (G @ u[t].ravel())[ch] + noise_gain * eps[t, ch]

where G = H_b A. Here A (73 x 60501, sparse) averages the PDE grid over each
source, H is the monopolar lead field (19 electrodes x 73 sources, in V per A*m),
H_b forms bipolar channels by subtracting two rows of H, and eps is S2 noise on
the PDE time grid (271 samples, 2 s apart, 0 to 540 s). Three of the 21 noise
channels (T7-FT9, FT9-FT10, FT10-T8) have no lead-field row, because FT9 and FT10
are not among the 19 electrodes. This leaves 18 usable bipolar channels. The E4
inversion uses the monopolar operator H_matrix.npz in this way, and it drops
further noisy channels per window (15 to 18 channels in the main case). The source
scale (1e-9, 1e-8, 1e-7 A*m) is illustrative and not calibrated to physiology.

## Limits of the data

- S1: three animals, one electrode, preliminary. The baseline comes from rats that later received zolpidem, so it is a pre-injection baseline and not a saline control.
- S2: three files from one subject. Noise windows are not independent samples of EEG noise in general.
- S3: adult template, pediatric noise. Source locations are not to be over-interpreted.
- No dataset links a CSD signal to scalp EEG in the same subject. The sensor data are constructed.

## Sources to cite

- S1: RepOD dataset "Cortical spreading depolarisation electrophysiology data", DOI 10.18150/HIOPYS, CC BY 4.0, and the associated paper (Pomierny et al., 2025) [CONFIRM: full reference].
- S2: Guttag J. CHB-MIT Scalp EEG Database, version 1.0.0. PhysioNet; 2010. DOI 10.13026/C2K01R. Also the original publication (Shoeb, 2009, MIT PhD thesis).
- S3: MNE-Python and FreeSurfer [FILL: references and the data license at the download source].