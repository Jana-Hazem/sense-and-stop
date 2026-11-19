# S1 data card: rat cortical spreading depolarization (CSD) recordings

Status: DRAFT. Items marked [CONFIRM] were filled from the source paper, the file headers, the extraction code or a look at the QC pictures; the person who has the recording notes should confirm them.

## 1. What this dataset is

- DC-potential and Doppler blood-flow recordings from 3 rats (R4, R5, R6), one electrode per animal.
- Used only to set the timescale of the PDE model (pulse FWHM, rise and decay times). It is not human migraine data.
- Recording dates in the file names: R4 and R5 on 26-03-2024, R6 on 27-03-2024.
- Species/strain, anesthesia, how CSD was triggered, electrode location, sampling rate: male Sprague-Dawley rats (240-280 g); isoflurane anaesthesia (4% induction, 1.5% maintenance in N2O/O2); CSD triggered by 1 M KCl (5 µL) on a cotton bud over the visual cortex, re-applied about every 15 min; one saline-filled glass DC electrode in motor cortex (M1, right hemisphere, 100-150 µm deep) and one laser Doppler probe 250-300 µm rostral to it. Sampling rates from the R4 file header: DC 10,000 Hz, Doppler 1,000 Hz, ECG 1,000 Hz; the R5 and R6 headers have not been inspected [CONFIRM]. Source: Pomierny et al. 2025, *Biomed. Pharmacother.* 189:118320, section 2.12 (the paper linked to RepOD record 10.18150/HIOPYS), plus the R4 header. All three recordings are from the zolpidem group ("Z" in the file names); the baseline waves are the ones before the injection.

## 2. Files

| File | Content |
|---|---|
| `26032024_Z_R4_v4_waves.csv` | 10 detected waves, rat R4 |
| `26032024_Z_R5_S2R_v4_waves.csv` | 12 detected waves, rat R5 |
| `27032024_Z_R6_v4_waves.csv` | 17 detected waves, rat R6 |
| `csd_per_animal_summary_clean.csv` | per-animal mean, SD, min and max (baseline and post-injection groups) |
| `s1_wave_selection.csv` | keep/exclude decision for each of the 26 pre-injection waves |
| `s1_manual_exclusions_reconstructed.csv` | the 14 excluded waves with reasons |

Raw traces are not stored in this repository. Location of raw traces: RepOD, doi:10.18150/HIOPYS (CC BY 4.0), files `26032024_Z_R4.smr`, `26032024_Z_R5.S2R.smr` and `27032024_Z_R6.smr`. A local copy is kept outside Git in `data/raw/S1-repod-csd/`; checksums are in `data/MANIFEST.md` (R4 verified; the R5 and R6 checksum results still have to be recorded there [CONFIRM]).

## 3. Columns of the wave files

`wave_id`, `phase` (pre_injection or post_injection), `dc_peak_time_s`, `dc_amp_V`, `rise_s`, `decay_s`, `fwhm_s`, `rise_10pct_s`, `decay_10pct_s`, `doppler_peak_time_s`, `dc_to_doppler_delay_s`, `doppler_amp_V`, `noise_sigma_V`, `prominence_threshold_V`, `flag`.

- All times are in seconds, amplitudes in volts.
- `fwhm_s`, `rise_10pct_s` and `decay_10pct_s` are the quantities used for the model. `rise_s` and `decay_s` are not used: `rise_s` equals exactly 100 s in many rows and `decay_s` reaches 350 s in two, which looks like a search-window limit. Confirm: confirmed from the extraction code. `rise_s` is the time from the lowest point in the 100 s before the peak to the peak, and `decay_s` the time from the peak to the lowest point in the 350 s after it (the windows of the source paper). A value of exactly 100 s or 350 s means the lowest point sat at the edge of its window, so no clear minimum was found. In the wave files `rise_s` is 100 s in 18 of 39 rows (R4: 2, R5: 5, R6: 11) and `decay_s` is 350 s in 2 rows (both R5).
- Detection used a prominence threshold per animal: 0.0305 V (R4), 0.030 V (R5), 0.0487 V (R6). For R4 and R6 this equals 6 times `noise_sigma_V`. For R5, 6 times sigma would be 0.023 V, so the used value is higher. Rule behind the R5 value: the threshold is max(6 x `noise_sigma_V`, 0.03 V). The 0.03 V is a fixed minimum in the extraction code (`prom_min_v`). For R5, sigma is 0.00386 V, so 6 sigma = 0.023 V falls below the minimum and the 0.03 V minimum was used.
- `flag = near_start_settling` marks waves at the start of a recording.
- What the injection was and when it happened: zolpidem, 1.5 mg/kg free base, intraperitoneal (dose from the source paper; it is not recorded in the files). The time is read from the hand-typed note that starts with "zol" in each file: R4 at 5579.6 s (note "zolp[iden..."), R5 at 6177.1 s ("zolpidemip"), R6 at 5433.3 s ("zolp[idem"). These typed times are approximate (seconds). A wave is labelled pre_injection if its peak is before this time.

## 4. Counts

| Stage | Waves |
|---|---|
| Detected | 39 (26 pre-injection, 13 post-injection) |
| Excluded | 14 (11 pre-injection, 3 post-injection) |
| Kept | 25 |
| of which baseline (pre-injection) used for model timing | 15 (R4: 5, R5: 4, R6: 6) |
| of which post-injection, kept as a separate group | 10 (R4: 4, R5: 2, R6: 4) |

## 5. How waves were selected

- Four waves were removed because of the `near_start_settling` flag (R4 wave 1, R5 wave 1, R6 waves 1 and 2).
- Ten waves were removed by eye (7 pre-injection, 3 post-injection). The reasons were reconstructed afterwards from the values in the wave files; they are not notes made at the time of selection.
- The selection was made after seeing the data.
- The 15 baseline waves are consistent with one rule: no settling flag, and the Doppler peak follows the DC peak by 20 to 120 s. The kept waves have delays of 28 to 97 s; the excluded unflagged waves have delays of -23, 16, 132, 192, 283, 286, 291 s. The 20 to 120 s window was chosen after the fact and describes the selection; it was not the original criterion. The rule does not use FWHM.
- Borderline decisions (R6 wave 5 pre-injection; R5 wave 12 post-injection): no reason was recorded at the time of selection [CONFIRM: edit this if you wrote one]. Reasons found on a later look at the QC pictures. R6 wave 5 (1656.5 s) is a 10 s wide, 0.034 V bump on a drifting baseline, 6 s before the first typed KCl note; the large Doppler response comes about 190 s later (delay 192 s), so the DC and Doppler peaks do not pair. A real wave may have occurred there, but the detected DC peak is probably a KCl-application artifact. R5 wave 12 (9019.5 s) is a broad (FWHM 129 s), low-amplitude (0.036 V) bump with no matching Doppler response (Doppler amplitude 0.032 V, delay 156 s).
- R6 wave 17 (post-injection) is excluded with no numeric anomaly; the reason is a suspected trace artifact that needs checking on the trace: on the QC panel (not the raw trace) the peak is at 9292.5 s, 105 s before the end-of-experiment note ("cullcull" at 9397 s). After the peak the DC falls to about -0.12 V and both channels go flat at the cull, so the wave lies in the end-of-experiment artifact zone. The raw trace itself has not been opened [CONFIRM on the trace].

## 6. Resulting summary values (baseline, n = 3 animals)

| Quantity | Mean of animal means | SD across animals |
|---|---|---|
| FWHM | 42.76 s | 11.32 s |
| Rise time 10% | 35.19 s | 10.43 s |
| Decay time 10% | 40.68 s | 10.98 s |

Pooled over the 15 waves the means are 44.07, 36.00 and 41.93 s. The FWHM mean changes to 53.5 s if only the settling flag is used for exclusion and to 54.9 s if no wave is excluded (see `results/E2_rat_comparison/tables/E2_selection_sensitivity.csv`).

## 7. Limits of use

- Three animals, one electrode each; preliminary.
- The selection rule was set after viewing the data.
- Propagation speed is not measured in this dataset (single electrode); it comes from the literature.
- Not a human or migraine recording.

## 8. Provenance

- Per-animal summary checked against the wave files: all means and SDs match.
- Scripts: `scripts/build_e2_tables_and_dotplot.py`.
