# E2 (S1 rat side): draft Experiments and Results text

Status: DRAFT for Member 2. Numbers come from `results/E2_rat_comparison/tables/`. Items in **[FILL]** or **[M3]** need information from outside the files checked so far. The model-versus-rat comparison itself belongs to Member 3.

---

## Experiments: E2, rat S1 data

**Question.** What pulse timescale do the rat recordings set for the model, and how well is that timescale known?

**Data.** DC-potential and Doppler blood-flow recordings from three rats (R4, R5, R6), one electrode each, recorded on 26 and 27 March 2024. The animals were male Sprague-Dawley rats under isoflurane anesthesia. CSD was triggered by applying 1 M KCl to the visual cortex at about 15-minute intervals, and the DC potential was recorded with one electrode in right motor cortex, with a laser Doppler probe 250 to 300 µm rostral to it (source: Pomierny et al. 2025; **[CONFIRM]** against the recording notes). All three recordings are from the zolpidem group: zolpidem (1.5 mg/kg, intraperitoneal) was given during the recording, and the waves before it ("pre-injection") are the baseline used here. The data are used only to set the model's timing. They are anesthetized, KCl-evoked rat recordings, not spontaneous human migraine data.

**Event detection.** Waves were detected on the DC trace with a per-animal prominence threshold: 0.0305 V (R4), 0.030 V (R5) and 0.0487 V (R6). For R4 and R6 this is 6 times the noise standard deviation of the trace; the threshold is the larger of 6 times the noise standard deviation and a fixed minimum of 0.03 V, which applies in R5. In total 39 waves were detected: 10 in R4, 12 in R5 and 17 in R6 (26 before the injection and 13 after it).

**Measured quantities.** For each wave: full width at half maximum of the DC pulse (FWHM), 10% rise time and 10% decay time, all in seconds, plus the delay between the DC peak and the Doppler peak. The 10% rise and decay times are used for the comparison; the 100%-based `rise_s` and `decay_s` columns are not, because they sit at the edge of their search windows (100 s before and 350 s after the peak) in 18 and 2 of the 39 waves, respectively, which means no clear minimum was found.

**Wave selection.** Four waves carry the `near_start_settling` flag (the first wave of R4 and R5, the first two of R6) and were removed. Ten more waves were removed by eye, 7 before the injection and 3 after it. The selection was made after the data had been seen, and the written reasons were reconstructed afterwards from the recorded values. The 15 pre-injection waves that remain (R4: 5, R5: 4, R6: 6) form the baseline set used for calibration. Their selection is consistent with a single description: no settling flag and a DC-to-Doppler delay between 20 and 120 s. The kept waves have delays of 28 to 97 s, and the excluded unflagged waves have delays of -23, 16, 132, 192, 283, 286 and 291 s. This window was chosen after the result was known, so it describes the selection and was not the original criterion. It does not use FWHM. Ten post-injection (after zolpidem) waves were kept as a separate group and are not part of the baseline. Two decisions were borderline and no reason was recorded at the time: R6 wave 5, whose DC peak does not pair with the Doppler response about 190 s later and may be a KCl-application artifact, and the post-injection R5 wave 12, a broad low-amplitude bump with no matching Doppler response.

**Aggregation.** Each rat is summarized by the mean of its retained waves. The calibration target is the mean of the three animal means, so each rat counts equally. Spread is reported as the SD of the three animal means (n = 3) and as the range; the within-animal SD is shown for information only. Individual waves are plotted as points (n = 15); no box plots are used, because n is small.

**Selection sensitivity.** The three targets were recomputed under two looser selection rules: excluding only the settling flag (22 waves) and excluding nothing (26 pre-injection waves).

**Scripts and data.** `scripts/build_e2_tables_and_dotplot.py`, inputs in `results/E2_rat_comparison/data/`, jitter seed 20261119 (affects only the horizontal position of the dots).

---

## Results: E2, rat S1 data

**Table 2 (per animal).** `E2_s1_per_animal.csv`. Retained baseline waves; mean (within-animal SD).

| Rat | n | FWHM (s) | Rise 10% (s) | Decay 10% (s) |
|---|---|---|---|---|
| R4 | 5 | 36.20 (2.95) | 25.40 (1.95) | 47.20 (7.29) |
| R5 | 4 | 36.25 (8.46) | 34.00 (15.68) | 28.00 (6.38) |
| R6 | 6 | 55.83 (12.24) | 46.17 (22.31) | 46.83 (17.01) |
| Mean of animal means (n = 3) | | 42.76 | 35.19 | 40.68 |
| SD across animals (n = 3) | | 11.32 | 10.43 | 10.98 |

**Text.** Across the three rats the mean FWHM of the DC wave is 42.76 s (SD across animals 11.32 s, range of animal means 36.20 to 55.83 s), the mean 10% rise time is 35.19 s and the mean 10% decay time is 40.68 s (Table 2, Figure 1a). Two rats have similar FWHM (R4 36.20 s, R5 36.25 s) and R6 is longer (55.83 s), so the mean is set largely by one animal. Within-animal spread is smallest in R4 (SD 2.95 s) and larger in R5 and R6 (8.46 s and 12.24 s). Pooled over the 15 waves, the means are 44.07 s (FWHM), 36.00 s (rise) and 41.93 s (decay), with a wave-to-wave FWHM SD of 13.05 s.

**Sensitivity to selection (Figure 1b, Table 3).**

| Selection rule | Waves | FWHM (s) | Rise 10% (s) | Decay 10% (s) |
|---|---|---|---|---|
| A. Baseline used for calibration | 15 | 42.76 | 35.19 | 40.68 |
| B. All pre-injection, excluding only the settling flag | 22 | 53.54 | 40.46 | 42.85 |
| C. All pre-injection waves | 26 | 54.85 | 37.17 | 45.54 |

The mean FWHM rises from 42.76 s to 53.54 s and 54.85 s when the selection is relaxed, because the removed waves include very long (101, 169, 178 and 201 s) and very short (9 and 10 s) pulses. The mean decay time changes less (40.68 to 45.54 s) and the mean rise time changes by 5.3 s at most.

**Link to the model.** **[M3: add the model-versus-rat comparison table here once regenerated. The model's 10% decay time (28.2 s) is shorter than the rat mean (40.68 s); report it as a mismatch.]**

**Scope.** These results come from three anesthetized rats with one electrode each, recorded under KCl-evoked CSD. The selection rule was set after viewing the data, so the targets are preliminary. Propagation speed was not measured and comes from the literature. The recordings are not human migraine data, and the results make no claim about migraine.

---

## Captions

**Figure 1.** FWHM of the DC wave for each detected wave, grouped by rat (3 rats, one electrode each). (a) The 15 retained baseline waves; each dot is one wave, the black bar is the mean of that rat's retained waves (R4 n = 5, R5 n = 4, R6 n = 6), and the dashed line is the mean of the three animal means (42.76 s). (b) All 26 pre-injection waves; open circles are the 11 excluded waves. Linear axes, seconds.

**Table 2.** Pulse timing of the retained baseline waves per rat, with the spread across rats (n = 3). Within-animal SD is shown for information only.

**Table 3.** Mean of animal means of FWHM, 10% rise and 10% decay time under three wave-selection rules.