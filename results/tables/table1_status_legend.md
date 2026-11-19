# Status legend for Table 1 and Table 1b

The `status_class` column assigns each row one of six classes. The original
`status` wording is kept next to it in each table.

| Class | Meaning |
|---|---|
| measured | Taken or recomputed directly from the recorded S1 data (the wave files). |
| literature | Taken from a published source. The source and range are given in the notes. |
| derived | Computed from other values in the tables by a stated formula or rule. |
| simulated | Produced by the model or a simulation run, not by an experiment. |
| assumed | A modelling choice or illustrative value, not measured or sourced. |
| unknown | Not known. It is left unknown and not filled with a guess. |

## Rules

- Unknown values stay unknown. They are written in the notes and are not replaced by a plausible number.
- Rows still marked "M3 / M4 / M5 to confirm" are open items owned by other members. They are not confirmed.
- A "measured" value from S1 is an animal-level summary (n = 3 rats). Spread across animals is preferred over wave-to-wave spread.
- Baseline values use the 15 pre-injection waves kept by the selection file. Post-injection waves are a separate condition.