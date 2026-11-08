# scDrugAudit - reproduction materials

Reproduction materials for a benchmark of pipeline choices in single-cell-derived
druggable target prioritisation, with an evidence-grading framework. This repository
contains the connectivity-calibration outputs and the R scripts that produced them.

## Battery #5: LINCS L1000 connectivity calibration

The connectivity battery was run against the real, full LINCS L1000 Level-5
reference: Bioconductor ExperimentHub resource EH3226 (`signatureSearchData` v0.1,
`lincs.h5`), 12,328 Entrez genes x 45,956 signature-level columns, derived from GEO
GSE92742. Three fixed specifications were calibrated for the naive and rigorous
pipelines (aggregation `absmax` as primary, `mean`, and `absmax` with `--min-ct 2`),
each with 2,000 permutations and Benjamini-Hochberg FDR control.

**Reporting boundary.** None of the six calibrations yielded an FDR < 0.05 hit. This
zero-hit outcome is a constructive consequence of the permutation resolution floor
(m = 45,956; p_min = 1/2001 = 4.9975e-4; BH rank-1 threshold = 1.0880e-6), not
evidence either for or against drug connectivity. It must not be reported as absence
of drug connectivity, nor as detection of a connectivity signal.

## Contents

| Path | Content |
| --- | --- |
| `reports/` | Full-library run summary, crop and parallel equivalence checks, tie diagnostics, EH3226 structure and export reports |
| `results/lincs_full_primary/` | Primary specification (`absmax`) calibration tables and summary |
| `results/lincs_full_mean/` | `mean` aggregation sensitivity run |
| `results/lincs_full_minct2/` | `absmax --min-ct 2` sensitivity run |
| `scripts/` | R scripts for connectivity scoring, calibration and the parallel runner |

Each calibration CSV has one row per signature-level column
(`compound, score, p_value, q_value, significant`). `score` is a signed value;
summaries report absolute scores. The `significant` flag is computed from the
unrounded q-value, while `q_value` is rounded to four decimals.

## Data availability

The Level-5 reference matrix and the 3,070-gene query projection are not
redistributed here. They can be regenerated from ExperimentHub EH3226; the
projection is query-specific and must be re-exported for a different query gene
set. The 2.46 GB HDF5 source and the exported matrix remain local analysis inputs.

## Licence

GPL-3.0 (see `LICENSE`).
