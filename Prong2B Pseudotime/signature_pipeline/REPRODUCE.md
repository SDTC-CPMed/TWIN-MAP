# Reproducing the published results

The presets in `presets/` rerun the published prostate, colorectal and lung signatures from their saved
trajectories. A run on a new dataset does not use this file.

The data are not in this repository. Each rerun needs two files from the published analysis, built by the
trajectory scripts in this repository:

| Preset | Trajectory built by | `cds` | graph_test results |
|---|---|---|---|
| `prostate.json` | `Prong2B Pseudotime/prostate/Normal to Prostate Cancer.R` | `cds_object.rds` | `pseudotime_significant_adj_p_val_status_OK.rds` |
| `colorectal.json` | `Prong2B Pseudotime/colorectal/Normal to Carcinoma.R` | `cds_object.rds` | `psuedotime_graph_test_results_crc.rds` |
| `lung.json` | `Prong2B Pseudotime/Lung/05_trajectory.R`, `Prong2B Pseudotime/Lung/14_graph_test.R` | `cds_main.rds` | `pseudotime_graph_test_results_luca.rds` |

Each preset holds the published choices, with the reason in each section's `why`, and the published signature in
`signature.expected`. Its paths are empty.

## Rerunning

1. Copy the preset into your own folder and set `cds`, `graph_test.saved_results` (the graph_test results file)
   and `out_dir` (where the outputs go) to paths on your machine.
2. Run the `signature` stage:

   ```
   Rscript "Prong2B Pseudotime/signature_pipeline/signature_pipeline.R" signature <your copy of the preset>.json
   ```

   It stops if the signature differs from `signature.expected`, and prints "Matches signature.expected" if not.
