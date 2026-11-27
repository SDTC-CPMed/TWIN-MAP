# Dummy UK Biobank data

UK Biobank data cannot be shared, so this folder holds **synthetic** stand-ins for the files `UKBB_model.R` reads.
Every value is randomly generated; nothing comes from UK Biobank. Only the file names, column names and answer codes
follow the UK Biobank layout. Use them to run and test the model code, not to draw any conclusion.

**A weak signal is planted on purpose.** In each organ's incident cancer cases, that organ's signature proteins
shift up or down and its known marker shifts up, more strongly closer to diagnosis. The models therefore give AUCs
above 0.5, falling with years before diagnosis. These numbers say nothing about real biology.

## Files

`input/` is written by `make_dummy_data.R` (seeded, so it rebuilds the same files):

| File | Contents |
|---|---|
| `Olink_proteomics_data_2ndPhase_transposed_decoded2UNIportID.txt` | 5,000 participants × 283 NPX-like protein values, tab-separated, `PID` = eid; about 3% missing |
| `ukb_clinic.rda` | `eid` and the clinical fields the model uses (dates, age, sex, ethnicity, centre, BMI, HbA1c, smoking, alcohol, diabetes) |
| `ukb_<disease>_subset_clinic.rda` | `eid` and `diagnosis_time` for each organ's precursor and cancer (6 files) |
| `ukb_healthy_controls.rda` | `eid` of 1,200 healthy controls |
| `ukb_AllCancers_subset_clinic.rda` | `eid` and `diagnosis_time` of every cancer case |

The `.rda` files are R objects saved with `saveRDS`, as the model reads them with `readRDS`.

Per organ: 260 incident cancers (170 of them after the precursor condition), 20 cancers diagnosed before the blood
draw, and 600 people with the precursor only. Prostate groups are men only. Eids start at 9000001.

The protein panel is not the real Olink panel: it holds every signature and marker gene in `config/` (under the
current symbols ACP3 and H2AZ1 for ACPP and H2AFZ; MALAT1 is non-coding, so it has no assay), 230 made-up proteins
`DUMMY001`–`DUMMY230` for the random panels (two of them over 25% missing, which the model filters out), and two
multi-gene assays.

## Running the model on it

From the repository root, with the R packages in `environment/R_packages.txt`:

```
EARLYT_UKBB_SERVER="$PWD/Prong3A UKBB models/dummy_data/input" UKBB_OUT_DIR=<an output folder> UKBB_NBOOT=10 \
  Rscript "Prong3A UKBB models/UKBB_model.R" colorectal
```

`UKBB_NBOOT` sets the bootstrap iterations (1,000 by default) and `UKBB_OUT_DIR` keeps the results out of this
folder. To rebuild `input/`: `Rscript "Prong3A UKBB models/dummy_data/make_dummy_data.R"`.
