# Signature pipeline: script reference

`signature_pipeline.R` takes a pseudotime trajectory (a saved Monocle 3 `cds`) and produces the trajectory genes
and a trajectory signature, the same way for every dataset. This file says what the script needs, what it computes
and what it writes. The steps of a run and the decisions with their options are in [GUIDE.md](GUIDE.md) (S1–S4);
rerunning the published analyses is in [REPRODUCE.md](REPRODUCE.md).

## Principles

1. **The script makes no analysis decisions.** Every decision is a required field in a JSON config file. The code has
   no default values for them; a missing field stops the run with a message naming the field and its options.
2. **Every decision has a reason.** Each config section has a `why` field.
3. **The config is saved with the results.** Every output folder contains the exact config used.

## Input

A saved Monocle 3 `cell_data_set` (`.rds`) with:

- raw integer counts in the `counts` assay, for all genes to be tested (not only variable genes);
- a learned principal graph on the `UMAP` reduction and pseudotime from `order_cells()`;
- gene symbols as row names (unique).

Cells with infinite pseudotime are allowed; they are counted and reported.

## Running

```
Rscript "Prong2B Pseudotime/signature_pipeline/signature_pipeline.R" <graph_test|signature> <config.json>
Rscript "Prong2B Pseudotime/signature_pipeline/signature_pipeline.R" diagnose <significance|signature> <config.json>
```

| Stage | What it does | Cost |
|---|---|---|
| `graph_test` | Stage A: Moran's I for every gene. Saved once; later stages read it. | minutes; the slow stage |
| `diagnose <significance\|signature>` | Prints the evidence for one decision. Changes nothing. | fast |
| `signature` | Stage B: significance filter, report, signature rule | fast |

Stage B reads stage A's results, or an earlier results file named in `graph_test.saved_results` (stage A is then
skipped), so stage B can be tuned without rerunning stage A.

## What each step computes

### Step 1: graph_test (stage A)

`graph_test(cds, neighbor_graph)` on all genes, in batches of `batch_size` genes with a checkpoint file per batch, so
an interrupted run resumes. After all batches, q values are recomputed with Benjamini–Hochberg over all genes with
status OK and a non-missing P value, so batching does not change q. Genes with a missing P value (for example genes
detected in a single cell) get q = 1.

### Step 2: significance filter

Keeps genes with status OK and `q_value < significance.q_max`: the **trajectory genes**, ranked by Moran's I
(largest first). `diagnose significance` prints the count at q < 0.01, 0.05 and 0.1.

### Step 3: report

Prints genes tested, converged (status OK) and significant, and the number of trajectory genes above each
`report.morans_I_cutoffs` value. Changes no gene list.

### Step 4: signature rule

Removes `exclude_genes` and genes matching `exclude_patterns` from the ranked trajectory genes, then applies
`signature.rule`:

- `top_n`: the first `n`;
- `morans_I_min`: Moran's I strictly greater than `morans_I_min`;
- `elbow`: Moran's I at or above the elbow, the point of the ranked curve furthest from the straight line joining
  its ends.

Each signature gene is flagged if it is mitochondrial, ribosomal, MALAT1/NEAT1, sex-linked, histone or
stress-response, and marked if it is in `platform_genes`. If `expected` is set and the signature differs, the run
stops.

`diagnose signature` prints the signature size under top 10/20/50, Moran's I > 0.2/0.3/0.5 and the elbow, and the
top 50 trajectory genes with their flags (and platform coverage, if set).

## Config file

A new run starts from [config_template.json](config_template.json). `cds`, `out_dir` and `saved_results` are paths
on the user's machine. They can be written in full, or as `${NAME}/...` for a folder set in a local `paths.env`
(see [paths.env.example](../../paths.env.example)), which `config.R` expands.

| Section | Fields | Decision |
|---|---|---|
| top level | `name`, `cds`, `out_dir` | run identity and input |
| `graph_test` | `neighbor_graph`, `why`; `cores`, `batch_size`, `saved_results` | S1 |
| `significance` | `q_max`, `why` | S2 |
| `report` | `morans_I_cutoffs` | S3 |
| `signature` | `rule`, `n`, `morans_I_min`, `exclude_patterns`, `exclude_genes`, `platform_genes`, `expected`, `why` | S4 |

## Outputs

```
<out_dir>/
  graph_test/            graph_test_results.rds, batch_NNN.rds checkpoints
  <name>/
    config.json          the config as run
    report.txt           graph_test counts, cut-off table, signature
    trajectory_genes.csv step 2
    signature.csv        step 4: rank, gene, Moran's I, q, flag, measured on platform
    sessionInfo.txt
```
