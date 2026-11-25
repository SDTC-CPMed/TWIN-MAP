# Pseudotime signature: guide

From a single-cell dataset to a trajectory signature: build a pseudotime trajectory, test every gene along it, and
select the signature. Follow it by hand, or with an AI coding agent (see the repository README), which takes the
same steps and asks you each decision. [PIPELINE.md](PIPELINE.md) describes the script used in steps 5–7. To rerun the
published analyses instead, see [REPRODUCE.md](REPRODUCE.md).

## Before you start

- R with monocle3, Seurat, harmony, dplyr, Matrix and jsonlite (versions used for the published analyses:
  `environment/R_packages.txt`).
- Your single-cell object, with raw counts and metadata columns for patient or donor, condition or stage, and
  cell type.
- A run folder on your machine (`out_dir`), outside this repository.

## Principles

1. **Every decision is recorded with its reason**: a `# why:` line in `trajectory.R`, or a `why` field in the
   config.
2. **Everything random is seeded and everything hand-picked is saved**: seed Harmony, UMAP, clustering and
   downsampling; save the root cells or nodes to an `.rds` next to the `cds`.
3. **Nothing saved is overwritten**: each run has its own `name` and folder, and a changed decision is a new name.
4. **The signature is fixed before any downstream test**: once accepted (step 7), a change to it is a new name.

## Steps

| # | Step | What you do | Result |
|---|---|---|---|
| 1 | Setup | Check that R loads the packages above. | |
| 2 | Intake | Load your object. Note the cells, genes, whether raw integer counts are present for all genes, the existing reductions, and every metadata column with its levels. Choose the patient, condition and cell-type columns, and the run's `name` and `out_dir`. | |
| 3 | Build the trajectory | Make decisions T1–T11 (below), in order. Write `<out_dir>/<name>/trajectory.R` with each choice as a named setting at the top, followed by its `# why:`, and run it. The published trajectory scripts (REPRODUCE.md) are worked examples. | `cds.rds` meeting PIPELINE.md "Input"; the root cells or nodes (`.rds`); UMAPs by condition, patient and pseudotime |
| 4 | Trajectory checkpoint | Look at the UMAPs. Accept the trajectory, or revisit a T decision under a new `name`. | |
| 5 | graph_test | Copy `config_template.json` to `<out_dir>/<name>/config.json`. Fill `name`, `cds` and `out_dir`, make decision S1, and set `cores` and `batch_size`. Run the `graph_test` stage. | `<out_dir>/graph_test/graph_test_results.rds` |
| 6 | Signature decisions | Make decisions S2–S4, running `diagnose <step>` where listed. | the config's `significance`, `report` and `signature` sections |
| 7 | Signature | Run the `signature` stage and review `report.txt` and `signature.csv`. Accept the signature, or revisit an S decision under a new `name`. | `trajectory_genes.csv`, `signature.csv`, `report.txt` |
| 8 | Record | Write `<out_dir>/<name>/DECISIONS.md`: each decision, the choice, the reason, and anything noted at the checkpoints. | `DECISIONS.md` |

## Decisions

For each decision, look at the evidence, weigh the options and what each one changes, and record your choice and
reason. The options listed are a starting set; any other choice the evidence supports is open.

### Building the trajectory (T1–T11)

Each choice is a named setting at the top of `trajectory.R`, followed by a `# why:` line.

| # | Decision | Evidence | Options and what they change |
|---|---|---|---|
| T1 | **Cells in the trajectory** | cell type × condition table; cells per patient × condition | Which cell types: only the lineage that runs from the control state to the disease state; each extra neighbouring state lengthens the path, another lineage adds a branch. Downsampling: a cap per patient × condition (seeded) stops a few large patients dominating the graph; none keeps every cell. |
| T2 | **Genes for the embedding** | genes detected per cell; the most expressed genes; share of counts from `MT-` and `RP[SL]` genes | All genes; a detection floor (for example ≥ 1% of cells); excluded classes (`MT-`, `RP[SL]`, MALAT1, sex-linked genes), which removes technical and sex axes from the embedding. Affects the embedding only, not which genes are tested (T10). |
| T3 | **Normalisation and variable genes** | mean-variance plot with the variable genes marked | LogNormalize (scale factor 1e4) or SCTransform; 2,000, 3,000 or 5,000 variable genes; selection method (Seurat vst, loess residual of variance). More genes keep weaker programmes in the embedding and add noise. |
| T4 | **PCs** | variance per PC (elbow); which PCs separate condition and which separate patient | PCs computed and PCs used; clipping of scaled values (for example at 10). Too few PCs drop condition-linked variation; too many add noise. |
| T5 | **Batch correction** | UMAP by patient and by sample, before and after; share of each cell's neighbours from the same patient and from the same condition | None; Harmony on patient; Harmony on sample (removes condition when each sample is one patient and one condition); theta and max_iter. Correction should lower same-patient neighbours and keep same-condition neighbours. Always seeded. |
| T6 | **2-D embedding** | UMAPs by condition and patient under each candidate setting | Dimensions in, neighbours, min_dist, metric; always seeded. The graph is learned on this embedding, so one connected body with few islands gives one connected trajectory. |
| T7 | **Clusters and partitions** | UMAP by cluster and by partition; condition mix of each cluster | Monocle `cluster_cells` (resolution) or Seurat clusters copied into the `cds`; Monocle's partitions or one partition. Each partition gets its own graph, and cells in a partition without a root get infinite pseudotime. The clusters also feed the root (T9). Always seeded. |
| T8 | **Graph learning** | graph overlaid on the UMAP: nodes, leaves, branch points, regions the graph does not reach | `use_partition`, `close_loop`, `learn_graph_control` (for example `minimal_branch_len`, `ncenter`). Longer minimal branches give a simpler graph that can miss edge regions. Seeded. |
| T9 | **Root** | UMAP by condition with candidate roots marked; condition mix of the candidate clusters | Cells of the clusters with the highest share of control cells; the graph node nearest the centroid of control cells; nodes picked by hand. Chosen from the control group before pseudotime is seen, and saved (cells or nodes) next to the `cds`. |
| T10 | **Genes in the saved `cds`** | gene count; genes with zero counts; duplicated symbols | All genes with raw counts (zero-count genes fail graph_test and are left out of q); a subset; how duplicated symbols are made unique. If the graph was learned on fewer genes, move it onto a `cds` with full counts. Sets which genes the signature can contain. |
| T11 | **Sensitivity fits** (optional) | the main and each sensitivity fit's UMAP side by side; pseudotime correlation with the main fit | None; another root, embedding, batch correction or cell set. Each fit is saved under its own name and does not replace the main fit. |

### Selecting the signature (S1–S4)

Each choice is a field in `<out_dir>/<name>/config.json`; each section has a `why`. `diagnose <step>` prints the
evidence and changes nothing. What the script computes at each step is in [PIPELINE.md](PIPELINE.md).

| # | Decision | Config field | Evidence | Options and what they change |
|---|---|---|---|---|
| S1 | **graph_test neighbourhood** | `graph_test.neighbor_graph` | the trajectory checkpoint UMAPs | `principal_graph`: Moran's I along the learned trajectory; `knn`: along each cell's nearest neighbours, regardless of the graph. |
| | technical | `graph_test.cores`, `batch_size`, `saved_results` | cores available; genes in the `cds` | `batch_size` null runs one batch; a number checkpoints each batch so a stopped run resumes. `saved_results` null runs graph_test; a path reuses an earlier result. Neither changes the result. |
| S2 | **Significance** | `significance.q_max` | `diagnose significance`: trajectory genes at q < 0.01, 0.05 and 0.1 | Any q cut-off. Sets the trajectory genes, from which the signature is taken. |
| S3 | **Report cut-offs** | `report.morans_I_cutoffs` | the Moran's I range among trajectory genes | Any list of Moran's I values. Reporting only: the count of trajectory genes above each; changes no gene list. |
| S4 | **Signature rule** | `signature.rule` with `n` or `morans_I_min` | `diagnose signature`: signature size under top 10/20/50, Moran's I > 0.2/0.3/0.5, and the elbow | `top_n`: a fixed number of genes; `morans_I_min`: every gene with Moran's I strictly above a value, so the size follows the data; `elbow`: every gene at or above the point of the ranked Moran's I curve furthest from the line joining its ends. |
| | exclusions | `signature.exclude_patterns`, `exclude_genes` | `diagnose signature`: candidates flagged as mitochondrial, ribosomal, MALAT1/NEAT1, sex-linked, histone or stress-response; any contaminating cell type found at the checkpoint | Regular expressions and gene symbols removed before the rule is applied; the next-ranked genes take their places. None keeps every trajectory gene eligible. |
| | platform | `signature.platform_genes` | `diagnose signature`: how many candidates are measured on the platform | A one-column file of the genes measurable downstream (for example an Olink panel), or null. Reported only; does not filter. |
| | expected | `signature.expected` | none | Null for a new dataset. A gene list makes the run stop if the signature differs (used to guard a rerun). |
