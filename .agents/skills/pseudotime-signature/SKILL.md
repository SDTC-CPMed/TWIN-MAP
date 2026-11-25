---
name: pseudotime-signature
description: Guided pseudotime signature run on a new single-cell dataset; builds the trajectory, tests genes along it and selects the signature, with every decision made by the user. Use only when the user explicitly asks for this skill.
---

# Pseudotime signature

You guide the user through [GUIDE.md](<../../../Prong2B Pseudotime/signature_pipeline/GUIDE.md>): its steps, in
order, and every decision in them (T1–T11, S1–S4), with its principles. The script used in steps 5–7 is described
in [PIPELINE.md](<../../../Prong2B Pseudotime/signature_pipeline/PIPELINE.md>). If the user wants to rerun the
published analyses, follow [REPRODUCE.md](<../../../Prong2B Pseudotime/signature_pipeline/REPRODUCE.md>) instead.
Run commands from the repository root.

## Rules

1. **Show, then ask.** Before every question or choice put to the user, display everything the answer depends
   on, in the message itself; the question comes last. For a decision that is the evidence GUIDE.md lists (a plot,
   a table, or `diagnose` output) and the options with what each one changes. For any other question (columns,
   paths, cores, accepting a result), it is the facts the answer rests on.
2. **Every decision is the user's.** Take the decisions one at a time. Give the options GUIDE.md lists, add any
   other option you judge relevant to this dataset and the evidence just shown, and always end with **your own
   choice**, so the user can give one that is not listed. Write their choice and their reason (`why`) into
   `trajectory.R` or the config. When the user asks what you would pick, say so and why, then still let them
   choose.
3. **A run writes only into its own folder.** Choices go into the run's `trajectory.R` and `config.json`;
   everything else the run records (checkpoint notes, findings, the hand-off) goes into new MD files in
   `<out_dir>/<name>/`. The repository's files (this skill, GUIDE.md, PIPELINE.md, REPRODUCE.md, the template,
   the presets and the scripts) stay unchanged during a run.
4. **Long runs go to the background** where your environment allows (the trajectory script on large data,
   graph_test). Report progress from the log; wait for the run to finish before reading its outputs.
5. **Every plot is a PNG, or HTML when it is interactive.** Draw it to a local temp file and copy it into
   `<out_dir>/<name>/plots/` (writing straight into a cloud-synced folder can fail midway).

## Steps

Take GUIDE.md's steps in order. What you add to each, and when it is done:

| # | You also | Done when |
|---|---|---|
| 1 | Check the packages in one `Rscript -e` call. | every package loads |
| 2 | Ask first for the path to the user's object and for `out_dir`, then report the object and ask for the columns and `name`. | the user has confirmed the columns, `name` and `out_dir` |
| 3 | For each T decision, produce the evidence plots and tables and save them in the run folder. Write and run `trajectory.R`. | `cds.rds` loads with raw integer counts for all genes; the root `.rds` and the three UMAP PNGs (`plots/umap_<colour>.png`) exist |
| 4 | Show the three UMAPs. | the user has accepted the trajectory in words |
| 5 | Ask for `cores` and `batch_size`; run `graph_test` in the background. | `graph_test_results.rds` has one row per gene in the `cds` |
| 6 | Run `diagnose <step>` before each decision that lists one. | `significance`, `report` and `signature` are set, `significance` and `signature` each with a `why` |
| 7 | Show `report.txt` and `signature.csv`: genes, Moran's I, flags, platform coverage if set. | the user has accepted the signature in words |
| 8 | Write `DECISIONS.md`; tell the user where the config, `trajectory.R`, outputs and this file are. | every decision in the config and `trajectory.R` appears in `DECISIONS.md` |
