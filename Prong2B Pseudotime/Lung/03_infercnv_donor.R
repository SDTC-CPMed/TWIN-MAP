# Stage 2 of the LuCA LUAD redo (LOCKED_PLAN.md): inferCNV for one donor, same method as
# inferCNV/inferCNV_LUAD.R. Run from this folder, one process per donor: Rscript 03_infercnv_donor.R <donor_id>
# Observations: AT2, transitional, club and tumour cells by sample, plus held-out immune cells from
# tumour samples as a negative control. Reference: the remaining immune cells and the stromal cells.
# Calls are made across donors in 04_cnv_calls.R.
source("../inferCNV/cnv_helpers.R")
options(scipen = 100)

donor <- commandArgs(TRUE)[1]
# Data and result folders: set in paths.env (see config.R)
if (!exists("earlyt_path")) local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); d <- if (length(a)) dirname(normalizePath(gsub("~+~", " ", a[1], fixed = TRUE))) else getwd()
  while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "config.R"))) stop("config.R not found: run this script from inside the EARLY-T code folder")
  Sys.setenv(EARLYT_CODE = d); source(file.path(d, "config.R")) })
out_dir <- file.path(earlyt_path("EARLYT_DATA"), "LuCa", "luad_10x_20261103")
cnv_dir <- file.path(out_dir, "cnv"); dir.create(cnv_dir, showWarnings = FALSE)
out_file <- file.path(cnv_dir, paste0(donor, ".rds"))
if (file.exists(out_file)) quit(save = "no")
t0 <- Sys.time()

cells <- readRDS(file.path(out_dir, "cells.rds"))$cells %>% filter(donor_id == !!donor)
genes <- readRDS(file.path(out_dir, "genes.rds"))
counts <- readRDS(file.path(out_dir, paste0("counts_", donor, ".rds")))
rownames(counts) <- genes$gene[match(rownames(counts), genes$gene_id)]

# Up to 100 immune cells per tumour sample are held out, and never over a third of the donor's immune cells
set.seed(14)
imm <- cells %>% filter(label == "Ref_immune")
ctrl <- imm %>% filter(origin == "tumor_primary") %>% group_by(sample) %>% group_modify(~ slice_sample(.x, n = min(nrow(.x), 100))) %>% ungroup() %>% pull(cell)
if (length(ctrl) > nrow(imm) / 3) ctrl <- sample(ctrl, floor(nrow(imm) / 3))

# Too few stromal cells to form their own reference group: merge them with the immune reference
stromal <- cells$cell[cells$label == "Ref_stromal"]
ref_ann <- c(setNames(rep("ref_immune", nrow(imm)), imm$cell), setNames(rep(if (length(stromal) >= 10) "ref_stromal" else "ref_immune", length(stromal)), stromal))
ref_ann <- ref_ann[!names(ref_ann) %in% ctrl]

epi <- cells %>% filter(!grepl("^Ref", label))
ann <- data.frame(group = c(paste0("epi_", epi$sample), paste0("ctrlimm_", cells$sample[match(ctrl, cells$cell)]), ref_ann),
                  row.names = c(epi$cell, ctrl, names(ref_ann)))
go <- gene_order_for(rownames(counts))
run_dir <- file.path(cnv_dir, paste0("run_", donor))
obj <- CreateInfercnvObject(raw_counts_matrix = counts[rownames(go), rownames(ann)], annotations_file = ann,
                            gene_order_file = go, ref_group_names = unique(ref_ann))
obj <- infercnv::run(obj, cutoff = 0.1, out_dir = run_dir, cluster_by_groups = TRUE, denoise = TRUE, HMM = FALSE,
                     num_threads = 2, no_plot = TRUE, no_prelim_plot = TRUE, save_rds = FALSE, write_expr_matrix = FALSE)

# ---- Per-cell scores ----
E <- obj@expr.data[, c(epi$cell, ctrl)]
score <- colMeans((E - 1)^2)
md <- cells %>% filter(cell %in% colnames(E)) %>% select(cell, donor_id, dataset, sample, origin, label)
rownames(md) <- md$cell; md <- md[colnames(E), ]
md$compartment <- ifelse(md$cell %in% ctrl, "immune_control", "epithelial")
md$cnv_score <- score

# Epithelial cells: split-half correlation with their own sample's CNV profile
e_cells <- md$cell[md$compartment == "epithelial"]
md[e_cells, "cnv_cor"] <- split_half_cor(E[, e_cells], md[e_cells, "sample"], score)
# Immune controls: correlation with the epithelial CNV profile of the same sample
for (s in unique(md[ctrl, "sample"])) {
  cc <- ctrl[md[ctrl, "sample"] == s]
  md[cc, "cnv_cor"] <- safe_cor(E[, cc, drop = FALSE] - 1, cnv_profile(E, e_cells[md[e_cells, "sample"] == s], score))
}
# Non-tumour-labelled cells: correlation with the profile of this donor's labelled tumour cells,
# which do not contribute to it. Tells whether bridge cells carry the same CNVs as the tumour.
tum <- md$cell[md$label == "Tumour" & md$origin == "tumor_primary"]
others <- setdiff(e_cells, tum)
md$cor_to_tumour <- NA_real_
md[others, "cor_to_tumour"] <- safe_cor(E[, others, drop = FALSE] - 1, cnv_profile(E, tum, score))

saveRDS(list(calls = md, n_genes = nrow(E), n_ref = length(ref_ann), n_ctrl = length(ctrl), mins = as.numeric(difftime(Sys.time(), t0, units = "mins"))), out_file)
unlink(run_dir, recursive = TRUE)
