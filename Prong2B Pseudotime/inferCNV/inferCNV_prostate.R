# inferCNV for the prostate normal -> cancer trajectory (support for pseudotime).
# dge_E.rds is epithelial-only, so the reference is normal epithelium: a random half of
# the author-labelled non-malignant cells in the four paired-normal samples (stratified by
# sample x cell type). The other half is held out as the null for calling.
source("cnv_helpers.R")  # run from this folder

# Data and result folders: set in paths.env (see config.R)
if (!exists("earlyt_path")) local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); d <- if (length(a)) dirname(normalizePath(gsub("~+~", " ", a[1], fixed = TRUE))) else getwd()
  while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "config.R"))) stop("config.R not found: run this script from inside the EARLY-T code folder")
  Sys.setenv(EARLYT_CODE = d); source(file.path(d, "config.R")) })
base_dir <- file.path(earlyt_path("EARLYT_DATA"), "Prostate Cancer")
out_dir  <- file.path(base_dir, "inferCNV_20261013"); dir.create(out_dir, showWarnings = FALSE)

o   <- readRDS(file.path(base_dir, "dge_E.rds"))
cds <- readRDS(file.path(base_dir, "cds_object.rds"))
m <- o@meta.data
m$sample <- as.character(m$orig.ident)

set.seed(1)
ref <- m %>% tibble::rownames_to_column("cell") %>%
  filter(sampletype == "Paired_N", malignancy == "Non-Malignant") %>%
  group_by(sample, ID) %>% slice_sample(prop = 0.5) %>% pull(cell)
obs <- setdiff(rownames(m), ref)

ann <- data.frame(group = c(paste0("obs_", m[obs, "sample"]), rep("ref_normal_epi", length(ref))),
                  row.names = c(obs, ref))
obj <- run_infercnv(LayerData(o, layer = "counts"), ann, "ref_normal_epi", file.path(out_dir, "run"))

# ---- Per-cell scores ----
E <- obj@expr.data[, obs]
score <- colMeans((E - 1)^2)
md <- data.frame(m[obs, c("sample", "sampletype", "ID", "malignancy")], cnv_score = score)
md$cnv_cor <- split_half_cor(E, md$sample, score)
md$pseudotime <- pseudotime(cds)[obs]

# Null = held-out author-non-malignant cells from paired-normal samples
is_null <- md$sampletype == "Paired_N" & md$malignancy == "Non-Malignant"
thr <- quantile(md$cnv_cor[is_null], 0.99)
md$cnv_call <- ifelse(md$cnv_cor > thr, "CNV+", "CNV-")

# Positive control: TMPRSS2-ERG interstitial deletion region (chr21 ~38.4-41.5 Mb, GRCh38)
go <- obj@gene_order
del_genes <- rownames(go)[go$chr == "chr21" & go$start > 38.3e6 & go$stop < 41.6e6]
md$chr21_ERG_TMPRSS2 <- colMeans(E[intersect(del_genes, rownames(E)), , drop = FALSE])
saveRDS(md, file.path(out_dir, "inferCNV_cell_calls_prostate.rds"))

# ---- Summaries ----
pt_ok <- md %>% filter(is.finite(pseudotime))
sink(file.path(out_dir, "inferCNV_summary_prostate.txt"))
cat("Threshold (99th pct, held-out non-malignant paired-normal epithelium):", signif(thr, 3), "\n")
cat("Reference cells:", length(ref), " observations:", length(obs), "\n\n")
cat("Agreement with author malignancy labels, by sample type\n")
print(with(md, ftable(sampletype, malignancy, cnv_call)))
cat("\nCNV+ share by sample\n")
print(as.data.frame(md %>% group_by(sampletype, sample) %>%
        summarise(n = n(), author_malignant = mean(malignancy == "Malignant"), cnv_pos = mean(cnv_call == "CNV+"),
                  cor_median = median(cnv_cor), .groups = "drop")), digits = 2)
cat("\nCNV+ share by author cell type (ID)\n")
print(as.data.frame(md %>% group_by(ID) %>% summarise(n = n(), cnv_pos = mean(cnv_call == "CNV+"), .groups = "drop")), digits = 2)
cat("\nPositive control: mean inferCNV signal over chr21 ERG-TMPRSS2 region (", length(del_genes), "genes; <1 = loss)\n")
print(as.data.frame(md %>% group_by(ID) %>% summarise(n = n(), median = median(chr21_ERG_TMPRSS2),
                                                      share_below_0.97 = mean(chr21_ERG_TMPRSS2 < 0.97), .groups = "drop")), digits = 3)
cat("\nCNV signal by pseudotime decile (finite pseudotime)\n"); print(as.data.frame(by_pseudotime(pt_ok)), digits = 2)
cat("\nCNV signal by pseudotime decile, per sample type\n"); print(as.data.frame(by_pseudotime(pt_ok, "sampletype")), digits = 2)
cat("\nSpearman(cnv_cor, pseudotime): overall", round(cor(pt_ok$cnv_cor, pt_ok$pseudotime, method = "spearman"), 2), "\n")
print(as.data.frame(pt_ok %>% group_by(sample) %>% summarise(n = n(), rho = cor(cnv_cor, pseudotime, method = "spearman"), .groups = "drop")), digits = 2)
sink()

ggsave(file.path(out_dir, "prostate_cnv_cor_vs_pseudotime.png"),
       ggplot(pt_ok, aes(pseudotime, cnv_cor, colour = malignancy)) + geom_point(size = 0.3, alpha = 0.4) +
         geom_hline(yintercept = thr, linetype = 2) + facet_wrap(~sample, ncol = 5) + theme_classic(base_size = 9) +
         scale_colour_manual(values = c("Malignant" = "#E8735A", "Non-Malignant" = "#4ABFB5"), name = "Author label") +
         labs(x = "Pseudotime", y = "Split-half CNV correlation (own sample)",
              title = "Prostate: CNV signal along pseudotime (dashed: 99th pct of held-out normal epithelium)"),
       width = 12, height = 7, dpi = 200, bg = "white")
