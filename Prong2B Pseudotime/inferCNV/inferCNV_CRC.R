# inferCNV for the CRC normal -> adenoma -> carcinoma trajectory (support for pseudotime).
# Reference: stromal + endothelial cells and up to 200 immune cells per sample.
# Observations: the 9,067 epithelial trajectory cells, plus held-out immune cells from
# adenoma/carcinoma samples as a negative control for sample-level (non-CNV) artefacts.
source("cnv_helpers.R")  # run from this folder

# Data and result folders: set in paths.env (see config.R)
if (!exists("earlyt_path")) local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); d <- if (length(a)) dirname(normalizePath(gsub("~+~", " ", a[1], fixed = TRUE))) else getwd()
  while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "config.R"))) stop("config.R not found: run this script from inside the EARLY-T code folder")
  Sys.setenv(EARLYT_CODE = d); source(file.path(d, "config.R")) })
base_dir <- file.path(earlyt_path("EARLYT_DATA"), "CRC")
out_dir  <- file.path(base_dir, "inferCNV_20261013"); dir.create(out_dir, showWarnings = FALSE)

o   <- readRDS(file.path(base_dir, "seurat_annotated_input.rds"))
cds <- readRDS(file.path(base_dir, "cds_object.rds"))
m <- o@meta.data
m$sample <- paste(m$patient, m$stage, sep = "_")
epi <- colnames(cds)
immune_types <- c("T cell", "B cell", "Macrophage", "Mast cell", "Regulatory T(Treg) cell")

set.seed(1)
imm <- m %>% tibble::rownames_to_column("cell") %>% filter(celltype %in% immune_types)
ctrl <- imm %>% filter(stage != "normal") %>% group_by(sample) %>% slice_sample(n = 200) %>% pull(cell)
ref_imm <- imm %>% filter(!cell %in% ctrl) %>% group_by(sample) %>% slice_sample(n = 200) %>% pull(cell)
ref <- c(rownames(m)[m$celltype %in% c("Stromal cell", "Endothelial cell")], ref_imm)

ann <- data.frame(group = c(paste0("epi_", m[epi, "sample"]), paste0("ctrlimm_", m[ctrl, "sample"]),
                            ifelse(m[ref, "celltype"] %in% immune_types, "ref_immune", "ref_stromal")),
                  row.names = c(epi, ctrl, ref))
obj <- run_infercnv(LayerData(o, assay = "RNA", layer = "counts"), ann, c("ref_immune", "ref_stromal"),
                    file.path(out_dir, "run"))

# ---- Per-cell scores ----
E <- obj@expr.data[, c(epi, ctrl)]
score <- colMeans((E - 1)^2)
md <- data.frame(m[c(epi, ctrl), c("patient", "stage", "sample", "celltype")],
                 compartment = rep(c("epithelial", "immune_control"), c(length(epi), length(ctrl))),
                 cnv_score = score)
md$pseudotime <- NA_real_; md[epi, "pseudotime"] <- pseudotime(cds)[epi]

# Epithelial cells: split-half correlation with their own sample's CNV profile
md[epi, "cnv_cor"] <- split_half_cor(E[, epi], md[epi, "sample"], score)
# Immune controls: correlation with the epithelial CNV profile of the same sample
for (s in unique(md[ctrl, "sample"])) {
  prof <- cnv_profile(E, epi[md[epi, "sample"] == s], score)
  cc <- ctrl[md[ctrl, "sample"] == s]
  md[cc, "cnv_cor"] <- safe_cor(E[, cc] - 1, prof)
}
# Secondary: correlation of every epithelial cell with its patient's carcinoma profile
md$cor_to_carcinoma <- NA_real_
for (p in unique(md$patient)) {
  prof <- cnv_profile(E, epi[md[epi, "patient"] == p & md[epi, "stage"] == "carcinoma"], score)
  cc <- epi[md[epi, "patient"] == p]
  md[cc, "cor_to_carcinoma"] <- safe_cor(E[, cc] - 1, prof)
}

# Null = epithelial cells from normal mucosa; CNV+ above its 99th percentile
null <- md$cnv_cor[md$compartment == "epithelial" & md$stage == "normal"]
thr  <- quantile(null, 0.99)
md$cnv_call <- ifelse(md$cnv_cor > thr, "CNV+", "CNV-")
saveRDS(md, file.path(out_dir, "inferCNV_cell_calls_CRC.rds"))

# ---- Summaries ----
e <- md %>% filter(compartment == "epithelial")
sink(file.path(out_dir, "inferCNV_summary_CRC.txt"))
cat("Threshold (99th pct, normal-mucosa epithelium, split-half correlation):", signif(thr, 3), "\n\n")
cat("CNV+ share by patient x stage (epithelial trajectory cells)\n")
print(as.data.frame(md %>% group_by(compartment, patient, stage) %>%
        summarise(n = n(), cnv_pos = mean(cnv_call == "CNV+"), cor_median = median(cnv_cor), .groups = "drop")), digits = 2)
cat("\nCNV+ share by cell type x stage\n")
print(as.data.frame(e %>% group_by(celltype, stage) %>% summarise(n = n(), cnv_pos = mean(cnv_call == "CNV+"), .groups = "drop")), digits = 2)
cat("\nCorrelation with the patient's carcinoma CNV profile (median) by patient x stage\n")
print(as.data.frame(e %>% group_by(patient, stage) %>% summarise(cor_to_carcinoma = median(cor_to_carcinoma), .groups = "drop")), digits = 2)
cat("\nCNV signal by pseudotime decile (all patients)\n"); print(as.data.frame(by_pseudotime(e)), digits = 2)
cat("\nCNV signal by pseudotime decile, per patient\n"); print(as.data.frame(by_pseudotime(e, "patient")), digits = 2)
cat("\nSpearman(cnv_cor, pseudotime): overall and within patient x stage\n")
cat("overall:", round(cor(e$cnv_cor, e$pseudotime, method = "spearman"), 2), "\n")
print(as.data.frame(e %>% group_by(patient, stage) %>% summarise(n = n(), rho = cor(cnv_cor, pseudotime, method = "spearman"), .groups = "drop")), digits = 2)
sink()

ggsave(file.path(out_dir, "CRC_cnv_cor_vs_pseudotime.png"),
       ggplot(e, aes(pseudotime, cnv_cor, colour = stage)) + geom_point(size = 0.4, alpha = 0.5) +
         geom_hline(yintercept = thr, linetype = 2) + facet_wrap(~patient) + theme_classic() +
         scale_colour_manual(values = c(normal = "#4ABFB5", adenoma = "#E6A817", carcinoma = "#E8735A")) +
         labs(x = "Pseudotime", y = "Split-half CNV correlation (own sample)",
              title = "CRC: CNV signal along pseudotime (dashed: 99th pct of normal mucosa)"),
       width = 11, height = 4, dpi = 200, bg = "white")
