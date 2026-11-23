library(SingleCellExperiment)
library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)
library(SeuratWrappers)
library(harmony)     # used: patient batch correction (RunHarmony)

# Data and result folders: set in paths.env (see config.R)
if (!exists("earlyt_path")) local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); d <- if (length(a)) dirname(normalizePath(gsub("~+~", " ", a[1], fixed = TRUE))) else getwd()
  while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "config.R"))) stop("config.R not found: run this script from inside the EARLY-T code folder")
  Sys.setenv(EARLYT_CODE = d); source(file.path(d, "config.R")) })
base_dir <- file.path(earlyt_path("EARLYT_DATA"), "CRC")

# Load the Seurat object
seurat_obj <- readRDS(file.path(base_dir, "seurat_annotated_input.rds"))

epi_cells <- subset(seurat_obj,
                    celltype %in% c("Epithelial cell", "Goblet cell", "Transit-amplifying cell"))

# --- Re-normalise and find variable features on the subset ---
DefaultAssay(epi_cells) <- "RNA"
epi_cells <- NormalizeData(epi_cells)
epi_cells <- FindVariableFeatures(epi_cells, nfeatures = 3000)
epi_cells <- ScaleData(epi_cells)

epi_cells <- RunPCA(epi_cells, npcs = 50)

# --- Harmony batch correction on patient ---
# NOTE: the run of record had no set.seed() before RunHarmony, so a full rerun gives a
# slightly different embedding and graph (15 of 20 signature genes recovered; Supp. Methods S5).
# Left unseeded so the saved objects stay the record; seed it for any new analysis.
epi_cells <- RunHarmony(epi_cells,
                        group.by.vars = "patient",
                        reduction = "pca",
                        dims.use = 1:20,
                        reduction.save = "harmony",
                        theta = 3,
                        max_iter = 20)

# --- Re-run UMAP on Harmony embedding ---
epi_cells <- RunUMAP(epi_cells,
                     reduction = "harmony",
                     dims = 1:20,
                     reduction.name = "umap_harmony")

# --- Clustering on Harmony embedding ---
epi_cells <- FindNeighbors(epi_cells, reduction = "harmony", dims = 1:20)
epi_cells <- FindClusters(epi_cells, resolution = 0.3)

epi_cells$is_root <- epi_cells$seurat_clusters %in% c("7", "8")
table(epi_cells$is_root)  # should be 678 TRUE (the root cells, see order_cells below)

# Set the Harmony UMAP as the default UMAP
# SeuratWrappers looks for a reduction named "umap" specifically
epi_cells[["umap"]] <- epi_cells[["umap_harmony"]]

# Now convert
cds <- as.cell_data_set(epi_cells)

reducedDim(cds, "PCA") <- epi_cells@reductions$harmony@cell.embeddings


recreate_partition <- function(cds, seurat_obj) {

  cluster_assignments <- setNames(
    factor(seurat_obj$seurat_clusters),
    colnames(seurat_obj)
  )

  partition_assignments <- setNames(
    rep(1, ncol(seurat_obj)),
    colnames(seurat_obj)
  )
  partition_assignments <- factor(partition_assignments)

  cds@clusters$UMAP$clusters <- cluster_assignments
  cds@clusters$UMAP$partitions <- partition_assignments

  return(cds)
}

cds <- recreate_partition(cds, epi_cells)

# --- Learn the trajectory graph ---
# use_partition = FALSE keeps all cells in one trajectory
# close_loop = FALSE appropriate for a linear progression
cds <- learn_graph(cds,
                   use_partition = FALSE,
                   close_loop = FALSE)

# Root: the 678 cells of Seurat clusters 7 and 8, the two clusters with the highest
# proportion of histologically normal cells (97.7% and 94.0%). Cluster numbers can
# change if Harmony/UMAP/Louvain are rerun, so the barcodes were saved to
# root_cells.rds and are read back here; this reproduces the saved pseudotime exactly.
# root_cells <- colnames(epi_cells)[epi_cells$seurat_clusters %in% c("7", "8")]   # original definition
root_cells <- readRDS(file.path(base_dir, "root_cells.rds"))
stopifnot(all(root_cells %in% colnames(cds)))
cds <- order_cells(cds, root_cells = root_cells)

# Figure versions without the root/branch labels
dir.create(file.path(base_dir, "figures"), showWarnings = FALSE)
p1_fig <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,
                     label_roots = FALSE, label_branch_points = FALSE) + ggtitle("Pseudotime")
p2_fig <- plot_cells(cds, color_cells_by = "stage", label_cell_groups = FALSE, label_leaves = FALSE,
                     label_roots = FALSE, label_branch_points = FALSE) + ggtitle("Stage")
ggsave(file.path(base_dir, "figures", "CRC_UMAP_pseudotime.png"), p1_fig, width = 6, height = 5, dpi = 300)
ggsave(file.path(base_dir, "figures", "CRC_UMAP_stage.png"), p2_fig, width = 6, height = 5, dpi = 300)

# Extract pseudotime and metadata
pt_df <- data.frame(
  pseudotime = pseudotime(cds),
  stage = epi_cells$stage,
  patient = epi_cells$patient,
  celltype = epi_cells$celltype,
  cluster = epi_cells$seurat_clusters
)

# Remove Inf pseudotime (disconnected cells)
pt_df <- pt_df %>% filter(is.finite(pseudotime))

# checking if normal cell expression being large ranged is a result of high pseudotime values
# on island clusters

# How many normal cells have unexpectedly high pseudotime?
pt_df %>%
  filter(stage == "normal") %>%
  summarise(
    total = n(),
    high_pt = sum(pseudotime > 30),
    pct_high = round(100 * high_pt / total, 1)
  )

# Which clusters are driving the high-pseudotime normal cells?
pt_df %>%
  filter(stage == "normal", pseudotime > 30) %>%
  dplyr::count(cluster) %>%
  arrange(desc(n))

# Is it patient-specific?
pt_df %>%
  filter(stage == "normal", pseudotime > 30) %>%
  dplyr::count(patient) %>%
  arrange(desc(n))

## DECISION : Normal cells are genuinely mixed in with adenoma and carcinoma cells,
## and not majorly contributed by one patient either. So we keep them as it is,
## do not exclude any clusters

library(ggpubr)
#### Statistical test to confirm
# Kruskal-Wallis + pairwise Wilcoxon to formally confirm ordering
ggplot(pt_df, aes(x = factor(stage, levels = c("normal", "adenoma", "carcinoma")),
                  y = pseudotime, fill = stage)) +
  geom_boxplot(outlier.size = 0.5) +
  stat_compare_means(comparisons = list(c("normal", "adenoma"),
                                        c("adenoma", "carcinoma"),
                                        c("normal", "carcinoma")),
                     method = "wilcox.test") +
  stat_compare_means(method = "kruskal.test", label.y = 60) +
  scale_x_discrete(limits = c("normal", "adenoma", "carcinoma")) +
  theme_classic() +
  ggtitle("Pseudotime distribution by stage")

# Exact P values for the Methods (the ggpubr labels print "< 2.2e-16"). Pairwise
# tests are BH-adjusted; the figure labels above are unadjusted, but every pair is
# below 1e-170 either way.
pt_df %>% group_by(stage) %>% summarise(n = n(), median_pt = median(pseudotime))
kruskal.test(pseudotime ~ stage, data = pt_df)$p.value   # 0 = below the smallest double
pairwise.wilcox.test(pt_df$pseudotime, pt_df$stage, p.adjust.method = "BH")$p.value

# --- Run graph_test ---
# cores = 4, adjust to your machine
pr_test_res <- graph_test(cds,
                          neighbor_graph = "principal_graph",
                          cores = 4)

pr_test_res <- pr_test_res %>% dplyr::mutate(gene_short_name = rownames(.))

saveRDS(pr_test_res, file = file.path(base_dir, "psuedotime_graph_test_results_crc.rds"))

# graph_test summary (numbers reported in the Methods): genes tested, converged,
# significant, and significant above each Moran's I cut-off (converged tests only)
r_ok <- pr_test_res %>% filter(status == "OK")
cat("graph_test -- tested:", nrow(pr_test_res),
    "| converged:", nrow(r_ok),
    "| q < 0.05:", sum(r_ok$q_value < 0.05, na.rm = TRUE), "\n")
data.frame(
  morans_I_cutoff = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30),
  n_genes = sapply(c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30), function(t) {
    sum(r_ok$q_value < 0.05 & r_ok$morans_I > t, na.rm = TRUE)
  })
)

# --- Filter significant genes ---
sig_genes <- pr_test_res %>%
  dplyr::filter(q_value < 0.05, status == 'OK') %>%
  dplyr::arrange(desc(morans_I))

nrow(sig_genes)

saveRDS(sig_genes, file = file.path(base_dir, "significant_genes_adjusted_pval_005_status_OK.rds"))


# --- Add pseudotime to metadata ---
epi_cells$pseudotime <- pseudotime(cds)

saveRDS(cds, file = file.path(base_dir, "cds_object.rds"))

# Saved to a subfolder: save_monocle_objects writes its own cds_object.rds, which in
# base_dir would overwrite the saveRDS() copy above that the Methods numbers are read from.
monocle3::save_monocle_objects(cds = cds, directory_path = file.path(base_dir, "monocle_cds"), hdf5_assays = FALSE)

saveRDS(epi_cells, file = file.path(base_dir, "epi_cells.rds"))

# --- DE between stages using the epithelial subset ---
# Normal vs Adenoma
Idents(epi_cells) <- "stage"
de_normal_vs_adenoma <- FindMarkers(epi_cells,
                                    ident.1 = "adenoma",
                                    ident.2 = "normal",
                                    test.use = "wilcox",
                                    min.pct = 0.1,
                                    logfc.threshold = 0.25)

# Adenoma vs Carcinoma
de_adenoma_vs_carcinoma <- FindMarkers(epi_cells,
                                       ident.1 = "carcinoma",
                                       ident.2 = "adenoma",
                                       test.use = "wilcox",
                                       min.pct = 0.1,
                                       logfc.threshold = 0.25)

# --- Find genes consistently upregulated across both transitions ---
# These are monotonically increasing normal -> adenoma -> carcinoma
sig_up_both <- intersect(
  rownames(de_normal_vs_adenoma)[de_normal_vs_adenoma$avg_log2FC > 0.5 &
                                   de_normal_vs_adenoma$p_val_adj < 0.05],
  rownames(de_adenoma_vs_carcinoma)[de_adenoma_vs_carcinoma$avg_log2FC > 0.5 &
                                      de_adenoma_vs_carcinoma$p_val_adj < 0.05]
)

# Consistently downregulated
sig_down_both <- intersect(
  rownames(de_normal_vs_adenoma)[de_normal_vs_adenoma$avg_log2FC < -0.5 &
                                   de_normal_vs_adenoma$p_val_adj < 0.05],
  rownames(de_adenoma_vs_carcinoma)[de_adenoma_vs_carcinoma$avg_log2FC < -0.5 &
                                      de_adenoma_vs_carcinoma$p_val_adj < 0.05]
)

length(sig_up_both)    # genes consistently rising across progression
length(sig_down_both)  # genes consistently falling

# Genes that are both pseudotime-correlated AND stage-discriminating
panel_up <- intersect(sig_up_both, sig_genes$gene_short_name)
panel_down <- intersect(sig_down_both, sig_genes$gene_short_name)

length(panel_up)
length(panel_down)

# Rank by Moran's I within the filtered set
panel_up_ranked <- sig_genes %>%
  dplyr::filter(gene_short_name %in% panel_up) %>%
  dplyr::arrange(desc(morans_I))

panel_down_ranked <- sig_genes %>%
  dplyr::filter(gene_short_name %in% panel_down) %>%
  dplyr::arrange(desc(morans_I))

head(panel_up_ranked, 20)
head(panel_down_ranked, 20)

genes_to_plot <- unique(c(panel_up_ranked$gene_short_name,panel_down_ranked$gene_short_name))

# Extract expression matrix for these genes
expr_df <- FetchData(epi_cells,
                     vars = c(genes_to_plot, "stage")) %>%
  tibble::rownames_to_column("cell") %>%
  mutate(pseudotime = pseudotime(cds)[cell]) %>%
  filter(is.finite(pseudotime)) %>%
  tidyr::pivot_longer(cols = all_of(genes_to_plot),
                      names_to = "gene",
                      values_to = "expression")

#expression curves — one panel per gene
ggplot(expr_df, aes(x = pseudotime, y = expression, color = stage)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("normal"       = "#E6A817",
                                "adenoma"   = "#E91E8C",
                                "carcinoma" = "#9B59B6"
  )) +
  facet_wrap(~gene, scales = "free_y", ncol = 3) +
  theme_classic() +
  ggtitle("Smoothed gene expression along pseudotime") +
  xlab("Pseudotime") + ylab("Expression")


# Plot upregulated genes only
p_up <- expr_df %>%
  filter(gene %in% panel_up_ranked$gene_short_name) %>%
  ggplot(aes(x = pseudotime, y = expression, color = stage)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("normal" = "#E6A817",
                                "adenoma" = "#E91E8C",
                                "carcinoma" = "#9B59B6")) +
  facet_wrap(~gene, scales = "free_y", ncol = 4) +
  theme_classic() +
  ggtitle("Upregulated genes along pseudotime") +
  xlab("Pseudotime") + ylab("Expression")

# Plot downregulated genes only
p_down <- expr_df %>%
  filter(gene %in% panel_down_ranked$gene_short_name) %>%
  ggplot(aes(x = pseudotime, y = expression, color = stage)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("normal" = "#E6A817",
                                "adenoma" = "#E91E8C",
                                "carcinoma" = "#9B59B6")) +
  facet_wrap(~gene, scales = "free_y", ncol = 4) +
  theme_classic() +
  ggtitle("Downregulated genes along pseudotime") +
  xlab("Pseudotime") + ylab("Expression")

ggsave(file.path(base_dir, "panel_up.pdf"), p_up, width = 16, height = 20)
ggsave(file.path(base_dir, "panel_down.pdf"), p_down, width = 16, height = 10)

# Colorectal signature = the 20 trajectory-dependent genes with the highest Moran's I
# (carried forward to Prong 3A). The hard-coded list is the one used in Prong 3A.
signature <- head(sig_genes$gene_short_name, 20)
stopifnot(setequal(signature,
                   c("ZG16", "GSN", "MUC2", "GUCA2A", "BEST4", "SLC26A3", "AQP8", "OTOP2", "FCGBP", "ANKRD36C",
                     "MALAT1", "HMGB2", "CA4", "TFF3", "PHGR1", "CA7", "TSPAN1", "FTH1", "CCL5", "H2AFZ")))


library(DESeq2)

# Pseudobulk by patient x stage
counts_mat <- AggregateExpression(
  epi_cells,
  group.by = c("patient", "stage"),
  assay = "RNA",
  slot = "counts"
)$RNA

# Build metadata
sample_meta <- data.frame(
  sample_id = colnames(counts_mat),
  stringsAsFactors = FALSE
) %>%
  mutate(
    stage = case_when(
      grepl("_normal$", sample_id)    ~ "normal",
      grepl("_adenoma$", sample_id)   ~ "adenoma",
      grepl("_carcinoma$", sample_id) ~ "carcinoma"
    ),
    patient = gsub("_normal$|_adenoma$|_carcinoma$", "", sample_id),
    stage = factor(stage, levels = c("normal", "adenoma", "carcinoma")),
    patient = factor(patient)
  )
rownames(sample_meta) <- sample_meta$sample_id

# Sanity check before proceeding
sample_meta[, c("patient", "stage")]

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = sample_meta,
  design    = ~ patient + stage
)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

dds <- DESeq(dds)

res_normal_v_adenoma   <- results(dds, contrast = c("stage", "adenoma", "normal"))
res_adenoma_v_carcinoma <- results(dds, contrast = c("stage", "carcinoma", "adenoma"))
res_normal_v_carcinoma  <- results(dds, contrast = c("stage", "carcinoma", "normal"))

res_df_all <- as.data.frame(res_normal_v_carcinoma) %>% filter(!is.na(padj))

saveRDS(res_df_all, file = file.path(base_dir, "DE_Results_all.rds"))   # carcinoma vs normal
# The other two contrasts reported in the Methods
saveRDS(as.data.frame(res_normal_v_adenoma) %>% filter(!is.na(padj)),
        file = file.path(base_dir, "DE_Results_adenoma_vs_normal.rds"))
saveRDS(as.data.frame(res_adenoma_v_carcinoma) %>% filter(!is.na(padj)),
        file = file.path(base_dir, "DE_Results_carcinoma_vs_adenoma.rds"))

# Intersect in BOTH directions (|log2FC| > 1), as for prostate. The colorectal
# trajectory's most autocorrelated genes are largely differentiated colonocyte and
# goblet genes LOST with progression (15 of the 20 signature genes are significantly
# down in carcinoma vs normal, none up); testing upregulation only would discard
# that half of the programme.
res_df_sig <- res_df_all %>% filter(abs(log2FoldChange) > 1, padj < 0.05) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
  arrange(desc(abs(log2FoldChange)))
table(res_df_sig$direction)

# Overlap with the trajectory-dependent genes (q < 0.05). Nearly every DE gene passes
# that filter alone (13,849 genes), so the overlap is also reported at Moran's I > 0.2.
res_df_sig_pseudotime <- res_df_sig %>%
  filter(rownames(.) %in% sig_genes$gene_short_name)
res_df_sig_pseudotime_I02 <- res_df_sig %>%
  filter(rownames(.) %in% sig_genes$gene_short_name[sig_genes$morans_I > 0.2])
table(res_df_sig_pseudotime$direction)       # 521 up / 527 down
table(res_df_sig_pseudotime_I02$direction)   # 75 up / 215 down

# Overlap with the FindMarkers progression panel: upregulated only (kept as originally run)
res_df_sig_progression_panel <- res_df_sig %>%
  filter(direction == "up", rownames(.) %in% panel_up_ranked$gene_short_name)
cat("Carcinoma vs normal, sig up overlapping the upregulated progression panel:",
    nrow(res_df_sig_progression_panel), "\n")   # 88

# ---- Numbers reported in the Methods / Supplementary Methods ----
# (graph_test counts are printed above, stage tests after the boxplot,
# late-pseudotime normal cells after pt_df)

# Signature genes in the pseudobulk DE results, carcinoma vs normal (fold change, adjusted P)
res_df_all %>% filter(rownames(.) %in% signature) %>%
  select(baseMean, log2FoldChange, padj) %>%
  tibble::rownames_to_column("gene") %>%
  mutate(call = case_when(padj < 0.05 & log2FoldChange >  1 ~ "up",
                          padj < 0.05 & log2FoldChange < -1 ~ "down",
                          TRUE ~ "ns")) %>%
  arrange(log2FoldChange)
setdiff(signature, rownames(res_df_all))   # signature genes dropped by the count filter
