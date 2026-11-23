library(SingleCellExperiment)
library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)
library(SeuratWrappers)

# Data and result folders: set in paths.env (see config.R)
if (!exists("earlyt_path")) local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); d <- if (length(a)) dirname(normalizePath(gsub("~+~", " ", a[1], fixed = TRUE))) else getwd()
  while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "config.R"))) stop("config.R not found: run this script from inside the EARLY-T code folder")
  Sys.setenv(EARLYT_CODE = d); source(file.path(d, "config.R")) })
base_dir <- file.path(earlyt_path("EARLYT_DATA"), "Prostate Cancer")

# Load the Seurat object
seurat_obj <- readRDS(file.path(base_dir, "dge_E.rds"))
sce <- NormalizeData(seurat_obj) 

DimPlot(sce, group.by = "malignancy")
DimPlot(sce, group.by = "orig.ident")

# plot UMAP
DimPlot(sce, reduction = "umap", group.by = "malignancy")

umap_embeddings <- Embeddings(sce, "umap")

cds <- as.cell_data_set(sce)

reducedDims(cds)[["UMAP"]] <- umap_embeddings

cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)

plot_cells(cds,
           color_cells_by = "malignancy",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

# Root: originally chosen interactively (order_cells(cds) opens a GUI where the
# root nodes are clicked). The 29 nodes selected that way (all in the
# non-malignant compartment) were saved to root_nodes.rds; reading them back
# makes the ordering reproducible without the GUI.
root_nodes <- readRDS(file.path(base_dir, "root_nodes.rds"))
stopifnot(all(root_nodes %in% igraph::V(principal_graph(cds)[["UMAP"]])$name))
cds <- order_cells(cds, root_pr_nodes = root_nodes)
# Saved to a subfolder so a rerun doesn't overwrite the saved cds_object.rds in base_dir;
# move it up deliberately if the rerun is to replace it.
monocle3::save_monocle_objects(cds, file.path(base_dir, "monocle_cds"), verbose = TRUE)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

pseudotime_res <- graph_test(cds, neighbor_graph="principal_graph")

saveRDS(pseudotime_res, file.path(base_dir, "pseudotime_res_normal_to_prostate_cancer.rds"))

sig_genes <- pseudotime_res %>% arrange(desc(abs(pseudotime_res$morans_I))) %>%
  filter(status == "OK") %>% 
  filter(q_value < 0.05) %>% filter(morans_I > 0.5)

saveRDS(sig_genes, file.path(base_dir, "sig_genes_normal_to_prostate_cancer.rds"))


plot_genes_in_pseudotime(cds[c("AR", "KLK3", "C3"),],
                         color_cells_by = "malignancy",label_by_short_name = F)


genes_to_plot <- unique(rownames(sig_genes))

# Extract expression matrix for these genes
expr_df <- FetchData(seurat_obj,
                     vars = c(genes_to_plot, "malignancy")) %>%
  tibble::rownames_to_column("cell") %>%
  mutate(pseudotime = pseudotime(cds)[cell]) %>%
  filter(is.finite(pseudotime)) %>%
  tidyr::pivot_longer(cols = all_of(genes_to_plot),
                      names_to = "gene",
                      values_to = "expression")

#expression curves — one panel per gene
ggplot(expr_df, aes(x = pseudotime, y = expression, color = malignancy)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("Malignant"       = "#E6A817",
                                "Non-Malignant"   = "#E91E8C"
                                )) +
  facet_wrap(~gene, scales = "free_y", ncol = 3) +
  theme_classic() +
  ggtitle("Smoothed gene expression along pseudotime") +
  xlab("Pseudotime") + ylab("Expression")


library(DESeq2)

# Pseudobulk by sample
counts_mat <- AggregateExpression(
  seurat_obj,
  group.by = c("orig.ident", "malignancy"),
  assay = "RNA",
  slot = "counts"
)$RNA

sig_genes <- readRDS(file.path(base_dir, "sig_genes_normal_to_prostate_cancer.rds"))
# Pseudobulk DE is run genome-wide (as in the saved DE_results_all.rds, 13,502
# genes). Restricting counts to the 11 trajectory genes before DESeq() would
# distort size-factor/dispersion estimation and the BH adjustment.

# Build metadata
sample_meta <- data.frame(
  sample_id = colnames(counts_mat),
  stringsAsFactors = FALSE
) %>%
  mutate(
    condition = case_when(
      grepl("Non-Malignant", sample_id) ~ "Non_Malignant",
      grepl("Malignant", sample_id)     ~ "Malignant"
    ),
    orig.ident = gsub("_Malignant$|_Non-Malignant$", "", sample_id),
    condition = factor(condition, levels = c("Non_Malignant", "Malignant"))
  )
rownames(sample_meta) <- sample_meta$sample_id

# Check cell counts per pseudobulk group before proceeding
table(seurat_obj@meta.data$orig.ident, seurat_obj@meta.data$malignancy)
sample_meta[, c("orig.ident", "condition")]

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = sample_meta,
  design    = ~ orig.ident + condition
)


keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Malignant", "Non_Malignant"))

summary(res)
plotMA(res)

res_df_all <- as.data.frame(res) %>% filter(!is.na(padj))

saveRDS(res_df_all, file.path(base_dir, "DE_results_all.rds"))
# Intersect in BOTH directions (|log2FC| > 1). The prostate trajectory runs from
# non-malignant to malignant epithelium, and its most autocorrelated genes are
# largely luminal/secretory and epithelial-identity genes LOST with malignancy
# (MSMB, OLFM4, KRT19, PIGR, KRT15, C3 are all significantly down). Testing
# upregulation only would discard that half of the programme and left the
# intersection empty (no signature gene is significantly up).
res_df_sig <- res_df_all %>% filter(abs(log2FoldChange) > 1, padj < 0.05) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
  arrange(desc(abs(log2FoldChange)))

sig_genes <- sig_genes %>% filter(status == 'OK',q_value < 0.05)
res_df_sig_and_pseudotime <- res_df_sig %>% filter(rownames(.) %in% rownames(sig_genes))
table(res_df_sig_and_pseudotime$direction)

# ---- Numbers reported in the Methods / Supplementary Methods ----

# graph_test summary: genes tested, converged, significant, and by Moran's I cut-off
r_ok <- pseudotime_res %>% filter(status == "OK")
cat("graph_test -- tested:", nrow(pseudotime_res),
    "| converged:", nrow(r_ok),
    "| q < 0.05:", sum(r_ok$q_value < 0.05, na.rm = TRUE), "\n")
data.frame(morans_I_cutoff = c(0.1, 0.2, 0.3, 0.5),
           n_genes = sapply(c(0.1, 0.2, 0.3, 0.5), function(t)
             sum(r_ok$q_value < 0.05 & r_ok$morans_I > t, na.rm = TRUE)))

# Pseudobulk DE totals
cat("DE genes (|log2FC| > 1, padj < 0.05):", nrow(res_df_sig),
    "| up:", sum(res_df_sig$direction == "up"),
    "| down:", sum(res_df_sig$direction == "down"), "\n")

# Signature genes in the pseudobulk DE results (fold change, adjusted P)
res_df_all %>% filter(rownames(.) %in% rownames(sig_genes)) %>%
  select(baseMean, log2FoldChange, padj) %>%
  tibble::rownames_to_column("gene")
setdiff(rownames(sig_genes), rownames(res_df_all))   # signature genes dropped by the count filter
