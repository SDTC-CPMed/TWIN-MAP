# Stage 7 of the LuCA LUAD redo (LOCKED_PLAN.md): Moran's I along the principal graph of the main
# trajectory, as in Normal to Lung Adenocarcinoma.R. Genes run in batches with checkpoints, and the
# q values are recomputed over all genes at the end. Cells are the unit, so P values ignore donor.
suppressPackageStartupMessages({ library(monocle3); library(dplyr); library(Matrix); library(igraph) })

# Data and result folders: set in paths.env (see config.R)
if (!exists("earlyt_path")) local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); d <- if (length(a)) dirname(normalizePath(gsub("~+~", " ", a[1], fixed = TRUE))) else getwd()
  while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "config.R"))) stop("config.R not found: run this script from inside the EARLY-T code folder")
  Sys.setenv(EARLYT_CODE = d); source(file.path(d, "config.R")) })
out_dir <- file.path(earlyt_path("EARLYT_DATA"), "LuCa", "luad_10x_20261103")
gt_dir <- file.path(out_dir, "graph_test"); dir.create(gt_dir, showWarnings = FALSE)
report <- file.path(gt_dir, "graph_test_summary.txt"); cat("", file = report)
say <- function(...) cat(..., "\n", file = report, append = TRUE)

tr <- readRDS(file.path(out_dir, "trajectory.rds")); m <- tr$cells; rownames(m) <- m$cell
genes <- readRDS(file.path(out_dir, "genes.rds"))
counts <- do.call(cbind, lapply(unique(m$donor_id), function(d) {
  co <- readRDS(file.path(out_dir, paste0("counts_", d, ".rds")))
  co[, colnames(co) %in% m$cell, drop = FALSE]
}))[, m$cell]
sym <- genes$gene[match(rownames(counts), genes$gene_id)]
keep <- !is.na(sym); counts <- counts[keep, ]; rownames(counts) <- make.unique(sym[keep])

# The same graph as 05_trajectory.R, fitted on the saved UMAP
cds <- new_cell_data_set(counts, cell_metadata = m, gene_metadata = data.frame(gene_short_name = rownames(counts), row.names = rownames(counts)))
reducedDim(cds, "UMAP") <- tr$umap[m$cell, ]
cds <- cluster_cells(cds, resolution = 1e-3, random_seed = 42)
cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE)
root_cells <- m$cell[m$group == "AT2 normal"]
nodes <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst)
root <- names(which.min(sqrt(rowSums(sweep(nodes, 2, colMeans(tr$umap[root_cells, ]))^2))))
cds <- order_cells(cds, root_pr_nodes = root)
pt <- pseudotime(cds)[m$cell]
say("=== Stage 7 graph_test ===")
say("Cells", ncol(cds), "; genes", nrow(cds), "; largest pseudotime difference from the saved trajectory", signif(max(abs(pt - tr$pt[m$cell])), 3))
stopifnot(max(abs(pt - tr$pt[m$cell])) < 1e-6)
saveRDS(cds, file.path(gt_dir, "cds_main.rds"))

batches <- split(rownames(cds), ceiling(seq_len(nrow(cds)) / 2000))
for (b in names(batches)) {
  f <- file.path(gt_dir, sprintf("batch_%03d.rds", as.integer(b)))
  if (file.exists(f)) next
  t0 <- Sys.time()
  res <- graph_test(cds[batches[[b]], ], neighbor_graph = "principal_graph", cores = 4)
  saveRDS(res, f)
  cat(sprintf("%s batch %s of %d: %.1f min\n", format(Sys.time()), b, length(batches), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

res <- bind_rows(lapply(names(batches), function(b) readRDS(file.path(gt_dir, sprintf("batch_%03d.rds", as.integer(b))))))
stopifnot(nrow(res) == nrow(cds), !anyDuplicated(res$gene_short_name))
res$q_value <- 1; ok <- res$status == "OK" & !is.na(res$p_value)
res$q_value[ok] <- p.adjust(res$p_value[ok], method = "BH")
saveRDS(res, file.path(gt_dir, "pseudotime_graph_test_results_luca.rds"))

sig <- res %>% filter(status == "OK", q_value < 0.05) %>% arrange(desc(morans_I))
write.csv(sig, file.path(gt_dir, "graph_test_sig_status_OK.csv"), row.names = FALSE)
write.csv(data.frame(gene = res$gene_short_name[ok]), file.path(gt_dir, "graph_test_tested_genes.csv"), row.names = FALSE)
say("Status OK", sum(ok), "; FAIL", sum(!ok), "; status OK and q < 0.05:", nrow(sig))
say("Moran's I in significant genes, quartiles", paste(round(quantile(sig$morans_I, c(.25, .5, .75)), 3), collapse = ", "))
say("Genes at Moran's I above 0.1, 0.2, 0.3, 0.5:", paste(sapply(c(.1, .2, .3, .5), function(t) sum(sig$morans_I > t)), collapse = ", "))
say("\nDONE")
