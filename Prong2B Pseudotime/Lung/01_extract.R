# Stage 0 of the LuCA LUAD redo (LOCKED_PLAN.md): raw counts per donor for the discovery cells,
# plus immune and stromal reference cells for inferCNV. raw/X is one large CSR matrix, so rows are
# read in blocks of nearby cells. Finished blocks are kept, so a rerun resumes.
suppressPackageStartupMessages({ library(rhdf5); library(Matrix); library(dplyr) })

# Data and result folders: set in paths.env (see config.R)
if (!exists("earlyt_path")) local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); d <- if (length(a)) dirname(normalizePath(gsub("~+~", " ", a[1], fixed = TRUE))) else getwd()
  while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "config.R"))) stop("config.R not found: run this script from inside the EARLY-T code folder")
  Sys.setenv(EARLYT_CODE = d); source(file.path(d, "config.R")) })
luca_h5 <- file.path(earlyt_path("EARLYT_DATA"), "LuCa", "eff99db0-14fb-4cea-bdda-61c700bb9a06.h5ad")
out_dir <- file.path(earlyt_path("EARLYT_DATA"), "LuCa", "luad_10x_20261103")
block_dir <- file.path(out_dir, "blocks"); dir.create(block_dir, recursive = TRUE, showWarnings = FALSE)
progress <- file.path(out_dir, "extract_progress.txt")
report <- file.path(out_dir, "extract_summary.txt"); cat("", file = report)
say <- function(...) cat(..., "\n", file = report, append = TRUE)
show <- function(x) capture.output(print(x)) %>% paste(collapse = "\n") %>% say()

h <- H5Fopen(luca_h5, flags = "H5F_ACC_RDONLY")
categorical <- function(path) {
  x <- h5read(h, path)
  k <- as.integer(x$codes); k[k < 0] <- NA
  as.character(x$categories)[k + 1]
}
obs_cols <- c("donor_id", "dataset", "sample", "assay", "disease", "origin", "cell_type_tumor", "cell_type_major",
              "tumor_stage", "uicc_stage", "doublet_status")
obs <- as.data.frame(lapply(setNames(obs_cols, obs_cols), function(n) categorical(paste0("/obs/", n))), stringsAsFactors = FALSE)
obs$cell <- as.character(h5read(h, "/obs/_index"))
obs$total_counts <- as.numeric(h5read(h, "/obs/total_counts"))
obs$row <- seq_len(nrow(obs))

# ---- Donors and cells ----
tumour_label <- "^Tumor cells (LUAD|NSCLC)"
v2 <- obs %>% filter(assay == "10x 3' v2", origin %in% c("tumor_primary", "normal_adjacent"))
luad_donors <- unique(obs$donor_id[obs$disease == "lung adenocarcinoma"])
donors <- v2 %>% filter(donor_id %in% luad_donors) %>% group_by(donor_id) %>%
  summarise(tumour = sum(grepl(tumour_label, cell_type_tumor) & origin == "tumor_primary"),
            at2_normal = sum(cell_type_tumor == "Alveolar cell type 2" & origin == "normal_adjacent"),
            dataset = first(dataset), uicc_stage = first(uicc_stage), tumor_stage = first(tumor_stage), .groups = "drop") %>%
  filter(tumour >= 50) %>% mutate(paired = at2_normal >= 20)
stopifnot(nrow(donors) == 32, sum(donors$paired) == 25)

cand <- v2 %>% filter(donor_id %in% donors$donor_id) %>%
  mutate(label = case_when(grepl(tumour_label, cell_type_tumor) ~ "Tumour", cell_type_tumor == "Alveolar cell type 2" ~ "AT2",
                           cell_type_tumor == "transitional club/AT2" ~ "Transitional", cell_type_tumor == "Club" ~ "Club",
                           grepl("^(T cell|NK cell|B cell|Plasma|Macrophage|Monocyte)", cell_type_major) ~ "Ref_immune",
                           grepl("^(Endothelial|Stromal)", cell_type_major) ~ "Ref_stromal", TRUE ~ NA_character_)) %>%
  filter(!is.na(label))
set.seed(14)
ref_cap <- c(Ref_immune = 700, Ref_stromal = 300)
keep_ref <- cand %>% filter(label %in% names(ref_cap)) %>% group_by(donor_id, label) %>%
  group_modify(~ slice_sample(.x, n = min(nrow(.x), ref_cap[[.y$label]]))) %>% ungroup() %>% pull(cell)
cells <- cand %>% filter(!label %in% names(ref_cap) | cell %in% keep_ref) %>% arrange(row)

scvi <- h5read(h, "/obsm/X_scVI"); scanvi <- h5read(h, "/obsm/X_scANVI")
stopifnot(ncol(scvi) == nrow(obs), ncol(scanvi) == nrow(obs))
scvi <- t(scvi[, cells$row, drop = FALSE]); scanvi <- t(scanvi[, cells$row, drop = FALSE])
rownames(scvi) <- rownames(scanvi) <- cells$cell
saveRDS(list(cells = cells, donors = donors, X_scVI = scvi, X_scANVI = scanvi), file.path(out_dir, "cells.rds"))

gene_id <- as.character(h5read(h, "/raw/var/_index"))
gene_name <- tryCatch(categorical("/raw/var/feature_name"), error = function(e) gene_id)
genes <- data.frame(gene_id = gene_id, gene = sub("_ENSG.*$", "", gene_name), stringsAsFactors = FALSE)
saveRDS(genes, file.path(out_dir, "genes.rds"))
rm(obs, v2, cand, scvi, scanvi); invisible(gc())

# ---- Raw counts, in blocks of nearby rows ----
indptr <- as.numeric(h5read(h, "/raw/X/indptr"))
max_nnz <- 1e7
block <- integer(nrow(cells)); b <- 1; first_ptr <- indptr[cells$row[1]]
for (i in seq_len(nrow(cells))) {
  if (indptr[cells$row[i] + 1] - first_ptr > max_nnz) { b <- b + 1; first_ptr <- indptr[cells$row[i]] }
  block[i] <- b
}
n_block <- max(block)
for (b in seq_len(n_block)) {
  f <- file.path(block_dir, sprintf("block_%04d.rds", b))
  if (file.exists(f)) next
  t0 <- Sys.time()
  idx <- which(block == b); rows <- cells$row[idx]
  from <- indptr[rows[1]]; to <- indptr[rows[length(rows)] + 1]
  x <- as.numeric(h5read(h, "/raw/X/data", start = from + 1, count = to - from))
  j <- as.integer(h5read(h, "/raw/X/indices", start = from + 1, count = to - from)) + 1L
  n <- indptr[rows + 1] - indptr[rows]; at <- indptr[rows] - from
  take <- unlist(lapply(seq_along(rows), function(k) if (n[k] > 0) seq(at[k] + 1, at[k] + n[k]) else integer(0)))
  saveRDS(list(i = j[take], j = rep(idx, n), x = x[take]), f, compress = FALSE)
  cat(sprintf("%s block %d of %d: %d cells, %.0f s\n", format(Sys.time(), "%H:%M"), b, n_block, length(idx),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))), file = progress, append = TRUE)
}
H5Fclose(h); h5closeAll()

parts <- lapply(seq_len(n_block), function(b) readRDS(file.path(block_dir, sprintf("block_%04d.rds", b))))
counts <- sparseMatrix(i = unlist(lapply(parts, `[[`, "i")), j = unlist(lapply(parts, `[[`, "j")), x = unlist(lapply(parts, `[[`, "x")),
                       dims = c(nrow(genes), nrow(cells)), dimnames = list(genes$gene_id, cells$cell))
rm(parts); invisible(gc())

# ---- Checks ----
say("=== Stage 0 extraction ===")
say("Donors:", nrow(donors), " paired:", sum(donors$paired), " cells:", ncol(counts), " genes:", nrow(counts))
say("Cells by label x origin"); show(table(cells$label, cells$origin))
say("Doublet status"); show(table(cells$doublet_status, useNA = "ifany"))
say("All values integer:", all(counts@x == round(counts@x)), " max:", max(counts@x))
lib <- Matrix::colSums(counts)
say(sprintf("Library size vs obs total_counts: Spearman %.4f, median ratio %.3f, cells differing by more than 1%%: %d",
            cor(lib, cells$total_counts, method = "spearman"), median(lib / cells$total_counts), sum(abs(lib / cells$total_counts - 1) > 0.01)))
say("Median library size by label"); show(tapply(lib, cells$label, median))
say("Median genes detected by label"); show(tapply(Matrix::colSums(counts > 0), cells$label, median))
say("Donors"); show(as.data.frame(donors))

for (d in donors$donor_id) saveRDS(counts[, cells$donor_id == d, drop = FALSE], file.path(out_dir, paste0("counts_", d, ".rds")))
say("\nDONE")
