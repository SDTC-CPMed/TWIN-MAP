# Shared inferCNV helpers for the CRC and prostate pseudotime support analyses (13 Oct 2026)
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(infercnv)  # EnsDb used via :: to avoid masking dplyr::filter
})

# Autosomal gene positions (Ensembl 86 / GRCh38) for the genes present in the data
gene_order_for <- function(genes) {
  g <- ensembldb::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, columns = "gene_name",
                        filter = AnnotationFilter::SeqNameFilter(as.character(1:22)))
  go <- data.frame(gene = g$gene_name, chr = as.integer(as.character(GenomicRanges::seqnames(g))),
                   start = BiocGenerics::start(g), stop = BiocGenerics::end(g)) %>%
    filter(gene %in% genes, gene != "") %>% distinct(gene, .keep_all = TRUE) %>% arrange(chr, start)
  data.frame(chr = paste0("chr", go$chr), start = go$start, stop = go$stop, row.names = go$gene)
}

run_infercnv <- function(counts, ann, ref_groups, out_dir) {
  if (file.exists(file.path(out_dir, "run.final.infercnv_obj")))
    return(readRDS(file.path(out_dir, "run.final.infercnv_obj")))  # reuse a finished run
  go <- gene_order_for(rownames(counts))
  obj <- CreateInfercnvObject(raw_counts_matrix = counts[rownames(go), rownames(ann)], annotations_file = ann,
                              gene_order_file = go, ref_group_names = ref_groups)
  infercnv::run(obj, cutoff = 0.1, out_dir = out_dir, cluster_by_groups = TRUE,
                denoise = TRUE, HMM = FALSE, num_threads = 8, output_format = "png",
                no_prelim_plot = TRUE, save_rds = FALSE)
}

# Cells whose denoised profile is flat (all 1) have undefined correlation: no CNV signal -> 0
safe_cor <- function(X, prof) { r <- suppressWarnings(cor(X, prof)[, 1]); r[is.na(r)] <- 0; r }

cnv_profile <- function(E, cells, score, top = 0.1) {
  s <- score[cells]
  rowMeans(E[, cells[s >= quantile(s, 1 - top)], drop = FALSE]) - 1
}

# Split-half correlation: each cell is correlated with the CNV profile (top-10% CNV-score
# cells) of the other half of its own sample, so no cell contributes to its own reference profile.
split_half_cor <- function(E, sample, score, seed = 1) {
  set.seed(seed)
  out <- setNames(rep(NA_real_, ncol(E)), colnames(E))
  for (s in unique(sample)) {
    cells <- colnames(E)[sample == s]
    half  <- sample(rep(1:2, length.out = length(cells)))
    for (h in 1:2) {
      prof <- cnv_profile(E, cells[half != h], score)
      out[cells[half == h]] <- safe_cor(E[, cells[half == h], drop = FALSE] - 1, prof)
    }
  }
  out
}

# CNV signal along pseudotime: share of CNV+ cells and median correlation per decile
by_pseudotime <- function(md, group = NULL) {
  md$pt_decile <- dplyr::ntile(md$pseudotime, 10)
  md %>% group_by(across(all_of(c(group, "pt_decile")))) %>%
    summarise(n = n(), pt_median = median(pseudotime), cnv_pos = mean(cnv_call == "CNV+"),
              cor_median = median(cnv_cor), .groups = "drop")
}
