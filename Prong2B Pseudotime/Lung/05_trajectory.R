# LuCA LUAD trajectory: Monocle3 pseudotime from normal AT2 cells through transitional cells
# to CNV+ tumour cells, with stop rules 1, 3, 4 and 5. Rule 2 (leave one donor out) is in 06_leave_one_out.R.
# Inputs: cells.rds, genes.rds and counts_<donor>.rds from 01_extract.R; cnv_calls.rds from 04_cnv_calls.R.
suppressPackageStartupMessages({ library(monocle3); library(dplyr); library(Matrix); library(irlba); library(harmony); library(uwot); library(RANN); library(igraph) })

# Data and result folders: set in paths.env (see config.R)
if (!exists("earlyt_path")) local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); d <- if (length(a)) dirname(normalizePath(gsub("~+~", " ", a[1], fixed = TRUE))) else getwd()
  while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "config.R"))) stop("config.R not found: run this script from inside the EARLY-T code folder")
  Sys.setenv(EARLYT_CODE = d); source(file.path(d, "config.R")) })
out_dir <- file.path(earlyt_path("EARLYT_DATA"), "LuCa", "luad_10x_20261103")
report <- file.path(out_dir, "trajectory_summary.txt"); cat("", file = report)
say <- function(...) cat(..., "\n", file = report, append = TRUE)
show <- function(x, digits = 3) capture.output(print(as.data.frame(x), digits = digits, row.names = FALSE)) %>% paste(collapse = "\n") %>% say()

inp <- readRDS(file.path(out_dir, "cells.rds")); genes <- readRDS(file.path(out_dir, "genes.rds"))
cn <- readRDS(file.path(out_dir, "cnv_calls.rds"))
donors <- inp$cells %>% distinct(donor_id) %>% left_join(inp$donors %>% select(donor_id, paired, tumor_stage, uicc_stage), by = "donor_id")
paired <- donors$donor_id[donors$paired]
e <- cn$calls %>% filter(compartment == "epithelial") %>% mutate(group = case_when(label == "Tumour" ~ "Tumour", tissue == "normal" ~ paste(label, "normal"), TRUE ~ paste(label, "tumour sample")))
is_null <- e$origin == "normal_adjacent" & e$label != "Tumour"
q_study <- tapply(e$cnv_cor[is_null], e$study[is_null], quantile, .99)
e$pos_pooled <- e$cnv_call == "CNV+"; e$pos_study <- e$cnv_cor > q_study[e$study]

counts_all <- do.call(cbind, lapply(donors$donor_id, function(d) readRDS(file.path(out_dir, paste0("counts_", d, ".rds")))))
sym <- genes$gene[match(rownames(counts_all), genes$gene_id)]

pick_cells <- function(pos, club = FALSE) {
  labs <- c("AT2", "Transitional", if (club) "Club")
  x <- e %>% filter((label == "Tumour" & tissue == "tumour" & .data[[pos]]) | label %in% labs)
  set.seed(14)
  x <- x %>% group_by(donor_id, label, origin) %>% group_modify(~ slice_sample(.x, n = min(nrow(.x), 500))) %>% ungroup() %>% as.data.frame()
  x$is_pos <- x[[pos]]; rownames(x) <- x$cell
  x
}

# Harmony embedding: genes in >= 1% of cells (mitochondrial, ribosomal, MALAT1, XIST, RPS4Y1 excluded),
# log-normalise, 3,000 variable genes, scale, 30 PCs, Harmony on sample (same recipe as 02_early_check.R)
harmony_embed <- function(m) {
  co <- counts_all[, m$cell]
  use <- rowMeans(co > 0) >= 0.01 & !grepl("^MT-|^RP[SL][0-9]|^MALAT1$|^XIST$|^RPS4Y1$", sym)
  ln <- co[use, ]; ln@x <- log1p(ln@x / rep(Matrix::colSums(co), diff(ln@p)) * 1e4)
  mu <- Matrix::rowMeans(ln); v <- Matrix::rowMeans(ln^2) - mu^2
  hvg <- names(sort(setNames(residuals(loess(v ~ mu, span = 0.3)), names(mu)), decreasing = TRUE))[1:3000]
  sc <- t(scale(t(as.matrix(ln[hvg, ])))); sc[sc > 10] <- 10
  set.seed(14); pc <- irlba(t(sc), nv = 30); pca <- pc$u %*% diag(pc$d); rownames(pca) <- m$cell
  set.seed(14); h <- RunHarmony(pca, m, "sample", verbose = FALSE); rownames(h) <- m$cell
  list(emb = h, hvg = hvg)
}

fit_graph <- function(emb, m, root_group, hvg, seed = 14) {
  set.seed(seed); um <- umap(emb[m$cell, ], n_neighbors = 30, min_dist = 0.1, metric = "euclidean", n_threads = 4, seed = seed, batch = TRUE)
  rownames(um) <- m$cell
  cds <- new_cell_data_set(counts_all[hvg[1:50], m$cell], cell_metadata = m, gene_metadata = data.frame(gene_short_name = hvg[1:50], row.names = hvg[1:50]))
  reducedDim(cds, "UMAP") <- um
  cds <- cluster_cells(cds, resolution = 1e-3, random_seed = 42)
  cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE)
  root_cells <- m$cell[m$group == root_group]
  nodes <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst)
  root <- names(which.min(sqrt(rowSums(sweep(nodes, 2, colMeans(um[root_cells, ]))^2))))
  cds <- order_cells(cds, root_pr_nodes = root)
  g <- principal_graph(cds)[["UMAP"]]
  list(pt = pseudotime(cds)[m$cell], umap = um, root = root, root_degree = igraph::degree(g, root), partitions = length(unique(partitions(cds))), n_nodes = vcount(g), leaves = sum(igraph::degree(g) == 1))
}

geodesic <- function(emb, m, root_group) {
  x <- emb[m$cell, ]; nn <- nn2(x, k = 31)
  g <- graph_from_edgelist(cbind(rep(seq_len(nrow(x)), each = 30), as.vector(t(nn$nn.idx[, -1]))), directed = FALSE)
  E(g)$weight <- pmax(as.vector(t(nn$nn.dists[, -1])), 1e-6); g <- simplify(g, edge.attr.comb = "min")
  rc <- which(m$group == root_group); root <- rc[which.min(rowSums(sweep(x[rc, ], 2, colMeans(x[rc, ]))^2))]
  setNames(as.numeric(distances(g, v = root)), m$cell)
}

group_order <- function(pt, m) {
  md <- function(i) median(pt[i], na.rm = TRUE)
  c(AT2_normal = md(m$group == "AT2 normal"), AT2_tumour_sample = md(m$group == "AT2 tumour sample"), Trans_normal = md(m$group == "Transitional normal"),
    Trans_tumour_CNVneg = md(m$group == "Transitional tumour sample" & !m$is_pos), Trans_tumour_CNVpos = md(m$group == "Transitional tumour sample" & m$is_pos),
    Transitional_all = md(m$label == "Transitional"), Tumour_CNVpos = md(m$group == "Tumour"))
}

stop_rules <- function(pt, m) {
  ok <- is.finite(pt); dec <- rep(NA_integer_, length(pt)); dec[ok] <- ntile(pt[ok], 10)
  d1 <- data.frame(dec, donor = m$donor_id)[ok, ] %>% count(dec, donor) %>% group_by(dec) %>% mutate(s = n / sum(n)) %>% slice_max(s, n = 1, with_ties = FALSE) %>% ungroup()
  r3 <- data.frame(pt, m) %>% filter(donor_id %in% paired, group %in% c("Tumour", "AT2 normal"), is.finite(pt)) %>% group_by(donor_id) %>%
    summarise(tum = median(pt[group == "Tumour"]), at2 = median(pt[group == "AT2 normal"]), n_tum = sum(group == "Tumour"), .groups = "drop") %>% mutate(ok = n_tum >= 10 & tum > at2)
  r4 <- data.frame(pt, m) %>% filter(is_pos, is.finite(pt)) %>% group_by(donor_id) %>% filter(n() >= 20) %>% summarise(n = n(), rho = suppressWarnings(cor(pt, cnv_score, method = "spearman")), rho_cor = suppressWarnings(cor(pt, cnv_cor, method = "spearman")), .groups = "drop")
  go <- group_order(pt, m)
  list(r1_max = max(d1$s), r1_pass = max(d1$s) <= .25, decile = d1, r3_n = sum(r3$ok), r3_of = nrow(r3), r3_pass = sum(r3$ok) >= 0.8 * nrow(r3), r3 = r3,
       r4_pos = sum(r4$rho > 0, na.rm = TRUE), r4_of = sum(!is.na(r4$rho)), r4_med = median(r4$rho, na.rm = TRUE), r4_cor_pos = sum(r4$rho_cor > 0, na.rm = TRUE), r4_pass = sum(r4$rho > 0, na.rm = TRUE) >= (2 / 3) * sum(!is.na(r4$rho)), r4 = r4,
       r5_order = unname(go["AT2_normal"] < go["Transitional_all"] & go["Transitional_all"] < go["Tumour_CNVpos"]), go = go)
}

describe <- function(name, f, m, r) {
  say("\n\n=================", name, "=================")
  say("Cells", nrow(m), " donors", n_distinct(m$donor_id), " finite pseudotime", round(mean(is.finite(f$pt)), 3), "; graph nodes", f$n_nodes, " leaves", f$leaves, " partitions", f$partitions, "; root", f$root, "degree", f$root_degree)
  say("Group medians"); say(paste(capture.output(print(round(r$go, 2))), collapse = "\n"))
  say(sprintf("Rule 1 (donor balance): largest donor share of any decile %.3f : %s", r$r1_max, ifelse(r$r1_pass, "PASS", "FAIL")))
  say("Deciles"); show(r$decile %>% select(dec, donor, cells = n, share = s) %>% arrange(dec), 2)
  say(sprintf("Rule 3 (direction): %d of %d paired donors : %s", r$r3_n, r$r3_of, ifelse(r$r3_pass, "PASS", "FAIL")))
  say(sprintf("Rule 4 (CNV): rho positive in %d of %d donors, median rho %.2f : %s. Reported only, using cnv_cor: %d of %d", r$r4_pos, r$r4_of, r$r4_med, ifelse(r$r4_pass, "PASS", "FAIL"), r$r4_cor_pos, r$r4_of))
  say("Normal AT2 < transitional < CNV+ tumour:", r$r5_order)
}

# ---- Main: Harmony, 32 donors, pooled threshold ----
m <- pick_cells("pos_pooled")
say("=== Stage 3 trajectory ===")
say("Main cells"); show(m %>% count(group, is_pos))
h <- harmony_embed(m)
main <- fit_graph(h$emb, m, "AT2 normal", h$hvg); r_main <- stop_rules(main$pt, m)
describe("Main: Harmony, 32 donors, pooled threshold", main, m, r_main)
say("\nRule 3 by paired donor"); show(r_main$r3 %>% arrange(tum - at2), 2)
say("\nRule 4 by donor"); show(r_main$r4 %>% arrange(rho), 2)

pt <- main$pt; m$pt <- pt
tum <- m %>% filter(group == "Tumour")
rank_all <- rank(pt)
say("\nVariance in pseudotime rank explained, CNV+ tumour cells: donor", round(summary(lm(rank_all[m$group == "Tumour"] ~ tum$donor_id))$r.squared, 3))
span <- tum %>% group_by(donor_id) %>% summarise(n = n(), q10 = quantile(pt, .1), median = median(pt), q90 = quantile(pt, .9), .groups = "drop") %>% mutate(span = q90 - q10)
say("Range of pseudotime", round(max(pt), 2), "; donor 10-90% span of CNV+ tumour: median", round(median(span$span), 2), ", under 13% of range:", sum(span$span < .13 * max(pt)), "of", nrow(span))
say("CNV+ tumour by donor"); show(span %>% arrange(median), 2)
st <- span %>% left_join(donors, by = "donor_id")
say("\nStage, donor medians of CNV+ tumour pseudotime (support only)"); show(st %>% group_by(tumor_stage) %>% summarise(donors = n(), median = median(median), .groups = "drop"), 3)
say("Wilcoxon early vs advanced P", signif(wilcox.test(median ~ tumor_stage, st)$p.value, 2), "; Spearman with UICC stage", round(cor(st$median, as.integer(factor(st$uicc_stage, c("I", "II", "III", "IV"))), method = "spearman"), 2))
say("Kim_Lee only:"); show(st %>% filter(grepl("^Kim_Lee", donor_id)) %>% group_by(tumor_stage) %>% summarise(donors = n(), median = median(median), .groups = "drop"), 3)
say("Spearman of pseudotime with CNV score, all cells:", round(cor(pt, m$cnv_score, method = "spearman"), 2))

# ---- Sensitivity ----
out <- list(main = list(name = "Main", r = r_main))
add <- function(key, name, f, mm, ref = TRUE) {
  r <- stop_rules(f$pt, mm); describe(name, f, mm, r)
  common <- intersect(names(f$pt), names(pt))
  say("Spearman with main pseudotime:", round(cor(f$pt[common], pt[common], method = "spearman", use = "complete.obs"), 3), "on", length(common), "cells")
  out[[key]] <<- list(name = name, r = r, pt = f$pt)
}
add("seed", "Second UMAP seed", fit_graph(h$emb, m, "AT2 normal", h$hvg, seed = 15), m)
geo <- geodesic(h$emb, m, "AT2 normal")
r_geo <- stop_rules(geo, m); say("\n\n================= Geodesic distance in full Harmony embedding =================")
say("Unreachable", sum(!is.finite(geo)), "; Spearman with main", round(cor(geo, pt, method = "spearman", use = "complete.obs"), 3)); say("Group medians"); say(paste(capture.output(print(round(r_geo$go, 2))), collapse = "\n"))
say(sprintf("Rule 1 %.3f; Rule 3 %d of %d; Rule 4 %d of %d, median %.2f; order %s", r_geo$r1_max, r_geo$r3_n, r_geo$r3_of, r_geo$r4_pos, r_geo$r4_of, r_geo$r4_med, r_geo$r5_order))
out$geo <- list(name = "Geodesic, Harmony", r = r_geo)

shared <- m$cell[complete.cases(inp$X_scVI[m$cell, ])]; ms <- m[shared, ]
add("scvi", "X_scVI, 23 donors", fit_graph(inp$X_scVI, ms, "AT2 normal", h$hvg), ms)
r_h23 <- stop_rules(pt[shared], ms)
say("\nMain Harmony fit restricted to the 23 shared donors"); say(paste(capture.output(print(round(r_h23$go, 2))), collapse = "\n"))
r5 <- r_h23$r5_order && out$scvi$r$r5_order
say(sprintf("Rule 5 (both embeddings, 23 shared donors): Harmony %s, X_scVI %s : %s", r_h23$r5_order, out$scvi$r$r5_order, ifelse(r5, "PASS", "FAIL")))

ms2 <- pick_cells("pos_study"); h2 <- harmony_embed(ms2)
add("study", "Study-specific CNV thresholds", fit_graph(h2$emb, ms2, "AT2 normal", h2$hvg), ms2)
mc <- pick_cells("pos_pooled", club = TRUE); h3 <- harmony_embed(mc)
fc <- fit_graph(h3$emb, mc, "Club normal", h3$hvg); add("club", "Club root (Club cells added)", fc, mc)
say("Club median pseudotime"); say(paste(capture.output(print(round(tapply(fc$pt, mc$group, median), 2))), collapse = "\n"))

# ---- Verdict ----
tab <- bind_rows(lapply(out, function(o) data.frame(run = o$name, rule1_max_share = round(o$r$r1_max, 3), rule3 = paste(o$r$r3_n, "of", o$r$r3_of), rule4 = paste(o$r$r4_pos, "of", o$r$r4_of), rule4_median_rho = round(o$r$r4_med, 2), order = o$r$r5_order,
                                                 AT2_normal = round(o$r$go["AT2_normal"], 2), Transitional = round(o$r$go["Transitional_all"], 2), Tumour_CNVpos = round(o$r$go["Tumour_CNVpos"], 2))))
say("\n\n================= Summary ================="); show(tab)
say("\nMain verdict: Rule 1", r_main$r1_pass, "; Rule 3", r_main$r3_pass, "; Rule 4", r_main$r4_pass, "; Rule 5", r5, "; Rule 2 needs stage 4")
saveRDS(list(cells = m, pt = pt, umap = main$umap, embedding = h$emb, hvg = h$hvg, out = out, geodesic = geo), file.path(out_dir, "trajectory.rds"))
say("\nDONE")
