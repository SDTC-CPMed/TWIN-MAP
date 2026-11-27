# Signature pipeline: from a pseudotime trajectory to the trajectory signature.
# Every analysis decision is read from a config file; see PIPELINE.md.
#   Rscript signature_pipeline.R <graph_test|signature> <config.json>
#   Rscript signature_pipeline.R diagnose <significance|signature> <config.json>
suppressPackageStartupMessages(library(dplyr))

# Data and result folders: set in paths.env (see config.R)
if (!exists("earlyt_path")) local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); d <- if (length(a)) dirname(normalizePath(gsub("~+~", " ", a[1], fixed = TRUE))) else getwd()
  while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "config.R"))) stop("config.R not found: run this script from inside the EARLY-T code folder")
  Sys.setenv(EARLYT_CODE = d); source(file.path(d, "config.R")) })

args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: Rscript signature_pipeline.R <graph_test|signature> <config.json>\n       Rscript signature_pipeline.R diagnose <significance|signature> <config.json>"
if (length(args) < 2) stop(usage, call. = FALSE)
stage <- args[1]
diag_step <- if (stage == "diagnose") args[2] else NULL
cfg_file <- normalizePath(args[length(args)], mustWork = TRUE)
cfg <- jsonlite::read_json(cfg_file, simplifyVector = FALSE)

# ---- Config access: no defaults for analysis decisions ----
need <- function(path, options = NULL) {
  x <- cfg
  for (p in strsplit(path, ".", fixed = TRUE)[[1]]) x <- if (is.list(x)) x[[p]] else NULL
  if (is.null(x)) stop(sprintf("`%s` is not set in %s%s (see PIPELINE.md)", path, basename(cfg_file),
                               if (length(options)) paste0(": choose ", paste(sprintf('"%s"', options), collapse = ", ")) else ""), call. = FALSE)
  if (length(options) && is.character(x) && !all(x %in% options))
    stop(sprintf("`%s` is \"%s\": must be one of %s", path, x, paste(options, collapse = ", ")), call. = FALSE)
  x
}
has <- function(path) { x <- cfg; for (p in strsplit(path, ".", fixed = TRUE)[[1]]) x <- if (is.list(x)) x[[p]] else NULL; !is.null(x) }
expand <- function(p) { if (is.null(p)) return(NULL)
  for (v in regmatches(p, gregexpr("\\$\\{[A-Z_]+\\}", p))[[1]]) p <- sub(v, earlyt_path(gsub("[${}]", "", v)), p, fixed = TRUE); p }
chr <- function(x) unlist(x, use.names = FALSE)

name <- need("name")
out_dir <- expand(need("out_dir")); run_dir <- file.path(out_dir, name); gt_dir <- file.path(out_dir, "graph_test")
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
report <- file.path(run_dir, "report.txt")
say <- function(...) { txt <- paste(..., collapse = " "); cat(txt, "\n"); if (is.null(diag_step)) cat(txt, "\n", file = report, append = TRUE) }
show <- function(x, digits = 3) say(paste(capture.output(print(as.data.frame(x), digits = digits, row.names = FALSE)), collapse = "\n"))

# ---- Shared helpers ----
load_cds <- function() { suppressPackageStartupMessages(library(monocle3)); readRDS(expand(need("cds"))) }

gt_file <- function() { s <- expand(cfg$graph_test$saved_results); if (!is.null(s)) s else file.path(gt_dir, "graph_test_results.rds") }
load_gt <- function() {
  f <- gt_file(); if (!file.exists(f)) stop("graph_test results not found: ", f, " (run the graph_test stage, or set graph_test.saved_results)", call. = FALSE)
  r <- readRDS(f); if (is.null(r$gene_short_name)) r$gene_short_name <- rownames(r)
  r %>% transmute(gene = gene_short_name, status = as.character(status), p_value, morans_I, q_value)
}
trajectory_genes <- function(gt) gt %>% filter(status == "OK", q_value < need("significance.q_max")) %>% arrange(desc(morans_I))

# The point of the ranked curve furthest from the line joining its ends
elbow <- function(v) {
  s <- sort(v, decreasing = TRUE); n <- length(s); a <- seq_len(n) / n; b <- (s - min(s)) / (max(s) - min(s))
  d <- abs((a[n] - a[1]) * (b[1] - b) - (a[1] - a) * (b[n] - b[1])) / sqrt((a[n] - a[1])^2 + (b[n] - b[1])^2)
  s[which.max(d)]
}

flag_genes <- function(g) case_when(
  grepl("^MT-", g) ~ "mitochondrial", grepl("^RP[SL][0-9]", g) ~ "ribosomal", g %in% c("MALAT1", "NEAT1") ~ "nuclear lncRNA",
  g %in% c("XIST", "TSIX", "RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "ZFY", "USP9Y") ~ "sex-linked",
  grepl("^HIST|^H[1-4]|^H2A|^H2B", g) ~ "histone", grepl("^FOS|^JUN|^EGR[1-3]$|^HSPA1|^DNAJB1$|^IER[0-9]", g) ~ "stress response",
  TRUE ~ "")

platform <- function() { f <- expand(cfg$signature$platform_genes); if (is.null(f)) NULL else read.csv(f)[[1]] }

select_signature <- function(tg) {
  rule <- need("signature.rule", c("top_n", "morans_I_min", "elbow"))
  pats <- chr(need("signature.exclude_patterns")); ex <- chr(need("signature.exclude_genes"))
  keep <- !(tg$gene %in% ex); for (p in pats) keep <- keep & !grepl(p, tg$gene)
  cand <- tg[keep, ]
  sel <- switch(rule,
    top_n        = head(cand, need("signature.n")),
    morans_I_min = cand %>% filter(morans_I > need("signature.morans_I_min")),
    elbow        = cand %>% filter(morans_I >= elbow(cand$morans_I)))
  list(sig = sel, excluded = tg$gene[!keep])
}

# ---- Stage A: graph_test ----
run_graph_test <- function() {
  if (!is.null(cfg$graph_test$saved_results)) stop("graph_test.saved_results is set: stage A is skipped and those results are used", call. = FALSE)
  cds <- load_cds(); dir.create(gt_dir, showWarnings = FALSE)
  bs <- cfg$graph_test$batch_size; genes <- rownames(cds)
  batches <- if (is.null(bs)) list(`1` = genes) else split(genes, ceiling(seq_along(genes) / bs))
  for (b in names(batches)) {
    f <- file.path(gt_dir, sprintf("batch_%03d.rds", as.integer(b))); if (file.exists(f)) next
    t0 <- Sys.time()
    saveRDS(graph_test(cds[batches[[b]], ], neighbor_graph = need("graph_test.neighbor_graph", c("principal_graph", "knn")),
                       cores = need("graph_test.cores")), f)
    cat(sprintf("batch %s of %d: %.1f min\n", b, length(batches), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  res <- bind_rows(lapply(names(batches), function(b) readRDS(file.path(gt_dir, sprintf("batch_%03d.rds", as.integer(b))))))
  if (is.null(res$gene_short_name)) res$gene_short_name <- rownames(res)
  stopifnot(nrow(res) == nrow(cds), !anyDuplicated(res$gene_short_name))
  res$q_value <- 1; ok <- res$status == "OK" & !is.na(res$p_value)
  res$q_value[ok] <- p.adjust(res$p_value[ok], method = "BH")   # over all genes, so batching does not change q
  saveRDS(res, file.path(gt_dir, "graph_test_results.rds"))
  say("graph_test: genes", nrow(res), "; status OK", sum(ok))
}

# ---- Steps 2-4: significance, report, signature ----
run_signature <- function() {
  gt <- load_gt(); tg <- trajectory_genes(gt)
  say("== Steps 2-4: trajectory genes and signature ==")
  say("graph_test results:", gt_file())
  say("graph_test -- tested:", nrow(gt), "| converged:", sum(gt$status == "OK"), "| q <", need("significance.q_max"), ":", nrow(tg))
  if (all(gt$q_value < need("significance.q_max"))) say("NOTE: the results file holds only significant genes, so 'tested' and 'converged' are not the full counts")
  cut <- unlist(need("report.morans_I_cutoffs"))
  show(data.frame(morans_I_cutoff = cut, n_genes = sapply(cut, function(t) sum(tg$morans_I > t))))
  write.csv(tg, file.path(run_dir, "trajectory_genes.csv"), row.names = FALSE)

  s <- select_signature(tg); sig <- s$sig
  if (length(s$excluded)) say("Excluded before the rule:", length(s$excluded), "genes")
  pg <- platform()
  out <- sig %>% mutate(rank = row_number(), flag = flag_genes(gene), on_platform = if (is.null(pg)) NA else gene %in% pg) %>%
    select(rank, gene, morans_I, q_value, flag, on_platform)
  write.csv(out, file.path(run_dir, "signature.csv"), row.names = FALSE)
  say(sprintf("Signature (%s): %d genes", need("signature.rule"), nrow(out)))
  say(paste(out$gene, collapse = ", "))
  if (any(out$flag != "")) say("Flagged:", paste(sprintf("%s (%s)", out$gene[out$flag != ""], out$flag[out$flag != ""]), collapse = ", "))
  if (!is.null(pg)) say("Measured on platform:", sum(out$on_platform), "of", nrow(out))
  exp <- chr(cfg$signature$expected)
  if (length(exp)) {
    if (!setequal(exp, out$gene)) stop("Signature differs from signature.expected. Missing: ", paste(setdiff(exp, out$gene), collapse = ", "),
                                       "; extra: ", paste(setdiff(out$gene, exp), collapse = ", "), call. = FALSE)
    say("Matches signature.expected")
  }
}

# ---- Diagnostics: evidence for each decision, changes nothing ----
run_diagnose <- function(step) {
  if (step == "significance") {
    gt <- load_gt(); ok <- gt %>% filter(status == "OK")
    show(data.frame(q_max = c(0.01, 0.05, 0.1), trajectory_genes = sapply(c(0.01, 0.05, 0.1), function(q) sum(ok$q_value < q))))
  } else if (step == "signature") {
    tg <- trajectory_genes(load_gt()); el <- elbow(tg$morans_I)
    show(data.frame(rule = c("top 10", "top 20", "top 50", "Moran's I > 0.2", "Moran's I > 0.3", "Moran's I > 0.5", sprintf("elbow (I >= %.3f)", el)),
                    genes = c(10, 20, 50, sum(tg$morans_I > .2), sum(tg$morans_I > .3), sum(tg$morans_I > .5), sum(tg$morans_I >= el))))
    top <- head(tg, 50) %>% mutate(rank = row_number(), flag = flag_genes(gene))
    pg <- platform(); if (!is.null(pg)) top$on_platform <- top$gene %in% pg
    show(top %>% select(rank, gene, morans_I, any_of(c("flag", "on_platform"))))
  } else stop("diagnose step must be significance or signature", call. = FALSE)
}

# ---- Run ----
if (is.null(diag_step)) {
  if (stage == "signature") cat("", file = report)
  cat(sprintf("\n# %s, stage %s, %s\n", name, stage, format(Sys.time())), file = report, append = TRUE)
  file.copy(cfg_file, file.path(run_dir, "config.json"), overwrite = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
}
switch(stage,
  graph_test = run_graph_test(),
  signature  = run_signature(),
  diagnose   = run_diagnose(diag_step),
  stop(usage, call. = FALSE))
