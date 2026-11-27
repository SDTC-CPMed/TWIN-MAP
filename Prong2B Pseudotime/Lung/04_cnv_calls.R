# Stage 2 of the LuCA LUAD redo (LOCKED_PLAN.md): CNV calls across the per-donor inferCNV runs.
# Null = non-tumour-labelled epithelial cells from adjacent normal tissue; CNV+ above its 99th percentile,
# the same rule as inferCNV/inferCNV_LUAD.R.
suppressPackageStartupMessages(library(dplyr))
# Data and result folders: set in paths.env (see config.R)
if (!exists("earlyt_path")) local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); d <- if (length(a)) dirname(normalizePath(gsub("~+~", " ", a[1], fixed = TRUE))) else getwd()
  while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "config.R"))) stop("config.R not found: run this script from inside the EARLY-T code folder")
  Sys.setenv(EARLYT_CODE = d); source(file.path(d, "config.R")) })
out_dir <- file.path(earlyt_path("EARLYT_DATA"), "LuCa", "luad_10x_20261103")
report <- file.path(out_dir, "cnv_summary.txt"); cat("", file = report)
say <- function(...) cat(..., "\n", file = report, append = TRUE)
show <- function(x, digits = 2) capture.output(print(as.data.frame(x), digits = digits, row.names = FALSE)) %>% paste(collapse = "\n") %>% say()

donors <- readRDS(file.path(out_dir, "cells.rds"))$donors
runs <- lapply(donors$donor_id, function(d) readRDS(file.path(out_dir, "cnv", paste0(d, ".rds"))))
names(runs) <- donors$donor_id
md <- bind_rows(lapply(runs, `[[`, "calls")) %>% mutate(study = sub("^(Leader_Merad|Lambrechts_Thienpont).*", "\\1", dataset), tissue = ifelse(origin == "normal_adjacent", "normal", "tumour"))
say("=== Stage 2 inferCNV ===")
say("Donors", length(runs), "; minutes per donor: median", round(median(sapply(runs, `[[`, "mins")), 2), ", max", round(max(sapply(runs, `[[`, "mins")), 2), "; genes kept: median", median(sapply(runs, `[[`, "n_genes")), ", range", paste(range(sapply(runs, `[[`, "n_genes")), collapse = "-"))

is_null <- md$compartment == "epithelial" & md$origin == "normal_adjacent" & md$label != "Tumour"
thr <- quantile(md$cnv_cor[is_null], 0.99)
md$cnv_call <- ifelse(md$cnv_cor > thr, "CNV+", "CNV-")
say("Null cells", sum(is_null), "; threshold (99th percentile of split-half correlation)", signif(thr, 3))
say("Null 99th percentile by study"); show(md[is_null, ] %>% group_by(study) %>% summarise(cells = n(), q99 = quantile(cnv_cor, .99), q50 = median(cnv_cor), .groups = "drop"))

say("\nImmune controls: CNV+ share"); ic <- md %>% filter(compartment == "immune_control")
say("Overall", round(mean(ic$cnv_call == "CNV+"), 4), "of", nrow(ic), "; donors above 0:", sum(tapply(ic$cnv_call == "CNV+", ic$donor_id, mean) > 0), "; largest donor share", round(max(tapply(ic$cnv_call == "CNV+", ic$donor_id, mean)), 3))

e <- md %>% filter(compartment == "epithelial")
say("\nCNV+ share by label and tissue"); show(e %>% group_by(label, tissue) %>% summarise(cells = n(), donors = n_distinct(donor_id), cnv_pos = mean(cnv_call == "CNV+"), cor_median = median(cnv_cor), .groups = "drop"), 3)

share <- function(l, t) e %>% filter(label == l, tissue == t) %>% group_by(donor_id) %>% summarise(n = n(), cnv_pos = mean(cnv_call == "CNV+"), .groups = "drop")
tbl <- donors %>% select(donor_id, paired, uicc_stage)
for (x in list(c("Tumour", "tumour", "tum"), c("Transitional", "tumour", "trans_tum"), c("AT2", "tumour", "at2_tum"), c("Transitional", "normal", "trans_nor"), c("AT2", "normal", "at2_nor"))) {
  s <- share(x[1], x[2]); tbl[[paste0(x[3], "_n")]] <- s$n[match(tbl$donor_id, s$donor_id)]; tbl[[paste0(x[3], "_pos")]] <- round(s$cnv_pos[match(tbl$donor_id, s$donor_id)], 2)
}
tbl$tum_pos_n <- sapply(tbl$donor_id, function(d) sum(e$donor_id == d & e$label == "Tumour" & e$tissue == "tumour" & e$cnv_call == "CNV+"))
say("\nBy donor"); show(tbl %>% arrange(desc(tum_pos)))

say("\nStop rule groundwork, counts only (no pseudotime yet)")
say("Paired donors with at least 10 CNV+ tumour cells:", sum(tbl$paired & tbl$tum_pos_n >= 10), "of", sum(tbl$paired))
say("Donors with at least 20 CNV+ tumour cells:", sum(tbl$tum_pos_n >= 20), "of", nrow(tbl))
say("Donors with tumour CNV+ share under 20%:", sum(tbl$tum_pos < .2, na.rm = TRUE))

# ---- The bridge: do transitional cells in tumour samples carry the donor's tumour CNVs? ----
say("\nCorrelation with the donor's labelled-tumour CNV profile, median by label and tissue")
show(e %>% filter(!is.na(cor_to_tumour)) %>% group_by(label, tissue) %>% summarise(cells = n(), median = median(cor_to_tumour), q90 = quantile(cor_to_tumour, .9), .groups = "drop"), 3)
tt <- e %>% filter(label == "Transitional", tissue == "tumour")
say("Transitional cells in tumour samples: CNV+", sum(tt$cnv_call == "CNV+"), "of", nrow(tt), "; donors with at least 20 such cells:", n_distinct(tt$donor_id[tt$donor_id %in% names(which(table(tt$donor_id) >= 20))]))
td <- tt %>% group_by(donor_id) %>% filter(n() >= 20) %>% summarise(n = n(), cnv_pos = mean(cnv_call == "CNV+"), cor_to_tumour = median(cor_to_tumour), .groups = "drop")
say("Spread of CNV+ share across those donors:"); say(paste(capture.output(print(round(quantile(td$cnv_pos), 2))), collapse = "\n"))
say("Spearman, donor-level, transitional CNV+ share vs tumour CNV+ share:", round(cor(td$cnv_pos, tbl$tum_pos[match(td$donor_id, tbl$donor_id)], method = "spearman", use = "complete.obs"), 2))

# ---- Cells available for stage 3 ----
set.seed(14)
traj <- e %>% filter((label == "Tumour" & tissue == "tumour" & cnv_call == "CNV+") | (label %in% c("AT2", "Transitional") & !(tissue == "normal" & cnv_call == "CNV+")))
say("\nCells for stage 3 if CNV+ is required for tumour cells, and normal-tissue CNV+ cells are removed"); show(traj %>% count(label, tissue, cnv_call))
saveRDS(list(calls = md, threshold = thr, donors = tbl), file.path(out_dir, "cnv_calls.rds"))
say("\nDONE")
