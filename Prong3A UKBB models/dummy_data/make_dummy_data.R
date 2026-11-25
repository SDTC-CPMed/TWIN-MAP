# Synthetic stand-in for the UK Biobank inputs of UKBB_model.R. Every value is randomly generated; nothing is taken
# from UK Biobank. Only the file names, column names and answer codes follow the UK Biobank layout the model reads.
# A weak signal is planted on purpose: each organ's signature proteins and known marker shift in that organ's
# cancer cases, more strongly closer to diagnosis. Results from these data mean nothing about real biology.
#   Rscript make_dummy_data.R [output folder]      (default: input/ next to this script)
suppressPackageStartupMessages(library(jsonlite))
set.seed(2026)

args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)))
out <- if (length(args)) args[1] else file.path(script_dir, "input")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
cfg <- lapply(c(prostate = "prostate", colorectal = "colorectal", lung = "lung"),
              function(o) read_json(file.path(script_dir, "..", "config", paste0(o, ".json")), simplifyVector = TRUE))

# ---- Settings ----
n_people     <- 5000
n_incident   <- 260   # incident cancer cases per organ, diagnosed 0-14 years after blood draw
n_progress   <- 170   # of those, also diagnosed with the organ's precursor condition first
n_prevalent  <- 20    # cancer diagnosed before blood draw (excluded by the model)
n_dis1_only  <- 600   # precursor condition, no cancer
n_other_canc <- 150   # other cancers (in AllCancers, not in any organ's cancer file)
n_hc         <- 1200  # healthy controls
n_dummy_prot <- 230   # proteins outside every signature (the random-panel pool)
shift_at_dx  <- 0.3   # planted shift in SD units at diagnosis ...
decay_years  <- 6     # ... falling by exp(-years before diagnosis / decay_years)
censor_date  <- as.Date("2022-12-30")

# ---- Participants and core clinical fields ----
eid <- as.character(9000000 + seq_len(n_people))
assess <- as.Date("2006-03-01") + sample(0:1700, n_people, replace = TRUE)
sex <- sample(c("Female", "Male"), n_people, replace = TRUE)
pick <- function(levels, p, n = n_people) sample(levels, n, replace = TRUE, prob = p)
clinic <- data.frame(
  eid = eid,
  date_of_attending_assessment_centre_f53_0_0 = assess,
  date_of_death_f40000_0_0 = as.Date(NA), date_of_death_f40000_1_0 = as.Date(NA),
  date_lost_to_followup_f191_0_0 = as.Date(NA),
  age_at_recruitment_f21022_0_0 = sample(40:70, n_people, replace = TRUE),
  sex_f31_0_0 = factor(sex),
  ethnic_background_f21000_0_0 = factor(pick(c("British", "Irish", "Any other white background", "Indian", "Pakistani",
                                               "Chinese", "Caribbean", "African", "Any other mixed background",
                                               "Other ethnic group", "Prefer not to answer"),
                                             c(.80, .03, .04, .02, .01, .01, .02, .02, .02, .02, .01))),
  uk_biobank_assessment_centre_f54_0_0 = factor(pick(paste("Centre", 1:10), rep(.1, 10))),
  body_mass_index_bmi_f21001_0_0 = round(rnorm(n_people, 27.4, 4.7), 2),
  glycated_haemoglobin_hba1c_f30750_0_0 = round(rnorm(n_people, 36, 5), 1),
  smoking_status_f20116_0_0 = factor(pick(c("Never", "Previous", "Current", "Prefer not to answer"), c(.55, .34, .10, .01))),
  alcohol_intake_frequency_f1558_0_0 = factor(pick(c("Daily or almost daily", "Three or four times a week", "Once or twice a week",
                                                     "One to three times a month", "Special occasions only", "Never",
                                                     "Prefer not to answer"), c(.20, .23, .26, .11, .12, .07, .01))),
  diabetes_diagnosed_by_doctor_f2443_0_0 = factor(pick(c("No", "Yes", "Do not know", "Prefer not to answer"), c(.94, .05, .005, .005))),
  stringsAsFactors = FALSE)
clinic$body_mass_index_bmi_f21001_0_0[runif(n_people) < .01] <- NA
clinic$glycated_haemoglobin_hba1c_f30750_0_0[runif(n_people) < .07] <- NA

# ---- Disease groups (disjoint across organs; prostate groups are men only) ----
free <- eid
take <- function(n, male_only = FALSE) {
  pool <- if (male_only) free[sex[match(free, eid)] == "Male"] else free
  x <- sample(pool, n); free <<- setdiff(free, x); x
}
years_after <- function(e, lo = .05) {           # a date between blood draw and the end of follow-up
  a <- assess[match(e, eid)]; span <- pmin(14, as.numeric(censor_date - a) / 365.25)
  a + round(runif(length(e), lo, span) * 365.25)
}
years_before <- function(e) assess[match(e, eid)] - round(runif(length(e), .1, 10) * 365.25)

groups <- list(); dx <- list()
for (o in names(cfg)) {
  male <- o == "prostate"
  inc <- take(n_incident, male); prev <- take(n_prevalent, male); d1only <- take(n_dis1_only, male)
  prog <- inc[seq_len(n_progress)]
  inc_dx <- years_after(inc)
  # Precursor diagnosis for progressors: before the cancer, at any time from 10 years before blood draw
  prog_d1 <- pmin(inc_dx[seq_len(n_progress)] - 30, years_before(prog) + round(runif(n_progress, 0, 10) * 365.25))
  d1only_dx <- ifelse(runif(n_dis1_only) < .5, years_before(d1only), years_after(d1only))
  dx[[cfg[[o]]$dis2]] <- data.frame(eid = c(inc, prev), diagnosis_time = c(inc_dx, years_before(prev)))
  dx[[cfg[[o]]$dis1]] <- data.frame(eid = c(prog, d1only), diagnosis_time = as.Date(c(prog_d1, d1only_dx)))
  groups[[o]] <- list(inc = inc, inc_dx = inc_dx)
}
other_canc <- take(n_other_canc)
hc <- take(n_hc)

# Deaths (5%) and loss to follow-up (1%), never before a participant's cancer diagnosis
case_dx <- do.call(rbind, dx[sapply(cfg, `[[`, "dis2")])
for (i in sample(n_people, .05 * n_people)) {
  last <- max(c(assess[i], case_dx$diagnosis_time[case_dx$eid == eid[i]]))
  if (last + 60 < censor_date) clinic$date_of_death_f40000_0_0[i] <- last + sample(60:as.numeric(censor_date - last), 1)
}
for (i in sample(setdiff(seq_len(n_people), which(!is.na(clinic$date_of_death_f40000_0_0))), .01 * n_people)) {
  last <- max(c(assess[i], case_dx$diagnosis_time[case_dx$eid == eid[i]]))
  if (last + 60 < censor_date) clinic$date_lost_to_followup_f191_0_0[i] <- last + sample(60:as.numeric(censor_date - last), 1)
}

# ---- Proteins (NPX-like values) ----
sig_genes <- unique(unlist(lapply(cfg, function(x) c(x$signature, x$known_marker))))
# Current symbols for two renamed genes (the model maps old to new); MALAT1 is non-coding, so it has no assay
renamed <- c(ACPP = "ACP3", H2AFZ = "H2AZ1")
assays <- unique(c(ifelse(sig_genes %in% names(renamed), renamed[sig_genes], sig_genes)))
assays <- setdiff(assays, "MALAT1")
assays <- c(assays, sprintf("DUMMY%03d", seq_len(n_dummy_prot)), "IL12A_IL12B", "LGALS7_LGALS7B")
prot <- matrix(rnorm(n_people * length(assays)), n_people, dimnames = list(eid, assays))
prot <- sweep(prot, 2, runif(length(assays), .4, 1.2), `*`)
prot <- sweep(prot, 2, rnorm(length(assays), 0, 1), `+`)

# Planted signal: per organ, each signature protein shifts up or down, and the known marker up, in that organ's
# incident cases
for (o in names(cfg)) {
  g <- c(cfg[[o]]$signature, cfg[[o]]$known_marker)
  g <- intersect(ifelse(g %in% names(renamed), renamed[g], g), assays)
  inc <- groups[[o]]$inc
  yrs <- as.numeric(groups[[o]]$inc_dx - assess[match(inc, eid)]) / 365.25
  w <- shift_at_dx * exp(-yrs / decay_years)
  for (p in g) {
    s <- (if (p %in% cfg[[o]]$known_marker) 1 else sample(c(-1, 1), 1)) * sd(prot[, p])
    prot[inc, p] <- prot[inc, p] + s * w * runif(length(inc), .5, 1.5)
  }
}
# Missing values: about 3% everywhere; two proteins above the model's 25% filter
prot[runif(length(prot)) < .03] <- NA
for (p in sprintf("DUMMY%03d", n_dummy_prot - 0:1)) prot[runif(n_people) < .4, p] <- NA

# ---- Write the files the model reads ----
pt <- data.frame(PID = eid, round(prot, 4), check.names = FALSE)
write.table(pt, file.path(out, "Olink_proteomics_data_2ndPhase_transposed_decoded2UNIportID.txt"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")
saveRDS(clinic, file.path(out, "ukb_clinic.rda"))
for (d in names(dx)) saveRDS(dx[[d]], file.path(out, paste0("ukb_", d, "_subset_clinic.rda")))
saveRDS(data.frame(eid = hc), file.path(out, "ukb_healthy_controls.rda"))
all_canc <- rbind(case_dx, data.frame(eid = other_canc, diagnosis_time = years_after(other_canc)))
saveRDS(all_canc, file.path(out, "ukb_AllCancers_subset_clinic.rda"))

cat(sprintf("Wrote %d files to %s\n  %d participants, %d proteins (%d signature or marker assays)\n",
            length(list.files(out)), out, n_people, length(assays), length(assays) - n_dummy_prot - 2))
for (o in names(cfg)) cat(sprintf("  %-10s %d incident cancers (%d after the precursor), %d prevalent, %d precursor only\n",
                                  o, n_incident, n_progress, n_prevalent, n_dis1_only))
cat(sprintf("  %d other cancers, %d healthy controls\n", n_other_canc, n_hc))
