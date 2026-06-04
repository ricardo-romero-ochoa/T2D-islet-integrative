stages <- c(
  "01_data_import.R",
  "02_metadata_curation.R",
  "03_preprocessing_gene_mapping.R",
  "04_within_study_differential_expression.R",
  "05_cross_study_meta_analysis.R",
  "06_pathway_scoring_consensus.R",
  "07_module_construction.R",
  "08_beta_cell_specificity_analysis.R",
  "09_anchor_cohort_depth_analysis.R",
  "11_sensitivity_analyses.R",
  "12_refined_module_scores.R",
  "13_make_figures_2_to_6_refined.R"
)

if (!dir.exists("logs")) dir.create("logs", recursive = TRUE, showWarnings = FALSE)
log_file <- file.path("logs", sprintf("first_pass_run_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))

log_msg <- function(...) {
  msg <- sprintf(...)
  cat(msg, "
")
  cat(msg, "
", file = log_file, append = TRUE)
}

for (stage in stages) {
  log_msg("[START] %s", stage)
  ok <- tryCatch({
    source(file.path("R", stage), local = new.env(parent = globalenv()))
    TRUE
  }, error = function(e) {
    log_msg("[ERROR] %s :: %s", stage, conditionMessage(e))
    FALSE
  })
  if (!ok) stop(sprintf("Pipeline stopped at %s; see %s", stage, log_file), call. = FALSE)
  log_msg("[DONE] %s", stage)
}

log_msg("Full manuscript-results pipeline run completed successfully.")
