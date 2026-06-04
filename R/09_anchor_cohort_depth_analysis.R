source("R/_shared.R")
load_or_install(c("limma"), bioc = TRUE)
load_or_install(c("broom"), bioc = FALSE)

acc <- "GSE164416"
meta <- safe_read_rds(file.path("data/processed/metadata", paste0(acc, "_metadata_analysis.rds")))
scores <- readr::read_csv(file.path("data/processed/modules", paste0(acc, "_module_scores.csv")), show_col_types = FALSE)

sec_meta <- meta %>%
  dplyr::filter(include_secondary_analysis, disease_group %in% c("ND", "IGT", "T2D", "T3cD"))
sec_scores <- scores %>% dplyr::filter(sample_id %in% sec_meta$sample_id)

preferred_meta <- sec_meta %>% dplyr::filter(include_anchor_preferred_subset %in% TRUE)
preferred_scores <- sec_scores %>% dplyr::filter(sample_id %in% preferred_meta$sample_id)

run_trend <- function(df, label) {
  trend_df <- df %>%
    dplyr::filter(disease_group %in% c("ND", "IGT", "T2D")) %>%
    dplyr::mutate(stage_num = dplyr::recode(disease_group, ND = 0L, IGT = 1L, T2D = 2L) %>% as.integer())

  purrr::map_dfr(c("StressInflammation", "IdentitySecretion", "IsletDysfunctionScore"), function(var) {
    fit <- lm(reformulate("stage_num", response = var), data = trend_df)
    broom::tidy(fit) %>%
      dplyr::filter(term == "stage_num") %>%
      dplyr::mutate(variable = var, cohort_slice = label)
  })
}

trend_results <- dplyr::bind_rows(
  run_trend(sec_scores, "full_secondary"),
  run_trend(preferred_scores, "ins_high_subset")
)

write_csv_safe(sec_scores, "data/processed/modules/GSE164416_secondary_state_module_scores.csv")
write_csv_safe(preferred_scores, "data/processed/modules/GSE164416_secondary_state_module_scores_ins_high_subset.csv")
write_csv_safe(trend_results, "results/tables/Table_6_anchor_cohort_trend_tests.csv")
message2("09_anchor_cohort_depth_analysis complete.")
