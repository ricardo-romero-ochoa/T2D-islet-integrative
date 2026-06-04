source("R/_shared.R")

modules_path <- file.path("data", "external", "refined_core_modules.csv")
if (!file.exists(modules_path)) stop("Missing refined module file: data/external/refined_core_modules.csv")

refined_modules_tbl <- readr::read_csv(modules_path, show_col_types = FALSE) %>%
  dplyr::mutate(
    module = as.character(module),
    gene = clean_gene_symbol(gene)
  ) %>%
  dplyr::filter(!is.na(gene), gene != "")

module_list <- split(refined_modules_tbl$gene, refined_modules_tbl$module)
module_list <- lapply(module_list, unique)

score_one_study_refined <- function(accession) {
  expr <- safe_read_rds(file.path("data/processed/expression", paste0(accession, "_expression_genelevel.rds")))
  meta <- safe_read_rds(file.path("data/processed/metadata", paste0(accession, "_metadata_analysis.rds")))
  common <- intersect(colnames(expr), meta$sample_id)
  expr <- expr[, common, drop = FALSE]
  meta <- meta %>% dplyr::filter(sample_id %in% common)
  meta <- meta[match(common, meta$sample_id), , drop = FALSE]

  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  rownames(expr) <- clean_gene_symbol(rownames(expr))
  expr <- expr[!is.na(rownames(expr)) & rownames(expr) != "", , drop = FALSE]
  expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]

  out <- tibble::tibble(
    sample_id = common,
    study = accession,
    ImmuneStress = compute_module_score(expr, module_list[["ImmuneStress"]], zscore = TRUE),
    BetaCellIdentitySecretion = compute_module_score(expr, module_list[["BetaCellIdentitySecretion"]], zscore = TRUE)
  ) %>%
    dplyr::mutate(IsletDysfunctionScore = ImmuneStress - BetaCellIdentitySecretion) %>%
    dplyr::left_join(
      meta %>% dplyr::select(dplyr::any_of(c(
        "sample_id", "sample_label", "disease_group", "include_primary_analysis",
        "tissue_class", "cell_context", "sex", "age", "bmi", "hba1c", "in_ins_filtered_data_subset"
      ))),
      by = "sample_id"
    )

  write_csv_safe(out, file.path("data/processed/modules_refined", paste0(accession, "_refined_module_scores.csv")))
  out
}

all_scores <- purrr::map_dfr(ACCESSIONS, score_one_study_refined)
write_csv_safe(all_scores, file.path("data/processed/modules_refined", "refined_module_scores_all_studies.csv"))

summary_by_group <- all_scores %>%
  dplyr::group_by(study, disease_group) %>%
  dplyr::summarise(
    n = dplyr::n(),
    ImmuneStress_mean = mean(ImmuneStress, na.rm = TRUE),
    ImmuneStress_sd = stats::sd(ImmuneStress, na.rm = TRUE),
    BetaCellIdentitySecretion_mean = mean(BetaCellIdentitySecretion, na.rm = TRUE),
    BetaCellIdentitySecretion_sd = stats::sd(BetaCellIdentitySecretion, na.rm = TRUE),
    IsletDysfunctionScore_mean = mean(IsletDysfunctionScore, na.rm = TRUE),
    IsletDysfunctionScore_sd = stats::sd(IsletDysfunctionScore, na.rm = TRUE),
    .groups = "drop"
  )
write_csv_safe(summary_by_group, file.path("results/tables", "Table_refined_module_summary_by_group.csv"))

cohens_d <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]
  g <- droplevels(as.factor(g[ok]))
  if (length(levels(g)) != 2) return(NA_real_)
  x1 <- x[g == levels(g)[1]]
  x2 <- x[g == levels(g)[2]]
  if (length(x1) < 2 || length(x2) < 2) return(NA_real_)
  s1 <- stats::sd(x1); s2 <- stats::sd(x2)
  sp <- sqrt(((length(x1)-1)*s1^2 + (length(x2)-1)*s2^2)/(length(x1)+length(x2)-2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (mean(x2) - mean(x1)) / sp
}

primary_nd_t2d <- all_scores %>%
  dplyr::mutate(
    include_primary_analysis = dplyr::case_when(
      is.logical(include_primary_analysis) ~ include_primary_analysis,
      is.character(include_primary_analysis) ~ tolower(include_primary_analysis) %in% c("true", "t", "1", "yes"),
      is.numeric(include_primary_analysis) ~ include_primary_analysis != 0,
      TRUE ~ FALSE
    ),
    disease_group = as.character(disease_group)
  ) %>%
  dplyr::filter(dplyr::coalesce(include_primary_analysis, FALSE), disease_group %in% c("ND", "T2D"))

write_csv_safe(primary_nd_t2d, file.path("results/tables", "DEBUG_refined_primary_nd_t2d_rows.csv"))

module_names <- c("ImmuneStress", "BetaCellIdentitySecretion", "IsletDysfunctionScore")

effect_sizes <- purrr::map_dfr(unique(primary_nd_t2d$study), function(acc) {
  dat <- primary_nd_t2d %>% dplyr::filter(study == acc)
  if (nrow(dat) == 0 || length(unique(dat$disease_group)) < 2) return(NULL)

  purrr::map_dfr(module_names, function(mn) {
    x <- dat[[mn]]
    g <- factor(dat$disease_group, levels = c("ND", "T2D"))
    tibble::tibble(
      study = acc,
      module = mn,
      n_ND = sum(g == "ND", na.rm = TRUE),
      n_T2D = sum(g == "T2D", na.rm = TRUE),
      mean_ND = mean(x[g == "ND"], na.rm = TRUE),
      mean_T2D = mean(x[g == "T2D"], na.rm = TRUE),
      cohens_d = cohens_d(x, g),
      ttest_pvalue = tryCatch(stats::t.test(x ~ g)$p.value, error = function(e) NA_real_),
      wilcox_pvalue = tryCatch(stats::wilcox.test(x ~ g)$p.value, error = function(e) NA_real_)
    )
  })
})

write_csv_safe(effect_sizes, file.path("results/tables", "Table_refined_module_effect_sizes.csv"))
message2("12_refined_module_scores complete.")
