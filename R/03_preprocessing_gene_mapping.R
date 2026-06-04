source("R/_shared.R")

preprocess_one_study <- function(accession, collapse_rule = "max_iqr") {
  expr <- safe_read_rds(file.path("data/raw", accession, paste0(accession, "_expression_raw.rds")))
  feat <- safe_read_rds(file.path("data/raw", accession, paste0(accession, "_feature_raw.rds")))
  meta <- readr::read_csv(
    file.path("data/processed/metadata", paste0(accession, "_metadata_harmonized.csv")),
    show_col_types = FALSE
  )

  mat <- coerce_expression_matrix(expr, feature_df = feat, sample_ids = meta$sample_id, accession = accession)
  if (ncol(mat) == nrow(meta) && !all(colnames(mat) %in% meta$sample_id)) {
    common_guess <- intersect(colnames(mat), meta$sample_id)
    if (length(common_guess) == 0 && ncol(mat) == length(meta$sample_id)) {
      colnames(mat) <- meta$sample_id
    }
  }

  finite_vals <- mat[is.finite(mat)]
  if (length(finite_vals) == 0) {
    stop(glue("{accession}: expression matrix has no finite numeric values after coercion. For GSE164416, rerun 01_data_import.R after switching to the supplementary HTSeq counts file."))
  }

  if (mean(grepl("^[0-9]+$", rownames(mat))) > 0.8) {
    message2("{accession}: raw matrix rownames are mostly numeric probe IDs; attempting feature-table collapse to gene level.")
  }

  if (is_count_matrix(mat)) {
    keep <- rowSums(mat >= 5) >= max(2, floor(0.1 * ncol(mat)))
    mat <- mat[keep, , drop = FALSE]
    gene_mat <- matrix_to_gene_level(mat, feature_df = feat, collapse_rule = collapse_rule)
    if (mean(grepl("^[0-9]+$", rownames(gene_mat))) > 0.8) {
      stop(glue::glue("{accession}: gene-level matrix still has mostly numeric rownames after collapse; check feature annotations."))
    }
    norm_mat <- log2(gene_mat + 1)
    save_rds_safe(gene_mat, file.path("data/processed/expression", paste0(accession, "_expression_genelevel_counts.rds")))
    analysis_type <- "rnaseq"
  } else {
    if (max(finite_vals, na.rm = TRUE) > 100 && min(finite_vals, na.rm = TRUE) >= 0) {
      message2("{accession}: values look unlogged; applying log2(x + 1).")
      mat <- log2(mat + 1)
    }
    gene_mat <- matrix_to_gene_level(mat, feature_df = feat, collapse_rule = collapse_rule)
    if (mean(grepl("^[0-9]+$", rownames(gene_mat))) > 0.8) {
      stop(glue::glue("{accession}: gene-level matrix still has mostly numeric rownames after collapse; check feature annotations."))
    }
    norm_mat <- gene_mat
    analysis_type <- "microarray"
  }

  shared_samples <- intersect(colnames(norm_mat), meta$sample_id)
  if (length(shared_samples) == 0) stop(glue("{accession}: no sample overlap between expression matrix and metadata."))

  norm_mat <- norm_mat[, shared_samples, drop = FALSE]
  meta2 <- meta %>% dplyr::filter(sample_id %in% shared_samples)
  meta2 <- meta2[match(shared_samples, meta2$sample_id), , drop = FALSE]

  if (file.exists(file.path("data/processed/expression", paste0(accession, "_expression_genelevel_counts.rds")))) {
    count_mat <- safe_read_rds(file.path("data/processed/expression", paste0(accession, "_expression_genelevel_counts.rds")))
    count_mat <- count_mat[, shared_samples, drop = FALSE]
    save_rds_safe(count_mat, file.path("data/processed/expression", paste0(accession, "_expression_genelevel_counts.rds")))
  }

  save_rds_safe(norm_mat, file.path("data/processed/expression", paste0(accession, "_expression_genelevel.rds")))
  save_rds_safe(meta2, file.path("data/processed/metadata", paste0(accession, "_metadata_analysis.rds")))

  qc <- tibble::tibble(
    study = accession,
    analysis_type = analysis_type,
    n_genes = nrow(norm_mat),
    n_samples = ncol(norm_mat),
    collapse_rule = collapse_rule,
    sample_overlap = length(shared_samples)
  )
  write_csv_safe(qc, file.path("results/supplementary", paste0(accession, "_preprocessing_qc.csv")))
}

for (acc in ACCESSIONS) preprocess_one_study(acc)
message2("03_preprocessing_gene_mapping complete.")
