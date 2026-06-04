suppressPackageStartupMessages({
  core_pkgs <- c("dplyr", "tidyr", "readr", "stringr", "purrr", "tibble", "ggplot2", "glue", "janitor", "data.table")
  for (p in core_pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
    library(p, character.only = TRUE)
  }
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(BiocManager, quietly = TRUE, warn.conflicts = FALSE))
})

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

message2 <- function(..., .envir = parent.frame()) cat(glue::glue(..., .envir = .envir), "\n")
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

load_or_install <- function(pkgs, bioc = FALSE) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (bioc) {
        BiocManager::install(p, ask = FALSE, update = FALSE)
      } else {
        install.packages(p, repos = "https://cloud.r-project.org")
      }
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}

dir_create_safe <- function(...) {
  paths <- c(...)
  for (p in paths) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  invisible(paths)
}

dir_create_safe(
  "data/raw", "data/external",
  "data/processed/expression", "data/processed/metadata", "data/processed/deg",
  "data/processed/meta", "data/processed/pathways", "data/processed/modules", "data/processed/modules_refined",
  "results/figures", "results/tables", "results/supplementary", "logs", "manuscript"
)

ACCESSIONS <- c("GSE25724", "GSE20966", "GSE38642", "GSE164416")

ACCESSION_INFO <- tibble::tribble(
  ~study,      ~platform_guess, ~cell_context,         ~primary_role,
  "GSE25724",  "microarray",    "islet_bulk",          "supporting_discovery",
  "GSE20966",  "microarray",    "beta_cell_enriched",  "specificity",
  "GSE38642",  "microarray",    "islet_bulk",          "major_legacy_replication",
  "GSE164416", "rnaseq",        "lcm_islet",           "anchor_depth"
)

safe_read_rds <- function(path) {
  if (!file.exists(path)) stop(glue::glue("Missing file: {path}"), call. = FALSE)
  readRDS(path)
}

write_csv_safe <- function(df, path) {
  dir_create_safe(dirname(path))
  readr::write_csv(df, path, na = "")
}

save_rds_safe <- function(obj, path) {
  dir_create_safe(dirname(path))
  saveRDS(obj, path)
}

clean_gene_symbol <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_trim() %>%
    dplyr::na_if("") %>%
    stringr::str_replace_all("///.*$", "") %>%
    stringr::str_replace_all(" // .*", "") %>%
    stringr::str_replace_all(";.*$", "") %>%
    stringr::str_replace_all("\\s+", "") %>%
    toupper()
}

clean_numeric_value <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "--", "NA", "N/A", "na", "null")] <- NA_character_
  readr::parse_number(x, na = c("", "--", "NA", "N/A", "na", "null"))
}

standardize_sex <- function(x) {
  x <- stringr::str_to_lower(as.character(x))
  dplyr::case_when(
    stringr::str_detect(x, "^m(ale)?$") ~ "Male",
    stringr::str_detect(x, "^f(emale)?$") ~ "Female",
    TRUE ~ NA_character_
  )
}

parse_field_from_text <- function(x, field) {
  x <- as.character(x)
  out <- stringr::str_match(
    x,
    paste0("(?:^|\\|)\\s*", field, "\\s*:?\\s*([^|]+)")
  )[, 2]
  stringr::str_trim(out)
}

parse_bool_field <- function(x, field) {
  val <- parse_field_from_text(x, field)
  dplyr::case_when(
    stringr::str_to_upper(val) %in% c("TRUE", "T", "YES", "Y", "1") ~ TRUE,
    stringr::str_to_upper(val) %in% c("FALSE", "F", "NO", "N", "0") ~ FALSE,
    TRUE ~ NA
  )
}

collapse_characteristics <- function(df) {
  df %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) %>%
    tidyr::unite("joined_text", dplyr::everything(), sep = " | ", remove = FALSE, na.rm = TRUE)
}

extract_characteristics_table <- function(pheno) {
  pheno <- as.data.frame(pheno)
  out <- tibble::rownames_to_column(pheno, "sample_rowname")
  char_cols <- grep("^characteristics", colnames(out), ignore.case = TRUE, value = TRUE)
  keep <- unique(c(
    "sample_rowname", "geo_accession", "title", "source_name_ch1",
    "description", char_cols
  ))
  keep <- intersect(keep, colnames(out))
  out[, keep, drop = FALSE] %>% janitor::clean_names()
}

is_count_matrix <- function(mat) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  vals <- as.numeric(mat[seq_len(min(length(mat), 50000))])
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(FALSE)
  prop_integer_like <- mean(abs(vals - round(vals)) < 1e-8)
  prop_nonnegative <- mean(vals >= 0)
  prop_integer_like > 0.95 && prop_nonnegative > 0.99 && max(vals, na.rm = TRUE) > 50
}

guess_symbol_column <- function(df) {
  candidates <- c(
    "Gene Symbol", "Gene symbol", "GENE_SYMBOL", "Symbol", "SYMBOL",
    "gene_symbol", "Gene.Symbol", "GENE SYMBOL", "gene_assignment", "Gene",
    "gene_name", "external_gene_name", "gene", "ILMN_Gene", "ILMN_GENE",
    "Reporter Symbol", "Reporter.Symbol", "Gene Name", "GeneName"
  )
  hits <- intersect(candidates, colnames(df))
  if (length(hits) > 0) return(hits[[1]])

  symbol_like <- names(df)[vapply(df, function(col) {
    if (!is.character(col)) return(FALSE)
    tmp <- na.omit(unique(col))
    tmp <- tmp[nchar(tmp) > 0 & nchar(tmp) < 25]
    if (length(tmp) < 10) return(FALSE)
    mean(grepl("^[A-Za-z0-9.-]+$", tmp)) > 0.7
  }, logical(1))]
  if (length(symbol_like) > 0) symbol_like[[1]] else NULL
}

guess_feature_id_col <- function(df) {
  candidates <- c("ID", "ID_REF", "ProbeID", "probe_id", "probeset_id", "feature_id")
  hits <- intersect(candidates, colnames(df))
  if (length(hits) > 0) return(hits[[1]])
  NULL
}

coerce_expression_matrix <- function(expr, feature_df = NULL, sample_ids = NULL, accession = NA_character_) {
  expr_df <- as.data.frame(expr, check.names = FALSE, stringsAsFactors = FALSE)

  if (!is.null(sample_ids)) {
    keep <- intersect(colnames(expr_df), sample_ids)
    if (length(keep) >= 2) {
      other_cols <- setdiff(colnames(expr_df), keep)
      id_col <- intersect(other_cols, c("ID", "ID_REF", "ProbeID", "probe_id", "probeset_id"))
      if (length(id_col) > 0 && (is.null(rownames(expr_df)) || all(rownames(expr_df) %in% as.character(seq_len(nrow(expr_df)))))) {
        rownames(expr_df) <- as.character(expr_df[[id_col[[1]]]])
      }
      expr_df <- expr_df[, keep, drop = FALSE]
    }
  }

  if (is.null(rownames(expr_df)) || all(rownames(expr_df) %in% as.character(seq_len(nrow(expr_df))))) {
    id_col <- guess_feature_id_col(expr_df)
    if (!is.null(id_col)) {
      rownames(expr_df) <- as.character(expr_df[[id_col]])
      expr_df[[id_col]] <- NULL
    }
  }

  expr_df[] <- lapply(expr_df, function(x) {
    if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
    readr::parse_number(as.character(x), na = c("", "NA", "N/A", "null", "NULL"))
  })

  mat <- as.matrix(expr_df)
  storage.mode(mat) <- "numeric"

  rn <- rownames(mat)
  if ((is.null(rn) || anyNA(rn) || all(rn == "") || all(rn %in% as.character(seq_len(nrow(mat))))) && !is.null(feature_df) && nrow(feature_df) == nrow(mat)) {
    fid_col <- guess_feature_id_col(feature_df)
    if (!is.null(fid_col)) rownames(mat) <- as.character(feature_df[[fid_col]])
  }

  if (!is.null(sample_ids) && length(intersect(colnames(mat), sample_ids)) >= 2) {
    mat <- mat[, intersect(colnames(mat), sample_ids), drop = FALSE]
  }

  mat
}

collapse_probes_to_gene <- function(expr_df, feature_df, feature_id_col = NULL, rule = c("max_iqr", "max_mean")) {
  rule <- match.arg(rule)
  expr_df <- as.data.frame(expr_df)
  feature_df <- as.data.frame(feature_df)

  if (is.null(feature_id_col)) {
    feature_id_col <- intersect(c("ID", "ID_REF", "ProbeID", "probe_id", "probeset_id"), colnames(feature_df))
    if (length(feature_id_col) == 0) {
      feature_df <- tibble::rownames_to_column(feature_df, "feature_id")
      feature_id_col <- "feature_id"
    } else {
      feature_id_col <- feature_id_col[[1]]
    }
  }

  symbol_col <- guess_symbol_column(feature_df)
  if (is.null(symbol_col)) stop("Could not infer gene symbol column from feature annotations.")

  map <- feature_df %>%
    dplyr::transmute(
      feature_id = as.character(.data[[feature_id_col]]),
      gene = clean_gene_symbol(.data[[symbol_col]])
    ) %>%
    dplyr::filter(!is.na(gene), gene != "")

  expr_df <- tibble::rownames_to_column(expr_df, "feature_id")
  expr_df$feature_id <- as.character(expr_df$feature_id)
  merged <- dplyr::inner_join(expr_df, map, by = "feature_id")
  if (nrow(merged) == 0) stop("No rows left after feature-to-gene merge.")

  value_cols <- setdiff(colnames(merged), c("feature_id", "gene"))
  metrics <- merged %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      .mean_expr = mean(dplyr::c_across(dplyr::all_of(value_cols)), na.rm = TRUE),
      .iqr_expr = IQR(dplyr::c_across(dplyr::all_of(value_cols)), na.rm = TRUE)
    ) %>%
    dplyr::ungroup()

  chosen <- if (rule == "max_iqr") {
    metrics %>% dplyr::group_by(gene) %>% dplyr::slice_max(order_by = .iqr_expr, n = 1, with_ties = FALSE)
  } else {
    metrics %>% dplyr::group_by(gene) %>% dplyr::slice_max(order_by = .mean_expr, n = 1, with_ties = FALSE)
  }

  out <- chosen %>%
    dplyr::select(gene, dplyr::all_of(value_cols)) %>%
    as.data.frame()
  rownames(out) <- out$gene
  out$gene <- NULL
  as.matrix(out)
}

map_ensembl_to_symbol <- function(ids) {
  load_or_install(c("AnnotationDbi", "org.Hs.eg.db"), bioc = TRUE)
  cleaned <- stringr::str_replace(as.character(ids), "\\..*$", "")
  mapped <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = cleaned,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  clean_gene_symbol(unname(mapped))
}

matrix_to_gene_level <- function(mat, feature_df = NULL, collapse_rule = "max_iqr") {
  mat <- as.matrix(mat)
  if (!is.null(feature_df) && nrow(feature_df) > 0) {
    try_feature <- tryCatch(
      collapse_probes_to_gene(mat, feature_df, rule = collapse_rule),
      error = function(e) NULL
    )
    if (!is.null(try_feature)) return(try_feature)
  }

  rn <- rownames(mat)
  if (is.null(rn)) stop("Expression matrix lacks rownames; cannot derive gene-level matrix.")
  if (mean(grepl("^[0-9]+$", rn)) > 0.8) {
    stop("Expression matrix rownames are mostly numeric probe IDs; feature-to-gene mapping failed.")
  }

  if (mean(grepl("^ENSG[0-9]+", rn)) > 0.8) {
    genes <- map_ensembl_to_symbol(rn)
    keep <- !is.na(genes) & genes != ""
    mat <- mat[keep, , drop = FALSE]
    genes <- genes[keep]
    rownames(mat) <- genes
    return(rowsum(mat, group = rownames(mat)))
  }

  genes <- clean_gene_symbol(rn)
  keep <- !is.na(genes) & genes != ""
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- genes[keep]
  rowsum(mat, group = rownames(mat))
}

harmonize_deg_columns <- function(df, study) {
  req <- c("gene", "logFC", "SE", "pvalue")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) stop(glue::glue("Missing columns in DEG table for {study}: {paste(miss, collapse=', ')}"))
  df %>%
    dplyr::mutate(
      study = study,
      gene = clean_gene_symbol(gene),
      direction = dplyr::case_when(logFC > 0 ~ "up", logFC < 0 ~ "down", TRUE ~ "flat")
    ) %>%
    dplyr::filter(!is.na(gene), gene != "")
}

compute_module_score <- function(expr_mat, genes, zscore = TRUE) {
  genes <- intersect(clean_gene_symbol(genes), rownames(expr_mat))
  if (length(genes) < 2) return(rep(NA_real_, ncol(expr_mat)))
  x <- expr_mat[genes, , drop = FALSE]
  if (zscore) x <- t(scale(t(x)))
  colMeans(x, na.rm = TRUE)
}

validate_metadata <- function(df, accession) {
  problems <- c()
  if (anyDuplicated(df$sample_id) > 0) problems <- c(problems, "duplicate_sample_id")
  if (sum(!is.na(df$disease_group)) == 0) problems <- c(problems, "all_disease_groups_missing")
  if (accession %in% c("GSE25724", "GSE20966", "GSE38642") && !all(sort(unique(stats::na.omit(df$disease_group))) %in% c("ND", "T2D"))) {
    problems <- c(problems, "unexpected_groups_in_primary_case_control")
  }
  tibble::tibble(
    study = accession,
    n_samples = nrow(df),
    n_primary = sum(df$include_primary_analysis %in% TRUE, na.rm = TRUE),
    n_secondary = sum(df$include_secondary_analysis %in% TRUE, na.rm = TRUE),
    n_missing_group = sum(is.na(df$disease_group)),
    problems = if (length(problems) == 0) "none" else paste(problems, collapse = ";")
  )
}

save_plot_pdf <- function(plot_obj, path, width = 8, height = 6) {
  dir_create_safe(dirname(path))
  ggplot2::ggsave(filename = path, plot = plot_obj, width = width, height = height, device = cairo_pdf)
}
