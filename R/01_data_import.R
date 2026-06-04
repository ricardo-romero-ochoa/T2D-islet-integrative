source("R/_shared.R")
load_or_install(c("GEOquery", "Biobase"), bioc = TRUE)

choose_best_eset <- function(gset_list) {
  if (length(gset_list) == 1) return(gset_list[[1]])
  sample_counts <- purrr::map_int(gset_list, ~ ncol(Biobase::exprs(.x)))
  gset_list[[which.max(sample_counts)]]
}

import_gse164416_counts <- function(accession, outdir, pheno) {
  message2("{accession}: downloading supplementary HTSeq count matrix.")
  GEOquery::getGEOSuppFiles(
    accession,
    makeDirectory = FALSE,
    baseDir = outdir,
    fetch_files = TRUE,
    filter_regex = "htseq_counts"
  )

  count_path <- list.files(
    outdir,
    pattern = "htseq_counts.*\\.(txt|tsv)(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(count_path) == 0) {
    stop(glue("{accession}: could not locate supplementary HTSeq counts file after download."), call. = FALSE)
  }
  count_path <- count_path[[1]]

  tbl <- readr::read_tsv(count_path, show_col_types = FALSE, progress = FALSE)
  if (ncol(tbl) < 2) {
    stop(glue("{accession}: supplementary HTSeq count table has fewer than two columns."), call. = FALSE)
  }

  sample_candidates <- unique(c(
    rownames(pheno),
    if ("geo_accession" %in% colnames(pheno)) pheno$geo_accession else NULL,
    if ("title" %in% colnames(pheno)) pheno$title else NULL
  ))
  sample_candidates <- sample_candidates[!is.na(sample_candidates) & nzchar(sample_candidates)]

  sample_cols <- intersect(colnames(tbl), sample_candidates)
  if (length(sample_cols) < 10) {
    numeric_like <- vapply(tbl, function(x) {
      if (is.numeric(x) || is.integer(x)) return(TRUE)
      vals <- suppressWarnings(readr::parse_number(as.character(x), na = c("", "NA", "N/A", "null", "NULL")))
      mean(is.finite(vals)) > 0.9
    }, logical(1))
    sample_cols <- setdiff(names(tbl)[numeric_like], names(tbl)[1])
  }
  if (length(sample_cols) < 10) {
    stop(glue("{accession}: could not identify sample columns in supplementary HTSeq counts table."), call. = FALSE)
  }

  non_sample_cols <- setdiff(colnames(tbl), sample_cols)
  feature_col <- non_sample_cols[[1]]

  expr <- as.data.frame(tbl[, sample_cols, drop = FALSE], check.names = FALSE)
  rownames(expr) <- as.character(tbl[[feature_col]])

  feat <- tibble::tibble(feature_id = rownames(expr))
  if (length(non_sample_cols) > 1) {
    extra_col <- non_sample_cols[[2]]
    extra_vals <- as.character(tbl[[extra_col]])
    if (mean(!is.na(extra_vals) & nzchar(extra_vals)) > 0.5) {
      feat$gene_symbol <- extra_vals
    }
  }

  list(
    expr = expr,
    feat = feat,
    count_path = count_path,
    import_source = "supplementary_htseq_counts"
  )
}

download_one_geo <- function(accession) {
  message2("Importing {accession} ...")
  outdir <- file.path("data/raw", accession)
  dir_create_safe(outdir)

  gset_list <- GEOquery::getGEO(accession, GSEMatrix = TRUE, AnnotGPL = TRUE)
  eset <- choose_best_eset(gset_list)

  pheno <- Biobase::pData(eset) %>% as.data.frame()
  feat <- Biobase::fData(eset) %>% as.data.frame()
  import_source <- "series_matrix"

  if (identical(accession, "GSE164416")) {
    supp <- import_gse164416_counts(accession, outdir, pheno)
    expr <- supp$expr
    feat <- supp$feat
    import_source <- supp$import_source
  } else {
    expr <- Biobase::exprs(eset) %>% as.data.frame()
    if (ncol(expr) == length(Biobase::sampleNames(eset))) {
      colnames(expr) <- Biobase::sampleNames(eset)
    }
    if ((is.null(rownames(expr)) || all(rownames(expr) %in% as.character(seq_len(nrow(expr))))) && nrow(feat) == nrow(expr)) {
      fid_col <- intersect(c("ID", "ID_REF", "ProbeID", "probe_id", "probeset_id"), colnames(feat))
      if (length(fid_col) > 0) rownames(expr) <- as.character(feat[[fid_col[[1]]]])
    }
  }

  save_rds_safe(eset, file.path(outdir, paste0(accession, "_eset.rds")))
  save_rds_safe(expr, file.path(outdir, paste0(accession, "_expression_raw.rds")))
  save_rds_safe(pheno, file.path(outdir, paste0(accession, "_pheno_raw.rds")))
  save_rds_safe(feat, file.path(outdir, paste0(accession, "_feature_raw.rds")))

  write_csv_safe(
    tibble::tibble(sample_id = rownames(pheno) %||% Biobase::sampleNames(eset)),
    file.path(outdir, paste0(accession, "_sample_ids.csv"))
  )

  summary_tbl <- tibble::tibble(
    accession = accession,
    n_features = nrow(expr),
    n_samples = ncol(expr),
    platform = Biobase::annotation(eset) %||% NA_character_,
    title = Biobase::experimentData(eset)@title %||% NA_character_,
    is_count_like = is_count_matrix(as.matrix(expr)),
    import_source = import_source
  )
  write_csv_safe(summary_tbl, file.path(outdir, paste0(accession, "_import_summary.csv")))
}

for (acc in ACCESSIONS) download_one_geo(acc)
message2("01_data_import complete.")
