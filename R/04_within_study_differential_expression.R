source("R/_shared.R")
load_or_install(c("limma", "edgeR", "DESeq2"), bioc = TRUE)

run_de_one_study <- function(accession) {
  expr_norm <- safe_read_rds(file.path("data/processed/expression", paste0(accession, "_expression_genelevel.rds")))
  count_path <- file.path("data/processed/expression", paste0(accession, "_expression_genelevel_counts.rds"))
  expr_counts <- if (file.exists(count_path)) safe_read_rds(count_path) else NULL

  meta <- safe_read_rds(file.path("data/processed/metadata", paste0(accession, "_metadata_analysis.rds"))) %>%
    dplyr::filter(include_primary_analysis, disease_group %in% c("ND", "T2D"))

  sample_keep <- intersect(colnames(expr_norm), meta$sample_id)
  expr_norm <- expr_norm[, sample_keep, drop = FALSE]
  if (!is.null(expr_counts)) expr_counts <- expr_counts[, sample_keep, drop = FALSE]

  meta <- meta %>% dplyr::filter(sample_id %in% sample_keep)
  meta <- meta[match(sample_keep, meta$sample_id), , drop = FALSE]

  if (ncol(expr_norm) < 4 || dplyr::n_distinct(meta$disease_group) < 2) {
    warning(glue("{accession}: insufficient samples/groups for DE; skipping."))
    return(NULL)
  }

  group <- factor(meta$disease_group, levels = c("ND", "T2D"))
  design <- model.matrix(~ group)

  if (!is.null(expr_counts)) {
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(expr_counts), colData = meta, design = ~ disease_group)
    keep <- rowSums(DESeq2::counts(dds) >= 5) >= max(2, floor(0.1 * ncol(dds)))
    dds <- dds[keep, ]
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    res <- as.data.frame(DESeq2::results(dds, contrast = c("disease_group", "T2D", "ND"))) %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::transmute(
        gene,
        logFC = log2FoldChange,
        SE = lfcSE,
        pvalue = pvalue,
        FDR = padj
      )
  } else {
    fit <- limma::lmFit(expr_norm, design)
    fit <- limma::eBayes(fit)
    tt <- limma::topTable(fit, coef = "groupT2D", number = Inf, sort.by = "none")
    se_vec <- if ("SE" %in% colnames(tt)) tt$SE else abs(tt$logFC / tt$t)
    res <- tt %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::transmute(
        gene,
        logFC = logFC,
        SE = se_vec,
        pvalue = P.Value,
        FDR = adj.P.Val
      )
  }

  res <- harmonize_deg_columns(res, accession) %>%
    dplyr::mutate(
      n_case = sum(meta$disease_group == "T2D"),
      n_control = sum(meta$disease_group == "ND")
    )

  write_csv_safe(res, file.path("data/processed/deg", paste0(accession, "_DE_full.csv")))
  write_csv_safe(dplyr::filter(res, FDR < 0.05), file.path("data/processed/deg", paste0(accession, "_DE_sig.csv")))
  save_rds_safe(res, file.path("data/processed/deg", paste0(accession, "_DE_full.rds")))
  invisible(res)
}

res_list <- purrr::map(ACCESSIONS, run_de_one_study)
names(res_list) <- ACCESSIONS
message2("04_within_study_differential_expression complete.")
