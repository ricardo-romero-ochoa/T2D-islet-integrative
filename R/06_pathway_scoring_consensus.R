source("R/_shared.R")
load_or_install(c("GSVA", "msigdbr", "limma"), bioc = TRUE)

msig_h <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::mutate(gene_symbol = clean_gene_symbol(gene_symbol)) %>%
  dplyr::filter(!is.na(gene_symbol), gene_symbol != "")

hallmark_sets <- split(msig_h$gene_symbol, msig_h$gs_name)
hallmark_sets <- lapply(hallmark_sets, unique)

run_ssgsea_compat <- function(expr_mat, gene_sets) {
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"

  rn <- rownames(expr_mat)
  rn <- as.character(rn)
  rn <- sub("\\.[0-9]+$", "", rn)
  rn <- sub(" ///.*$", "", rn)
  rn <- trimws(rn)
  rn <- clean_gene_symbol(rn)
  rownames(expr_mat) <- rn

  keep_rn <- !is.na(rownames(expr_mat)) & rownames(expr_mat) != ""
  expr_mat <- expr_mat[keep_rn, , drop = FALSE]
  expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]

  keep_rows <- apply(expr_mat, 1, function(x) any(is.finite(x)))
  keep_cols <- apply(expr_mat, 2, function(x) any(is.finite(x)))
  expr_mat <- expr_mat[keep_rows, keep_cols, drop = FALSE]

  if (nrow(expr_mat) < 100) stop("Expression matrix has too few valid genes after filtering.")
  if (ncol(expr_mat) < 2) stop("Expression matrix has too few valid samples after filtering.")

  gene_sets2 <- lapply(gene_sets, function(gs) intersect(gs, rownames(expr_mat)))
  gene_sets2 <- gene_sets2[lengths(gene_sets2) >= 5]

  if (length(gene_sets2) == 0) {
    print(head(rownames(expr_mat), 30))
    stop("No gene sets retained after matching to expression matrix.")
  }

  if (utils::packageVersion("GSVA") >= "1.50.0") {
    ssgsea_par <- GSVA::ssgseaParam(exprData = expr_mat, geneSets = gene_sets2)
    scores <- GSVA::gsva(ssgsea_par, verbose = FALSE)
  } else {
    scores <- GSVA::gsva(
      expr = expr_mat,
      gset.idx.list = gene_sets2,
      method = "ssgsea",
      kcdf = "Gaussian",
      verbose = FALSE
    )
  }

  as.matrix(scores)
}

score_one_study <- function(accession) {
  expr <- safe_read_rds(file.path("data/processed/expression", paste0(accession, "_expression_genelevel.rds")))
  meta <- safe_read_rds(file.path("data/processed/metadata", paste0(accession, "_metadata_analysis.rds")))

  common <- intersect(colnames(expr), meta$sample_id)
  expr <- expr[, common, drop = FALSE]
  meta <- meta %>% dplyr::filter(sample_id %in% common)
  meta <- meta[match(common, meta$sample_id), , drop = FALSE]

  expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
  message2("{accession}: n_genes={nrow(expr)}, n_samples={ncol(expr)}")
  message2("{accession}: example gene IDs: {paste(head(rownames(expr), 10), collapse = ', ')}")

  scores <- run_ssgsea_compat(expr, hallmark_sets)
  save_rds_safe(scores, file.path("data/processed/pathways", paste0(accession, "_hallmark_scores.rds")))

  meta2 <- meta %>% dplyr::filter(include_primary_analysis, disease_group %in% c("ND", "T2D"))
  if (nrow(meta2) < 4 || length(unique(meta2$disease_group)) < 2) {
    stop(paste0(accession, ": insufficient ND/T2D samples for pathway differential analysis."))
  }

  scores2 <- scores[, meta2$sample_id, drop = FALSE]
  design <- model.matrix(~ factor(meta2$disease_group, levels = c("ND", "T2D")))
  fit <- limma::lmFit(scores2, design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("pathway") %>%
    dplyr::transmute(pathway, logFC = logFC, pvalue = P.Value, FDR = adj.P.Val, study = accession)
  write_csv_safe(tt, file.path("data/processed/pathways", paste0(accession, "_hallmark_differential.csv")))
  tt
}

res <- purrr::map_dfr(ACCESSIONS, function(acc) {
  message2("Scoring pathways: {acc}")
  score_one_study(acc)
})

consensus <- res %>%
  dplyr::group_by(pathway) %>%
  dplyr::summarise(
    k = dplyr::n(),
    mean_logFC = mean(logFC, na.rm = TRUE),
    direction_consistency = mean(sign(logFC) == sign(mean(logFC)), na.rm = TRUE),
    n_sig = sum(FDR < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(n_sig), desc(abs(mean_logFC)))

write_csv_safe(res, "data/processed/pathways/hallmark_all_studies_long.csv")
write_csv_safe(consensus, "data/processed/pathways/hallmark_consensus_summary.csv")
write_csv_safe(head(consensus, 50), "results/tables/Table_3_consensus_pathways.csv")

message2("06_pathway_scoring_consensus complete.")
