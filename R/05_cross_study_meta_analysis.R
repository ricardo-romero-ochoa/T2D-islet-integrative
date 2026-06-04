source("R/_shared.R")
load_or_install(c("metafor"), bioc = TRUE)

deg_files <- list.files("data/processed/deg", pattern = "_DE_full\\.csv$", full.names = TRUE)
deg_list <- purrr::map(deg_files, readr::read_csv, show_col_types = FALSE)
names(deg_list) <- basename(deg_files) %>% stringr::str_remove("_DE_full\\.csv$")
all_deg <- dplyr::bind_rows(deg_list)

meta_one_gene <- function(df) {
  df <- df %>% dplyr::filter(is.finite(logFC), is.finite(SE), SE > 0)
  if (nrow(df) < 2) return(NULL)
  fit <- tryCatch(
    metafor::rma.uni(yi = df$logFC, sei = df$SE, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  tibble::tibble(
    gene = df$gene[[1]],
    k = nrow(df),
    meta_logFC = as.numeric(fit$b),
    meta_SE = fit$se,
    meta_z = fit$zval,
    meta_pvalue = fit$pval,
    tau2 = fit$tau2,
    I2 = fit$I2,
    QEp = fit$QEp,
    direction_consistency = mean(sign(df$logFC) == sign(mean(df$logFC))),
    studies = paste(df$study, collapse = ";")
  )
}

meta_one_gene_loo <- function(df) {
  df <- df %>% dplyr::filter(is.finite(logFC), is.finite(SE), SE > 0)
  if (nrow(df) < 2) return(NULL)
  fit <- tryCatch(
    metafor::rma.uni(yi = df$logFC, sei = df$SE, method = "DL"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  tibble::tibble(
    gene = df$gene[[1]],
    k = nrow(df),
    meta_logFC = as.numeric(fit$b),
    meta_SE = fit$se,
    meta_z = fit$zval,
    meta_pvalue = fit$pval,
    tau2 = fit$tau2,
    I2 = fit$I2,
    QEp = fit$QEp,
    direction_consistency = mean(sign(df$logFC) == sign(mean(df$logFC))),
    studies = paste(df$study, collapse = ";")
  )
}

meta_tbl <- all_deg %>%
  split(.$gene) %>%
  purrr::map_dfr(meta_one_gene) %>%
  dplyr::mutate(meta_FDR = p.adjust(meta_pvalue, method = "BH")) %>%
  dplyr::arrange(meta_FDR, desc(abs(meta_logFC)))

high_conf <- meta_tbl %>%
  dplyr::filter(
    meta_FDR < 0.05,
    k >= 3,
    direction_consistency >= 0.75,
    is.na(I2) | I2 < 70
  )

write_csv_safe(meta_tbl, "data/processed/meta/meta_all_genes.csv")
write_csv_safe(high_conf, "data/processed/meta/meta_high_confidence_genes.csv")
write_csv_safe(dplyr::slice_head(high_conf, n = 100), "results/tables/Table_2_meta_top_genes.csv")

study_ids <- unique(all_deg$study)
gene_split <- split(all_deg, all_deg$gene)
message2("Starting leave-one-out meta-analysis for {length(study_ids)} studies and {length(gene_split)} genes.")

loo_list <- vector("list", length(study_ids))
names(loo_list) <- study_ids

for (i in seq_along(study_ids)) {
  study_drop <- study_ids[[i]]
  message2("LOO {i}/{length(study_ids)}: leaving out {study_drop}")

  res_i <- purrr::map_dfr(gene_split, function(df_gene) {
    df_sub <- df_gene %>% dplyr::filter(study != study_drop)
    out <- meta_one_gene_loo(df_sub)
    if (is.null(out)) return(NULL)
    out
  })

  if (nrow(res_i) > 0) {
    res_i <- res_i %>%
      dplyr::mutate(
        meta_FDR = p.adjust(meta_pvalue, method = "BH"),
        left_out = study_drop
      )
  }

  loo_list[[i]] <- res_i
  write_csv_safe(res_i, file.path("data/processed/meta", paste0("meta_leave_one_out_", study_drop, ".csv")))
  message2("LOO {i}/{length(study_ids)} complete: {study_drop}, rows={nrow(res_i)}")
}

loo_res <- dplyr::bind_rows(loo_list)
write_csv_safe(loo_res, "data/processed/meta/meta_leave_one_out_summary.csv")
message2("05_cross_study_meta_analysis complete.")
