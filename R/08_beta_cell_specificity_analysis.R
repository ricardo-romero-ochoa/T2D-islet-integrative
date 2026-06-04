source("R/_shared.R")

acc <- "GSE20966"
scores <- readr::read_csv(file.path("data/processed/modules", paste0(acc, "_module_scores.csv")), show_col_types = FALSE)
meta_genes <- readr::read_csv("data/processed/meta/meta_high_confidence_genes.csv", show_col_types = FALSE)
deg <- readr::read_csv(file.path("data/processed/deg", paste0(acc, "_DE_full.csv")), show_col_types = FALSE)

preservation <- meta_genes %>%
  dplyr::select(gene, meta_logFC, meta_FDR) %>%
  dplyr::left_join(deg %>% dplyr::select(gene, logFC, FDR), by = "gene") %>%
  dplyr::mutate(
    concordant_direction = sign(meta_logFC) == sign(logFC),
    detected_in_beta = !is.na(logFC)
  )

write_csv_safe(preservation, "data/processed/modules/GSE20966_preservation_table.csv")

summary_tbl <- preservation %>%
  dplyr::summarise(
    n_meta_genes = dplyr::n(),
    n_detected_in_beta = sum(detected_in_beta, na.rm = TRUE),
    n_concordant = sum(concordant_direction, na.rm = TRUE),
    prop_concordant = mean(concordant_direction, na.rm = TRUE)
  )

write_csv_safe(summary_tbl, "results/tables/Table_5_beta_cell_specificity_summary.csv")
message2("08_beta_cell_specificity_analysis complete.")
