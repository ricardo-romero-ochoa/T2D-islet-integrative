source("R/_shared.R")

meta_all <- readr::read_csv("data/processed/meta/meta_all_genes.csv", show_col_types = FALSE)
meta_loo <- readr::read_csv("data/processed/meta/meta_leave_one_out_summary.csv", show_col_types = FALSE)
modules <- readr::read_csv("data/processed/modules/final_modules_long.csv", show_col_types = FALSE)

sens1 <- meta_loo %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    mean_abs_logFC = mean(abs(meta_logFC), na.rm = TRUE),
    n_sig_retained = sum(meta_FDR < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(n_sig_retained), desc(mean_abs_logFC))

write_csv_safe(sens1, "results/supplementary/Sensitivity_leave_one_out_gene_stability.csv")

sens2 <- modules %>%
  dplyr::count(module, name = "n_genes") %>%
  dplyr::mutate(note = "Re-run 03 + 04 + 05 + 07 with alternative collapse rule if needed.")
write_csv_safe(sens2, "results/supplementary/Sensitivity_module_size_summary.csv")

message2("11_sensitivity_analyses complete.")
