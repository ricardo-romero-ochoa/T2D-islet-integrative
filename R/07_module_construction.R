source("R/_shared.R")

meta_tbl <- readr::read_csv("data/processed/meta/meta_high_confidence_genes.csv", show_col_types = FALSE)
path_tbl <- readr::read_csv("data/processed/pathways/hallmark_consensus_summary.csv", show_col_types = FALSE)

# Simple starter logic:
# - Upregulated genes = stress/inflammatory candidate pool
# - Downregulated genes = identity/secretory candidate pool
# You should refine these after inspecting pathway membership and literature context.

up_genes <- meta_tbl %>%
  dplyr::filter(meta_logFC > 0) %>%
  dplyr::slice_max(order_by = abs(meta_logFC), n = 40, with_ties = FALSE) %>%
  dplyr::pull(gene)

down_genes <- meta_tbl %>%
  dplyr::filter(meta_logFC < 0) %>%
  dplyr::slice_max(order_by = abs(meta_logFC), n = 40, with_ties = FALSE) %>%
  dplyr::pull(gene)

module_tbl <- tibble(
  module = c(rep("StressInflammation", length(up_genes)), rep("IdentitySecretion", length(down_genes))),
  gene = c(up_genes, down_genes)
)

write_csv_safe(module_tbl, "data/processed/modules/final_modules_long.csv")
write_csv_safe(module_tbl %>% dplyr::count(module, name = "n_genes"), "results/tables/Table_4_final_modules.csv")

for (acc in ACCESSIONS) {
  expr <- safe_read_rds(file.path("data/processed/expression", paste0(acc, "_expression_genelevel.rds")))
  meta <- safe_read_rds(file.path("data/processed/metadata", paste0(acc, "_metadata_analysis.rds")))

  scores <- tibble(
    sample_id = colnames(expr),
    study = acc,
    StressInflammation = compute_module_score(expr, up_genes),
    IdentitySecretion = compute_module_score(expr, down_genes)
  ) %>%
    dplyr::mutate(IsletDysfunctionScore = StressInflammation - IdentitySecretion) %>%
    dplyr::left_join(meta %>% dplyr::select(sample_id, disease_group, cell_context), by = "sample_id")

  write_csv_safe(scores, file.path("data/processed/modules", paste0(acc, "_module_scores.csv")))
}

message2("07_module_construction complete. Refine module composition after inspecting consensus pathways.")
