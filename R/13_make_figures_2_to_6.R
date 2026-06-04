find_repo_root <- function(start = getwd()) {
  cur <- normalizePath(start, winslash = "/", mustWork = FALSE)
  while (TRUE) {
    if (file.exists(file.path(cur, "R", "_shared.R"))) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

repo_root <- find_repo_root()
setwd(repo_root)
source(file.path(repo_root, "R", "_shared.R"))


suppressPackageStartupMessages({
  pkgs <- c("readr", "dplyr", "tidyr", "ggplot2", "stringr", "purrr", "tibble", "patchwork", "forcats", "scales")
  to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(to_install) > 0) install.packages(to_install, repos = "https://cloud.r-project.org")
  lapply(pkgs, library, character.only = TRUE)
})

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------
ensure_dir <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) stop("None of these paths exist:\n", paste(paths, collapse = "\n"))
  hit[[1]]
}

resolve_repo_file <- function(rel_path = NULL, filename = NULL, search_roots = c(repo_root, getwd()), must_exist = TRUE) {
  candidates <- character()
  if (!is.null(rel_path)) {
    candidates <- c(candidates,
      file.path(repo_root, rel_path),
      file.path(getwd(), rel_path),
      rel_path
    )
  }
  candidates <- unique(candidates)
  hit <- candidates[file.exists(candidates)]
  if (length(hit) > 0) return(normalizePath(hit[[1]], winslash = "/", mustWork = FALSE))

  if (!is.null(filename) && nzchar(filename)) {
    found <- character()
    for (root in unique(search_roots)) {
      if (dir.exists(root)) {
        f <- list.files(root, recursive = TRUE, full.names = TRUE)
        f <- f[basename(f) == filename]
        if (length(f) > 0) found <- c(found, f)
      }
    }
    found <- unique(found[file.exists(found)])
    if (length(found) > 0) return(normalizePath(found[[1]], winslash = "/", mustWork = FALSE))
  }

  if (must_exist) {
    msg <- c()
    if (!is.null(rel_path)) msg <- c(msg, paste0("expected relative path: ", rel_path))
    if (!is.null(filename)) msg <- c(msg, paste0("searched recursively for filename: ", filename))
    stop("Could not locate required file. ", paste(msg, collapse = "; "))
  }
  NA_character_
}

safe_read_csv <- function(paths = NULL, filename = NULL) {
  path <- if (!is.null(paths)) first_existing(paths) else resolve_repo_file(filename = filename)
  readr::read_csv(path, show_col_types = FALSE)
}

safe_read_rds <- function(paths = NULL, filename = NULL) {
  path <- if (!is.null(paths)) first_existing(paths) else resolve_repo_file(filename = filename)
  readRDS(path)
}

save_plot <- function(plot_obj, filename, width, height, dpi = 320) {
  ensure_dir(dirname(filename))
  ggsave(filename, plot_obj, width = width, height = height, dpi = dpi, bg = "white")
}

standardize_gene <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("\\.[0-9]+$", "", x)
  x <- sub(" ///.*$", "", x)
  x <- sub(" //.*$", "", x)
  x <- toupper(x)
  x[x %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  x
}

cohens_d <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  nx <- length(x); ny <- length(y)
  if (nx < 2 || ny < 2) return(NA_real_)
  sx <- stats::sd(x); sy <- stats::sd(y)
  sp <- sqrt(((nx - 1) * sx^2 + (ny - 1) * sy^2) / (nx + ny - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (mean(y) - mean(x)) / sp
}

se_d <- function(d, n1, n2) {
  if (any(!is.finite(c(d, n1, n2))) || n1 < 2 || n2 < 2) return(NA_real_)
  sqrt((n1 + n2) / (n1 * n2) + (d^2) / (2 * (n1 + n2 - 2)))
}

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
study_order <- c("GSE25724", "GSE20966", "GSE38642", "GSE164416")
fig_dir <- "results/figures"
tbl_dir <- "results/tables"
ensure_dir(fig_dir)
ensure_dir(tbl_dir)

refined_modules_file_candidates <- c(
  file.path(repo_root, "data", "external", "refined_core_modules.csv"),
  file.path("data", "external", "refined_core_modules.csv"),
  "refined_core_modules.csv",
  file.path(tbl_dir, "refined_core_modules.csv"),
  "/mnt/data/refined_core_modules.csv"
)
refined_modules_file <- refined_modules_file_candidates[file.exists(refined_modules_file_candidates)]
refined_modules_file <- if (length(refined_modules_file) > 0) refined_modules_file[[1]] else NA_character_

meta_high_conf_file <- resolve_repo_file(rel_path = "data/processed/meta/meta_high_confidence_genes.csv", filename = "meta_high_confidence_genes.csv")

effect_sizes_file <- resolve_repo_file(rel_path = file.path(tbl_dir, "Table_refined_module_effect_sizes.csv"), filename = "Table_refined_module_effect_sizes.csv")

summary_group_file <- resolve_repo_file(rel_path = file.path(tbl_dir, "Table_refined_module_summary_by_group.csv"), filename = "Table_refined_module_summary_by_group.csv")

# ------------------------------------------------------------
# Load refined modules
# ------------------------------------------------------------
default_refined_modules <- tibble::tribble(
  ~module, ~gene,
  "ImmuneStress", "MICB",
  "ImmuneStress", "HLA-DRA",
  "ImmuneStress", "HLA-DPA1",
  "ImmuneStress", "IL1R2",
  "ImmuneStress", "IL1RL1",
  "ImmuneStress", "IDO1",
  "ImmuneStress", "SERPING1",
  "ImmuneStress", "FPR3",
  "ImmuneStress", "LTB4R",
  "ImmuneStress", "GBP2",
  "ImmuneStress", "TNFRSF10A",
  "ImmuneStress", "CFH",
  "ImmuneStress", "ADORA3",
  "ImmuneStress", "APOL1",
  "BetaCellIdentitySecretion", "RASGRP1",
  "BetaCellIdentitySecretion", "PPP1R1A",
  "BetaCellIdentitySecretion", "ENTPD3",
  "BetaCellIdentitySecretion", "ADCYAP1",
  "BetaCellIdentitySecretion", "FFAR4",
  "BetaCellIdentitySecretion", "TMED6",
  "BetaCellIdentitySecretion", "PLCXD3",
  "BetaCellIdentitySecretion", "PDE8B",
  "BetaCellIdentitySecretion", "CASR",
  "BetaCellIdentitySecretion", "PFKFB2",
  "BetaCellIdentitySecretion", "ACLY",
  "BetaCellIdentitySecretion", "TGFBR3",
  "BetaCellIdentitySecretion", "ASB9",
  "BetaCellIdentitySecretion", "PPM1K"
)

refined_modules <- if (is.na(refined_modules_file)) {
  message("refined_core_modules.csv not found; using embedded refined module definitions.")
  default_refined_modules
} else {
  readr::read_csv(refined_modules_file, show_col_types = FALSE)
}

refined_modules <- refined_modules %>%
  mutate(
    module = as.character(module),
    gene = standardize_gene(gene)
  ) %>%
  filter(!is.na(gene))

module_list <- split(refined_modules$gene, refined_modules$module)
module_list <- lapply(module_list, unique)

# ------------------------------------------------------------
# Recompute refined sample-level scores from expression matrices
# ------------------------------------------------------------
score_one_study <- function(accession) {
  expr <- safe_read_rds(filename = paste0(accession, "_expression_genelevel.rds"))
  meta <- safe_read_rds(filename = paste0(accession, "_metadata_analysis.rds"))

  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  rownames(expr) <- standardize_gene(rownames(expr))
  expr <- expr[!is.na(rownames(expr)) & rownames(expr) != "", , drop = FALSE]
  expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]

  common <- intersect(colnames(expr), meta$sample_id)
  expr <- expr[, common, drop = FALSE]
  meta <- meta %>% filter(sample_id %in% common)
  meta <- meta[match(common, meta$sample_id), , drop = FALSE]

  score_module <- function(genes) {
    genes2 <- intersect(genes, rownames(expr))
    if (length(genes2) == 0) return(rep(NA_real_, ncol(expr)))
    colMeans(expr[genes2, , drop = FALSE], na.rm = TRUE)
  }

  out <- tibble(
    sample_id = common,
    study = accession,
    ImmuneStress = score_module(module_list[["ImmuneStress"]]),
    BetaCellIdentitySecretion = score_module(module_list[["BetaCellIdentitySecretion"]]),
    IsletDysfunctionScore = score_module(module_list[["ImmuneStress"]]) - score_module(module_list[["BetaCellIdentitySecretion"]]),
    disease_group = as.character(meta$disease_group),
    cell_context = as.character(meta$cell_context)
  )
  out
}

sample_scores <- purrr::map_dfr(study_order, score_one_study)
readr::write_csv(sample_scores, file.path(tbl_dir, "DEBUG_refined_sample_scores_recomputed.csv"))

# ------------------------------------------------------------
# Figure 2
# ------------------------------------------------------------

# 2A. DEG counts
deg_summary <- purrr::map_dfr(study_order, function(acc) {
  x <- safe_read_csv(paths = c(
    file.path("data/processed/deg", paste0(acc, "_DE_full.csv")),
    file.path("data/processed/deg", paste0(acc, "_deg_full.csv"))
  ))
  tibble(
    study = acc,
    direction = c("Up", "Down"),
    n_genes = c(
      sum(x$FDR < 0.05 & x$logFC > 0, na.rm = TRUE),
      -sum(x$FDR < 0.05 & x$logFC < 0, na.rm = TRUE)
    )
  )
}) %>% mutate(study = factor(study, levels = rev(study_order)))

p2a <- ggplot(deg_summary, aes(x = study, y = n_genes, fill = direction)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  coord_flip() +
  labs(title = "A. Per-study differential expression overview", x = NULL, y = "Significant genes (FDR < 0.05)") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")

# 2B. Heatmap of representative high-confidence genes
meta_hc <- readr::read_csv(meta_high_conf_file, show_col_types = FALSE)
rep_genes <- meta_hc %>%
  arrange(meta_FDR, desc(abs(meta_logFC))) %>%
  slice_head(n = 24) %>%
  pull(gene) %>% standardize_gene() %>% unique()

gene_heat <- purrr::map_dfr(study_order, function(acc) {
  x <- safe_read_csv(c(
    file.path("data/processed/deg", paste0(acc, "_DE_full.csv")),
    file.path("data/processed/deg", paste0(acc, "_deg_full.csv")),
    file.path("/mnt/data/data/processed/deg", paste0(acc, "_DE_full.csv")),
    file.path("/mnt/data/data/processed/deg", paste0(acc, "_deg_full.csv"))
  )) %>%
    mutate(gene = standardize_gene(gene)) %>%
    filter(gene %in% rep_genes) %>%
    transmute(study = acc, gene, logFC)
  x
})

gene_order <- meta_hc %>% mutate(gene = standardize_gene(gene)) %>%
  filter(gene %in% rep_genes) %>%
  arrange(meta_FDR, desc(abs(meta_logFC))) %>% pull(gene) %>% unique()

gene_heat <- gene_heat %>%
  group_by(gene) %>%
  mutate(z = as.numeric(scale(logFC))) %>%
  ungroup() %>%
  mutate(
    study = factor(study, levels = study_order),
    gene = factor(gene, levels = rev(gene_order))
  )

p2b <- ggplot(gene_heat, aes(x = study, y = gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.2) +
  labs(title = "B. Cross-cohort concordance of representative high-confidence genes", x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank())

# 2C. Forest plots for selected genes
forest_genes <- c("MICB", "HLA-DRA", "RASGRP1", "PPP1R1A", "ENTPD3", "ACLY")

gene_forest <- purrr::map_dfr(study_order, function(acc) {
  x <- safe_read_csv(c(
    file.path("data/processed/deg", paste0(acc, "_DE_full.csv")),
    file.path("data/processed/deg", paste0(acc, "_deg_full.csv")),
    file.path("/mnt/data/data/processed/deg", paste0(acc, "_DE_full.csv")),
    file.path("/mnt/data/data/processed/deg", paste0(acc, "_deg_full.csv"))
  )) %>%
    mutate(gene = standardize_gene(gene))
  se_col <- names(x)[tolower(names(x)) %in% c("se", "stderr", "logfc_se")]
  if (length(se_col) == 0) {
    x$SE <- NA_real_
  } else {
    x$SE <- x[[se_col[[1]]]]
  }
  x %>% filter(gene %in% forest_genes) %>%
    transmute(study = acc, gene, logFC, SE,
              lower = logFC - 1.96 * SE,
              upper = logFC + 1.96 * SE)
}) %>%
  mutate(study = factor(study, levels = rev(study_order)),
         gene = factor(gene, levels = forest_genes))

p2c <- ggplot(gene_forest, aes(x = logFC, y = study)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.15, orientation = "y", na.rm = TRUE) +
  geom_point(size = 1.9) +
  facet_wrap(~ gene, scales = "free_x", ncol = 2) +
  labs(title = "C. Forest plots for selected representative genes", x = "logFC (T2D vs ND)", y = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold"))

# 2D. Ranked meta-effects
meta_rank <- meta_hc %>% arrange(meta_FDR, desc(abs(meta_logFC))) %>% mutate(rank = row_number())
label_df <- meta_rank %>% slice_head(n = 15)
p2d <- ggplot(meta_rank, aes(x = rank, y = meta_logFC)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_text(data = label_df, aes(label = gene), size = 2.4, vjust = -0.35, check_overlap = TRUE) +
  labs(title = "D. Ranked summary of meta-analysis effect sizes", x = "Gene rank", y = "Meta-analysis logFC") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

fig2 <- (p2a | p2b) / (p2c | p2d)
save_plot(fig2, file.path(fig_dir, "Figure2_cross_study_core.png"), width = 14, height = 11)

# ------------------------------------------------------------
# Figure 3
# ------------------------------------------------------------

# 3A. Module schematic as gene tiles
module_labels <- refined_modules %>%
  filter(status == "core" | is.na(status)) %>%
  group_by(module) %>%
  mutate(idx = row_number(),
         row = ifelse(idx <= 7, 1, 2),
         col = idx - ifelse(idx <= 7, 0, 7),
         label = gene) %>%
  ungroup() %>%
  mutate(module = factor(module, levels = c("ImmuneStress", "BetaCellIdentitySecretion")))

p3a <- ggplot(module_labels, aes(x = col, y = row, fill = module)) +
  geom_tile(color = "white", linewidth = 0.4, width = 0.95, height = 0.95) +
  geom_text(aes(label = label), size = 3) +
  facet_wrap(~ module, ncol = 1) +
  scale_y_reverse() +
  labs(title = "A. Refined module architecture and gene composition", x = NULL, y = NULL) +
  theme_void(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "none")

# 3B. Heatmap of standardized module genes across cohorts
module_gene_heat <- purrr::map_dfr(study_order, function(acc) {
  x <- safe_read_csv(c(
    file.path("data/processed/deg", paste0(acc, "_DE_full.csv")),
    file.path("data/processed/deg", paste0(acc, "_deg_full.csv")),
    file.path("/mnt/data/data/processed/deg", paste0(acc, "_DE_full.csv")),
    file.path("/mnt/data/data/processed/deg", paste0(acc, "_deg_full.csv"))
  )) %>%
    mutate(gene = standardize_gene(gene)) %>%
    filter(gene %in% refined_modules$gene) %>%
    transmute(study = acc, gene, logFC)
  x
}) %>%
  left_join(refined_modules %>% dplyr::select(module, gene), by = "gene") %>%
  group_by(gene) %>%
  mutate(z = as.numeric(scale(logFC))) %>%
  ungroup()

gene_levels_mod <- refined_modules %>%
  mutate(module = factor(module, levels = c("ImmuneStress", "BetaCellIdentitySecretion"))) %>%
  arrange(module, gene) %>% pull(gene)

p3b <- ggplot(module_gene_heat %>%
                mutate(
                  study = factor(study, levels = study_order),
                  gene = factor(gene, levels = rev(gene_levels_mod))
                ),
              aes(x = study, y = gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  labs(title = "B. Compact heatmap of standardized module genes across cohorts", x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank())

# 3C. Formula panel
formula_df <- tibble(
  x = 0.5,
  y = 0.5,
  label = "IsletDysfunctionScore = ImmuneStress \u2212 BetaCellIdentitySecretion"
)
p3c <- ggplot(formula_df, aes(x, y, label = label)) +
  geom_text(size = 5.2, fontface = "bold") +
  xlim(0, 1) + ylim(0, 1) +
  labs(title = "C. Composite score formula") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0))

fig3 <- (p3a | p3b) / p3c +
  plot_layout(heights = c(3, 1.2))
save_plot(fig3, file.path(fig_dir, "Figure3_refined_modules.png"), width = 14, height = 11)

# ------------------------------------------------------------
# Figure 4
# ------------------------------------------------------------
score_long_nd_t2d <- sample_scores %>%
  filter(disease_group %in% c("ND", "T2D")) %>%
  dplyr::select(sample_id, study, disease_group, ImmuneStress, BetaCellIdentitySecretion, IsletDysfunctionScore) %>%
  pivot_longer(cols = c("ImmuneStress", "BetaCellIdentitySecretion", "IsletDysfunctionScore"),
               names_to = "module", values_to = "score") %>%
  mutate(
    study = factor(study, levels = study_order),
    disease_group = factor(disease_group, levels = c("ND", "T2D")),
    module = factor(module, levels = c("ImmuneStress", "BetaCellIdentitySecretion", "IsletDysfunctionScore"))
  )

p4 <- ggplot(score_long_nd_t2d, aes(x = disease_group, y = score, fill = disease_group)) +
  geom_violin(trim = FALSE, alpha = 0.35, color = NA) +
  geom_boxplot(width = 0.22, outlier.size = 0.8, alpha = 0.8) +
  facet_grid(module ~ study, scales = "free_y") +
  labs(x = NULL, y = "Module score") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "none")
save_plot(p4, file.path(fig_dir, "Figure4_cross_cohort_refined_module_scores.png"), width = 14, height = 9)

# ------------------------------------------------------------
# Figure 5
# ------------------------------------------------------------
# Panel A: GSE20966 beta-cell validation
p5a_dat <- sample_scores %>%
  filter(study == "GSE20966", disease_group %in% c("ND", "T2D")) %>%
  dplyr::select(sample_id, disease_group, ImmuneStress, BetaCellIdentitySecretion, IsletDysfunctionScore) %>%
  pivot_longer(cols = c("ImmuneStress", "BetaCellIdentitySecretion", "IsletDysfunctionScore"),
               names_to = "module", values_to = "score") %>%
  mutate(disease_group = factor(disease_group, levels = c("ND", "T2D")),
         module = factor(module, levels = c("ImmuneStress", "BetaCellIdentitySecretion", "IsletDysfunctionScore")))

p5a <- ggplot(p5a_dat, aes(x = disease_group, y = score, fill = disease_group)) +
  geom_violin(trim = FALSE, alpha = 0.35, color = NA) +
  geom_boxplot(width = 0.22, outlier.size = 0.8, alpha = 0.85) +
  facet_wrap(~ module, scales = "free_y", nrow = 1) +
  labs(title = "A. GSE20966 beta-cell–enriched cohort", x = NULL, y = "Module score") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "none")

# Panel B: GSE164416 anchor progression
p5b_dat <- sample_scores %>%
  filter(study == "GSE164416", disease_group %in% c("ND", "IGT", "T2D", "T3cD")) %>%
  dplyr::select(sample_id, disease_group, ImmuneStress, BetaCellIdentitySecretion, IsletDysfunctionScore) %>%
  pivot_longer(cols = c("ImmuneStress", "BetaCellIdentitySecretion", "IsletDysfunctionScore"),
               names_to = "module", values_to = "score") %>%
  mutate(disease_group = factor(disease_group, levels = c("ND", "IGT", "T2D", "T3cD")),
         module = factor(module, levels = c("ImmuneStress", "BetaCellIdentitySecretion", "IsletDysfunctionScore")))

p5b <- ggplot(p5b_dat, aes(x = disease_group, y = score, fill = disease_group)) +
  geom_violin(trim = FALSE, alpha = 0.35, color = NA) +
  geom_boxplot(width = 0.22, outlier.size = 0.6, alpha = 0.85) +
  facet_wrap(~ module, scales = "free_y", nrow = 1) +
  labs(title = "B. GSE164416 anchor cohort across disease categories", x = NULL, y = "Module score") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "none")

fig5 <- p5a / p5b
save_plot(fig5, file.path(fig_dir, "Figure5_beta_cell_and_anchor_progression.png"), width = 14, height = 9)

# ------------------------------------------------------------
# Figure 6
# ------------------------------------------------------------
fx <- readr::read_csv(effect_sizes_file, show_col_types = FALSE) %>%
  mutate(
    study = factor(study, levels = rev(study_order)),
    module = factor(module, levels = c("ImmuneStress", "BetaCellIdentitySecretion", "IsletDysfunctionScore")),
    se = purrr::pmap_dbl(list(cohens_d, n_ND, n_T2D), se_d),
    lower = cohens_d - 1.96 * se,
    upper = cohens_d + 1.96 * se
  )

p6 <- ggplot(fx, aes(x = cohens_d, y = study)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.15, orientation = "y") +
  geom_point(size = 2) +
  facet_wrap(~ module, scales = "free_x", ncol = 1) +
  labs(x = "Effect size (Cohen's d)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))
save_plot(p6, file.path(fig_dir, "Figure6_module_effect_size_forests.png"), width = 7.5, height = 9)

# ------------------------------------------------------------
# Optional exports for manual assembly or review
# ------------------------------------------------------------
write_csv(deg_summary, file.path(tbl_dir, "Figure2A_deg_summary_data.csv"))
write_csv(gene_heat, file.path(tbl_dir, "Figure2B_gene_heatmap_data.csv"))
write_csv(gene_forest, file.path(tbl_dir, "Figure2C_gene_forest_data.csv"))
write_csv(module_gene_heat, file.path(tbl_dir, "Figure3B_module_gene_heatmap_data.csv"))
write_csv(score_long_nd_t2d, file.path(tbl_dir, "Figure4_module_score_long.csv"))
write_csv(p5a_dat, file.path(tbl_dir, "Figure5A_beta_cell_long.csv"))
write_csv(p5b_dat, file.path(tbl_dir, "Figure5B_anchor_long.csv"))
write_csv(fx, file.path(tbl_dir, "Figure6_effect_sizes_with_ci.csv"))

message("Figures written to: ", normalizePath(fig_dir))
message("Supporting tables written to: ", normalizePath(tbl_dir))
