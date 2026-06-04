source("R/_shared.R")

ensure_raw_cols <- function(df, cols = c("geo_accession", "title", "source_name_ch1")) {
  for (nm in cols) if (!nm %in% colnames(df)) df[[nm]] <- NA_character_
  df
}

curate_gse25724 <- function(raw) {
  tmp <- raw %>% ensure_raw_cols() %>% collapse_characteristics()
  tmp %>%
    dplyr::transmute(
      study = "GSE25724",
      sample_id = dplyr::coalesce(geo_accession, sample_rowname, title),
      sample_label = title,
      donor_id = stringr::str_extract(stringr::str_to_lower(title), "rep\\d+") %>% stringr::str_to_upper(),
      disease_group = dplyr::case_when(
        stringr::str_detect(stringr::str_to_lower(title), "non[- ]?diabetic") ~ "ND",
        stringr::str_detect(stringr::str_to_lower(title), "type 2 diabetic|t2d|diabetic") ~ "T2D",
        TRUE ~ NA_character_
      ),
      tissue_class = "pancreatic_islet",
      cell_context = "islet_bulk",
      source_name = source_name_ch1,
      sex = NA_character_,
      age = NA_real_,
      bmi = NA_real_,
      hba1c = NA_real_,
      in_ins_filtered_data_subset = NA,
      parser_note = "Parsed from sample title: Non-diabetic islets vs Type 2 diabetic islets",
      joined_text = joined_text
    )
}

curate_gse20966 <- function(raw) {
  tmp <- raw %>% ensure_raw_cols() %>% collapse_characteristics()
  tmp %>%
    dplyr::transmute(
      study = "GSE20966",
      sample_id = dplyr::coalesce(geo_accession, sample_rowname, title),
      sample_label = title,
      donor_id = stringr::str_extract(stringr::str_to_lower(title), "donor\\d+") %>% stringr::str_to_upper(),
      disease_group = dplyr::case_when(
        stringr::str_detect(stringr::str_to_lower(title), "non[- ]?diabetic") ~ "ND",
        stringr::str_detect(stringr::str_to_lower(title), "(^|_)diabetic") ~ "T2D",
        TRUE ~ NA_character_
      ),
      tissue_class = "pancreatic_islet",
      cell_context = "beta_cell_enriched",
      source_name = source_name_ch1,
      sex = NA_character_,
      age = NA_real_,
      bmi = NA_real_,
      hba1c = NA_real_,
      in_ins_filtered_data_subset = NA,
      parser_note = "Parsed from title pattern beta-cells_non-diabetic / beta-cells_diabetic",
      joined_text = joined_text
    )
}

curate_gse38642 <- function(raw) {
  tmp <- raw %>% ensure_raw_cols() %>% collapse_characteristics()

  status <- parse_field_from_text(tmp$joined_text, "status")
  sex    <- parse_field_from_text(tmp$joined_text, "sex")
  age    <- parse_field_from_text(tmp$joined_text, "age")
  bmi    <- parse_field_from_text(tmp$joined_text, "bmi")
  hba1c  <- parse_field_from_text(tmp$joined_text, "hba1c")

  safe_chr <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  }

  txt <- stringr::str_to_lower(paste(
    safe_chr(tmp$title),
    safe_chr(tmp$source_name_ch1),
    safe_chr(tmp$joined_text),
    safe_chr(status),
    sep = " | "
  ))

  disease_group <- dplyr::case_when(
    stringr::str_detect(txt, "non[- ]?diab|non diabetic|nondiabetic|control|normal glucose|nd\\b") ~ "ND",
    stringr::str_detect(txt, "type ?2|t2d|diab|hypergly") ~ "T2D",
    TRUE ~ NA_character_
  )

  tmp %>%
    dplyr::transmute(
      study = "GSE38642",
      sample_id = dplyr::coalesce(geo_accession, sample_rowname, title),
      sample_label = title,
      donor_id = dplyr::coalesce(
        stringr::str_extract(title, "ID\\d+"),
        stringr::str_extract(source_name_ch1, "ID\\d+")
      ),
      disease_group = disease_group,
      tissue_class = "pancreatic_islet",
      cell_context = "islet_bulk",
      source_name = source_name_ch1,
      sex = standardize_sex(sex),
      age = clean_numeric_value(age),
      bmi = clean_numeric_value(bmi),
      hba1c = clean_numeric_value(hba1c),
      in_ins_filtered_data_subset = NA,
      parser_note = "Parsed from title/source/characteristics with broad ND/T2D matching",
      joined_text = joined_text
    )
}

curate_gse164416 <- function(raw) {
  tmp <- raw %>% ensure_raw_cols() %>% collapse_characteristics()
  tmp %>%
    dplyr::transmute(
      study = "GSE164416",
      sample_id = dplyr::coalesce(geo_accession, sample_rowname, title),
      sample_label = title,
      donor_id = stringr::str_extract(title, "DP\\d+"),
      disease_group = dplyr::case_when(
        stringr::str_detect(title, "_ND$") ~ "ND",
        stringr::str_detect(title, "_IFG$") ~ "IFG",
        stringr::str_detect(title, "_IGT$") ~ "IGT",
        stringr::str_detect(title, "_T3cD$") ~ "T3cD",
        stringr::str_detect(title, "_T2D$") ~ "T2D",
        TRUE ~ NA_character_
      ),
      tissue_class = "pancreatic_islet",
      cell_context = "lcm_islet",
      source_name = source_name_ch1,
      sex = NA_character_,
      age = NA_real_,
      bmi = NA_real_,
      hba1c = NA_real_,
      in_ins_filtered_data_subset = parse_bool_field(joined_text, "in_ins_filtered_data_subset"),
      parser_note = "Parsed from title suffix and characteristic in_ins_filtered_data_subset",
      joined_text = joined_text
    )
}

curate_one_study <- function(accession) {
  pheno <- safe_read_rds(file.path("data/raw", accession, paste0(accession, "_pheno_raw.rds")))
  raw <- extract_characteristics_table(pheno)

  out <- switch(
    accession,
    GSE25724 = curate_gse25724(raw),
    GSE20966 = curate_gse20966(raw),
    GSE38642 = curate_gse38642(raw),
    GSE164416 = curate_gse164416(raw),
    stop(glue::glue("No parser implemented for {accession}"))
  ) %>%
    dplyr::mutate(
      include_primary_analysis = disease_group %in% c("ND", "T2D"),
      include_secondary_analysis = accession == "GSE164416" & disease_group %in% c("ND", "IFG", "IGT", "T2D", "T3cD"),
      include_anchor_preferred_subset = accession == "GSE164416" & include_secondary_analysis & in_ins_filtered_data_subset %in% TRUE
    ) %>%
    dplyr::select(
      study, sample_id, sample_label, donor_id, disease_group,
      tissue_class, cell_context, sex, age, bmi, hba1c,
      include_primary_analysis, include_secondary_analysis,
      include_anchor_preferred_subset, in_ins_filtered_data_subset,
      source_name, parser_note, joined_text
    )

  out
}

all_meta <- purrr::map_dfr(ACCESSIONS, curate_one_study)
validation_tbl <- purrr::map_dfr(ACCESSIONS, ~ validate_metadata(dplyr::filter(all_meta, study == .x), .x))

for (acc in ACCESSIONS) {
  tmp <- dplyr::filter(all_meta, study == acc)
  write_csv_safe(tmp, file.path("data/processed/metadata", paste0(acc, "_metadata_harmonized.csv")))

  summary_tbl <- tmp %>%
    dplyr::count(disease_group, include_primary_analysis, include_secondary_analysis, name = "n") %>%
    dplyr::mutate(study = acc, .before = 1)
  write_csv_safe(summary_tbl, file.path("results/supplementary", paste0(acc, "_metadata_group_summary.csv")))
}

write_csv_safe(all_meta, "data/processed/metadata/all_metadata_harmonized.csv")
write_csv_safe(validation_tbl, "results/supplementary/metadata_parser_validation.csv")

cohort_tbl <- all_meta %>%
  dplyr::count(study, disease_group, cell_context, include_primary_analysis, name = "n")
write_csv_safe(cohort_tbl, "results/tables/Table_1_cohort_characteristics.csv")

message2("02_metadata_curation complete. Accession-specific parsers applied for all four studies.")
