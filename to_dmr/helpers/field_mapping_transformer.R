#' Output ctcc_dmr_fields.csv in a tidy format
#' The working directory must be the same as this script's directory
#' in order for this script to execute properly.

library(dplyr)
library(tidyr)
library(stringr)

format_input <- function(ctcc_dmr_fields) {
  ctcc_dmr_fields %>%
    select(dmr_variable, at_home_pd_baseline, at_home_pd_m12, at_home_pd_m24,
           at_home_pd_m36, super_pd_baseline, super_pd_on, super_pd_off) %>%
    pivot_longer(!dmr_variable, names_to="cohort_visit", values_to="clinical_variable") %>%
    mutate(cohort = ifelse(
              str_detect(cohort_visit, "at_home_pd"), "at-home-pd", "super-pd"),
           visit = str_extract(cohort_visit, "[^_]+$"),
           visit = case_when(
             visit == "baseline" ~ "Baseline",
             visit == "m12" ~ "12 months",
             visit == "m24" ~ "24 months",
             visit == "m36" ~ "36 months",
             visit == "on" ~ "Physician_ON",
             visit == "off" ~ "Physician_OFF")) %>%
    select(dmr_variable, cohort, visit, clinical_variable)
}

main <- function() {
  ctcc_dmr_fields <- readr::read_csv("resources/ctcc_dmr_fields.csv")
  dmr_dic <- readr::read_csv("resources/pdbp_complete_data_dictionary.csv") %>%
    select(dmr_variable = field_name, form_name)
  formatted_ctcc_dmr_fields <- format_input(ctcc_dmr_fields) %>%
    inner_join(dmr_dic, by="dmr_variable") %>%
    select(dmr_variable, form_name, dplyr::everything())
  fname <- "../resources/dmr_to_clinical_field_mapping.csv"
  readr::write_csv(formatted_ctcc_dmr_fields, fname)
}

main()
