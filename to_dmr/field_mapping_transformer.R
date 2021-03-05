#' Output ctcc_dmr_fields.csv in a tidy format

library(dplyr)
library(tidyr)
library(stringr)

read_args <- function() {
  parser <- argparse::ArgumentParser()
  parser$add_argument("--input", default="syn23593083")
  parser$add_argument("--output", default="syn16809549")
  args <- parser$parse_args()
  return(args)
}

format_input <- function(ctcc_dmr_fields) {
  ctcc_dmr_fields %>%
    select(dmr_variable, at_home_pd_baseline, at_home_pd_m12, at_home_pd_m24,
           super_pd_baseline, super_pd_on, super_pd_off) %>%
    pivot_longer(!dmr_variable, names_to="cohort_visit", values_to="clinical_variable") %>%
    mutate(cohort = ifelse(
              str_detect(cohort_visit, "at_home_pd"), "at-home-pd", "super-pd"),
           visit = str_extract(cohort_visit, "[^_]+$"),
           visit = case_when(
             visit == "baseline" ~ "Baseline",
             visit == "m12" ~ "12 months",
             visit == "m24" ~ "24 months",
             visit == "on" ~ "Physician_ON",
             visit == "off" ~ "Physician_OFF")) %>%
    select(dmr_variable, cohort, visit, clinical_variable)
}

read_from_synapse <- function(synapse_id) {
  f <- synapser::synGet(synapse_id)
  df <- readr::read_csv(f$path)
  return(df)
}

store_to_synapse <- function(formatted_ctcc_dmr_fields, synapse_parent) {
  fname <- "tidy_ctcc_dmr_fields.csv"
  readr::write_csv(formatted_ctcc_dmr_fields, fname)
  f <- synapser::File(fname, synapse_parent)
  synapser::synStore(f)
  unlink(fname)
}

main <- function() {
  args <- read_args()
  synapser::synLogin()
  ctcc_dmr_fields <- read_from_synapse(args$input)
  formatted_ctcc_dmr_fields <- format_input(ctcc_dmr_fields)
  store_to_synapse(
    formatted_ctcc_dmr_fields = formatted_ctcc_dmr_fields,
    synapse_parent = args$output)
}

main()
