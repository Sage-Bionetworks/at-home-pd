#' This script outputs a CSV to OUTPUT_PARENT containing a mapping of GUIDs to
#' their cohort. The possible cohorts are "at-home-pd" and "super-pd".

library(tidyverse)
library(synapser)

CLINICAL_DATA <- "syn17051543"
OUTPUT_PARENT <- "syn16809549"

get_clinical_data <- function(syn_id) {
  f <- synapser::synGet(syn_id)
  clinical <- readr::read_csv(f$path)
  return(clinical)
}

create_mapping <- function(clinical) {
  mapping <- clinical %>%
    select(guid, redcap_event_name) %>%
    mutate(cohort = case_when(
      str_detect(redcap_event_name, "Arm 1") ~ "at-home-pd",
      str_detect(redcap_event_name, "Arm 2") ~ "super-pd")) %>%
    distinct(guid, cohort)
  return(mapping)
}

store_to_synapse <- function(mapping, parent) {
  fname <- "guid_to_cohort_mapping.csv"
  readr::write_csv(mapping, fname)
  f <- synapser::File(fname, parent = parent)
  synapser::synStore(f)
  unlink(fname)
}

main <- function() {
  synapser::synLogin()
  clinical <- get_clinical_data(CLINICAL_DATA)
  mapping <- create_mapping(clinical)
  store_to_synapse(mapping, OUTPUT_PARENT)
}

main()
