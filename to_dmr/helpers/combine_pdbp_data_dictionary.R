#' Combine PDBP/DMR data dictionaries into a single data dictionary file which
#' contains all possible fields exactly once and store to `OUTPUT_PARENT`.

library(tidyverse)
library(synapser)

DICTIONARY_PARENT <- "syn24171951"
OUTPUT_PARENT <- "syn24171950"

fetch_dictionaries <- function(dictionary_parent) {
  dictionary_entities <- synGetChildren(dictionary_parent,
                                        includeTypes = list("file"))
  dictionary_entities <- dictionary_entities$asList()
  combined_dictionaries <- purrr::map_dfr(dictionary_entities, function(dict) {
    f <- synGet(dict$id)
    form_name <- f$get("name") %>%
      stringr::str_remove("PDBP_") %>%
      stringr::str_remove("_dataElementExport.csv")
    non_form_specific_fields <- c(
          "GUID", "VisitDate", "SiteName", "AgeVal",
          "VisitTypPDBP", "AgeYrs", "AgeRemaindrMonths")
    df <- read_csv(f$path) %>%
      rename(field_name = `variable name`) %>%
      mutate(form_name = case_when(
          field_name %in% non_form_specific_fields ~ "",
          TRUE ~ form_name)) %>%
    return(df)
  })
  return(combined_dictionaries)
}

store_to_synapse <- function(dictionary, parent) {
  fname <- "pdbp_complete_data_dictionary.csv"
  write_csv(dictionary, fname)
  f <- synapser::File(fname, parent = parent)
  synStore(f)
  unlink(fname)
}

main <- function() {
  synLogin()
  dictionaries <- fetch_dictionaries(dictionary_parent = DICTIONARY_PARENT)
  sorted_fields <- dictionaries %>%
    count(field_name)
  sorted_dictionary <- sorted_fields %>%
    inner_join(dictionaries) %>%
    arrange(desc(n)) %>%
    distinct(field_name, form_name, .keep_all = TRUE) %>%
    select(field_name, form_name, dplyr::everything()) %>%
    select(-n)
  store_to_synapse(sorted_dictionary, parent = OUTPUT_PARENT)
}

main()
