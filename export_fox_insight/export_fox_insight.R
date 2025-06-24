library(foxden)
library(tidyverse)
library(synapser)

OUTPUT_PARENT <- "syn21670408"

get_fox_surveys <- function(download_path = ".") {
  all_tables <- list_tables()
  table_paths <- names(all_tables) %>%
    purrr::map(download_table, download_to = download_path)
  identifier_mapping <- get_identifier_mapping() 
  filtered_tables <- purrr::map(table_paths, function(p) {
    survey <- read_csv(p, guess_max = 100000) %>% 
      inner_join(identifier_mapping, by = "fox_insight_id") %>% 
      select(guid, fox_insight_id, dplyr::everything())
    return(survey)
  })
  names(filtered_tables) <- names(all_tables)
  return(filtered_tables)
}

get_identifier_mapping <- function() {
  identifier_mapping_path <- download_table(table = "at_home_pd_user_ids",
                                            database = "substudies")
  identifier_mapping <- read_csv(identifier_mapping_path) %>% 
    select(fox_insight_id, guid = ancillary_study_identifier) %>% 
    drop_na()
  return(identifier_mapping)
}

store_fox_surveys <- function(fox_surveys, write_path = ".") {
  purrr::map2(fox_surveys, names(fox_surveys), function(df, survey_name) {
    if (nrow(df)) {
      fname <- file.path(write_path, paste0(survey_name, ".csv"))
      write_csv(df, fname)
      f <- synapser::File(fname, parent = OUTPUT_PARENT)
      synStore(f)
    }
  })
}

main <- function() {
  fox_login(Sys.getenv("foxEmail"), Sys.getenv("foxPassword"))
  synLogin(authToken=Sys.getenv("synapseAccessToken"))
  fox_surveys <- get_fox_surveys(download_path = "/root")
  store_fox_surveys(fox_surveys, write_path = "/root")
}

main()
