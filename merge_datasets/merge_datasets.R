library(synapser)
library(tidyverse)

ROCHESTER <- "syn18637150"
ROCHESTER_DATA_DICTIONARY <- "syn17051559"
MJFF_PARENT <- "syn18637133"
BRIDGE_SUMMARY <- "syn18681888"
BRIDGE_MAPPING <- list(
    "syn17015960" = "syn18681888",
    "syn17015065" = "syn18681890",
    "syn17014786" = "syn18681891",
    "syn17014785" = "syn18681893",
    "syn17014784" = "syn18681894",
    "syn17014783" = "syn18681898",
    "syn17014781" = "syn18681899",
    "syn17014780" = "syn18681900",
    "syn17014779" = "syn18681901",
    "syn17014778" = "syn18681902",
    "syn17014777" = "syn18681903",
    "syn17014776" = "syn18681904",
    "syn17014775" = "syn18681905",
    "syn17014782" = "syn18683083",
    "syn20712373" = "syn20715163")
TABLE_OUTPUT <- "syn18693245"

read_syn_csv <- function(syn_id, encoding = "UTF-8") {
  f <- synapser::synGet(syn_id)
  df <- readr::read_csv(f$path,
                        locale = readr::locale(encoding = encoding))
  return(df)
}

read_syn_table <- function(syn_id) {
  q <- synapser::synTableQuery(paste("select * from", syn_id))
  table <- q$asDataFrame() %>%
    as_tibble()
  return(table)
}

summarize_rochester <- function() {
  dataset <- read_syn_csv(ROCHESTER)
  visit_dates <- dplyr::distinct(dataset, guid, visstatdttm)
  data_dictionary <- read_syn_csv(ROCHESTER_DATA_DICTIONARY)
  forms <- unique(data_dictionary[["Form Name"]])
  summarized_dataset <- merge(visit_dates, forms, all = TRUE) %>%
    rename(activity = y, createdOn = visstatdttm) %>%
    mutate(activity = as.character(activity),
           createdOn = lubridate::as_datetime(createdOn),
           source = "ROCHESTER") %>%
    as_tibble()
  summarized_dataset <- summarized_dataset %>%
mutate(hash_key = if_else(is.na(createdOn),
                      str_c(guid, "NA"),
                              str_c(guid, createdOn)),
           recordId = unlist(purrr::map(
             hash_key, ~ digest::digest(., algo = "md5")))) %>%
    select(-hash_key)
  return(summarized_dataset)
}

summarize_mjff <- function() {
  files <- synapser::synGetChildren(MJFF_PARENT)$asList()
  names(files) <- lapply(files, function(f) f$id)
  datasets <- purrr::map(files, ~ read_syn_csv(.$id))
  summarized_dataset <- purrr::map2_dfr(files, datasets, function(f, df) {
    activity <- stringr::str_match(
        f$name, "deidentified_((\\w|-)+)(_\\d+)?\\.csv")[, 2]
    if (str_ends(activity, "_\\d\\d")) {
        stringr::str_sub(activity, nchar(activity) - 2, nchar(activity)) <- ""
    }
    if (hasName(df, "study_date")) {
      summarized_dataset <- df %>%
        select(guid, createdOn = study_date) %>%
        mutate(createdOn = lubridate::as_datetime(createdOn),
               activity = activity,
               source = "MJFF")
    } else {
      summarized_dataset <- df %>%
        select(guid) %>%
        mutate(createdOn = NA, activity = activity, source = "MJFF")
    }
  })
  summarized_dataset <- summarized_dataset %>%
    mutate(hash_key = str_c(guid, createdOn),
           recordId = unlist(purrr::map(
             hash_key, ~ digest::digest(., algo = "md5")))) %>%
    select(-hash_key)
  return(summarized_dataset)
}

summarize_bridge <- function() {
  original_tables <- read_syn_table(BRIDGE_SUMMARY) %>%
    select(recordId, guid = externalId, createdOn, appVersion,
           phoneInfo, dataGroups, activity = originalTable) %>%
    filter(activity != "sms-messages-sent-from-bridge-v1",
           activity != "StudyBurstReminder-v1") %>%
    mutate(source = "MPOWER",
           createdOn = lubridate::as_datetime(
             createdOn, tz = "America/Los_Angeles"))
  return(original_tables)
}

mutate_participant_week_day <- function(summarized_all) {
  first_activity <- summarized_all %>%
    filter(source == "MPOWER") %>%
    group_by(guid) %>%
    summarise(first_activity_time = min(createdOn, na.rm = T))
  summarized_all <- dplyr::full_join(summarized_all, first_activity)
  summarized_all_mpower <- summarized_all %>%
    filter(source == "MPOWER") %>%
    mutate(createdOnDate = lubridate::as_date(
        createdOn, tz = "America/Los_Angeles"),
        dayInStudy = (createdOnDate - lubridate::as_date(
             first_activity_time, tz = "America/Los_Angeles")) + 1)
  summarized_all <- summarized_all %>%
    anti_join(summarized_all_mpower, by = "recordId") %>%
    bind_rows(summarized_all_mpower) %>%
    select(-first_activity_time, -createdOnDate)
  return(summarized_all)
}

update_store_merged_datasets <- function(summarized_all) {
  preexisting_summary <- synTableQuery(
    paste("select * from", TABLE_OUTPUT))$asDataFrame()
  new_records <- summarized_all %>%
    anti_join(preexisting_summary, by = "recordId")
  synapser::synStore(synapser::Table(TABLE_OUTPUT, new_records),
                     executed = list(paste0(
        "https://github.com/Sage-Bionetworks/at-home-pd/blob/master/",
        "merge_datasets/merge_datasets.R")))
}

main <- function() {
  synapser::synLogin(
    Sys.getenv("synapseUsername"), Sys.getenv("synapsePassword"))
  summarized_mjff <- summarize_mjff()
  summarized_rochester <- summarize_rochester()
  summarized_bridge <- summarize_bridge()
  summarized_all <- bind_rows(
    summarized_mjff, summarized_rochester, summarized_bridge) %>%
    mutate_participant_week_day() %>%
    select(recordId, guid, source, activity, dplyr::everything())
  update_store_merged_datasets(summarized_all)
}

main()

