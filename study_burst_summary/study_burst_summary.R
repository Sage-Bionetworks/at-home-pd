#' This script produces a table with columns:
#'
#' * guid (str)
#' * guid_prefix (str)
#' * study_burst (str)
#' * study_burst_start_date (str)
#' * study_burst_end_date (str)
#' * days_completed (int)
#' * study_burst_successful (bool)
#'
#' And stores this result to TABLE_OUTPUT (global var below).
#'
#' This is more granular than the table produced by compliance_overview.R,
#' which summarizes study burst compliance at the study burst level (Y1,Q1, etc.)
library(synapser)
library(tidyverse)
library(bridgeclient)

HEALTH_DATA_SUMMARY_TABLE <- Sys.getenv("inputTable")
TABLE_OUTPUT <- Sys.getenv("outputTable")
GUID_PREFIX_LENGTH <- Sys.getenv("guidPrefixLength")

read_syn_table <- function(syn_id) {
  q <- synapser::synTableQuery(paste("select * from", syn_id))
  table <- q$asDataFrame() %>%
    as_tibble()
  return(table)
}

get_timezone_as_integer <- function(createdOnTimeZone) {
  # If there is no timezone information we make the conservative
  # (for a US user) estimate that the time zone is Pacific
  if (is.na(createdOnTimeZone)) {
    return(-8)
  } else {
    cotz_integer <- as.integer(as.integer(createdOnTimeZone) / 100)
    return(cotz_integer)
  }
}

fetch_mpower <- function(health_data_summary_table) {
  mpower <- read_syn_table(health_data_summary_table)
  mpower$createdOnTimeZoneInteger <- unlist(purrr::map(mpower$createdOnTimeZone,
                                                       get_timezone_as_integer))
  mpower <- mpower %>%
    mutate(createdOnLocalTime = createdOn + lubridate::hours(createdOnTimeZoneInteger)) %>%
    rename(guid = externalId)
  return(mpower)
}

get_schedule <- function(guids) {
    schedule <- purrr::map_dfr(unique(guids), function(guid) {
      participant_info <- bridgeclient::get_participant(
              external_id = guid)
      activity_events <- bridgeclient::get_activity_events(
              user_id = participant_info$id) %>%
        bind_rows() %>%
        filter(str_detect(eventId, "custom:activityBurst") |
                 eventId == "study_start_date") %>%
        mutate(study_burst = as.character(as.integer(
          str_extract(.$eventId, "\\d+")) - 1),
               study_burst = ifelse(
                 is.na(study_burst), "0", study_burst),
               study_burst_start_date = lubridate::as_date(
                 lubridate::as_datetime(timestamp, tz = "US/Pacific")),
               study_burst_end_date = study_burst_start_date + lubridate::days(19),
               guid = guid) %>%
        distinct(guid, study_burst, study_burst_start_date, study_burst_end_date) %>%
        arrange(guid, study_burst)
    })
    schedule$study_burst <- dplyr::recode(
      schedule$study_burst, "0" = "Y1,Q1", "1" = "Y1,Q2", "2" = "Y1,Q3",
      "3" = "Y1,Q4", "4" = "Y2,Q1", "5" = "Y2,Q2", "6" = "Y2,Q3",
      "7" = "Y2,Q4", "8" = "Y3,Q1", "9" = "Y3,Q2", "10" = "Y3,Q3",
      "11" = "Y3,Q4", "12" = "Y4,Q1", "13" = "Y4,Q2", "14" = "Y4,Q3",
      "15" = "Y4,Q4", "16" = "Y5,Q1", "17" = "Y5,Q2", "18" = "Y5,Q3",
      "19" = "Y5,Q4", "20" = "Y6,Q1", "21" = "Y6,Q2", "22" = "Y6,Q3",
      "23" = "Y6,Q4", "24" = "Y7,Q1", "25" = "Y7,Q2", "26" = "Y7,Q3",
      "27" = "Y7,Q4", "28" = "Y8,Q1", "29" = "Y8,Q2", "30" = "Y8,Q3",
      "31" = "Y8,Q4", "32" = "Y9,Q1", "33" = "Y9,Q2", "34" = "Y9,Q3",
      "35" = "Y9,Q4", "36" = "Y10,Q1", "37" = "Y10,Q2", "38" = "Y10,Q3",
      "39" = "Y10,Q4", "40" = "Y11,Q1", "41" = "Y11,Q2", "42" = "Y11,Q3",
      "43" = "Y11,Q4", "44" = "Y12,Q1", "45" = "Y12,Q2", "46" = "Y12,Q3",
      "47" = "Y12,Q4")
    schedule <- schedule %>%
      mutate(study_burst_start_date = as.character(study_burst_start_date),
             study_burst_end_date = as.character(study_burst_end_date))
    return(schedule)
}

build_study_burst_schedule <- function(mpower, previous_schedule_id) {
  previous_schedule <- read_syn_table(previous_schedule_id)
  new_participants <- mpower %>%
    filter(!is.na(guid)) %>%
    anti_join(previous_schedule, by = "guid") %>%
    distinct(guid)
  all_participant_guids <- c(new_participants$guid, previous_schedule$guid)
  bridge_login(study = "sage-mpower-2",
               email = Sys.getenv("bridgeUsername"),
               password = Sys.getenv("bridgePassword"))
  current_schedule <- get_schedule(all_participant_guids)
  schedule_diff <- anti_join(
    current_schedule,
    previous_schedule,
    by = c("guid", "study_burst")
  )
  if (nrow(schedule_diff) > 0) {
    new_schedule <- previous_schedule %>%
      bind_rows(schedule_diff)
  } else {
    new_schedule <- previous_schedule
  }
  new_schedule <- new_schedule %>%
      distinct(guid, study_burst, study_burst_start_date, study_burst_end_date) %>%
      arrange(guid, study_burst)
  return(new_schedule)
}

build_study_burst_summary <- function(mpower, study_burst_schedule, guid_prefix_length) {
  #' We make the conservative (for a US user) assumption that the current day is
  #' relative to Pacific time.
  days_completed <- purrr::pmap_dfr(study_burst_schedule,
    function(guid_, study_burst, study_burst_start_date, study_burst_end_date) {
      study_burst_start_date <- lubridate::as_date(study_burst_start_date)
      study_burst_end_date <- lubridate::as_date(study_burst_end_date)
      relevant_activities <- mpower %>%
        filter(guid == guid_,
               createdOnLocalTime >= study_burst_start_date,
               createdOnLocalTime <= study_burst_end_date + lubridate::days(1))
      days_completed <- relevant_activities %>%
        mutate(createdOnDate = lubridate::as_date(createdOnLocalTime)) %>%
        group_by(createdOnDate) %>%
          summarize(success = ("Tremor-v3" %in% originalTable &
                               "WalkAndBalance-v1" %in% originalTable &
                               "Tapping-v4" %in% originalTable))
      days_completed_this_burst <- sum(days_completed$success)
      # If the participant has not yet finished this study burst, store NA for days completed
      if (study_burst_end_date >= lubridate::today()) {
        days_completed_this_burst <- NA
      }
      result <- tibble(
        guid = guid_,
        study_burst = study_burst,
        study_burst_start_date = study_burst_start_date,
        study_burst_end_date = study_burst_end_date,
        days_completed = days_completed_this_burst,
        study_burst_successful = days_completed >= 10)
      return(result)
  })
  study_burst_summary <- days_completed %>%
    mutate(study_burst_start_date = as.character(study_burst_start_date),
           study_burst_end_date = as.character(study_burst_end_date),
           guid_prefix = str_extract(
              guid, glue::glue("^.{{{guid_prefix_length}}}"))) %>%
    arrange(guid, study_burst) %>%
    select(guid, guid_prefix, study_burst, study_burst_start_date,
           study_burst_end_date, days_completed, study_burst_successful)
  return(study_burst_summary)
}

store_to_synapse <- function(study_burst_summary, table_output) {
  q <- synapser::synTableQuery(paste("select * from", table_output))
  synapser::synDelete(q) # Remove preexisting rows
  t <- synapser::synTable(table_output, study_burst_summary)
  synStore(t, executed=list(
    paste0("https://github.com/Sage-Bionetworks/at-home-pd/",
    "blob/master/study_burst_summary/study_burst_summary.R")))
}

main <- function() {
  synapser::synLogin(Sys.getenv("synapseUsername"), Sys.getenv("synapsePassword"))
  mpower <- fetch_mpower(HEALTH_DATA_SUMMARY_TABLE)
  study_burst_schedule <- build_study_burst_schedule(mpower, TABLE_OUTPUT)
  study_burst_summary <- build_study_burst_summary(
      mpower, study_burst_schedule, GUID_PREFIX_LENGTH)
  store_to_synapse(study_burst_summary, TABLE_OUTPUT)
}

main()
