#' This script produces a table derived from STUDY_BURST_SUMMARY (below) with columns:
#'
#' * study_burst (str)
#' * complete (int)
#' * partially_complete (int)
#' * to_complete (int)
#' * percent_complete (double)
#' * percent_no_activity (double)
#' * percent_any_activity (double)
#'
#' And stores the output to TABLE_OUTPUT (global var below). The sum of
#' percent_no_activity and percent_any_activity may be greater than 100
#' because percent_no_activity is proportion who did not complete any
#' study burst days, whereas percent_any_activity is proportion who
#' completed any active task during their study burst period.
library(synapser)
library(tidyverse)

STUDY_BURST_SUMMARY <- Sys.getenv("studyBurstSummaryTable")
TABLE_OUTPUT <- Sys.getenv("outputTable")
HEALTH_DATA_SUMMARY_TABLE <- Sys.getenv("studyBurstSummaryTable")

read_syn_table <- function(syn_id) {
  q <- synTableQuery(paste("select * from", syn_id))
  table <- q$asDataFrame() %>%
    as_tibble()
  return(table)
}

store_to_synapse <- function(compliance_overview) {
  q <- synTableQuery(paste("select * from", TABLE_OUTPUT))
  synDelete(q) # Remove preexisting rows
  t <- synapser::Table(TABLE_OUTPUT, compliance_overview)
  synStore(t)
}

build_compliance_overview <- function(study_burst_summary, any_activity_counts) {
  num_no_activity <- study_burst_summary %>%
    filter(lubridate::as_date(study_burst_start_date) <= lubridate::today(),
           days_completed == 0) %>%
    count(study_burst, name = "no_activity")
  num_any_activity <- any_activity_counts %>%
    group_by(study_burst) %>%
    summarize(any_activity = sum(any_activity > 0, na.rm=TRUE))
  compliance_overview <- study_burst_summary %>%
    group_by(study_burst) %>%
    summarize(complete = sum(study_burst_successful, na.rm = T),
              partially_complete = sum(!study_burst_successful, na.rm = T),
              to_complete = sum(is.na(study_burst_successful)),
              percent_compliant =  round(complete / (complete + partially_complete), 2)) %>%
    left_join(num_no_activity) %>%
    inner_join(num_any_activity) %>%
    mutate(no_activity = replace_na(no_activity, 0),
           partially_complete = partially_complete - no_activity,
           percent_no_activity = round(
              no_activity / (no_activity + complete + partially_complete), 2),
           percent_any_activity = round(
              any_activity / (no_activity + complete + partially_complete), 2)) %>%
    select(study_burst, complete, partially_complete, no_activity,
           to_complete, percent_compliant, percent_no_activity, percent_any_activity) %>%
    arrange(study_burst)
  return(compliance_overview)
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

fetch_health_data_summary_table <- function(health_data_summary_table) {
  mpower <- read_syn_table(health_data_summary_table)
  mpower$createdOnTimeZoneInteger <- unlist(purrr::map(mpower$createdOnTimeZone,
                                                       get_timezone_as_integer))
  mpower <- mpower %>%
    mutate(createdOnLocalTime = createdOn + lubridate::hours(createdOnTimeZoneInteger)) %>%
    rename(guid = externalId)
  return(mpower)
}

get_any_activity_counts <- function(health_data_summary_table, study_burst_summary) {
  study_burst_schedule <- study_burst_summary %>%
    select(guid, study_burst, study_burst_start_date, study_burst_end_date)
  days_any_activity <- purrr::pmap_dfr(study_burst_schedule,
    function(guid_, study_burst, study_burst_start_date, study_burst_end_date) {
      study_burst_start_date <- lubridate::as_date(study_burst_start_date)
      study_burst_end_date <- lubridate::as_date(study_burst_end_date)
      relevant_activities <- health_data_summary_table %>%
        filter(guid == guid_,
               createdOnLocalTime >= study_burst_start_date,
               createdOnLocalTime <= study_burst_end_date + lubridate::days(1))
      days_any_activity <- relevant_activities %>%
        mutate(createdOnDate = lubridate::as_date(createdOnLocalTime)) %>%
        group_by(createdOnDate) %>%
        summarize(any_activity = ("Tremor-v3" %in% originalTable |
                               "WalkAndBalance-v1" %in% originalTable |
                               "Tapping-v4" %in% originalTable))
      any_activity_this_burst <- sum(days_any_activity$any_activity)
      # If the participant has not yet finished this study burst, store NA for days completed
      if (study_burst_end_date >= lubridate::today()) {
        any_activity_this_burst <- NA
      }
      result <- tibble(
        guid = guid_,
        study_burst = study_burst,
        any_activity = any_activity_this_burst)
      return(result)
  })
  return(days_any_activity)
}

main <- function() {
  synLogin(Sys.getenv("synapseUsername"), Sys.getenv("synapsePassword"))
  study_burst_summary <- read_syn_table(STUDY_BURST_SUMMARY)
  health_data_summary_table <- fetch_health_data_summary_table(
      health_data_summary_table = HEALTH_DATA_SUMMARY_TABLE)
  any_activity_counts <- get_any_activity_counts(
      health_data_summary_table = health_data_summary_table,
      study_burst_summary = study_burst_summary)
  compliance_overview <- build_compliance_overview(
     study_burst_summary = study_burst_summary,
     any_activity_counts = any_activity_counts)
  store_to_synapse(compliance_overview)
}

main()
