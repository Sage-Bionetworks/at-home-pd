#' This script produces a table with columns:
#' 
#' * guid (str)
#' * televisits (int)
#' * completed_annual_televisit (bool)
#' * completed_active_tasks_annual_visit (bool)
#' 
#' And stores the output to TABLE_OUTPUT (global var below)
library(tidyverse)
library(synapser)

HEALTH_DATA_SUMMARY_TABLE <- "syn17015960"
TABLE_OUTPUT <- "syn22911752"
CLINICAL_DATA <- "syn17051543"

read_syn_table <- function(syn_id) {
  q <- synTableQuery(paste("select * from", syn_id))
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

# replace once https://www.pivotaltracker.com/story/show/174677990 is resolved
fetch_health_data_summary <- function() {
  mpower <- read_syn_table(HEALTH_DATA_SUMMARY_TABLE)
  mpower$createdOnTimeZoneInteger <- unlist(purrr::map(mpower$createdOnTimeZone,
                                                       get_timezone_as_integer))
  mpower <- mpower %>% 
    mutate(createdOnLocalTime = createdOn + lubridate::hours(createdOnTimeZoneInteger))
  first_activity <- mpower %>% 
    filter(originalTable != "sms-messages-sent-from-bridge-v1",
           !is.na(createdOnLocalTime)) %>% 
    group_by(externalId) %>%
    summarize(first_activity = lubridate::as_date(min(createdOnLocalTime, na.rm = T)),
              currentDayInStudy = as.integer(
                lubridate::today(tz="America/Los_Angeles") - first_activity),
              currentlyInStudyBurst = currentDayInStudy %% 90 < 20)
  mpower <- left_join(mpower, first_activity) %>% 
    mutate(dayInStudy = as.integer(lubridate::as_date(createdOnLocalTime) - first_activity))
  return(mpower)
}

count_active_tasks <- function(health_data_summary) {
  active_tasks <- c("Tapping-v4", "Tremor-v3", "WalkAndBalance-v1")
  all_possible_dates <- health_data_summary %>% 
    distinct(externalId, first_activity) %>% 
    filter(!is.na(first_activity)) %>% 
    purrr::pmap_dfr(function(externalId, first_activity) {
      tibble(externalId = externalId,
             date = seq.Date(first_activity, lubridate::today(), by = "day"))
    })
  active_task_counter <- health_data_summary %>% 
    filter(originalTable %in% active_tasks) %>% 
    mutate(activity = case_when(
      originalTable == "Tapping-v4" ~ "tapping",
      originalTable == "Tremor-v3" ~ "tremor",
      originalTable == "WalkAndBalance-v1" ~ "walk"),
      date = lubridate::as_date(createdOnLocalTime)) %>% 
    group_by(externalId, date) %>%
    summarize(completed_motor_activities = ("tremor" %in% activity &
                                            "tapping" %in% activity &
                                            "walk" %in% activity)) %>% 
    full_join(all_possible_dates) %>% 
    mutate(completed_motor_activities = replace_na(completed_motor_activities, FALSE)) %>% 
    arrange(externalId, date)
  return(active_task_counter)   
}

# for answering how many participants did their activities w.r.t. tele-visit
clinical_data <- function() {
  clinical <- readr::read_csv(synGet(CLINICAL_DATA)$path) %>% 
    filter(redcap_event_name %in% c("Month 12 (Arm 1: Arm 1)",
                                    "Month 24 (Arm 1: Arm 1)")) %>% 
    select(externalId = guid, redcap_event_name, visstatdttm) %>% 
    filter(!is.na(visstatdttm)) %>%
    mutate(date = lubridate::as_date(visstatdttm))
  return(clinical)
}

store_to_synapse <- function(results) {
  q <- synTableQuery(paste("select * from", TABLE_OUTPUT))
  synDelete(q) # Remove preexisting rows
  t <- synTable(TABLE_OUTPUT, results)  
  synStore(t)
}

main <- function() {
  synLogin(authToken=Sys.getenv("synapseAccessToken"))
  health_data_summary <- fetch_health_data_summary()
  active_task_counter <- count_active_tasks(health_data_summary)
  clinical <- clinical_data()
  results <- clinical %>% 
    inner_join(active_task_counter, by = c("externalId", "date")) %>% 
    mutate( 
      completed_active_tasks_year_one = (
        redcap_event_name == "Month 12 (Arm 1: Arm 1)" & completed_motor_activities),
      completed_active_tasks_year_two = (
        redcap_event_name == "Month 24 (Arm 1: Arm 1)" & completed_motor_activities)) %>% 
    group_by(externalId) %>% 
    summarize(televisits = n(),
              completed_active_tasks_year_one = any(completed_active_tasks_year_one),
              completed_active_tasks_year_two = any(completed_active_tasks_year_two),
              completed_year_one_televisit = any(redcap_event_name == "Month 12 (Arm 1: Arm 1)"),
              completed_year_two_televisit = any(redcap_event_name == "Month 24 (Arm 1: Arm 1)")) %>% 
    select(guid = externalId, televisits,
           completed_year_one_televisit, completed_year_two_televisit,
           completed_active_tasks_year_one, completed_active_tasks_year_two)
  store_to_synapse(results) 
}

main()
