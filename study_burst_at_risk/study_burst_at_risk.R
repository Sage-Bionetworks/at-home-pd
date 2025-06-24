library(synapser)
library(tidyverse)

MPOWER_SUMMARY <- "syn18693245"
USER_OFFSET <- "syn18637903"
TABLE_OUTPUT <- "syn20709661"

read_syn_table <- function(syn_id) {
  q <- synTableQuery(paste("select * from", syn_id))
  table <- q$asDataFrame() %>% 
    as_tibble() 
  return(table)
}

fetch_mpower <- function() {
  user_offset <- read_syn_table(USER_OFFSET) %>% 
    select(-ROW_ID, -ROW_VERSION)
  mpower <- read_syn_table(MPOWER_SUMMARY) %>% 
    filter(source == "MPOWER") %>% 
    left_join(user_offset, by = c("guid", "source")) %>% 
    mutate(createdOn = lubridate::with_tz(createdOn, "America/Los_Angeles"),
           createdOn = createdOn - lubridate::days(day_offset)) %>% 
    select(-day_offset)
  return(mpower)
}

identify_at_risk_users <- function(mpower) {
  mpower <- mpower %>% 
    group_by(guid) %>%
    mutate(first_activity = lubridate::as_date(min(createdOn)),
           currentDayInStudy = as.integer(
             lubridate::today(tz="America/Los_Angeles") - first_activity),
           currentlyInStudyBurst = currentDayInStudy %% 90 < 20,
           currentDayInStudyBurst = {ifelse(currentlyInStudyBurst, currentDayInStudy %% 90, NA)},
           currentStudyBurstNumber = {ifelse(currentlyInStudyBurst, currentDayInStudy %/% 90, NA)},
           studyBurstStart = currentStudyBurstNumber * 90,
           studyBurstEnd = studyBurstStart + 20,
           daysRemainingInStudyBurst = studyBurstEnd - currentDayInStudy) %>% 
    filter(str_detect(activity, "StudyBurst"))
  current_study_burst_users <- mpower %>% 
    filter(currentlyInStudyBurst) %>%
    summarize(daysRemainingInStudyBurst = median(daysRemainingInStudyBurst),
              currentStudyBurst = median(currentStudyBurstNumber))
  at_risk_users <- mpower %>% 
    filter(currentlyInStudyBurst,
           dayInStudy > studyBurstStart & dayInStudy < studyBurstEnd) %>% 
    summarize(daysCompletedInStudyBurst = n_distinct(dayInStudy)) %>% 
    full_join(current_study_burst_users, by = "guid") %>% 
    mutate(daysCompletedInStudyBurst = {ifelse(
      is.na(daysCompletedInStudyBurst), 0, daysCompletedInStudyBurst)}) %>% 
    filter((daysCompletedInStudyBurst == 0 & daysRemainingInStudyBurst <= 18) | 
            daysRemainingInStudyBurst + daysCompletedInStudyBurst <= 15,
           daysCompletedInStudyBurst < 10,
           daysCompletedInStudyBurst + daysRemainingInStudyBurst >= 5)
  at_risk_past_activity <- mpower %>% 
    ungroup() %>% 
    semi_join(at_risk_users, by = "guid") %>% 
    mutate(previousStudyBurstStart = studyBurstStart - 90,
           previousStudyBurstEnd = studyBurstEnd - 90) %>% 
    filter(previousStudyBurstStart >= 0,
           dayInStudy >= previousStudyBurstStart,
           dayInStudy <= previousStudyBurstEnd) %>% 
    group_by(guid) %>% 
    summarize(daysCompletedPreviousStudyBurst = n_distinct(dayInStudy))
  at_risk_users <- at_risk_users %>% 
    left_join(at_risk_past_activity, by = "guid") %>% 
    mutate(daysCompletedPreviousStudyBurst = {ifelse(
      is.na(daysCompletedPreviousStudyBurst) & currentStudyBurst > 0, # bursts are 0-indexed
      0, daysCompletedPreviousStudyBurst)}) %>% 
    mutate(currentStudyBurst = currentStudyBurst + 1) # bursts are now 1-indexed
  return(at_risk_users) 
}

store_to_synapse <- function(at_risk_users) {
  q <- synTableQuery(paste("select * from", TABLE_OUTPUT))
  synDelete(q) # Remove preexisting rows
  t <- synTable(TABLE_OUTPUT, at_risk_users)  
  synStore(t)
}

main <- function() {
  synLogin(authToken=Sys.getenv("synapseAccessToken"))
  mpower <- fetch_mpower()
  at_risk_users <- identify_at_risk_users(mpower)
  store_to_synapse(at_risk_users)
}

main()
