#################################################
#' Script to generate AT-HOME-PD UPDRS scores
#' to synapse. This script should be run from the repo root.
#' @author: Aryton Tediarjo
##################################################
library(tidyverse)
library(tidyr)
library(data.table)
library(synapser)
library(glue)
source("compute_updrs_score/scoring_utils.R")
synLogin()

#' Synapse reference
CLINICAL_DATA_SYN_ID <- "syn18565500"
FOX_DATA_SYN_ID <- "syn21670565"
EXCEL_LOOKUP_SYN_ID <- "syn25050918"
SYNAPSE_PARENT_ID <- "syn16809549"
SUPER_PARENT_ID <- "syn25715411"
FOX_PARENT_ID <- "syn16809550"
OUTPUT_FILENAME <- "computed_ahpd_updrs_scores.csv"
FOX_OUTPUT_FILENAME <- glue("fox_{OUTPUT_FILENAME}")

#' instantiate github
GIT_URL <- file.path(
  "https://github.com/Sage-Bionetworks/at-home-pd/blob/master",
  "compute_updrs_score",
  "compute_updrs_score.R")

#' Final scoring metrics
SCORES <- c(
  "UPDRS1", "UPDRS2", "UPDRS3",
  "UPDRS3R", "UPDRS3", "UPDRSAMB",
  "UPDRSAMR", "UPDRSPRO", "MAX_RTTR",
  "MAX_PKTR", "UPDRSTOT", "TRM_TYPE")

#' function to get AT-HOME-PD clincal data
#' and clean the data by removing testing
#' and filtering based on redcap events
get_ahpd_clinical_data <- function(){
  clinical_records <- fread(synGet(CLINICAL_DATA_SYN_ID)$path) %>%
    tibble::as_tibble(.)
  section_1_3_4 <- clinical_records %>%
    dplyr::filter(!stringr::str_detect(guid, "TEST"),
                  redcap_event_name %in% c(
                    "baseline_arm_1", "month_12_arm_1", "month_24_arm_1"),
                  redcap_repeat_instrument != "reportable_event") %>%
    dplyr::mutate(
      visit = case_when(
          redcap_event_name == "baseline_arm_1" ~ "Baseline",
          redcap_event_name == "month_12_arm_1" ~ "Month 12",
          redcap_event_name == "month_24_arm_1" ~ "Month 24")) %>%
    dplyr::select(guid,
                  redcap_event_name,
                  visit,
                  visstatdttm,
                  viscompltyn,
                  c("mdsupdrs_dttm":"painfloffstatdystn"),
                  matches("neck|postinst_instr"))
  section_1_2 <- clinical_records %>%
    dplyr::filter(!stringr::str_detect(guid, "TEST"),
                  redcap_event_name %in% c(
                    "baseline_pre_visit_arm_1",
                    "month_12_pre_visit_arm_1",
                    "month_24_pre_visit_arm_1")) %>%
    dplyr::mutate(visit = case_when(
      redcap_event_name == "baseline_pre_visit_arm_1" ~ "Baseline",
      redcap_event_name == "month_12_pre_visit_arm_1" ~ "Month 12",
      redcap_event_name == "month_24_pre_visit_arm_1" ~ "Month 24")) %>%
    dplyr::select(guid,
                  visit,
                  qstnnreinfoprovdrt:mdsupdrsfreezingscore,
                  qstnnreinfoprovdrt_m12:mdsupdrsfreezingscore_m12) %>%
    pivot_longer(!guid:visit) %>%
    drop_na(value) %>%
    mutate(name = if_else(
      str_detect(name, "_m12"),
      str_remove(name, "_m12"),
      name),
      name = if_else(
        name %in% c("sleepprobscore", "daytmsleepscore",
                    "painothrsensscore", "urnryprobscore",
                    "constipprobscore", "slivadroolscore",
                    "chwngswllwngscore", "eatingtskscore",
                    "handwritingscore", "gttngoutbedscore",
                    "wlkngbalancescore"),
        stringr::str_c("mdsupdrs", name),
        name)) %>%
    pivot_wider(id_cols = guid:visit)
  section_1_2_annual_survey <- clinical_records %>%
    dplyr::filter(!stringr::str_detect(guid, "TEST"),
                  assessdate_fall_m12_as != "",
                  redcap_event_name %in% c(
                    "month_36_arm_1",
                    "month_48_arm_1",
                    "month_60_arm_1")) %>%
    dplyr::mutate(visit = case_when(
      redcap_event_name == "month_36_arm_1" ~ "Month 36",
      redcap_event_name == "month_48_arm_1" ~ "Month 48",
      redcap_event_name == "month_60_arm_1" ~ "Month 60")) %>%
    dplyr::select(guid,
                  redcap_event_name,
                  visit,
                  visstatdttm = assessdate_fall_m12_as,
                  qstnnreinfoprovdrt_m12_as:mdsupdrsfreezingscore_m12_as) %>%
    pivot_longer(!guid:visstatdttm) %>%
    drop_na(value) %>%
    mutate(name = if_else(
      str_detect(name, "_m12_as"),
      str_remove(name, "_m12_as"),
      name),
      name = if_else(
        name %in% c("sleepprobscore", "daytmsleepscore",
                    "painothrsensscore", "urnryprobscore",
                    "constipprobscore", "slivadroolscore",
                    "chwngswllwngscore", "eatingtskscore",
                    "handwritingscore", "gttngoutbedscore",
                    "wlkngbalancescore"),
        stringr::str_c("mdsupdrs", name),
        name)) %>%
    pivot_wider(id_cols = guid:visstatdttm)
  scores <- section_1_3_4 %>%
    dplyr::left_join(section_1_2, by = c("guid", "visit")) %>%
    dplyr::bind_rows(section_1_2_annual_survey) %>%
    dplyr::select(createdOn = visstatdttm, everything()) %>%
    tibble::as_tibble(.)
  return(scores)
}

#' Format input for SUPER-PD physician visit
#'
#' @param field_mapping A list mapping DMR to clinical field names
#' @param physician_visit One of "Physician_ON" or "Physician_OFF"
get_super_physician_scores <- function(field_mapping, physician_visit) {
  clinical_records <- fread(synGet(CLINICAL_DATA_SYN_ID)$path) %>%
    tibble::as_tibble(.)
  relevant_fields <- field_mapping %>%
    filter(form_name == "MDS-UPDRS",
           cohort == "super-pd",
           visit == physician_visit) %>%
    drop_na()
  date_col <- case_when(
    physician_visit == "Physician_OFF" ~ "time_mdsupdrs",
    physician_visit == "Physician_ON" ~ "time_mdsupdrs_off")
  scores <- clinical_records %>%
    dplyr::filter(guid != "TESTING",
                  !is.na(mdsupdrsoffon)) %>%
    dplyr::select(guid, createdOn = {{ date_col }}, relevant_fields$clinical_variable)
  col_map <- readxl::read_excel(synGet(EXCEL_LOOKUP_SYN_ID)$path) %>%
    dplyr::select(
      field_name = `Variable / Field Name`,
      form_name = `Form Name`,
      ctcc_name = `CTCC Name`)
  scores_with_ctcc_names <- scores %>%
    pivot_longer(-c("guid", "createdOn"), names_to="field_name") %>%
    inner_join(col_map, by = "field_name") %>%
    select(guid, createdOn, ctcc_name, value) %>%
    drop_na() %>%
    distinct(guid, ctcc_name, .keep_all=T) %>%
    mutate(ctcc_name = unlist(purrr::map(ctcc_name, ~ paste0("C_", .)))) %>%
    pivot_wider(id_cols=c("guid", "createdOn"),
                names_from="ctcc_name",
                values_from="value")
  updrs_scores <- scores_with_ctcc_names %>%
    compute_updrs_total_scores(join_cols=c("guid", "createdOn")) %>%
    select(guid, createdOn, UPDRS1, UPDRS2, UPDRS3, UPDRS3R,
           UPDRS4, UPDRSAMB, UPDRSAMR, UPDRSPRO)
  return(updrs_scores)
}

get_fox_data <- function(visit_date_mapping) {
  movement_survey <- as_tibble(fread(synGet(FOX_DATA_SYN_ID)$path))
  section_2 <- movement_survey %>%
    dplyr::mutate(
      visit = get_associated_visit_type(guid, study_date, visit_date_mapping),
      MoveWho = dplyr::case_match(MoveWho, 1 ~ 2, 2 ~ 1, 3 ~ 3),
      study_date = as.character(study_date)) %>%
    dplyr::select(
      guid = guid,
      visit = visit,
      createdOn = study_date,
      mdsupdrschwngswllwngscore = MoveChew,
      mdsupdrsdressingscore = MoveDress,
      mdsupdrseatingtskscore = MoveEat,
      mdsupdrsfreezingscore = MoveFreeze,
      hobbieothractscore = MoveHobby,
      mdsupdrshygienescore = MoveHygiene,
      mdsupdrsslivadroolscore = MoveSaliva,
      mdsupdrsturngbedscore = MoveSleep,
      mdsupdrsspeechscore = MoveSpeech,
      mdsupdrstremorscore = MoveTremor,
      mdsupdrsgttngoutbedscore = MoveUp,
      mdsupdrswlkngbalancescore = MoveWalk,
      qstnnreinfoprovdrt = MoveWho,
      mdsupdrshandwritingscore = MoveWrite)
  return(section_2)
}

get_associated_visit_type <- function(guids, study_dates, visit_date_mapping) {
  visit_type <- purrr::map2(
    guids,
    study_dates,
    function(fox_guid, study_date) {
      visits <- visit_date_mapping %>%
        filter(guid == fox_guid,
               redcap_event_name %in% list(
                 "baseline_arm_1",
                 "month_12_arm_1",
                 "month_24_arm_1",
                 "month_36_arm_1",
                 "month_48_arm_1",
                 "month_60_arm_1"
               )) %>%
        arrange(visstatdttm)
      visit_type <- purrr::map2(visits$visstatdttm,
                                visits$VisitTypPDBP,
                                function(visit_date, visit_type) {
          #' This map returns a list containing visits which took place after
          #' this survey had been completed. Since we are checking each visit
          #' in chronological order, the first non-NA element in this list is
          #' the visit at which this survey was completed.
          if (is.na(study_date)) {
            # If no start date provided, assume earliest visit type
            return(visit_type)
          } else if (study_date <= visit_date) {
            return(visit_type)
          }
          return(NA_character_)
        }) %>%
        purrr::discard(is.na) %>%
        dplyr::first()
      if (is.null(visit_type)) {
        # In this case the survey date was later than the
        # most recent visit. We associate this with the next visit, or the last possible visit.
        last_available_visit_type <- dplyr::last(visits$VisitTypPDBP)
        visit_type <- case_when(
          last_available_visit_type == "Baseline" ~ "Month 12",
          last_available_visit_type == "Month 12" ~ "Month 24",
          last_available_visit_type == "Month 24" ~ "Month 36",
          last_available_visit_type == "Month 36" ~ "Month 48",
          last_available_visit_type == "Month 48" ~ "Month 60",
          last_available_visit_type == "Month 60" ~ "Month 60",
        )
      }
      return(visit_type)
  })
  return(unlist(visit_type))
}

#' Get VisitTypPDBP field from cohort and event name
get_visit_type <- function(cohort, redcap_event_name) {
  visit_type <- case_when(
    cohort == "at-home-pd" && str_starts(redcap_event_name, "baseline") ~ "Baseline",
    cohort == "at-home-pd" && str_starts(redcap_event_name, "screening") ~ "Baseline",
    cohort == "at-home-pd" && str_starts(redcap_event_name, "month_12") ~ "Month 12",
    cohort == "at-home-pd" && str_starts(redcap_event_name, "month_24") ~ "Month 24",
    cohort == "at-home-pd" && str_starts(redcap_event_name, "month_36") ~ "Month 36",
    cohort == "at-home-pd" && str_starts(redcap_event_name, "month_48") ~ "Month 48",
    cohort == "at-home-pd" && str_starts(redcap_event_name, "month_60") ~ "Month 60",
    cohort == "super-pd" ~ "Baseline",
    TRUE ~ "Baseline") # Logs and Premature Withdrawal
  return(visit_type)
}

#' Build a mapping of GUID to visit type and date
build_visit_date_mapping <- function() {
  cohorts  <- readr::read_csv(synGet("syn24173690")$path) %>%
    rename(study_cohort = cohort)
  clinical <- readr::read_csv(synGet(CLINICAL_DATA_SYN_ID)$path) %>%
    inner_join(cohorts, by = "guid") %>%
    filter(!str_detect(guid, "TEST"))
  visit_date_mapping <- clinical %>%
    filter(!is.na(visstatdttm) | !is.na(annual_survey_timestamp),
           redcap_event_name != "Screening (Arm 1: Arm 1)") %>%
    mutate(visstatdttm = case_when(
        is.na(visstatdttm) ~ annual_survey_timestamp,
        !is.na(visstatdttm) ~ visstatdttm
    )) %>%
    filter(!is.na(visstatdttm)) %>%
    distinct(guid, redcap_event_name, visstatdttm, study_cohort)
  visit_date_mapping$VisitTypPDBP <- unlist(
        purrr::map2(
          visit_date_mapping$study_cohort,
          visit_date_mapping$redcap_event_name,
          get_visit_type))
  return(visit_date_mapping)
}

store_to_synapse <- function(df, fname, parent, used) {
  readr::write_csv(df, fname)
  f <- synapser::File(fname, parent)
  synapser::synStore(f, activity = Activity(
    "compute updrs scores",
    used = used,
    executed = GIT_URL))
  unlink(fname)
}

main <- function(){
  field_mapping <- readr::read_csv("to_dmr/resources/dmr_to_clinical_field_mapping.csv")
  visit_date_mapping <- build_visit_date_mapping()
  #' compute updrs score and write to .csv
  super_physician_on <- get_super_physician_scores(field_mapping, "Physician_ON")
  store_to_synapse(df = super_physician_on,
                   fname = "mdsupdrs_physician_on_med_scores.csv",
                   parent = SUPER_PARENT_ID,
                   used = c(CLINICAL_DATA_SYN_ID, EXCEL_LOOKUP_SYN_ID))
  super_physician_off <- get_super_physician_scores(field_mapping, "Physician_OFF")
  store_to_synapse(df = super_physician_off,
                   fname = "mdsupdrs_physician_off_med_scores.csv",
                   parent = SUPER_PARENT_ID,
                   used = c(CLINICAL_DATA_SYN_ID, EXCEL_LOOKUP_SYN_ID))
  clinical_data <- get_ahpd_clinical_data() %>%
    select(
      -mdsupdrsttlhrawkdysknum,
      -mdsupdrsttlhrdysknum,
      -mdsupdrsprcntdyskval,
      -ttlhrawkoffstatenu,
      -ttlhroffdemndystni,
      -ttlhroffwdystnianu,
      -prcntoffdystniaval,
      -mdsupdrsttlhroffnum,
      -mdsupdrsprcntoffval)
  fox_data <- get_fox_data(visit_date_mapping)
  medication <- clinical_data %>%
      dplyr::select(guid, visit,
                    C_PDSTATE = ptclinstateprknsnm,
                    C_ONLDOPA = mdsupdrsptntuseldopaind,
                    C_LSTDSMIN = mdsupdrslstldopadosetm)
  scores <- clinical_data %>%
      map_column_names("ahpd") %>%
      run_updrs_scoring(
        join_cols = c("guid", "visit", "createdOn")) %>%
      distinct(guid, visit, createdOn, .keep_all=TRUE)
  fox_scores <- fox_data %>%
      map_column_names("ahpd") %>%
      run_updrs_scoring(
        join_cols = c("guid", "visit", "createdOn")) %>%
      distinct(guid, visit, createdOn, .keep_all=TRUE)

  combined_output <- list(
    medication = medication,
    scores = scores) %>%
    purrr::reduce(dplyr::inner_join, by = c("guid", "visit"))

  #' store result to synapse
  store_to_synapse(combined_output, OUTPUT_FILENAME, SYNAPSE_PARENT_ID, used = CLINICAL_DATA_SYN_ID)
  store_to_synapse(fox_scores, FOX_OUTPUT_FILENAME, FOX_PARENT_ID, used = FOX_DATA_SYN_ID)
}

main()

