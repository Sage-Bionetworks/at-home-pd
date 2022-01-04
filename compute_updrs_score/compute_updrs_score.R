#################################################
#' Script to generate AT-HOME-PD UPDRS scores
#' to synapse
#' @author: Aryton Tediarjo
##################################################
library(tidyverse)
library(tidyr)
library(data.table)
library(synapser)
source("compute_updrs_score/scoring_utils.R")
synLogin()

#' Synapse reference
CLINICAL_DATA_SYN_ID <- "syn18565500"
EXCEL_LOOKUP_SYN_ID <- "syn25050918"
SYNAPSE_PARENT_ID <- "syn16809549"
SUPER_PARENT_ID <- "syn25715411"
OUTPUT_FILENAME <- "computed_ahpd_updrs_scores.tsv"
DMR_FIELD_MAPPING <- "syn25056102"

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
  visit_info <- clinical_records %>%
    dplyr::filter(!stringr::str_detect(guid, "TEST"),
                  redcap_event_name %in% c(
                    "baseline_arm_1", "month_12_arm_1", "month_24_arm_1"),
                  redcap_repeat_instrument != "reportable_event") %>%
    dplyr::mutate(visit = case_when(
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
  score_info <- clinical_records %>%
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
    pivot_wider(guid:visit)
  scores <- visit_info %>%
    dplyr::left_join(score_info, by = c("guid", "visit")) %>%
    dplyr::select(createdOn = visstatdttm, everything()) %>%
    tibble::as_tibble(.)
  return(scores)
}

#' Format input for SUPER-PD physician visit
#'
#' @param field_mapping The file referenced above by DMR_FIELD_MAPPING
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

store_to_synapse <- function(df, fname, parent) {
  readr::write_csv(df, fname)
  f <- synapser::File(fname, parent)
  synapser::synStore(f, activity = Activity(
    "compute updrs scores",
    used = c(CLINICAL_DATA_SYN_ID,
             EXCEL_LOOKUP_SYN_ID),
    executed = GIT_URL))
  unlink(fname)
}

main <- function(){
  field_mapping <- readr::read_csv(synapser::synGet(DMR_FIELD_MAPPING)$path)
  #' compute updrs score and write to .tsv
  super_physician_on <- get_super_physician_scores(field_mapping, "Physician_ON")
  store_to_synapse(df = super_physician_on,
                   fname = "mdsupdrs_physician_on_med_scores.csv",
                   parent = SUPER_PARENT_ID)
  super_physician_off <- get_super_physician_scores(field_mapping, "Physician_OFF")
  store_to_synapse(df = super_physician_off,
                   fname = "mdsupdrs_physician_off_med_scores.csv",
                   parent = SUPER_PARENT_ID)
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
  list(
    medication = clinical_data %>% 
      dplyr::select(guid, visit, 
                    C_PDSTATE = ptclinstateprknsnm,
                    C_ONLDOPA = mdsupdrsptntuseldopaind,
                    C_LSTDSMIN = mdsupdrslstldopadosetm),
    scores = clinical_data %>%
      map_column_names("ahpd") %>%
      run_updrs_scoring(
        join_cols = c("guid", "visit", "createdOn")) %>%
      distinct(guid, visit, createdOn, .keep_all=TRUE)) %>%
    purrr::reduce(dplyr::inner_join, by = c("guid", "visit")) %>%
    readr::write_tsv(OUTPUT_FILENAME)
  
  #' store result to synapse
  f <- File(OUTPUT_FILENAME, SYNAPSE_PARENT_ID)
  synStore(f, activity = Activity(
    "compute updrs scores",
    used = c(CLINICAL_DATA_SYN_ID,
             EXCEL_LOOKUP_SYN_ID),
    executed = GIT_URL))
  unlink(OUTPUT_FILENAME)
}

main()

