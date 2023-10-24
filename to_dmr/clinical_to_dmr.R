#' Given AT-HOME PD clinical data from REDCap, translate each record into a
#' format which conforms to the appropriate PDBP DMR schema. This also
#' incorporates Fox Insight data.

library(synapser)
library(tidyverse)
library(lubridate)
library(glue)
library(logger)

AHPD_MDSUPDRS_SCORES <- "syn25050919"
SUPER_OFF_MDSUPDRS_SCORES <- "syn25165548"
SUPER_ON_MDSUPDRS_SCORES <- "syn25165546"
FOX_MDSUPDRS_SCORES <- "syn51104164"
CLINICAL_DATA_DICTIONARY <- "syn21740194"
CLINICAL_DATA <- "syn17051543"
FIELD_MAPPING <- "syn25056102"
VALUE_MAPPING <- "syn25155671"
FORM_TO_DATETIME_MAPPING <- "syn25575806"
FORM_TO_FORM_MAPPING <- "syn25758571"
FORM_TO_HEADER_MAPPING <- "syn26118733"
AHPD_CLINICAL_FORMS <- list(
  "inclusion_exclusion", "participant_demographics", "moca",
  "modified_schwab_and_england_adl", "concomitant_medication_log",
  "reportable_event", "conclusion")
AHPD_MDSUPDRS_FORMS <- list(
  "prebaseline_survey", "previsit_survey", "mdsupdrs", "annual_survey")
# Not included in the SUPER forms is the participant_mdsupdrs_survey, since
# we parse it as part of the mdsupdrs_physician_exam.
SUPER_CLINICAL_FORMS <- list(
  "concomitant_medications", "inclusion_exclusion_spd", "substudy_moca",
  "substudy_mdsupdrs_part_iii", "moca_spd", "demographics_spd",
  "mdsupdrs_physician_exam", "pdq39", "reportable_event", "conclusion")
OUTPUT_PARENT <- "syn25759357"

#' Read a CSV file from Synapse as a tibble
read_synapse_csv <- function(synapse_id) {
  f <- synapser::synGet(synapse_id)
  df <- readr::read_csv(f$path)
  return(df)
}

read_synapse_json <- function(synapse_id) {
  f <- synapser::synGet(synapse_id)
  j <- jsonlite::read_json(f$path)
  return(j)
}

#' Determine if this clinical record conforms to a clinical form
#'
#' Most forms have a required field where the form collection date is recorded.
#' If this field is empty, we know this record does not contain the form info.
#' The exception to this case is the clinical form `concomitant_medication_log`
#' @param record A one-row dataframe from the clinical data containing a single record
#' @param visit_date_col The mandatory field where the form collection date
#' is recorded
#' @return boolean
has_form_info <- function(record, form_name, visit_date_col) {
  if (is.null(visit_date_col)) {
    repeat_instrument  <- record[["redcap_repeat_instrument"]]
    if (!is.na(repeat_instrument) &&
        form_name == "concomitant_medication_log" &&
        repeat_instrument == "Concomitant Medication Log") {
      return(TRUE)
    } else if (!is.na(repeat_instrument) &&
        form_name == "reportable_event" &&
        repeat_instrument == "Reportable Event") {
      return(TRUE)
    } else { # This isn't a form relevant to this study
      return(FALSE)
    }
  } else if (!is.na(record[[visit_date_col]])) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Parse mandatory DMR schema fields from a clinical record
get_universal_fields <- function(record, visit_date_col, dob_mapping,
                                 cohort, redcap_event_name) {
  if (is.null(visit_date_col)) { # concomitant_medication_log
    visit_date <- NA_character_
    age_in_months <- NA_character_
  } else {
    visit_date <- lubridate::format_ISO8601(lubridate::as_datetime(record[[visit_date_col]]))
    age_in_months <- get_age_in_months(
          current_date = visit_date,
          dob_mapping = dob_mapping,
          participant_id = record[["guid"]])
  }
  universal_fields <- tibble::tibble(
      GUID = record[["guid"]],
      VisitDate = visit_date,
      SiteName = "AT-HOME-PD_University of Rochester",
      AgeVal = age_in_months,
      VisitTypPDBP = get_visit_type(cohort, redcap_event_name)) %>%
  mutate_all(as.character)
  universal_fields[["AgeYrs"]] <- as.character(as.integer(
    universal_fields[["AgeVal"]]) %/% 12)
  universal_fields[["AgeRemaindrMonths"]] <- as.character(as.integer(
    universal_fields[["AgeVal"]]) %% 12)
  return(universal_fields)
}

#' Compute current age in months
#'
#' The following participants do not provide the necessary fields to compute
#' their current age:
#' 1 NIHGE434YJLLA
#' 2 NIHMR963TPLWF
#' 3 PDGA-859-GTZ
#' 4 PDRJ-711-MKF
#'
#' @param current_date The date this record's data was collected
#' @param dob_mapping A mapping created by `build_dob_mapping`
#' @param participant_id The GUID of this participant
#' @return The participant's current age in months (as character type)
get_age_in_months <- function(current_date, dob_mapping, participant_id) {
  if (is.na(current_date)) {
      return(NA_character_)
  }
  dob_record <- dob_mapping %>%
    filter(guid == participant_id)
  if (nrow(dob_record) == 0) {
    return(NA_character_)
  } else {
    dob <- dob_record[["dob"]]
  }
  lifespan <- dob %--% current_date
  age_in_months <- lubridate::time_length(lifespan, unit = "month")
  rounded_age_in_months <- round(age_in_months, digits = 1)
  return(as.character(rounded_age_in_months))
}

#' Get VisitTypPDBP field from cohort and event name
get_visit_type <- function(cohort, redcap_event_name) {
  visit_type <- case_when(
    cohort == "at-home-pd" && str_detect(redcap_event_name, "Baseline") ~ "Baseline",
    cohort == "at-home-pd" && str_detect(redcap_event_name, "Screening") ~ "Baseline",
    cohort == "at-home-pd" && str_detect(redcap_event_name, "Month 12") ~ "12 months",
    cohort == "at-home-pd" && str_detect(redcap_event_name, "Month 24") ~ "24 months",
    cohort == "at-home-pd" && str_detect(redcap_event_name, "Month 36") ~ "36 months",
    cohort == "at-home-pd" && str_detect(redcap_event_name, "Month 48") ~ "48 months",
    cohort == "at-home-pd" && str_detect(redcap_event_name, "Month 60") ~ "60 months",
    cohort == "at-home-pd" && str_detect(redcap_event_name, "Reconsent") ~ "36 months",
    cohort == "at-home-pd" && str_detect(redcap_event_name, "Premature") ~ "60 months",
    cohort == "super-pd" ~ "Baseline",
    TRUE ~ "Baseline") # Logs and Premature Withdrawal
  return(visit_type)
}

#' Build a mapping of GUID to date of birth
build_dob_mapping <- function(clinical) {
  ahpd_dob <- clinical %>%
      filter(!is.na(dob)) %>%
      select(guid, dob)
  super_dob <- clinical %>%
      filter(!is.na(demo_age), !is.na(visitdate)) %>%
      mutate(dob = (visitdate - lubridate::dyears(demo_age))) %>%
      select(guid, dob)
  all_dob <- bind_rows(ahpd_dob, super_dob)
  return(all_dob)
}

#' Build a mapping of GUID to visit type and date
build_visit_date_mapping <- function(clinical) {
  visit_date_mapping <- clinical %>%
    filter(!is.na(visstatdttm) | !is.na(annual_survey_timestamp),
           redcap_event_name != "Screening (Arm 1: Arm 1)") %>%
    mutate(visstatdttm = case_when(
        is.na(visstatdttm) ~ annual_survey_timestamp,
        !is.na(visstatdttm) ~ visstatdttm
    )) %>%
    filter(!is.na(visstatdttm)) %>%
    distinct(guid, redcap_event_name, visstatdttm, study_cohort)
  visit_date_mapping <- visit_date_mapping %>%
    mutate(VisitTypPDBP = unlist(
             purrr::map2(
               visit_date_mapping$study_cohort,
               visit_date_mapping$redcap_event_name,
               get_visit_type)))
  return(visit_date_mapping)
}

#' Parse concomitant medication record for DMR, AT-HOME PD Cohort
#'
#' Participants can have multiple concomitant medication records for a single
#' visit, hence this type of record needs to be handled specially.
#' This function extracts fields MedctnPriorConcomRteTyp, MedctnPriorConcomName,
#' MedctnPriorConcomDoseMsrTxt, MedctnPriorConcomDoseUoM, MedctnPriorConcomFreqTxt,
#' and MedctnPriorConcomPD. Other fields specific to this form are set to NA_character_
#' Super PD data is stored in a different format, see
#' parse_concomitant_medication_record_spd
#' @param record A one-row dataframe from the clinical data containing a single record
#' @param visit_date_mapping A mapping of GUID to visit type and date
#' @param value_mapping The value mapping. A list with heirarchy form > field identifier.
#' @return A tibble with fields specific to concomitant medications
parse_concomitant_medication_record_ahpd <- function(record, visit_date_mapping,
                                                     value_mapping, dob_mapping) {
  med_map <- value_mapping[["concomitant_medication_log"]]
  pd_meds <- value_map(med_map, "pd_meds", record$pd_meds)
  pd_med_other <- value_map(med_map, "pd_med_other", record$pd_med_other)
  is_pd_med <- !is.na(record$pd_med_yn) && record$pd_med_yn == "Yes"
  recognized_pd_med <- is_pd_med && !(is.na(pd_meds) && is.na(pd_med_other))
  if (is_pd_med && !recognized_pd_med) {
    log_warn(glue("Neither {record$pd_meds} nor {record$pd_med_other}",
                  " is a recognized PD medication. We will include this medication",
                  " under the field `Other Medications.MedctnPriorConcomName`"))
  }
  is_other_medication <- !recognized_pd_med
  route <- value_map(med_map, "route", record$route)
  route_oth <- value_map(med_map, "route_oth", record$route_oth)
  freq <- value_map(med_map, "freq", record$freq)
  freq_oth <- value_map(med_map, "freq_oth", record$freq_oth, other_specify=TRUE)
  units  <- value_map(med_map, "units", record$units)
  units_oth  <- value_map(med_map, "units", record$units_oth, other_specify=TRUE)
  indication <- value_map(med_map, "indication", record$indication, as_is=TRUE)
  visits <- visit_date_mapping %>%
    filter(guid == record$guid,
           redcap_event_name %in% list(
             "Baseline (Arm 1: Arm 1)",
             "Month 12 (Arm 1: Arm 1)",
             "Month 24 (Arm 1: Arm 1)",
             "Month 36 (Arm 1: Arm 1)",
             "Month 48 (Arm 1: Arm 1)",
             "Month 60 (Arm 1: Arm 1)"
           )) %>%
    arrange(visstatdttm)
  visit_type <- purrr::map2(visits$visstatdttm,
                            visits$VisitTypPDBP,
                            function(visit_date, visit_type) {
      #' This map returns a list containing visits which took place after
      #' this medication had been started. Since we are checking each visit
      #' in chronological order, the first non-NA element in this list is
      #' the visit at which this medication was reported.
      if (is.na(record$start_dt)) {
        # If no start date provided, assume earliest visit type
        log_warn(glue("No medication start date provided for medication",
                      " either {pd_meds} or {pd_med_other} or {record$non_pd_med}.",
                      " Assuming VisitTypPDBP = 'Baseline'."))
        return(visit_type)
      } else if (record$start_dt <= visit_date) {
        return(visit_type)
      }
      return(NA_character_)
    }) %>%
    purrr::discard(is.na) %>%
    dplyr::first()
  if (is.null(visit_type)) {
    # In this case the start date of the medication was later than the
    # most recent visit. We assume this was reported at the most recent visit.
    visit_type <- dplyr::last(visits$VisitTypPDBP)
    log_warn(glue("Medication start date for medication",
                  " either {pd_meds} or {pd_med_other} or {record$non_pd_med}",
                  " is later than most recent visit.",
                  " Assuming VisitTypPDBP = '{visit_type}'."))
  }
  visit_date <- visits[visits$VisitTypPDBP == visit_type,][["visstatdttm"]]
  age_in_months <- get_age_in_months(
        current_date = visit_date,
        dob_mapping = dob_mapping,
        participant_id = record[["guid"]])
  dmr_record <- tibble(
      VisitTypPDBP = visit_type,
      VisitDate = as.character(visit_date),
      AgeVal = as.character(age_in_months),
      `Parkinson's Disease Medications.MedctnPriorConcomRteTyp` = case_when(
        recognized_pd_med && !is.na(route) ~ route,
        recognized_pd_med && !is.na(route_oth) ~ route_oth,
        TRUE ~ NA_character_),
      `Other Medications.MedctnPriorConcomRteTyp` = case_when(
        is_other_medication && !is.na(route) ~ route,
        is_other_medication && !is.na(route_oth) ~ route_oth,
        TRUE ~ NA_character_),
      `Parkinson's Disease Medications.MedctnPriorConcomPD` = case_when(
        recognized_pd_med && !is.na(pd_meds) ~ pd_meds,
        recognized_pd_med && !is.na(pd_med_other) ~ pd_med_other,
      TRUE ~ NA_character_),
      `Other Medications.MedctnPriorConcomName` = case_when(
        is_other_medication && !is_pd_med ~ record$non_pd_med,
        is_other_medication && record$pd_meds != "Other" ~ record$pd_meds,
        is_other_medication && !is.na(record$pd_med_other) ~ record$pd_med_other,
        TRUE ~ NA_character_),
      `Other Medications.MedctnPriorConcomIndTxt` = tolower(indication),
      `Parkinson's Disease Medications.MedctnPriorConcomFreqTxt` = case_when(
        recognized_pd_med && !is.na(freq) ~ freq,
        recognized_pd_med && !is.na(freq_oth) ~ freq_oth,
        TRUE ~ NA_character_),
      `Other Medications.MedctnPriorConcomFreqTxt` = case_when(
        is_other_medication && !is.na(freq) ~ freq,
        is_other_medication && !is.na(freq_oth) ~ freq_oth,
        TRUE ~ NA_character_),
      `Parkinson's Disease Medications.MedctnPriorConcomDoseMsrTxt` =
        case_when(
          recognized_pd_med ~ record$dose,
          TRUE ~ NA_character_),
      `Other Medications.MedctnPriorConcomDoseMsrTxt` =
        case_when(
          is_other_medication ~ record$dose,
          TRUE ~ NA_character_),
      `Parkinson's Disease Medications.MedctnPriorConcomDoseUoM` = case_when(
        recognized_pd_med && record$units == "other" ~ ifelse(!is.na(units_oth), units_oth, NA_character_),
        recognized_pd_med && !is.na(units) ~ units,
        recognized_pd_med && !is.na(units_oth) ~ units_oth,
        TRUE ~ NA_character_),
      `Other Medications.MedctnPriorConcomDoseUoM` = case_when(
        is_other_medication && record$units == "other" ~ ifelse(!is.na(units_oth), units_oth, NA_character_),
        is_other_medication && !is.na(units) ~ units,
        is_other_medication && !is.na(units_oth) ~ units_oth,
        TRUE ~ NA_character_),
      `Parkinson's Disease Medications.MedctnPriorConcomMinsLstDose` = NA_character_,
      `Parkinson's Disease Medications.MedctnPriorConcomHrsLstDose` = NA_character_)
  dmr_record[["AgeYrs"]] <- as.character(as.integer(
    dmr_record[["AgeVal"]]) %/% 12)
  dmr_record[["AgeRemaindrMonths"]] <- as.character(as.integer(
    dmr_record[["AgeVal"]]) %% 12)
  return(dmr_record)
}

#' Remove duplicates from DMR-formatted concomitant medications
#'
#' For medications which do not have a start date (AHPD cohort) and were all
#' given VisitTypPDBP = 'Baseline', or for repeat medications (SUPER-PD cohort),
#' keep only unique medication names per each visit, keeping the most recently
#' reported medication.
remove_duplicate_concomitant_medication_record <- function(concomitant_medications) {
  # reverse order so that we keep most recent entry
  concomitant_medications_pd <- concomitant_medications %>%
    filter(!is.na(`Parkinson's Disease Medications.MedctnPriorConcomPD`)) %>%
    arrange(desc(row_number())) %>%
    distinct(GUID,
             VisitTypPDBP,
             `Parkinson's Disease Medications.MedctnPriorConcomPD`,
             .keep_all=TRUE)
  concomitant_medications_other <- concomitant_medications %>%
    filter(!is.na(`Other Medications.MedctnPriorConcomName`)) %>%
    arrange(desc(row_number())) %>%
    distinct(GUID,
             VisitTypPDBP,
             `Other Medications.MedctnPriorConcomName`,
             .keep_all=TRUE)
  concomitant_medications_clean <- bind_rows(
     concomitant_medications_pd,
     concomitant_medications_other)
  return(concomitant_medications_clean)
}

#' Parse concomitant medication record for DMR, SUPER PD Cohort
#'
#' Participants can have multiple concomitant medication records for a single
#' visit, hence this type of record needs to be handled specially.
#' This function extracts fields MedctnPriorConcomRteTyp, MedctnPriorConcomName,
#' MedctnPriorConcomDoseMsrTxt, MedctnPriorConcomDoseUoM, MedctnPriorConcomFreqTxt,
#' and MedctnPriorConcomPD. Other fields specific to this form are set to NA_character_.
#' In contrast to AT-HOME PD, SUPER-PD medications are recorded as a single record
#' for each participant. Fields for a specific medication can be grouped by their
#' numeric suffix, e.g., conmed_1, conmed_dose_amt_1, conmed_dose_unit_1, etc.
#' AT-HOME PD data is stored in a different format, see
#' parse_concomitant_medication_record_ahpd
#' @param record A one-row dataframe from the clinical data containing a single record
#' @param value_mapping The value mapping. A list with heirarchy form > field identifier.
#' @return A tibble with fields specific to concomitant medications
parse_concomitant_medication_record_spd <- function(record, value_mapping) {
  med_map <- value_mapping[["concomitant_medications"]]
  num_medications <- as.integer(record$conmed_num)
  dmr_records <- purrr::map_dfr(1:num_medications, function(n) {
    if (!hasName(record, glue("conmed_{n}")) || is.na(record[[glue("conmed_{n}")]])) {
      return(tibble())
    }
    conmed <- value_map(med_map, "conmed", record[[glue("conmed_{n}")]])
    is_pd_med <- !is.na(conmed) || record[[glue("conmed_pd_{n}")]] == "Yes"
    recognized_pd_med <- is_pd_med && !is.na(conmed)
    if (is_pd_med && !recognized_pd_med) {
      conmed_raw <- record[[glue("conmed_{n}")]]
      log_warn(glue("{conmed_raw} is not a recognized PD medication.",
                    " We will include this medication",
                    " under the field `Other Medications.MedctnPriorConcomName`"))
    }
    is_other_medication <- !recognized_pd_med
    conmed_dose_amt <- as.character(record[[glue("conmed_dose_amt_{n}")]])
    conmed_dose_unit <- value_map(
      med_map, "conmed_dose_unit", record[[glue("conmed_dose_unit_{n}")]])
    conmed_dose_unit_other <- value_map(
      med_map, "conmed_dose_unit_other",
      record[[glue("conmed_dose_unit_other_{n}")]], other_specify=TRUE)
    conmed_dose_frequency <- value_map(
      med_map, "conmed_dose_frequency", record[[glue("conmed_dose_frequency_{n}")]])
    conmed_dose_frequency_other <- value_map(
      med_map, "conmed_dose_frequency_other",
      record[[glue("conmed_dose_frequency_other_{n}")]], other_specify=TRUE)
    conmed_dose_route <- value_map(
      med_map, "conmed_dose_route", record[[glue("conmed_dose_route_{n}")]])
    conmed_dose_route_other <- value_map(
      med_map, "conmed_dose_route_other", record[[glue("conmed_dose_route_other_{n}")]])
    conmed_indication <- value_map(
      med_map, "conmed_indication", record[[glue("conmed_indication_{n}")]], as_is=TRUE)
    dmr_record <- tibble(
      `Parkinson's Disease Medications.MedctnPriorConcomRteTyp` = case_when(
        recognized_pd_med && !is.na(conmed_dose_route) ~ conmed_dose_route,
        recognized_pd_med && !is.na(conmed_dose_route_other) ~ conmed_dose_route_other,
        TRUE ~ NA_character_),
      `Other Medications.MedctnPriorConcomRteTyp` = case_when(
        is_other_medication && !is.na(conmed_dose_route) ~ conmed_dose_route,
        is_other_medication && !is.na(conmed_dose_route_other) ~ conmed_dose_route_other,
        TRUE ~ NA_character_),
      `Parkinson's Disease Medications.MedctnPriorConcomPD` = case_when(
        recognized_pd_med ~ conmed,
        TRUE ~ NA_character_),
      `Other Medications.MedctnPriorConcomName` = case_when(
        is_other_medication ~ record[[glue("conmed_{n}")]],
        TRUE ~ NA_character_),
      `Other Medications.MedctnPriorConcomIndTxt` = tolower(conmed_indication),
      `Parkinson's Disease Medications.MedctnPriorConcomFreqTxt` = case_when(
        recognized_pd_med && !is.na(conmed_dose_frequency) ~ conmed_dose_frequency,
        recognized_pd_med && !is.na(conmed_dose_frequency_other) ~ conmed_dose_frequency_other,
        TRUE ~ NA_character_),
      `Other Medications.MedctnPriorConcomFreqTxt` = case_when(
        is_other_medication && !is.na(conmed_dose_frequency) ~ conmed_dose_frequency,
        is_other_medication && !is.na(conmed_dose_frequency_other) ~ conmed_dose_frequency_other,
        TRUE ~ NA_character_),
      `Parkinson's Disease Medications.MedctnPriorConcomDoseMsrTxt` =
        case_when(
          recognized_pd_med ~ conmed_dose_amt,
          TRUE ~ NA_character_),
      `Other Medications.MedctnPriorConcomDoseMsrTxt` =
        case_when(
          is_other_medication ~ conmed_dose_amt,
          TRUE ~ NA_character_),
      `Parkinson's Disease Medications.MedctnPriorConcomDoseUoM` = case_when(
        recognized_pd_med && record[[glue("conmed_dose_unit_{n}")]] == "other" ~ ifelse(
              !is.na(conmed_dose_unit_other), conmed_dose_unit_other, NA_character_),
        recognized_pd_med && !is.na(conmed_dose_unit) ~ conmed_dose_unit,
        recognized_pd_med && !is.na(conmed_dose_unit_other) ~ conmed_dose_unit_other,
        TRUE ~ NA_character_),
      `Other Medications.MedctnPriorConcomDoseUoM` = case_when(
        is_other_medication && record[[glue("conmed_dose_unit_{n}")]] == "other" ~ ifelse(
              !is.na(conmed_dose_unit_other), conmed_dose_unit_other, NA_character_),
        is_other_medication && !is.na(conmed_dose_unit) ~ conmed_dose_unit,
        is_other_medication && !is.na(conmed_dose_unit_other) ~ conmed_dose_unit_other,
        TRUE ~ NA_character_),
      `Parkinson's Disease Medications.MedctnPriorConcomMinsLstDose` = NA_character_,
      `Parkinson's Disease Medications.MedctnPriorConcomHrsLstDose` = NA_character_)
    return(dmr_record)
  })
  return(dmr_records)
}

#' Parse MDS-UPDRS record for AT-HOME PD cohort
#'
#' @param record A one-row dataframe containing a single record.
#' This record consists of the combined fields of the clinical forms `mdsupdrs` and
#' either `prebaseline_survey' or `previsit_survey`
#' @param field_mapping The DMR to clinical field mapping.
#' @param value_mapping The value mapping. A list with heirarchy (form) > (field identifier).
#' @param scores A dataframe containing MDS-UPDRS section scores
#' @param date_field The field in `record` which contains the date of the survey
#' @return A tibble with fields specific to the PDBP_MDS-UPDRS form.
parse_mdsupdrs_ahpd <- function(record, field_mapping, value_mapping,
                                scores=NULL, date_field=NULL) {
  this_visit <- case_when(
      record$visittyppdbpmds %in% c("Baseline", "Screening") ~ "Baseline",
      record$visittyppdbpmds == "Month 12" ~ "12 months",
      record$visittyppdbpmds == "Month 24" ~ "24 months",
      record$visittyppdbpmds == "Month 36" ~ "36 months",
      record$visittyppdbpmds == "Month 48" ~ "36 months",
      record$visittyppdbpmds == "Month 60" ~ "36 months")
  this_field_mapping <- field_mapping %>%
    filter(form_name == "MDS-UPDRS",
           cohort == "at-home-pd",
           visit == this_visit)
  dmr_record <- purrr::map2_dfr(
      this_field_mapping$dmr_variable,
      this_field_mapping$clinical_variable,
      function(dmr_variable, clinical_variable) {
        value <- tibble(
          name = dmr_variable,
          value = value_map_updrs(value_mapping, record, clinical_variable))
        return(value)
      })
  if (nrow(dmr_record) == 0) {
    return(tibble())
  }
  if (!is.null(scores) && !is.null(date_field)) {
    section_scores <- mdsupdrs_section_scores(
        scores = scores,
        participant_id = record$guid,
        date_of_exam = as.Date(record[[date_field]]))
  } else {
    log_warn(glue("Did not find MDS-UPDRS section scores for",
                  " {record$guid} at visit {this_visit} and",
                  " date {record[[date_field]]}"))
    section_scores <- tibble(
        name = character(),
        value = character())
  }
  dmr_record <- dmr_record %>%
    anti_join(section_scores, by = "name") %>%
    bind_rows(section_scores) %>%
    pivot_wider(names_from="name", values_from="value")
  return(dmr_record)
}

#' Parse MDS-UPDRS record for SUPER PD cohort
#'
#' SUPER-PD cohort participants completed a full MDS-UPDRS examination in both OFF
#' and ON states at a physician visit (form participant_mdsupdrs_survey sections 1&2,
#' form mdsupdrs_physician_exam sections 1&3&4) and a separate, stand-alone section 3
#' (form substudy_mdsupdrs_part_iii). The physician visit data is contained in a single
#' record and the stand-alone section 3 data is in another record.
#'
#' @param record A one-row dataframe from the clinical data containing a single record
#' @param field_mapping The DMR to clinical field mapping.
#' @param value_mapping The value mapping. A list with
#' heirarchy (form) > (field identifier) > (values).
#' @param scores A list of dataframes containing MDS-UPDRS section scores for the
#' ON and OFF medication exam. The list names should be
#' "Physician_ON" and "Physician_OFF".
#' @return A tibble with two records, if this is a physician administered exam (usually --
#' one participant only completed an ON exam), or a tibble with a single record if this is the
#' stand-alone section 3 exam. All fields are specific to the DMR's MDS-UPDRS form.
parse_mdsupdrs_spd <- function(record, field_mapping, value_mapping, scores) {
  event <- case_when(
      !is.na(record$time_mdsupdrs) ~ "physician",
      !is.na(record$mdsupdrs_sub_dttm) ~ "substudy")
  date_of_exam <- case_when(
      !is.na(record$time_mdsupdrs) ~ as.Date(record$time_mdsupdrs),
      !is.na(record$mdsupdrs_sub_dttm) ~ as.Date(record$mdsupdrs_sub_dttm))
  if (event == "physician") {
    # All MDS-UPDRS sections, for both OFF and ON med states
    dmr_records <- purrr::map_dfr(c("No", "Yes"), function(med_status) {
        # This one specific user had their OFF and ON exams on separate days
        date_of_exam <- case_when(
            record$guid == "NIHRH896FYACQ" && med_status == "No" ~ as.Date(record$time_mdsupdrs),
            record$guid == "NIHRH896FYACQ" && med_status == "Yes" ~ as.Date(record$time_mdsupdrs_off),
            TRUE ~ date_of_exam)
        this_visit <- case_when(med_status == "No" ~ "Physician_OFF",
                                med_status == "Yes" ~ "Physician_ON")
        this_field_mapping <- field_mapping %>%
          filter(form_name == "MDS-UPDRS",
                 cohort == "super-pd",
                 visit == this_visit)
        dmr_record <- purrr::map2_dfr(
            this_field_mapping$dmr_variable,
            this_field_mapping$clinical_variable,
            function(dmr_variable, clinical_variable) {
              value <- tibble(
                name = dmr_variable,
                value = value_map_updrs(value_mapping, record, clinical_variable))
            })
        section_scores <- mdsupdrs_section_scores(
            scores = scores[[this_visit]],
            participant_id = record$guid,
            date_of_exam = date_of_exam)
        dmr_record <- dmr_record %>%
          anti_join(section_scores, by = "name") %>%
          bind_rows(section_scores) %>%
          pivot_wider(names_from="name", values_from="value")
        return(dmr_record)
      })
    return(dmr_records)
  } else if (event == "substudy") {
    # Section 3
    this_field_mapping <- field_mapping %>%
      filter(form_name == "MDS-UPDRS",
             cohort == "super-pd",
             visit == "Baseline")
    dmr_record <- purrr::map2_dfr(
        this_field_mapping$dmr_variable,
        this_field_mapping$clinical_variable,
        function(dmr_variable, clinical_variable) {
          value <- tibble(
            name = dmr_variable,
            value = value_map_updrs(value_mapping, record, clinical_variable))
          return(value)
        })
    dmr_record_score <- dmr_record %>%
      mutate(value = as.integer(value)) %>%
      filter(!is.na(value),
             value < 10)
    mean_value <- round(mean(dmr_record_score$value))
    # 4 missing rigidity questions and 1 missing postural stability question
    dmr_record_score <- dmr_record_score %>%
      summarize(score = sum(value) + 5*mean_value)
    dmr_record <- dmr_record %>%
      pivot_wider(names_from="name", values_from="value") %>%
      mutate(MDSUPDRS_PartIIIScore = dmr_record_score$score) %>%
      mutate_all(as.character)
    return(dmr_record)
  }
}

#' Parse MDS-UPDRS record for Fox Movement survey
#'
#' @param record A one-row dataframe containing a single record.
#' This record consists of the combined fields of the clinical forms `mdsupdrs` and
#' either `prebaseline_survey' or `previsit_survey`
#' @param field_mapping The DMR to clinical field mapping.
#' @param scores A dataframe containing MDS-UPDRS section scores.
#' @param date_field The field in `record` which contains the date of the survey.
#' @return A tibble with fields specific to the PDBP_MDS-UPDRS form.
parse_mdsupdrs_fox <- function(record, field_mapping, scores=NULL, date_field=NULL) {
  this_field_mapping <- field_mapping %>%
    filter(form_name == "MDS-UPDRS",
           cohort == "at-home-pd",
           visit == "Baseline")
  dmr_record <- purrr::map2_dfr(
      this_field_mapping$dmr_variable,
      this_field_mapping$clinical_variable,
      function(dmr_variable, clinical_variable) {
        value <- tibble(
          name = dmr_variable,
          value = ifelse(
              hasName(record, clinical_variable),
              as.character(record[[clinical_variable]]),
              NA_character_))
        return(value)
  })
  if (nrow(dmr_record) == 0) {
    return(tibble())
  }
  section_scores <- mdsupdrs_section_scores(
      scores = scores,
      participant_id = record[["guid"]],
      date_of_exam = as.Date(record[[date_field]]))
  dmr_record <- dmr_record %>%
    anti_join(section_scores, by = "name") %>%
    bind_rows(section_scores) %>%
    pivot_wider(names_from="name", values_from="value")
  return(dmr_record)
}

#' Parse MOCA scores for the AT-HOME PD cohort
#'
#' MOCA exam was administered to AT-HOME PD cohort at baseline, 12, and 24
#' month visits. The same form was used at each visit, hence the fields
#' are the same across visits.
#'
#' @param record A one-row dataframe from the clinical data containing
#' a single record
#' @param field_mapping The DMR to clinical field mapping.
#' @param value_mapping The value mapping. A list with
#' heirarchy (form) > (field identifier) > (values).
#' @return A tibble with fields specific to the DMR MoCA form
parse_moca_ahpd <- function(record, field_mapping, value_mapping) {
  this_field_mapping <- field_mapping %>%
    filter(form_name == "MoCA",
           cohort == "at-home-pd",
           visit == "Baseline")
  this_value_mapping <- value_mapping[["moca"]]
  dmr_record <- purrr::map2_dfr(
    this_field_mapping$dmr_variable,
    this_field_mapping$clinical_variable,
    function(dmr_variable, clinical_variable) {
      this_key <- ifelse(is.na(clinical_variable),
                         NA_character_,
                         as.character(record[[clinical_variable]]))
      value <- tibble(
        name = dmr_variable,
        value = value_map(
          mapping = this_value_mapping,
          field = clinical_variable,
          key = this_key,
          as_is = TRUE))
      return(value)
  })
  dmr_record  <- dmr_record %>%
    pivot_wider(names_from = name, values_from = value)
  return(dmr_record)
}

#' Parse MOCA scores for the SUPER PD cohort
#'
#' MOCA exam was administered to SUPER PD cohort at the physician visit,
#' (form moca_spd) as well as to the substudy cohort
#' (Arm 2: Sub-study, form substudy_moca).
#'
#' @param record A one-row dataframe from the clinical data containing
#' a single record
#' @param field_mapping The DMR to clinical field mapping.
#' @param value_mapping The value mapping. A list with
#' heirarchy (form) > (field identifier) > (values).
#' @return A tibble with fields specific to the DMR MoCA form
parse_moca_spd <- function(record, field_mapping, value_mapping) {
  this_event <- case_when(
      !is.na(record$moca_dttm_v2) ~ "Baseline",
      !is.na(record$visitdate) ~ "Physician_ON")
  if (this_event == "Physician_ON") { # form moca_spd
    this_field_mapping <- field_mapping %>%
      filter(form_name == "MoCA",
             cohort == "super-pd",
             visit == "Physician_ON")
    this_value_mapping <- value_mapping[["moca_spd"]]
    dmr_record <- purrr::map2_dfr(
        this_field_mapping$dmr_variable,
        this_field_mapping$clinical_variable,
        function(dmr_variable, clinical_variable) {
          this_key <- ifelse(is.na(clinical_variable),
                             NA_character_,
                             as.character(record[[clinical_variable]]))
          value <- tibble(
            name = dmr_variable,
            value = value_map(
              mapping = this_value_mapping,
              field = clinical_variable,
              key = this_key,
              as_is = TRUE))
        })
    dmr_record <- dmr_record %>%
      pivot_wider(names_from = name, values_from = value) %>%
      mutate(
        MOCA_VisuospatialExec = as.character(sum(
          record$moca_1_spd, record$moca_2_spd, record$moca_3a_spd,
          record$moca_3b_spd, record$moca_3c_spd)),
        MOCA_Naming = as.character(sum(
          record$moca_4a_spd, record$moca_4b_spd, record$moca_4c_spd)),
        MOCA_DelydRecall = as.character(sum(
          record$moca_9a_spd, record$moca_9b_spd, record$moca_9c_spd,
          record$moca_9d_spd, record$moca_9e_spd)),
        MOCA_Orient = as.character(sum(
          record$moca_10a_spd, record$moca_10b_spd, record$moca_10c_spd,
          record$moca_10d_spd, record$moca_10e_spd, record$moca_10f_spd)),
        MOCA_Digits = as.character(sum(
          record$moca_5a_spd, record$moca_5b_spd)))
  } else if (this_event == "Baseline") { # form substudy_moca
    this_field_mapping <- field_mapping %>%
      filter(form_name == "MoCA",
             cohort == "super-pd",
             visit == "Baseline")
    this_value_mapping <- value_mapping[["substudy_moca"]]
    dmr_record <- purrr::map2_dfr(
      this_field_mapping$dmr_variable,
      this_field_mapping$clinical_variable,
      function(dmr_variable, clinical_variable) {
        this_key <- ifelse(is.na(clinical_variable),
                           NA_character_,
                           as.character(record[[clinical_variable]]))
        value <- tibble(
          name = dmr_variable,
          value = value_map(
            mapping = this_value_mapping,
            field = clinical_variable,
            key = this_key,
            as_is = TRUE))
        return(value)
    })
    dmr_record  <- dmr_record %>%
      pivot_wider(names_from = name, values_from = value)
  }
  return(dmr_record)
}

#' Parse PDQ-39 form for SUPER-PD cohort
#'
#' The SUPER-PD cohort took a single PDQ-39 exam at the physician visit.
#' Neither AT-HOME PD cohort or the sub-study cohort filled out this form.
#'
#' @param
#' @param record A one-row dataframe from the clinical data containing
#' a single record
#' @param field_mapping The DMR to clinical field mapping.
#' @return A tibble with fields specific to the DMR PDQ-39 form
parse_pdq_39 <- function(record, field_mapping) {
  this_field_mapping <- field_mapping %>%
    filter(form_name == "PDQ-39",
           cohort == "super-pd",
           visit == "Baseline")
  dmr_record <- purrr::map2_dfr(
    this_field_mapping$dmr_variable,
    this_field_mapping$clinical_variable,
    function(dmr_variable, clinical_variable) {
      if (dmr_variable == "PDQ_39_LackOfSuprtPrtnr" && record[["spouse_check"]] == "No") {
        this_value <- "No spouse or partner"
      } else {
        this_value <- ifelse(is.na(clinical_variable),
                             NA_character_,
                             str_extract(record[[clinical_variable]], "\\d"))
      }
      value <- tibble(
        name = dmr_variable,
        value = this_value)
      return(value)
  })
  sections <- tibble(
    clinical_variable = c("leisure", "housework", "bags", "mile", "yards",
                          "home", "public", "accompany", "fall", "confined",
                          "wash", "dress", "shoes", "writing", "cut", "spill",
                          "depressed", "isolated", "weepy", "angry", "anxious",
                          "future", "conceal", "avoid", "embarrassed", "worried",
                          "relationships", "spouse", "friends", "sleep",
                          "concentration", "memory_e3a8f9", "dreams", "speech",
                          "communicate", "ignored", "cramps", "aches", "temp"),
    section = c("mobility", "mobility", "mobility", "mobility", "mobility",
                "mobility", "mobility", "mobility", "mobility", "mobility",
                "adl", "adl", "adl", "adl", "adl", "adl",
                "emotional", "emotional", "emotional", "emotional",
                "emotional", "emotional",
                "stigma", "stigma", "stigma", "stigma",
                "social", "social", "social",
                "cognition", "cognition", "cognition", "cognition",
                "communication", "communication", "communication",
                "discomfort", "discomfort", "discomfort"),
  max_score = 4)
  dmr_sections <- sections %>%
    inner_join(this_field_mapping) %>%
    suppressMessages() %>%
    select(dmr_variable, section)
  section_scores <- sections %>%
    inner_join(this_field_mapping) %>%
    suppressMessages() %>%
    inner_join(dmr_record, by = c("dmr_variable" = "name")) %>%
    suppressMessages() %>%
    filter(value != "No spouse or partner") %>% # don't count this question
    group_by(section) %>%
    summarize(section_score = sum(as.integer(value)),
              max_section_score = sum(max_score),
              dimension_score = section_score / max_section_score * 100) %>%
    select(name = section, value = dimension_score) %>%
    mutate(name = case_when(
      name == "adl" ~ "PDQ_39_TotalScore_ADL",
      name == "cognition" ~ "PDQ_39_TotalScore_CogImpairmnt",
      name == "communication" ~ "PDQ_39_TotalScore_Communcation",
      name == "discomfort" ~ "PDQ_39_TotalScore_BodDiscomfrt",
      name == "emotional" ~ "PDQ_39_TotalScore_Emotional",
      name == "mobility" ~ "PDQ_39_TotalScore_Mobility",
      name == "social" ~ "PDQ_39_TotalScore_SocialSuprt",
      name == "stigma" ~ "PDQ_39_TotalScore_Stigma"),
      value = as.character(signif(value, 3)))
  dmr_record <- dmr_record %>%
    mutate(value = case_when(
      value == "0" ~ "Never",
      value == "1" ~ "Occasionally",
      value == "2" ~ "Sometimes",
      value == "3" ~ "Often",
      value == "4" ~ "Always or cannot do at all")) %>%
    anti_join(section_scores, by = "name") %>%
    suppressMessages() %>%
    bind_rows(section_scores)
  dmr_record  <- dmr_record %>%
    pivot_wider(names_from = name, values_from = value)
  return(dmr_record)
}

#' Parse Demographic info for the AT-HOME PD cohort
#'
#' Get data for DMR fields EmplmtStatus, EduLvlUSATypPDBP,
#' and EthnUSACat. These are the fields which are specific to the Demographics
#' DMR form. RaceExpndCatPDBP field info is not included in the clinical data.
#'
#' @param record A one-row dataframe from the clinical data containing
#' a single record
#' @param field_mapping The DMR to clinical field mapping.
#' @param value_mapping The value mapping. A list with
#' heirarchy (form) > (field identifier) > (values).
#' @return A tibble with fields specific to the DMR Demographics form
parse_demographics_ahpd <- function(record, field_mapping, value_mapping) {
  this_field_mapping <- field_mapping %>%
    filter(form_name == "Demographics",
           cohort == "at-home-pd",
           visit == "Baseline")
  this_value_mapping <- value_mapping[["participant_demographics"]]
  dmr_record <- purrr::map2_dfr(
      this_field_mapping$dmr_variable,
      this_field_mapping$clinical_variable,
      function(dmr_variable, clinical_variable) {
        this_key <- ifelse(is.na(clinical_variable),
                           NA_character_,
                           as.character(record[[clinical_variable]]))
        value <- tibble(
          name = dmr_variable,
          value = value_map(
            mapping = this_value_mapping,
            field = clinical_variable,
            key = this_key,
            as_is = TRUE))
      })
  dmr_record <- dmr_record %>%
    pivot_wider(names_from = name, values_from = value)
  return(dmr_record)
}

#' Parse Demographic info for the SUPER PD cohort
#'
#' Get fields unique to DMR Demographics form.
#'
#' @param record A one-row dataframe from the clinical data containing
#' a single record
#' @param field_mapping The DMR to clinical field mapping.
#' @param value_mapping The value mapping. A list with
#' heirarchy (form) > (field identifier) > (values).
#' @return A tibble with fields specific to the DMR Demographics form
parse_demographics_spd <- function(record, field_mapping, value_mapping) {
  this_field_mapping <- field_mapping %>%
    filter(form_name == "Demographics",
           cohort == "super-pd",
           visit == "Baseline")
  this_value_mapping <- value_mapping[["demographics_spd"]]
  dmr_record <- purrr::map2_dfr(
      this_field_mapping$dmr_variable,
      this_field_mapping$clinical_variable,
      function(dmr_variable, clinical_variable) {
        this_key <- ifelse(is.na(clinical_variable),
                           NA_character_,
                           as.character(record[[clinical_variable]]))
        value <- tibble(
          name = dmr_variable,
          value = value_map(
            mapping = this_value_mapping,
            field = clinical_variable,
            key = this_key,
            as_is = TRUE))
      })
  dmr_record <- dmr_record %>%
    pivot_wider(names_from = name, values_from = value)
  return(dmr_record)
}

#' Parse reportable events for either AT-HOME PD or SUPER PD
#'
#' This function parses responses to the clinical `reportable_event` form.
#' These events include adverse events, which map to the DMR's AdverseEvents
#' form.
#'
#' @param record A one-row dataframe from the clinical data containing
#' a single record
#' @param field_mapping The DMR to clinical field mapping.
#' @param value_mapping The value mapping. A list with
#' heirarchy (form) > (field identifier) > (values).
#' @return A tibble with fields specific to the DMR AdverseEvents form
parse_reportable_event <- function(record, field_mapping,
                                   value_mapping, visit_date_mapping, dob_mapping) {
  this_field_mapping <- field_mapping %>%
    filter(form_name == "AdverseEvents",
           cohort == "at-home-pd",
           visit == "Baseline")
  this_value_mapping <- value_mapping[["reportable_event"]]
  dmr_record <- purrr::map2_dfr(
      this_field_mapping$dmr_variable,
      this_field_mapping$clinical_variable,
      function(dmr_variable, clinical_variable) {
        this_key <- ifelse(is.na(clinical_variable),
                           NA_character_,
                           as.character(record[[clinical_variable]]))
        value <- tibble(
          name = dmr_variable,
          value = value_map(
            mapping = this_value_mapping,
            field = clinical_variable,
            key = this_key,
            as_is = TRUE))
      })
  dmr_record <- dmr_record %>%
    pivot_wider(names_from = name, values_from = value)
  dmr_record$AdvrsEvntDuringStudyInd <- "Yes"
  serious_event_codes <- c(
      record$evntcode___4, record$evntcode___5, record$evntcode___6,
      record$evntcode___8, record$evntcode___9, record$evntcode___10,
      record$evntcode___11)
  if (!all(is.na(serious_event_codes)) && any(serious_event_codes == "Checked")) {
    dmr_record$SeriousAdvrsEvntInd <- "Yes"
  } else {
    dmr_record$SeriousAdvrsEvntInd <- "No"
  }
  # Fatal/Death Life-Threatening/Disabling Severe Moderate Mild
  # case_when short-circuits -- so more severe events take precedence
  # In most cases, we don't know the exact severity grade. The five severity
  # grades are from the Common Terminology Criteria for Adverse Events v4.0 (CTCAE)
  dmr_record$AdvrsEvntSeverScale <- case_when(
   record$evntcode___8  == "Checked" ~ "Fatal/Death",
   record$evntcode___4  == "Checked" ~ "Life-Threatening/Disabling",
   record$evntcode___5  == "Checked" ~ "Severe",
   record$evntcode___6  == "Checked" ~ "Severe",
   dmr_record$SeriousAdvrsEvntInd == "Yes" ~ NA_character_)
  # Set Visit Date and Age fields based on visit type
  visits <- visit_date_mapping %>%
    filter(guid == record$guid)
  visit_date <- visits[visits$redcap_event_name == record$redcap_event_name,][["visstatdttm"]]
  visit_type <- visits[visits$redcap_event_name == record$redcap_event_name,][["VisitTypPDBP"]]
  age_in_months <- get_age_in_months(
        current_date = visit_date,
        dob_mapping = dob_mapping,
        participant_id = record[["guid"]])
  dmr_record[["VisitTypPDBP"]] <- visit_type
  dmr_record[["VisitDate"]] <- as.character(visit_date)
  dmr_record[["AgeVal"]] <- as.character(age_in_months)
  dmr_record[["AgeYrs"]] <- as.character(as.integer(
    dmr_record[["AgeVal"]]) %/% 12)
  dmr_record[["AgeRemaindrMonths"]] <- as.character(as.integer(
    dmr_record[["AgeVal"]]) %% 12)
  return(dmr_record)
}

#' Parse conclusion form for either AT-HOME PD or SUPER PD
#'
#' This function parses responses to the clinical `conclusion` form.
#' These can map to the DMR's EarlyTerminationQuest form.
#'
#' @param record A one-row dataframe from the clinical data containing
#' a single record
#' @param field_mapping The DMR to clinical field mapping.
#' @param value_mapping The value mapping. A list with
#' heirarchy (form) > (field identifier) > (values).
#' @return A tibble with fields specific to the DMR EarlyTerminationQuest form
parse_conclusion <- function(record, field_mapping, value_mapping) {
  if (record$subj_status == glue("Subject discontinued participation ",
                                "before the planned study conclusion")) {
    this_field_mapping <- field_mapping %>%
      filter(form_name == "EarlyTerminationQuest",
             cohort == "at-home-pd", # same form for both AHPD and SUPER
             visit == "Baseline")
    this_value_mapping <- value_mapping[["conclusion"]]
    dmr_record <- purrr::map2_dfr(
        this_field_mapping$dmr_variable,
        this_field_mapping$clinical_variable,
        function(dmr_variable, clinical_variable) {
          this_key <- ifelse(is.na(clinical_variable),
                             NA_character_,
                             as.character(record[[clinical_variable]]))
          value <- tibble(
            name = dmr_variable,
            value = value_map(
              mapping = this_value_mapping,
              field = clinical_variable,
              key = this_key,
              as_is = FALSE))
        })
    dmr_record <- dmr_record %>%
      pivot_wider(names_from = name, values_from = value)
  } else {
    return(tibble())
  }
  return(dmr_record)
}

#' Parse inclusion_exclusion form for AT-HOME PD
#'
#' This function parses responses to the clinical `inclusion_exclusion` form.
#' These map to the DMR's InclExclCriteria form.
#'
#' @param record A one-row dataframe from the clinical data containing
#' a single record
#' @param field_mapping The DMR to clinical field mapping.
#' @param value_mapping The value mapping. A list with
#' heirarchy (form) > (field identifier) > (values).
#' @return A tibble with fields specific to the DMR AdverseEvents form
parse_inclusion_exclusion_ahpd <- function(record, field_mapping, value_mapping) {
  this_field_mapping <- field_mapping %>%
    filter(form_name == "InclExclCriteria",
           cohort == "at-home-pd",
           visit == "Baseline")
  this_value_mapping <- value_mapping[["inclusion_exclusion"]]
  dmr_record <- purrr::map2_dfr(
      this_field_mapping$dmr_variable,
      this_field_mapping$clinical_variable,
      function(dmr_variable, clinical_variable) {
        this_key <- ifelse(is.na(clinical_variable),
                           NA_character_,
                           as.character(record[[clinical_variable]]))
        value <- tibble(
          name = dmr_variable,
          value = value_map(
            mapping = this_value_mapping,
            field = clinical_variable,
            key = this_key,
            as_is = TRUE))
      })
  dmr_record <- dmr_record %>%
    pivot_wider(names_from = name, values_from = value)
  dmr_record[["SubjectCaseInd"]] = "Yes"
  dmr_record[["SubjectCntrlInd"]] = "No"
  dmr_record[["PDBPInclusnXclusn_InclusnCase"]] <- glue(
      "Clinically diagnosed with Parkinson's Disease (or other ",
      "Neurodegenerative disease appropriate for study protocol)")
  return(dmr_record)
}

#' Parse inclusion_exclusion form for SUPER PD
#'
#' This function parses responses to the clinical `inclusion_exclusion_spd` form.
#' These can map to the DMR's InclExclCriteria form.
#'
#' @param record A one-row dataframe from the clinical data containing
#' a single record
#' @param field_mapping The DMR to clinical field mapping.
#' @param value_mapping The value mapping. A list with
#' heirarchy (form) > (field identifier) > (values).
#' @return A tibble with fields specific to the DMR AdverseEvents form
parse_inclusion_exclusion_spd <- function(record, field_mapping, value_mapping) {
  if (record$enrollconfirm == "No" && record$guid == "NIHGE434YJLLA") {
    # One participant in SUPER PD completed consent form but then went AWOL
    # There are no other forms collected from this participant.
    return(tibble())
  }
  this_field_mapping <- field_mapping %>%
    filter(form_name == "InclExclCriteria",
           cohort == "super-pd",
           visit == "Baseline")
  this_value_mapping <- value_mapping[["inclusion_exclusion_spd"]]
  dmr_record <- purrr::map2_dfr(
      this_field_mapping$dmr_variable,
      this_field_mapping$clinical_variable,
      function(dmr_variable, clinical_variable) {
        this_key <- ifelse(is.na(clinical_variable),
                           NA_character_,
                           as.character(record[[clinical_variable]]))
        value <- tibble(
          name = dmr_variable,
          value = value_map(
            mapping = this_value_mapping,
            field = clinical_variable,
            key = this_key,
            as_is = TRUE))
      })
  dmr_record <- dmr_record %>%
    pivot_wider(names_from = name, values_from = value)
  dmr_record[["SubjectCaseInd"]] = "Yes"
  dmr_record[["SubjectCntrlInd"]] = "No"
  dmr_record[["PDBPInclusnXclusn_InclusnCase"]] <- glue(
      "Clinically diagnosed with Parkinson's Disease (or other ",
      "Neurodegenerative disease appropriate for study protocol)")
  return(dmr_record)
}

parse_informed_consent_ahpd <- function(inclex_record, enroll_record) {
  # TODO only first informed consent record for each participant
  # goes to DMR
  dmr_record <- list()
  dmr_record[["EnrldStdyDateTime"]] <- inclex_record[["inexdttm"]]
  dmr_record[["InfConsntObtInd"]] <- inclex_record[["inex5"]]
  dmr_record[["InformConsntObtnDateTime"]] <- inclex_record[["inexdttm"]]
  dmr_record[["EnrldStdyInd"]] <- enroll_record[["enrollconfirm"]]
  dmr_record[["RandomizedDateTime"]] <- NA_character_
  dmr_record[["RandomizedInd"]] <- NA_character_
  dmr_record <- tibble::as_tibble(dmr_record)
  return(dmr_record)
}

parse_informed_consent_super <- function(inclex_record, enroll_record) {
  dmr_record <- list()
  dmr_record[["EnrldStdyDateTime"]] <- inclex_record[["inexdttm_spd"]]
  dmr_record[["InfConsntObtInd"]] <- inclex_record[["inex5_spd"]]
  dmr_record[["InformConsntObtnDateTime"]] <- inclex_record[["inexdttm_spd"]]
  dmr_record[["EnrldStdyInd"]] <- enroll_record[["enrollconfirm"]]
  dmr_record[["RandomizedDateTime"]] <- NA_character_
  dmr_record[["RandomizedInd"]] <- NA_character_
  dmr_record <- tibble::as_tibble(dmr_record)
  return(dmr_record)
}

#' Parse modified_schwab_and_england_adl form
#'
#' This function parses responses to the clinical
#' `modified_schwab_and_england_adl` form. Only AT-HOME PD participants
#' completed this form. These map to the DMR's ModSchwabAndEnglandScale form.
#'
#' @param record A one-row dataframe from the clinical data containing
#' a single record
#' @param field_mapping The DMR to clinical field mapping.
#' @param value_mapping The value mapping. A list with
#' heirarchy (form) > (field identifier) > (values).
#' @return A tibble with fields specific to the DMR ModSchwabAndEnglandScale form
parse_mod_schwab_and_england <- function(record, field_mapping, value_mapping) {
  this_field_mapping <- field_mapping %>%
    filter(form_name == "ModSchwabAndEnglandScale",
           cohort == "at-home-pd",
           visit == "Baseline")
  this_value_mapping <- value_mapping[["modified_schwab_and_england_adl"]]
  dmr_record <- purrr::map2_dfr(
      this_field_mapping$dmr_variable,
      this_field_mapping$clinical_variable,
      function(dmr_variable, clinical_variable) {
        this_key <- ifelse(is.na(clinical_variable),
                           NA_character_,
                           as.character(record[[clinical_variable]]))
        value <- tibble(
          name = dmr_variable,
          value = value_map(
            mapping = this_value_mapping,
            field = clinical_variable,
            key = this_key,
            as_is = TRUE))
      })
  dmr_record <- dmr_record %>%
    pivot_wider(names_from = name, values_from = value)
  return(dmr_record)
}

parse_fox_about_you <- function(dob_mapping, field_mapping, visit_date_mapping) {
  about_you <- read_synapse_csv("syn21670545")
  dmr_records <- purrr::pmap_dfr(about_you, function(...) {
    record <- list(...)
    universal_fields <- get_universal_fields(
        record = record,
        visit_date_col = "study_date",
        dob_mapping = dob_mapping,
        cohort = "at-home-pd",
        redcap_event_name = "Baseline")
    vital_signs <- tibble(
          WgtMeasr = record[["WeightKgs"]],
          HgtMeasr = record[["HeightCm"]])
    vital_signs <- bind_cols(universal_fields, vital_signs)
    return(vital_signs)
  })
  dmr_records <- dmr_records %>%
    mutate(VisitTypPDBP = get_associated_visit_type(GUID, VisitDate, visit_date_mapping))
  all_names <- field_mapping %>%
    filter(form_name == "Vital Signs") %>%
    distinct(dmr_variable)
  for (n in all_names$dmr_variable) {
    if (!hasName(dmr_records, n)) {
      dmr_records[[n]] <- NA_character_
    }
  }
  return(dmr_records)
}

parse_fox_environmental_exposure <- function(dob_mapping, field_mapping, visit_date_mapping) {
  alcohol <- read_synapse_csv("syn21670561")
  smoking <- read_synapse_csv("syn21670551")
  dmr_alcohol <- purrr::pmap_dfr(alcohol, function(...) {
    record <- list(...)
    universal_fields <- get_universal_fields(
        record = record,
        visit_date_col = "study_date",
        dob_mapping = dob_mapping,
        cohort = "at-home-pd",
        redcap_event_name = "Baseline")
    dmr_alcohol <- tibble(
      EverUsedAlcoholInd = eeq_map(record[["alq1"]]),
      AlcUseStrtAgeVal = as.character(record[["al3a_age"]]),
      AlcUseStopAgeVal = as.character(record[["al4_age"]]))
    dmr_alcohol <- bind_cols(universal_fields, dmr_alcohol)
    return(dmr_alcohol)
  })
  dmr_smoking <- purrr::pmap_dfr(smoking, function(...) {
    record <- list(...)
    universal_fields <- get_universal_fields(
        record = record,
        visit_date_col = "study_date",
        dob_mapping = dob_mapping,
        cohort = "at-home-pd",
        redcap_event_name = "Baseline")
    dmr_smoking <- tibble(
      EverUsedTobaccoInd = eeq_map(record[["sm1"]]),
      TobcoUseStrtAgeVal = as.character(record[["sm5astop1"]]),
      TobcoUseStopAgeVal = as.character(record[["sm5astart1"]]))
    dmr_smoking <- bind_cols(universal_fields, dmr_smoking)
    return(dmr_smoking)
  })
  dmr_behavior <- full_join(
      dmr_alcohol, dmr_smoking,
      by = c("GUID", "VisitDate", "SiteName", "AgeVal",
             "VisitTypPDBP", "AgeYrs", "AgeRemaindrMonths")) %>%
    mutate(VisitTypPDBP = get_associated_visit_type(GUID, VisitDate, visit_date_mapping))
  all_names <- field_mapping %>%
    filter(form_name == "BehavioralHistory") %>%
    distinct(dmr_variable)
  for (n in all_names$dmr_variable) {
    if (!hasName(dmr_behavior, n)) {
      dmr_behavior[[n]] <- NA_character_
    }
  }
  return(dmr_behavior)
}

eeq_map <- function(val) {
  return_val <- case_when(
    is.na(val) ~ NA_character_,
    val == 1 ~ "Yes",
    val == 2 ~ "No",
    val == 3 ~ "Unknown")
  return(return_val)
}

get_associated_visit_type <- function(guids, study_dates, visit_date_mapping) {
  visit_type <- purrr::map2(
    guids,
    study_dates,
    function(fox_guid, study_date) {
      visits <- visit_date_mapping %>%
        filter(guid == fox_guid,
               redcap_event_name %in% list(
                 "Baseline (Arm 1: Arm 1)",
                 "Month 12 (Arm 1: Arm 1)",
                 "Month 24 (Arm 1: Arm 1)",
                 "Month 36 (Arm 1: Arm 1)",
                 "Month 48 (Arm 1: Arm 1)",
                 "Month 60 (Arm 1: Arm 1)"
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
          last_available_visit_type == "Baseline" ~ "12 months",
          last_available_visit_type == "12 months" ~ "24 months",
          last_available_visit_type == "24 months" ~ "36 months",
          last_available_visit_type == "36 months" ~ "48 months",
          last_available_visit_type == "48 months" ~ "60 months",
          last_available_visit_type == "60 months" ~ "60 months",
        )
      }
      return(visit_type)
  })
  return(unlist(visit_type))
}

parse_fox_movement <- function(visit_date_mapping) {
  movement_survey <- read_csv(synGet("syn21670565")$path)
  section_2 <- movement_survey %>%
    dplyr::mutate(
      visittyppdbpmds = "Baseline", # Set later in `main`
      MoveWho = dplyr::case_when(
          MoveWho == 1 ~ 2,
          MoveWho == 2 ~ 1,
          MoveWho == 3 ~ 3),
      study_date = as.character(study_date)) %>%
    dplyr::select(
      guid = guid,
      visittyppdbpmds,
      mdsupdrs_dttm = study_date,
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

parse_fox_family_history <- function(dob_mapping, field_mapping, visit_date_mapping) {
  family_history <- read_synapse_csv("syn21670518")
  conditions <- list(
      Parkinson = "Parkinson's disease", Alzheimer = "Alzheimer's disease",
      ALS = "Amytrophic lateral sclerosis", Autism = "Autism",
      Dystonia = "Dystonia", Epilepsy = "Epilepsy", MS = "Multiple sclerosis",
      Stroke = "Stroke", Depression = "Depression",
      Suicide = "Suicide or suicide attempt", Bipolar = "Bi-polar disorder")
  relatives <- list(
      Moth = "Mother", Fath = "Father",
      Child = "Child", Grand = "Grandchild",
      Great = "Great-grandchild", Sib = "Sibling",
      HalfSib = "Half sibling", MatGrMoth = "Maternal grandmother",
      MatGrFath = "Maternal grandfather", MatAunt = "Maternal aunt",
      MatUnc = "Maternal Uncle", MatCous = "Maternal cousin",
      MatNieNep = "Maternal niece/nephew", PatGrMoth = "Paternal grandmother",
      PatGrFath = "Paternal grandfather", PatAunt = "Paternal aunt",
      PatUnc = "Paternal uncle", PatCous = "Paternal cousin",
      PatNieNep = "Paternal niece/nephew", Oth = "Other, specify")
  dmr_family_history <- purrr::pmap_dfr(family_history, function(...) {
    record <- list(...)
    universal_fields <- get_universal_fields(
        record = record,
        visit_date_col = "study_date",
        dob_mapping = dob_mapping,
        cohort = "at-home-pd",
        redcap_event_name = "Baseline")
    condition_records <- purrr::map_dfr(names(conditions), function(cond) {
      indicator_field <- str_c("Fam", cond, "Hx")
      these_relatives <- purrr::map(names(relatives), function(relat) {
        relation_field <- str_c("Fam", cond, relat)
        if (!is.na(record[[relation_field]]) && record[[relation_field]] == 1) {
            return(relatives[[relat]])
          } else {
            return(NULL)
          }
      })
      these_relatives <- purrr::compact(these_relatives)
        if (length(these_relatives) == 0) {
          these_relatives_str <- NA_character_
        } else {
          these_relatives_str <- str_c(these_relatives, collapse=",")
        }
      condition_records <- list()
        condition_records[[glue("{conditions[[cond]]}.FamHistMedclCondTyp")]] <- conditions[[cond]]
        condition_records[[glue("{conditions[[cond]]}.FamHistMedclCondInd")]] <- case_when(
              record[[indicator_field]] == 0 ~ "No",
              record[[indicator_field]] == 1 ~ "Yes",
              record[[indicator_field]] == 2 ~ "Unknown")
        condition_records[[glue("{conditions[[cond]]}.FamHistMedclCondReltvTyp")]] <- these_relatives_str
      return(condition_records)
    })
    dmr_family_history <- bind_cols(universal_fields, condition_records)
    return(dmr_family_history)
  })
  dmr_family_history <- dmr_family_history %>%
    mutate(VisitTypPDBP = get_associated_visit_type(GUID, VisitDate, visit_date_mapping))
  irrelevant_conditions <- c(
    "Ataxia", "Brain aneurysm", "Cancer", "Dementia", "Depression",
    "Diabetes mellitus", "Heart disease", "Hypertension",
    "Memory loss", "Migraines", "Muscle disease", "Schizophrenia",
    "Tourette syndrome", "Additional conditions")
  for (condition in irrelevant_conditions) {
    dmr_family_history[[glue("{condition}.FamHistMedclCondTyp")]] <- NA_character_
    dmr_family_history[[glue("{condition}.FamHistMedclCondInd")]] <- NA_character_
    dmr_family_history[[glue("{condition}.FamHistMedclCondReltvTyp")]] <- NA_character_
  }
  return(dmr_family_history)
}

#' Map a value from clinical to DMR
#'
#' This function can be conservative in the sense that it
#' returns NA in two cases:
#' 1. The value to be mapped is NA (key is NA)
#' 2. There are no mappings for this clinical value (key not in names of mapping[[field]])
#'
#' To instead return the value as-is if there are no value mappings for it, set
#' as_is to TRUE. This function is type-safe -- it always returns character type.
#'
#' @param mapping The value mapping for a single clinical form and its fields.
#' @param field The clinical field to map the values of.
#' @param key The value from the clinical data
#' @param as_is Whether to return the value as is if there are no mappings found.
#' @return NA_character_ or the DMR-compliant value
value_map <- function(mapping, field, key, as_is = FALSE, other_specify=FALSE) {
  if (is.null(key) || is.na(key) || is.null(field) || is.na(field)) {
    return(NA_character_)
  } else if (is.null(mapping) || !(key %in% names(mapping[[field]]))) {
      # if mapping doesn't exist or field has no mappings or key not mapped
      if (as_is) {
        return(as.character(key))
      } else if (other_specify) {
        return("Other, specify")
      } else {
        return(NA_character_)
      }
  }
  value <- mapping[[field]][[key]]
  return(value)
}

#' Map clinical MDS-UPDRS values to DMR permissible values
#'
#' @param mapping The value mapping. A list with heirarchy (form) > (field identifier) > (values).
#' @param record A one-row data frame or named vector
#' @param field The name of the clinical variable in `record` to map
#' @return An integer MDS-UPDRS score
value_map_updrs <- function(mapping, record, field) {
  # IF (
  # This DMR field has no mapping to the clinical data
  #   OR
  # The clinical field is empty in the clinical data)
  if (is.na(field) || is.null(record[[field]])) {
    return(NA_character_)
  }
  dmr_value <- str_extract(record[[field]], "\\d")
  if (is.na(dmr_value)) {
    # This value does not record a numeric response.
    # If this value needs mapping, it will be mapped under one of
    # the clinical MDS-UPDRS forms. Otherwise, we can use the value as-is.
    mdsupdrs_map <- value_map(mapping[["mdsupdrs"]], field, record[[field]])
    substudy_mdsupdrs_part_iii_map <- value_map(
        mapping[["substudy_mdsupdrs_part_iii"]], field, record[[field]])
    participant_mdsupdrs_survey_map <- value_map(
        mapping[["participant_mdsupdrs_survey"]], field, record[[field]])
    mdsupdrs_physician_exam_map <- value_map(
        mapping[["mdsupdrs_physician_exam"]], field, record[[field]])
    dmr_value <- case_when(
        !is.na(mdsupdrs_map) ~ mdsupdrs_map,
        !is.na(substudy_mdsupdrs_part_iii_map) ~ substudy_mdsupdrs_part_iii_map,
        !is.na(participant_mdsupdrs_survey_map) ~ participant_mdsupdrs_survey_map,
        !is.na(mdsupdrs_physician_exam_map) ~ mdsupdrs_physician_exam_map,
        !is.na(record[[field]]) ~ as.character(record[[field]]),
        TRUE ~ NA_character_)
  }
  return(dmr_value)
}

#' Get section scores for a specific participant and date
#'
#' @param scores A dataframe with MDS-UPDRS scores indexed by guid and createdOn
#' @param participant_id The participant's GUID
#' @param date_of_exam The date or datetime of the exam
#' @return row-wise section scores for convenient anti-join/row-bind use
mdsupdrs_section_scores <- function(scores, participant_id, date_of_exam) {
  score <- scores %>%
    mutate(createdOnDate = as.Date(createdOn)) %>%
    filter(guid == participant_id,
           createdOnDate == as.Date(date_of_exam),
           !is.na(UPDRS1) | !is.na(UPDRS2) | !is.na(UPDRS3) | !is.na(UPDRS4)) %>%
    select(all_of(c("UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4")))
  if (nrow(score) == 0) {
    section_scores <- tibble(
        name = character(),
        value = character())
    return(section_scores)
  }
  total_score <- {
    total_score <- score %>%
      pivot_longer(dplyr::everything()) %>%
      summarize(total_score = sum(value))
    as.character(total_score$total_score)
  }
  score <- score %>%
    mutate_all(as.character)
  section_scores <- tibble(
      name = c("MDSUPDRS_PartIScore", "MDSUPDRS_PartIIScore", "MDSUPDRS_PartIIIScore",
              "MDSUPDRS_PartIVScore", "MDSUPDRS_TotalScore"),
      value = c(score[["UPDRS1"]], score[["UPDRS2"]], score[["UPDRS3"]],
                 score[["UPDRS4"]], total_score)) %>%
    mutate_all(as.character)
  return(section_scores)
}

#' Get fields from potentially different events and forms
#'
#' This function is useful when data that would otherwise end up in a single DMR
#' record is spread across different forms and/or events. For example, the
#' MDS-UPDRS data for AHPD participants was partially collected as part of a
#' pre-baseline survey and the rest of the fields were collected at the baseline
#' visit. This function will combine fields from all events and forms
#' into a single record for each participant. This assumes that the two forms
#' do not share any fields (otherwise there would be a conflict).
#'
#' @param clinical The clinical data
#' @param clinical_dic The clinical data dictionary (for looking up form fields)
#' @param event_name The redcap event names to filter upon
#' @param form_name The clinical form names to filter upon
#' @param additional_fields Fields to include in the result which are neither
#' the participant identifier (guid), form specific fields, nor the visit_date_col.
#' @return A dataframe with one row per participant. Each column will be of
#' type character().
get_event_and_form_fields <- function(clinical, clinical_dic, event_name,
                                      form_name, visit_date_col,
                                      additional_fields=NULL, ignore_fields=NULL) {
  form_fields <- clinical_dic %>%
    filter(`Form Name` %in% form_name,
           `Field Type` != "descriptive") %>%
    distinct(`Variable / Field Name`)
  if (!is.null(ignore_fields)) {
    form_fields <- form_fields %>%
      filter(!(`Variable / Field Name` %in% ignore_fields))
  }
  event_records <- clinical %>%
    filter(redcap_event_name %in% event_name,
           !is.na(!!sym(visit_date_col))) %>%
    select(guid,
           {visit_date_col},
           all_of(form_fields[["Variable / Field Name"]]),
           all_of(additional_fields)) %>%
    mutate(across(.fns=as.character)) %>%
    pivot_longer(!guid) %>%
    drop_na() %>%
    pivot_wider(guid)
  return(event_records)
}

#' rename and reorder headers so that they match the templates
modify_header <- function(records, form_name, form_to_header_mapping) {
  header_df <- synGet(form_to_header_mapping[[form_name]])$path %>%
    read_csv(skip=1) %>%
    select(-record)
  actual_headers <- names(header_df)
  replacement_names <- purrr::map(names(records), function(field_name) {
    field_matches <- stringr::str_detect(actual_headers, field_name)
    total_matches <- sum(field_matches)
    if (total_matches == 1) {
      return(actual_headers[field_matches])
    } else if (total_matches == 0) {
      stop(glue("Did not find matching header in form structure template ",
                "for field name {field_name}"))
    } else { # return shortest/simplest match
      possible_fields <- actual_headers[field_matches]
      field_lengths <- unlist(lapply(possible_fields, nchar))
      shorter_match_length <- min(field_lengths)
      shorter_match_index <- match(shorter_match_length, field_lengths)
      return(possible_fields[shorter_match_index])
    }
  })
  names(records) <- unlist(replacement_names)
  return(records[actual_headers])
}

configure_logger <- function(guid, redcap_event_name) {
  identifier <- glue("{guid}/{redcap_event_name}")
  logger_layout <- logger::layout_glue_generator(
      format=paste0("{level} [", identifier, "/{fn}]: {msg}"))
  logger::log_layout(logger_layout)
}

parse_form <- function(record, form, field_mapping, value_mapping, ...) {
  kwargs <- list(...)
  parsed_form <- tibble()
  if (form == "inclusion_exclusion") {
    parsed_form <- parse_inclusion_exclusion_ahpd(
        record = record,
        field_mapping = field_mapping,
        value_mapping = value_mapping)
  } else if (form == "participant_demographics") {
    parsed_form <- parse_demographics_ahpd(
        record = record,
        field_mapping = field_mapping,
        value_mapping = value_mapping)
  } else if (form == "moca") {
    parsed_form <- parse_moca_ahpd(
        record = record,
        field_mapping = field_mapping,
        value_mapping = value_mapping)
  } else if (form == "modified_schwab_and_england_adl") {
    parsed_form <- parse_mod_schwab_and_england(
        record = record,
        field_mapping = field_mapping,
        value_mapping = value_mapping)
  } else if (form == "concomitant_medication_log") {
    parsed_form <- parse_concomitant_medication_record_ahpd(
        record = record,
        visit_date_mapping = kwargs$visit_date_mapping,
        value_mapping = value_mapping,
        dob_mapping = kwargs$dob_mapping)
  } else if (form == "concomitant_medications") {
    parsed_form <- parse_concomitant_medication_record_spd(
        record = record,
        value_mapping = value_mapping)
  } else if (form == "inclusion_exclusion_spd") {
    parsed_form <- parse_inclusion_exclusion_spd(
        record = record,
        field_mapping = field_mapping,
        value_mapping = value_mapping)
  } else if (form == "demographics_spd") {
    parsed_form <- parse_demographics_spd(
        record = record,
        field_mapping = field_mapping,
        value_mapping = value_mapping)
  } else if (form == "substudy_moca") {
    parsed_form <- parse_moca_spd(
        record = record,
        field_mapping = field_mapping,
        value_mapping = value_mapping)
  } else if (form == "moca_spd") {
    parsed_form <- parse_moca_spd(
        record = record,
        field_mapping = field_mapping,
        value_mapping = value_mapping)
  } else if (form == "substudy_mdsupdrs_part_iii") {
    parsed_form <- parse_mdsupdrs_spd(
        record = record,
        field_mapping = field_mapping,
        value_mapping = value_mapping,
        scores = kwargs$scores)
  } else if (form == "mdsupdrs_physician_exam") {
    parsed_form <- parse_mdsupdrs_spd(
        record = record,
        field_mapping = field_mapping,
        value_mapping = value_mapping,
        scores = kwargs$scores)
  } else if (form == "pdq39") {
    parsed_form <- parse_pdq_39(
        record = record,
        field_mapping = field_mapping)
  } else if (form == "reportable_event" &&
             !is.na(record$reportableeventyn) &&
             record$reportableeventyn == "Yes") {
    parsed_form <- parse_reportable_event(
        record = record,
        field_mapping = field_mapping,
        value_mapping = value_mapping,
        visit_date_mapping = kwargs$visit_date_mapping,
        dob_mapping = kwargs$dob_mapping)
  } else if (form == "conclusion") {
    parsed_form <- parse_conclusion(
        record = record,
        field_mapping = field_mapping,
        value_mapping = value_mapping)
  }
  return(parsed_form)
}

#' Conform clinical data with the schemas required by PDBP DMR
#'
#' Clinical data comes from two different cohorts, at-home-pd and super-pd.
#' SUPER PD participants have one screening event, one baseline event, and
#' one physician visit. AT-HOME PD participants have a screening event,
#' a baseline visit, a pre 12 month event, a 12 month event, a pre 24 month
#' event, a 24 month event, and month 36, 48, and 60 events. Concomitant
#' medications are not associated with an event in the redcap data.
main <- function() {
  synapser::synLogin()
  # Load all reference material and clinical data
  cohorts  <- readr::read_csv("resources/guid_to_cohort_mapping.csv")
  clinical_data_dictionary <- readr::read_csv("resources/clinical_data_dictionary.csv")
  field_mapping <- readr::read_csv("resources/dmr_to_clinical_field_mapping.csv")
  value_mapping <- jsonlite::read_json("resources/value_mapping.json")
  form_to_datetime_mapping <- jsonlite::read_json("resources/clinical_form_to_datetime_mapping.json")
  form_to_form_mapping <- jsonlite::read_json("resources/dmr_form_to_clinical_form_mapping.json")
  form_to_header_mapping <- jsonlite::read_json("resources/dmr_form_to_header_mapping.json")
  clinical <- read_synapse_csv(CLINICAL_DATA) %>%
    inner_join(cohorts, by = "guid") %>%
    filter(!str_detect(guid, "TEST"))
  dob_mapping <- build_dob_mapping(clinical)
  visit_date_mapping <- build_visit_date_mapping(clinical)
  ahpd_mdsupdrs_scores <- read_synapse_csv(AHPD_MDSUPDRS_SCORES)
  super_mdsupdrs_scores <- list(
    "Physician_ON" = read_synapse_csv(SUPER_ON_MDSUPDRS_SCORES),
    "Physician_OFF" = read_synapse_csv(SUPER_OFF_MDSUPDRS_SCORES))
  fox_mdsupdrs_scores <- read_synapse_csv(FOX_MDSUPDRS_SCORES)

  # Row-wise map each record conditional on its contents
  dmr_records <- purrr::pmap(clinical, function(...) {
    clinical_record <- list(...)
    configure_logger(guid = clinical_record[["guid"]],
                     redcap_event_name = clinical_record[["redcap_event_name"]])
    cohort <- clinical_record$study_cohort
    # Iterate through list of potential clinical forms we can parse from
    # this record. Given the clinical form, we can extract the universal fields
    # and pass this record into the appropriate `parse_*` function.
    if (cohort == "at-home-pd") {
      dmr_records <- purrr::map(AHPD_CLINICAL_FORMS, function(clinical_form) {
        visit_date_col <- form_to_datetime_mapping[[clinical_form]]
        if (!has_form_info(record = clinical_record,
                           form_name = clinical_form,
                           visit_date_col = visit_date_col)) {
          log_debug(glue("Missing data for clinical form {clinical_form}"))
          return(tibble())
        }
        universal_fields <- get_universal_fields(
          record = clinical_record,
          visit_date_col = visit_date_col,
          dob_mapping = dob_mapping,
          cohort = cohort,
          redcap_event_name = clinical_record[["redcap_event_name"]])
        form_specific_fields <- parse_form(
            record = clinical_record,
            form = clinical_form,
            field_mapping = field_mapping,
            value_mapping = value_mapping,
            visit_date_mapping = visit_date_mapping,
            dob_mapping = dob_mapping)
        if (clinical_form %in% c("concomitant_medication_log", "reportable_event")) {
          universal_fields <- universal_fields %>%
            select(-VisitTypPDBP, -VisitDate, -AgeYrs, -AgeVal, -AgeRemaindrMonths)
        }
        dmr_record <- bind_cols(universal_fields, form_specific_fields)
        return(dmr_record)
      })
      names(dmr_records) <- AHPD_CLINICAL_FORMS
    } else if (cohort == "super-pd") {
      dmr_records <- purrr::map(SUPER_CLINICAL_FORMS, function(clinical_form) {
        visit_date_col <- form_to_datetime_mapping[[clinical_form]]
        if (!has_form_info(record = clinical_record,
                           form_name = clinical_form,
                           visit_date_col = visit_date_col)) {
          log_debug(glue("Missing data for clinical form {clinical_form}"))
          return(tibble())
        }
        universal_fields <- get_universal_fields(
          record = clinical_record,
          visit_date_col = visit_date_col,
          dob_mapping = dob_mapping,
          cohort = cohort,
          redcap_event_name = clinical_record[["redcap_event_name"]])
        form_specific_fields <- parse_form(
            record = clinical_record,
            form = clinical_form,
            field_mapping = field_mapping,
            value_mapping = value_mapping,
            scores = super_mdsupdrs_scores)
        dmr_record <- bind_cols(universal_fields, form_specific_fields)
        return(dmr_record)
      })
      names(dmr_records) <- SUPER_CLINICAL_FORMS
    } else {
      log_error(glue("Unknown cohort encounted for GUID { clinical_record$guid }"))
      stop(glue("Unknown cohort encounted for GUID { clinical_record$guid }"))
    }
    return(dmr_records)
  })

  # Collate all forms of the same type
  all_form_names <- unique(c(AHPD_CLINICAL_FORMS, SUPER_CLINICAL_FORMS))
  clinical_forms <- purrr::map(all_form_names, function(form_name) {
    all_forms_of_this_type <- purrr::map_dfr(dmr_records, function(record) {
      if (has_name(record, form_name)) {
        return(record[[form_name]])
      } else {
        return(tibble())
      }
    })
    return(all_forms_of_this_type)
  })
  names(clinical_forms) <- all_form_names

  # Clean up same-visit repeat concomitant medication records
  clinical_forms[["concomitant_medication_log"]] <-
    remove_duplicate_concomitant_medication_record(
      concomitant_medications = clinical_forms[["concomitant_medication_log"]])
  clinical_forms[["concomitant_medications"]] <-
    remove_duplicate_concomitant_medication_record(
      concomitant_medications = clinical_forms[["concomitant_medications"]])

  #' The DMR's InformedConsent form has its info spread across two of the
  #' clinical forms. Fortunately, it is easy to stitch together.
  #' There is a redcap_repeat_instrument == "Informed Consent Log" for AHPD
  #' participants, but none of these fields are relevant to the DMR's
  #' InformedConsent form.
  ahpd_informed_consent <- {
    inclex_records <- clinical %>%
      filter(!is.na(inexdttm), study_cohort == "at-home-pd") %>%
      distinct(guid, .keep_all=TRUE)
    enroll_records <- clinical %>%
      filter(!is.na(enrollconfirm), study_cohort == "at-home-pd") %>%
      distinct(guid, .keep_all=TRUE)
    purrr::map_dfr(inclex_records[["guid"]], function(guid_) {
      inclex_record <- inclex_records %>%
        filter(guid == guid_)
      enroll_record <- enroll_records %>%
        filter(guid == guid_)
      universal_fields <- get_universal_fields(
        record = inclex_record,
        visit_date_col = "inexdttm",
        dob_mapping = dob_mapping,
        cohort = "at-home-pd",
        redcap_event_name = inclex_record[["redcap_event_name"]])
      dmr_record <- parse_informed_consent_ahpd(
          inclex_record = as.list(inclex_record),
          enroll_record = as.list(enroll_record))
      dmr_record <- bind_cols(universal_fields, dmr_record)
      return(dmr_record)
    })
  }
  super_informed_consent <- {
    inclex_records <- clinical %>%
      filter(!is.na(inexdttm_spd), study_cohort == "super-pd") %>%
      distinct(guid, .keep_all=TRUE)
    enroll_records <- clinical %>%
      filter(!is.na(enrollconfirm), study_cohort == "super-pd") %>%
      distinct(guid, .keep_all=TRUE)
    purrr::map_dfr(inclex_records[["guid"]], function(guid_) {
      inclex_record <- inclex_records %>%
        filter(guid == guid_)
      enroll_record <- enroll_records %>%
        filter(guid == guid_)
      universal_fields <- get_universal_fields(
        record = inclex_record,
        visit_date_col = "inexdttm_spd",
        dob_mapping = dob_mapping,
        cohort = "super-pd",
        redcap_event_name = inclex_record[["redcap_event_name"]])
      dmr_record <- parse_informed_consent_super(
          inclex_record = as.list(inclex_record),
          enroll_record = as.list(enroll_record))
      dmr_record <- bind_cols(universal_fields, dmr_record)
      return(dmr_record)
    })
  }
  # We will add the data we just curated after creating our `dmr_forms`
  # variable so that we can avoid doing something hacky in form_to_form_mapping.json

  # AHPD MDS-UPDRS sections are split up across different visits and forms.
  # We need to combine each into a single record with `get_event_and_form_fields`
  # before parsing.
  # The form pairings are prebaseline_survey/mdsupdrs for the baseline visit,
  # previsit_survey/mdsupdrs for the month 12 and 24 visits, and annual_survey
  # for the month 36, 48, and 60 visits. For these latter visits, only part of
  # Part 1 and all of Part 2 of MDSUPDRS is included in the annual_survey form.
  ahpd_pre_baseline_mdsupdrs <- get_event_and_form_fields(
      clinical = clinical,
      clinical_dic = clinical_data_dictionary,
      event_name = "Baseline Pre Visit Survey (Arm 1: Arm 1)",
      form_name = "prebaseline_survey",
      visit_date_col = "assessdate_fall")
  ahpd_actual_baseline_mdsupdrs <- get_event_and_form_fields(
      clinical = clinical,
      clinical_dic = clinical_data_dictionary,
      event_name = "Baseline (Arm 1: Arm 1)",
      form_name = "mdsupdrs",
      visit_date_col = "mdsupdrs_dttm")
  ahpd_baseline_mdsupdrs <- full_join(
      ahpd_pre_baseline_mdsupdrs, ahpd_actual_baseline_mdsupdrs, by="guid") %>%
      mutate(mdsupdrs_dttm = case_when(
                 is.na(mdsupdrs_dttm) ~ assessdate_fall,
                 !is.na(mdsupdrs_dttm) ~ mdsupdrs_dttm),
             visittyppdbpmds = "Baseline") %>%
      select(-assessdate_fall)
  ahpd_pre_month_12_mdsupdrs <- get_event_and_form_fields(
      clinical = clinical,
      clinical_dic = clinical_data_dictionary,
      event_name = "Month 12 Pre Visit Survey (Arm 1: Arm 1)",
      form_name = "previsit_survey",
      visit_date_col = "assessdate_fall_m12")
  ahpd_actual_month_12_mdsupdrs <- get_event_and_form_fields(
      clinical = clinical,
      clinical_dic = clinical_data_dictionary,
      event_name = "Month 12 (Arm 1: Arm 1)",
      form_name = "mdsupdrs",
      visit_date_col = "mdsupdrs_dttm")
  ahpd_month_12_mdsupdrs <- full_join(
      ahpd_pre_month_12_mdsupdrs, ahpd_actual_month_12_mdsupdrs, by="guid") %>%
      mutate(mdsupdrs_dttm = case_when(
                 is.na(mdsupdrs_dttm) ~ assessdate_fall_m12,
                 !is.na(mdsupdrs_dttm) ~ mdsupdrs_dttm),
             visittyppdbpmds = "Month 12") %>%
      select(-assessdate_fall_m12)
  ahpd_pre_month_24_mdsupdrs <- get_event_and_form_fields(
      clinical = clinical,
      clinical_dic = clinical_data_dictionary,
      event_name = "Month 24 Pre Visit Survey (Arm 1: Arm 1)",
      form_name = "previsit_survey",
      visit_date_col = "assessdate_fall_m12")
  ahpd_actual_month_24_mdsupdrs <- get_event_and_form_fields(
      clinical = clinical,
      clinical_dic = clinical_data_dictionary,
      event_name = "Month 24 (Arm 1: Arm 1)",
      form_name = "mdsupdrs",
      visit_date_col = "mdsupdrs_dttm")
  ahpd_month_24_mdsupdrs <- full_join(
      ahpd_pre_month_24_mdsupdrs, ahpd_actual_month_24_mdsupdrs, by="guid") %>%
      mutate(mdsupdrs_dttm = case_when(
                 is.na(mdsupdrs_dttm) ~ assessdate_fall_m12,
                 !is.na(mdsupdrs_dttm) ~ mdsupdrs_dttm),
             visittyppdbpmds = "Month 24") %>%
      select(-assessdate_fall_m12)
  ahpd_mdsupdrs_baseline_12_24 <- bind_rows(
      ahpd_baseline_mdsupdrs,
      ahpd_month_12_mdsupdrs,
      ahpd_month_24_mdsupdrs)
  dmr_records_ahpd_mdsupdrs_baseline_12_24 <- purrr::pmap_dfr(
        ahpd_mdsupdrs_baseline_12_24, function(...) {
    mdsupdrs_record <- list(...)
    universal_fields <- get_universal_fields(
      record = mdsupdrs_record,
      visit_date_col = "mdsupdrs_dttm",
      dob_mapping = dob_mapping,
      cohort = "at-home-pd",
      redcap_event_name = mdsupdrs_record[["visittyppdbpmds"]])
    form_specific_fields <- parse_mdsupdrs_ahpd(
        record = mdsupdrs_record,
        field_mapping = field_mapping,
        value_mapping = value_mapping,
        scores = ahpd_mdsupdrs_scores,
        date_field = "mdsupdrs_dttm")
    dmr_record <- bind_cols(universal_fields, form_specific_fields)
    return(dmr_record)
  })

  # The processing for 36, 48, and 60 months is a bit different
  superflous_annual_survey_fields <- c(
    "phone_as", "num_type_as", "emer_cont_as", "emerphone_as", "email_as",
    "street1_as", "street2_as", "city_as", "state_as", "zip_as",
    "pgi_ill_m12_as", "pgichange_as", "falllastmonth_m12_as")
  ahpd_month_36_mdsupdrs <- get_event_and_form_fields(
      clinical = clinical,
      clinical_dic = clinical_data_dictionary,
      event_name = "Month 36 (Arm 1: Arm 1)",
      form_name = "annual_survey",
      visit_date_col = "assessdate_fall_m12_as",
      additional_fields = "redcap_event_name",
      ignore_fields = superflous_annual_survey_fields)
  ahpd_month_48_mdsupdrs <- get_event_and_form_fields(
      clinical = clinical,
      clinical_dic = clinical_data_dictionary,
      event_name = "Month 48 (Arm 1: Arm 1)",
      form_name = "annual_survey",
      visit_date_col = "assessdate_fall_m12_as",
      additional_fields = "redcap_event_name",
      ignore_fields = superflous_annual_survey_fields)
  ahpd_month_60_mdsupdrs <- get_event_and_form_fields(
      clinical = clinical,
      clinical_dic = clinical_data_dictionary,
      event_name = "Month 60 (Arm 1: Arm 1)",
      form_name = "annual_survey",
      visit_date_col = "assessdate_fall_m12_as",
      additional_fields = "redcap_event_name",
      ignore_fields = superflous_annual_survey_fields)
  ahpd_mdsupdrs_36_48_60 <- bind_rows(
      ahpd_month_36_mdsupdrs,
      ahpd_month_48_mdsupdrs,
      ahpd_month_60_mdsupdrs) %>%
    mutate(visittyppdbpmds = case_when(
      startsWith(redcap_event_name, "Month 36") ~ "Month 36",
      startsWith(redcap_event_name, "Month 48") ~ "Month 48",
      startsWith(redcap_event_name, "Month 60") ~ "Month 60")) %>%
    select(-redcap_event_name)
  dmr_records_ahpd_mdsupdrs_36_48_60 <- purrr::pmap_dfr(
        ahpd_mdsupdrs_36_48_60, function(...) {
    mdsupdrs_record <- list(...)
    universal_fields <- get_universal_fields(
      record = mdsupdrs_record,
      visit_date_col = "assessdate_fall_m12_as",
      dob_mapping = dob_mapping,
      cohort = "at-home-pd",
      redcap_event_name = mdsupdrs_record[["visittyppdbpmds"]])
    form_specific_fields <- parse_mdsupdrs_ahpd(
        record = mdsupdrs_record,
        field_mapping = field_mapping,
        value_mapping = value_mapping,
        scores = ahpd_mdsupdrs_scores,
        date_field = "assessdate_fall_m12_as")
    dmr_record <- bind_cols(universal_fields, form_specific_fields)
    return(dmr_record)
  })

  # And again, Fox Movement survey has an entirely different format
  fox_movement <- parse_fox_movement(visit_date_mapping = visit_date_mapping)
  dmr_records_fox <- purrr::pmap_dfr(fox_movement, function(...) {
    mdsupdrs_record <- list(...)
    universal_fields <- get_universal_fields(
      record = mdsupdrs_record,
      visit_date_col = "mdsupdrs_dttm",
      dob_mapping = dob_mapping,
      cohort = "at-home-pd",
      redcap_event_name = mdsupdrs_record[["visittyppdbpmds"]]) %>%
      mutate(VisitTypPDBP = get_associated_visit_type(
                GUID, VisitDate, visit_date_mapping))
    form_specific_fields <- parse_mdsupdrs_fox(
        record = mdsupdrs_record,
        field_mapping = field_mapping,
        scores = fox_mdsupdrs_scores,
        date_field = "mdsupdrs_dttm")
    dmr_record <- bind_cols(universal_fields, form_specific_fields)
    return(dmr_record)
  })

  dmr_records_ahpd_mdsupdrs <- bind_rows(
      dmr_records_ahpd_mdsupdrs_baseline_12_24,
      dmr_records_ahpd_mdsupdrs_36_48_60,
      dmr_records_fox)
  clinical_forms[["mdsupdrs"]] <- dmr_records_ahpd_mdsupdrs

  # Combine clinical form data which map to the same DMR form
  dmr_forms <- purrr::map(names(form_to_form_mapping), function(dmr_form_name) {
    matching_forms <- form_to_form_mapping[[dmr_form_name]]
    this_dmr_form <- purrr::map_dfr(matching_forms, ~ clinical_forms[[.]])
    return(this_dmr_form)
  })
  names(dmr_forms) <- names(form_to_form_mapping)

  # Add informed consent data from earlier
  dmr_forms[["InformedConsent"]] <- bind_rows(
      ahpd_informed_consent, super_informed_consent)

  # Add Vital Signs, Behavioriol History and Family History DMR forms,
  # which only have representative fields from Fox
  dmr_forms[["VitalSigns"]] <- parse_fox_about_you(
      dob_mapping = dob_mapping,
      field_mapping = field_mapping,
      visit_date_mapping = visit_date_mapping)
  dmr_forms[["BehavioralHistory"]] <- parse_fox_environmental_exposure(
      dob_mapping = dob_mapping,
      field_mapping = field_mapping,
      visit_date_mapping = visit_date_mapping)
  dmr_forms[["FamilyHistory"]] <- parse_fox_family_history(
      dob_mapping = dob_mapping,
      field_mapping = field_mapping,
      visit_date_mapping = visit_date_mapping)

  # Store to Synapse
  purrr::map(names(dmr_forms), function(dmr_form_name) {
    form_data <- dmr_forms[[dmr_form_name]] %>%
      mutate(GUID = stringr::str_replace_all(GUID, "-", ""))
    form_data <- modify_header(
        records = form_data,
        form_name = dmr_form_name,
        form_to_header_mapping = form_to_header_mapping) %>%
      mutate(record = "x") %>%
      select(record, dplyr::everything())
    dir.create("forms", showWarnings = FALSE)
    fname <- glue("forms/{dmr_form_name}.csv")
    write_csv(form_data, fname)
    system(glue("sed -i '' '1s/^/{dmr_form_name} \\n/' {fname}"))
    f <- synapser::File(fname, parent = OUTPUT_PARENT)
    synStore(f, used = CLINICAL_DATA)
    unlink(fname)
  })
}

main()
