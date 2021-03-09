#' Given AT-HOME PD clinical data from REDCap, translate each record into a
#' format which conforms to the appropriate PDBP DMR schema.

library(dplyr)
library(glue)

#' Parse mandatory DMR schema fields from a clinical record
get_universal_fields <- function(record, visit_date_col, dob_mapping,
                                 cohort, redcap_event_name) {
  universal_fields <- tibble::tibble(
      GUID = record$guid,
      VisitDate = record[[visit_date_col]],
      SiteName = "AT-HOME-PD_University of Rochester",
      AgeVal = get_age_in_months(
          current_date = record[[visit_date_col]],
          dob = dob_mapping[[record$guid]]),
      VisitTypPDBP = get_visit_type(cohort, redcap_event_name))
  universal_fields[["AgeYrs"]] <- universal_fields$AgeVal %/% 12
  universal_fields[["AgeRemaindrMonths"]] <- universal_fields$AgeVal %% 12
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
get_age_in_months <- function(current_date, dob) {
  if (anyNA(c(current_date, dob))) {
      return(NA)
  }
  age_in_months <- lubridate::months(current_date - dob)
  return(age_in_months)
}

#' Get VisitTypPDBP field from cohort and event name
get_visit_type <- function(cohort, redcap_event_name) {
  visit_type <- case_when(
    cohort == "at-home-pd" && str_detect(redcap_event_name, "Baseline") ~ "Baseline",
    cohort == "at-home-pd" && str_detect(redcap_event_name, "Screening") ~ "Baseline",
    cohort == "at-home-pd" && str_detect(redcap_event_name, "Month 12") ~ "Month 12",
    cohort == "at-home-pd" && str_detect(redcap_event_name, "Month 24") ~ "Month 24",
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

#' Parse concomitant medication record for DMR, AT-HOME PD Cohort
#'
#' Participants can have multiple concomitant medication records for a single
#' visit, hence this type of record needs to be handled specially.
#' This function extracts fields MedctnPriorConcomRteTyp, MedctnPriorConcomName,
#' MedctnPriorConcomDoseMsrTxt, MedctnPriorConcomDoseUoM, MedctnPriorConcomFreqTxt,
#' and MedctnPriorConcomPD. Other fields specific to this form are set to NA.
#' Super PD data is stored in a different format, see
#' parse_concomitant_medication_record_spd
#' @param record A one-row dataframe from the clinical data containing a single record
#' @param mapping The value mapping. A list with heirarchy form > field identifier.
#' @return A tibble with fields specific to concomitant medications
parse_concomitant_medication_record_ahpd <- function(record, mapping) {
  med_map <- mapping[["concomitant_medication_log"]]
  is_pd_med <- record$pd_med_yn == "Yes"
  route <- value_map(med_map, "route", record$route)
  route_oth <- value_map(med_map, "route_oth", record$route_oth)
  pd_meds <- value_map(med_map, "pd_meds", record$pd_meds)
  pd_med_other <- value_map(med_map, "pd_med_other", record$pd_med_other)
  freq <- value_map(med_map, "freq", record$freq)
  freq_oth <- value_map(med_map, "freq_oth", record$freq_oth)
  units  <- value_map(med_map, "units", record$units)
  units_oth  <- value_map(med_map, "units", record$units_oth)
  dmr_record <- tibble(
      MedctnPriorConcomRteTyp = case_when(
        !is.na(route) ~ route,
        !is.na(route_oth) ~ route_oth,
        TRUE ~ NA),
      MedctnPriorConcomPD = case_when(
        is_pd_med && !is.na(pd_meds) ~ pd_meds,
        is_pd_med && !is.na(pd_med_other) ~ pd_med_other,
      TRUE ~ NA),
      MedctnPriorConcomName = case_when(
        is_pd_med && is.na(pd_med_other) ~ record$pd_med_other,
        !is_pd_med ~ record$non_pd_med),
      MedctnPriorConcomIndTxt = record$indication,
      MedctnPriorConcomFreqTxt = case_when(
        !is.na(freq) ~ freq,
        !is.na(freq_oth) ~ freq_oth,
        TRUE ~ NA),
      MedctnPriorConcomDoseMsrTxt = record$dose,
      MedctnPriorConcomDoseUoM = case_when(
        record$units == "other" ~ ifelse(!is.na(units_oth), units_oth, NA),
        !is.na(units) ~ units,
        !is.na(units_oth) ~ units_oth,
        TRUE ~ NA),
      MedctnPriorConcomMinsLstDose = NA_integer_,
      MedctnPriorConcomHrsLstDose = NA_integer_)
  return(dmr_record)
}


#' Parse concomitant medication record for DMR, SUPER PD Cohort
#'
#' Participants can have multiple concomitant medication records for a single
#' visit, hence this type of record needs to be handled specially.
#' This function extracts fields MedctnPriorConcomRteTyp, MedctnPriorConcomName,
#' MedctnPriorConcomDoseMsrTxt, MedctnPriorConcomDoseUoM, MedctnPriorConcomFreqTxt,
#' and MedctnPriorConcomPD. Other fields specific to this form are set to NA.
#' In contrast to AT-HOME PD, SUPER-PD medications are recorded as a single record
#' for each participant. Fields for a specific medication can be grouped by their
#' numeric suffix, e.g., conmed_1, conmed_dose_amt_1, conmed_dose_unit_1, etc.
#' AT-HOME PD data is stored in a different format, see
#' parse_concomitant_medication_record_ahpd
#' @param record A one-row dataframe from the clinical data containing a single record
#' @param mapping The value mapping. A list with heirarchy form > field identifier.
#' @return A tibble with fields specific to concomitant medications
parse_concomitant_medication_record_spd <- function(record, mapping) {
  med_map <- mapping[["concomitant_medications"]]
  num_medications <- as.integer(record$conmed_num)
  dmr_records <- purrr::map_dfr(1:num_medications, function(n) {
    is_pd_med <- record[[glue("conmed_pd_{n}")]] == "Yes"
    conmed <- value_map(med_map, "conmed", record[[glue("conmed_{n}")]])
    conmed_dose_amt <- record[[glue("conmed_dose_amt_{n}")]]
    conmed_dose_unit <- value_map(
      med_map, "conmed_dose_unit", record[[glue("conmed_dose_unit_{n}")]])
    conmed_dose_unit_other <- value_map(
      med_map, "conmed_dose_unit_other", record[[glue("conmed_dose_unit_other_{n}")]])
    conmed_dose_frequency <- value_map(
      med_map, "conmed_dose_frequency", record[[glue("conmed_dose_frequency_{n}")]])
    conmed_dose_frequency_other <- value_map(
      med_map, "conmed_dose_frequency_other",
      record[[glue("conmed_dose_frequency_other_{n}")]])
    conmed_dose_route <- value_map(
      med_map, "conmed_dose_route", record[[glue("conmed_dose_route_{n}")]])
    conmed_dose_route_other <- value_map(
      med_map, "conmed_dose_route_other", record[[glue("conmed_dose_route_other_{n}")]])
    conmed_indication <- record[[glue("conmed_indication_{n}")]]
    dmr_record <- tibble(
        MedctnPriorConcomRteTyp = case_when(
          !is.na(conmed_dose_route) ~ conmed_dose_route,
          !is.na(conmed_dose_route_other) ~ conmed_dose_route_other,
          TRUE ~ NA_character_),
        MedctnPriorConcomPD = case_when(
          !is.na(conmed) ~ conmed,
          TRUE ~ NA_character_),
        MedctnPriorConcomName = case_when(
          is.na(conmed) ~ record[[glue("conmed_{n}")]],
          TRUE ~ conmed),
        MedctnPriorConcomIndTxt = conmed_indication,
        MedctnPriorConcomFreqTxt = case_when(
          !is.na(conmed_dose_frequency) ~ conmed_dose_frequency,
          !is.na(conmed_dose_frequency_other) ~ conmed_dose_frequency_other,
          TRUE ~ NA_character_),
        MedctnPriorConcomDoseMsrTxt = conmed_dose_amt,
        MedctnPriorConcomDoseUoM = case_when(
          record[[glue("conmed_dose_unit_{n}")]] == "other" ~ ifelse(!is.na(conmed_dose_unit_other), conmed_dose_unit_other, NA_character_),
          !is.na(conmed_dose_unit) ~ conmed_dose_unit,
          !is.na(conmed_dose_unit_other) ~ conmed_dose_unit_other,
          TRUE ~ NA_character_),
        MedctnPriorConcomMinsLstDose = NA_integer_,
        MedctnPriorConcomHrsLstDose = NA_integer_)
    return(dmr_record)
  })
  return(dmr_records)
}

#' Parse MDS-UPDRS record for AT-HOME PD cohort
#'
#' @param record A one-row dataframe from the clinical data containing a single record
#' @param field_mapping
#' @param value_mapping The value mapping. A list with heirarchy (form) > (field identifier).
#' @return A tibble with fields specific to the PDBP_MDS-UPDRS form.
parse_mdsupdrs_ahpd <- function(record, field_mapping, value_mapping) {
  visit <- record$visittyppdbpmds


}

#' Parse MDS-UPDRS record for SUPER PD cohort
#'
#' @param record A one-row dataframe from the clinical data containing a single record
#' @param field_mapping
#' @param value_mapping The value mapping. A list with heirarchy (form) > (field identifier) > (values).
#' @return A tibble with two records, if this is a physician administered exam (usually --
#' one participant only completed an ON exam), or one record, if this is the
#' substudy participant survey. All fields are specific to the DMR's MDS-UPDRS form.
parse_mdsupdrs_spd <- function(record, field_mapping, value_mapping) {
  event <- case_when(
      !is.na(visitdate) ~ "physician",
      !is.na(mdsupdrs_sub_dttm) ~ "substudy")
  if (event == "physician") {
    dmr_records <- purrr::map_dfr(c("No", "Yes"), function(med_status) {
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
              value <- tibble(name = dmr_variable,
                              value = value_map_updrs(value_mapping, record, clinical_variable))
            })
        dmr_record <- dmr_record %>%
          pivot_wider(names_from="name", values_from="value")
        # TODO: Fill in UPDRS scores from file
        return(dmr_record)
      })

  } else if (event == "substudy") {

  }
}

value_map <- function(mapping, field, key) {
  if (is.na(key) || !(key %in% names(mapping[[field]]))) {
    return(NA_character_)
  }
  value <- mapping[[field]][[key]]
  return(value)
}

# Map clinical MDS-UPDRS values to DMR permissible values
#'
#' @param mapping The value mapping. A list with heirarchy (form) > (field identifier) > (values).
#' @param record A one-row data frame or named vector
#' @param field The name of the clinical variable in `record` to map
#' @return An integer MDS-UPDRS score
value_map_updrs <- function(mapping, record, field) {
  # IF (This DMR field has no mapping to the clinical data ||
  #     The clinical field is empty in the clinical data)
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

field_map <- function(mapping, field, cohort, visit) {

}

main <- function() {

}
