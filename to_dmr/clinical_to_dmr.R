#' Given AT-HOME PD clinical data from REDCap, translate each record into a
#' format which conforms to the appropriate PDBP DMR schema.

library(synapser)
library(dplyr)
library(glue)

AHPD_MDSUPDRS_SCORES <- "syn25050919"
SUPER_OFF_MDSUPDRS_SCORES <- "syn25165548"
SUPER_ON_MDSUPDRS_SCORES <- "syn25165546"
CLINICAL_DATA_DICTIONARY <- "syn21740194"

#' Read a CSV file from Synapse as a tibble
read_synapse_csv <- function(synapse_id) {
  f <- synapser::synGet(synapse_id)
  df <- readr::read_csv(f$path)
  return(df)
}

#' Read a TSV file from Synapse as a tibble
read_synapse_tsv <- function(synapse_id) {
  f <- synapser::synGet(synapse_id)
  df <- readr::read_tsv(f$path)
  return(df)
}

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
#' @param record A one-row dataframe containing a single record.
#' This record consists of the combined fields of the clinical forms `mdsupdrs` and
#' either `prebaseline_survey' or `previsit_survey`
#' @param field_mapping The DMR to clinical field mapping.
#' @param value_mapping The value mapping. A list with heirarchy (form) > (field identifier).
#' @return A tibble with fields specific to the PDBP_MDS-UPDRS form.
parse_mdsupdrs_ahpd <- function(record, field_mapping, value_mapping, scores) {
  # The same clinical forms were administered at 12 and 24 months
  this_visit <- case_when(
      has_name(record, "assessdate_fall") ~ "Baseline",
      has_name(record, "assessdate_fall_m12") ~ "12 Months")
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
      !is.na(record$visitdate) ~ "physician",
      !is.na(record$mdsupdrs_sub_dttm) ~ "substudy")
  date_of_exam <- case_when(
      !is.na(record$visitdate) ~ as.Date(record$visitdate),
      !is.na(record$mdsupdrs_sub_dttm) ~ as.Date(record$mdsupdrs_sub_dttm))
  if (event == "physician") {
    # All MDS-UPDRS sections, for both OFF and ON med states
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
        })
    # TODO Include part 3 section score if applicable (rigidity scores are missing)
    dmr_record <- dmr_record %>%
      pivot_wider(names_from="name", value_from="value")
    return(dmr_record)
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

#' Get section scores for a specific participant and date
#'
#' @param scores A dataframe with MDS-UPDRS scores indexed by guid and createdOn
#' @param participant_id The participant's GUID
#' @param date_of_exam The date or datetime of the exam
#' @return row-wise section scores for convenient anti-join/row-bind use
mdsupdrs_section_scores <- function(scores, participant_id, date_of_exam) {
  score <- scores %>%
    filter(guid == participant_id,
           as.Date(createdOn) == as.Date(date_of_exam)) %>%
    select(all_of(c("UPDRS1", "UPDRS2", "UPDRS3", "UPDRS4")))
  total_score <- {
    total_score <- score %>%
      pivot_longer(dplyr::everything()) %>%
      summarize(total_score = sum(value))
    as.character(total_score$total_score)
  }
  score <- score %>%
    mutate_all(as.character)
  section_scores <- tribble(
        ~name, ~value,
        "MDSUPDRS_PartIScore", score$UPDRS1,
        "MDSUPDRS_PartIIScore", score$UPDRS2,
        "MDSUPDRS_PartIIIScore", score$UPDRS3,
        "MDSUPDRS_PartIVScore", score$UPDRS4,
        "MDSUPDRS_TotalScore", total_score)
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
#' @return A dataframe with one row per participant. Each column will be of
#' type character().
get_event_and_form_fields <- function(clinical, clinical_dic, event_name, form_name) {
  form_fields <- clinical_dic %>%
    filter(`Form Name` %in% form_name,
           `Field Type` != "descriptive") %>%
    distinct(`Variable / Field Name`)
  event_records <- clinical %>%
    filter(redcap_event_name %in% event_name) %>%
    select(guid, all_of(form_fields[["Variable / Field Name"]])) %>%
    mutate(across(.fns=as.character)) %>%
    pivot_longer(!guid) %>%
    drop_na() %>%
    pivot_wider(guid)
  return(event_records)
}

main <- function() {
  synLogin()
  ahpd_mdsupdrs_scores <- read_synapse_tsv(AHPD_MDSUPDRS_SCORES)
  super_mdsupdrs_scores <- list(
    "Physician_ON" = read_synapse_csv(SUPER_ON_MDSUPDRS_SCORES),
    "Physician_OFF" = read_synapse_csv(SUPER_OFF_MDSUPDRS_SCORES))
  clinical_data_dictionary <- read_synapse_csv(CLINICAL_DATA_DICTIONARY)
}
