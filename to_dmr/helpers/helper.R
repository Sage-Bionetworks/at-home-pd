library(tidyverse)
library(synapser)

#' filter rows from a specific redcap_event_name or redcap_repeat_instrument
#' and drop columns which contain only NA
select_survey_and_columns <- function(df, re_name, cohort_name, ri_name=NULL) {
  if (is.null(ri_name)) {
    df %>%
      filter(redcap_event_name == re_name,
             study_cohort == cohort_name) %>%
      select_if(~ any(!is.na(.)))
  } else {
    df %>%
      filter(redcap_repeat_instrument == ri_name,
             study_cohort == cohort_name) %>%
      select_if(~ any(!is.na(.)))
  }
}

select_na <- function(df) {
  df %>%
    select_if(~ any(is.na(.)))
}

select_not_na <- function(df) {
  df %>%
    select_if(~ any(!is.na(.)))
}

which_cohort <- function(clinical, column) {
    relevant_records <- clinical %>%
        filter(!is.na({{ column }}))
    cohort_count <- relevant_records %>%
        count(study_cohort)
    survey_count <- relevant_records %>%
        count(redcap_event_name)
    print(cohort_count)
    print(survey_count)
}

get_dttm_fields <- function(clinical) {
  clinical %>%
    select_if(lubridate::is.Date)
}

examples <- function(clinical, col) {
  examples <- clinical %>%
    drop_na({{col}}) %>%
    distinct({{col}}) %>%
    arrange({{col}})
  return(examples)
}

permissible_values <- function(dmr_dic, field) {
  this_field <- dmr_dic %>%
    filter(field_name == field)
  all_values <- stringr::str_split(this_field[["permissible values"]], ";")
  for (v in unlist(all_values)) {
    print(v)
  }
}

values_to_map <- function(clinical, col) {
  ex <- examples(clinical, {{ col }}) %>%
    purrr::flatten() %>%
    unlist()
  ex_str <- paste0('"', stringr::str_c(ex, collapse='": "",\n"'), '": ""')
  write(ex_str, "value_mapping.json", append=TRUE)
}

spd_conmed <- function(clinical, col) {
  clinical %>%
    filter(!is.na(conmed_1)) %>%
    select(matches(glue::glue("{col}_\\d"))) %>%
    purrr::map(as.character) %>%
    as_tibble() %>%
    pivot_longer(dplyr::everything()) %>%
    distinct(value)
}

select_primary_columns <- function(clinical, add_col = NULL) {
    clinical %>%
        select(guid, study_cohort, redcap_event_name, visstatdttm,
               mdsupdrs_sub_dttm, {{ add_col }})
}

mdsupdrs_helper <- function(clinical, mapping) {
    mds_mapping <- mapping %>%
        filter(stringr::str_detect(dmr_variable, "MDSUPDRS")) %>%
        mutate(notes = lapply(notes, parse_mds_updrs_notes))
    redcap_events <- purrr::pmap(list(
            dmr_var = mds_mapping$dmr_variable,
            primary_var = mds_mapping$ctcc_variable,
            secondary_vars = mds_mapping$notes),
    function(dmr_var, primary_var, secondary_vars) {
        events_list <- list()
        every_var <- na.omit(c(primary_var, unlist(secondary_vars)))
        for (v in every_var) {
            this_var_events <- clinical %>%
                filter(!is.na(.[[v]])) %>%
                distinct(redcap_event_name)
            events_list <- c(events_list,
                             list(this_var_events$redcap_event_name))
        }
        names(events_list) <- every_var
        return(events_list)
    })
    names(redcap_events) <- mds_mapping$dmr_variable
    purrr::map2(names(redcap_events), redcap_events, function(.x, .y) {
        print(glue::glue("##### { .x } #####"))
    })
}

parse_mds_updrs_notes <- function(s) {
    s_short <- stringr::str_sub(s, 6) # remove "also "
    s_split <- stringr::str_split(s_short, " and ")
    return(s_split)
}

select_survey <- function(clinical, survey_name, .clinical_dic = clinical_dic) {
  survey_fields <- .clinical_dic %>%
    dplyr::filter(`Form Name` == survey_name) %>%
    dplyr::select(`Variable / Field Name`) %>%
    unlist() %>%
    as.vector()
  clinical_data <- clinical %>%
    dplyr::select(guid, redcap_event_name, tidyselect::any_of(survey_fields), study_cohort) # checkbox fields not included
  return(clinical_data)
}

#' Return the clinical forms this record contains
has_forms <- function(record) {
  forms <- purrr::map(names(form_to_datetime_mapping), function(f) {
      if (!is.null(form_to_datetime_mapping[[f]])) {
        if (!is.na(record[[form_to_datetime_mapping[[f]]]])) {
          return(f)
        }
      }
    })
  forms <- purrr::compact(forms)
  return(forms)
}

synapser::synLogin()
clinical <- synhelper::synGetFile("syn17051543")
mapping <- synhelper::synGetFile("syn23593083")
field_mapping <- synhelper::synGetFile("syn25056102")
value_mapping <- {
  tryCatch({
    jsonlite::read_json("value_mapping.json")
  }, error = function(e) {
      warning("using Synapse value mapping")
      f <- synGet("syn25155671")
      jsonlite::read_json(f$path)
  })
}
ahpd_mdsupdrs_scores <- synhelper::synGetFile("syn25050919")
super_on_mdsupdrs_scores <- synhelper::synGetFile("syn25165546")
super_off_mdsupdrs_scores <- synhelper::synGetFile("syn25165548")
dmr_dic <- synhelper::synGetFile("syn24171997")
clinical_dic  <- synhelper::synGetFile("syn21740194")
cohorts  <- synhelper::synGetFile("syn24173690") %>%
    rename(study_cohort = cohort)
clinical <- clinical %>%
    inner_join(cohorts, by = "guid")
all_possible_forms <- clinical %>%
    distinct(redcap_event_name)
# ahpd_surveys <- purrr::map(
#                            all_possible_forms$redcap_event_name, function(ren) {
#                                this_survey <- clinical %>%
#            select_survey_and_columns(ren, "at-home-pd")
#            return(names(this_survey))
#        })
# names(ahpd_surveys) <- all_possible_forms$redcap_event_name
# super_surveys <- purrr::map(
#        all_possible_forms$redcap_event_name, function(ren) {
#            this_survey <- clinical %>%
#                select_survey_and_columns(ren, "super-pd")
#            return(names(this_survey))
#        })
# names(super_surveys) <- all_possible_forms$redcap_event_name
