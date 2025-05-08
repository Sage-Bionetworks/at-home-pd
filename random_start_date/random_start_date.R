library(synapser)
library(tidyverse)

MJFF_USERS <- "syn21670519"
AT_HOME_PD_REDCAP <- "syn17051543"
AT_HOME_PD_USERS <- "syn16786935"
AT_HOME_PD_2_USERS <- "syn52789119"
AT_HOME_PD_2_REDCAP <- "syn66250726"
ROCHESTER_PARENT <- "syn18637131"
AT_HOME_PD_2_PARENT <- "syn66496678"
MJFF_PARENT <- "syn18637133"
BRIDGE_PARENT <- "syn12617210"
DEIDENTIFIED_DATA_OFFSET <- "syn18637903"
FOX_INSIGHT_EXPORTED_SURVEYS <- "syn21670408"
BRIDGE_MAPPING <- list(
    "syn17015960" = "syn18681888",
    "syn17015065" = "syn18681890",
    "syn17014786" = "syn18681891",
    "syn17014785" = "syn18681893",
    "syn17014784" = "syn18681894",
    "syn17014783" = "syn18681898",
    "syn17014781" = "syn18681899",
    "syn17014780" = "syn18681900",
    "syn17014779" = "syn18681901",
    "syn17014778" = "syn18681902",
    "syn17014777" = "syn18681903",
    "syn17014776" = "syn18681904",
    "syn17014775" = "syn18681905",
    "syn17014782" = "syn18683083",
    "syn20712373" = "syn20715163")

read_syn_csv <- function(syn_id, encoding = "UTF-8") {
  f <- synGet(syn_id)
  df <- read_csv(f$path,
                 locale = locale(encoding = encoding),
                 guess_max=10000)
  return(df)
}

read_syn_table <- function(syn_id) {
  q <- synTableQuery(paste("select * from", syn_id))
  table <- q$asDataFrame() %>%
    as_tibble() %>%
    select(-ROW_ID, -ROW_VERSION)
  return(table)
}

curate_user_list <- function() {
  mjff_users <- read_syn_csv(MJFF_USERS) %>%
    distinct(guid) %>%
    mutate(source = "MJFF")
  at_home_pd_redcap <- read_syn_csv(AT_HOME_PD_REDCAP) %>%
    distinct(guid) %>%
    mutate(source = "ROCHESTER")
  at_home_pd_one_users <- read_syn_table(AT_HOME_PD_USERS) %>%
    distinct(guid) %>%
    mutate(source = "MPOWER")
  at_home_pd_two_users <- read_syn_table(AT_HOME_PD_2_USERS) %>%
    distinct(guid) %>%
    mutate(source = "MPOWER")
  at_home_pd_two_redcap <- read_syn_csv(AT_HOME_PD_2_REDCAP) %>%
    mutate(guid = record_id) %>%
    distinct(guid) %>%
    mutate(source = "ROCHESTER")
  bridge_users <- bind_rows(at_home_pd_one_users, at_home_pd_two_users) %>%
    distinct(guid, .keep_all = TRUE)
  redcap_users <- bind_rows(at_home_pd_redcap, at_home_pd_two_redcap) %>%
    distinct(guid, .keep_all = TRUE)
  users <- bind_rows(mjff_users, redcap_users, bridge_users)
  users <- users %>%
    group_by(guid) %>%
    mutate(day_offset = round(runif(1, -10, 10))) %>%
    ungroup()
  return(users)
}

update_user_list <- function(users) {
  existing_users <- read_syn_table(DEIDENTIFIED_DATA_OFFSET) %>%
    mutate(guid = as.character(guid),
           source = as.character(source))
  brand_new_users <- anti_join(users, existing_users, by = "guid")
  new_sources <- anti_join(users, existing_users, by = c("guid", "source")) %>%
    anti_join(brand_new_users, by = "guid") %>%
    select(-day_offset)
  new_sources <- new_sources %>%
    bind_rows(existing_users) %>%
    group_by(guid) %>%
    mutate(day_offset = replace_na(day_offset, median(day_offset, na.rm = T))) %>%
    semi_join(new_sources, by = c("guid", "source"))
  new_users <- bind_rows(brand_new_users, new_sources)
  if (nrow(new_users) > 0) {
    updated_table <- synStore(Table(DEIDENTIFIED_DATA_OFFSET, new_users))
    all_users <- read_syn_table(updated_table$tableId)
  } else {
    all_users <- existing_users
  }
  return(all_users)
}

perturb_mjff_dates <- function(users) {
  mjff <- synGetChildren(FOX_INSIGHT_EXPORTED_SURVEYS)$asList()
  names(mjff) <- purrr::map(mjff, ~ .$id)
  mjff <- purrr::map(mjff, ~ read_syn_csv(.$id))
  mjff$syn21670519 <- NULL # users.csv, we don't want to export this (?)
  mjff$syn21670549 <- NULL # General.csv, we don't want to export this (?)
  mjff$syn21670554 <- NULL # Genetic.csv, we don't want to export this (?)
  date_cols <- c("study_date", "registration_date")
  mjff_dates <- mjff %>%
    purrr::map(function(df) {
      df <- df %>%
        select_if(names(df) %in% date_cols) %>%
        purrr::map_dfc(lubridate::as_datetime)
      return(df)
    })
  mjff <- purrr::map2(mjff, mjff_dates, function(df, df_dates) {
      if (nrow(df_dates) > 0) {
        df <- df %>%
          select_if(!(names(df) %in% date_cols)) %>%
          select(-fox_insight_id) %>%
          bind_cols(df_dates)
      }
      return(df)
    })
  mjff_perturbed <- purrr::map(mjff, function(df) {
    col_order <- names(df)
    if (any(date_cols %in% col_order)) {
      df <- perturb_dates(
        df = df,
        users = users,
        source_name = "MJFF",
        guid = "guid",
        date_cols = date_cols)
      }
    return(df[col_order])
  })
}

perturb_at_home_pd_dates <- function(users) {
  rochester <- read_syn_csv(AT_HOME_PD_REDCAP)
  date_cols <- c("assessdate_fall", "assessdate_fall_m12", "bl_partburdenblvisitdt",
                 "determinefall_visitdt", "compliance_dttm", "concldttm", "wddt",
                 "irbapprovedt", "subjsigdt", "start_dt", "inexdttm", "demo_dttm",
                 "dob", "screenorientdttm", "phoneorientdttm", "onsetdt", "mdsupdrs_dttm",
                 "moca_dttm", "mseadl", "cgi_dttm", "determinefall_visitdt",
                 "bl_partburdenblvisitdt", "compliance_dttm", "visstatdttm",
                 "stop_dt", "reslvdt",
                 "partburdenv1visitdt", "enrollment_confirmation_timestamp",
                 "prebaseline_survey_timestamp", "preference_and_burden_bl_timestamp",
                 "preference_and_burden_visit_timestamp", "previsit_survey_timestamp",
                 "redcap_survey_identifier", "reportable_event_timestamp",
                 "study_burst_reminders_timestamp", "visitdate", "inexdttm_spd",
                 "time_mdsupdrs", "date_medhx", "phoneorientdttmcmp",
                 "mdsupdrs_sub_dttm", "moca_dttm_v2")
  rochester <- rochester %>%
    mutate(preference_and_burden_bl_timestamp = replace(
      preference_and_burden_bl_timestamp,
      preference_and_burden_bl_timestamp == "[not completed]", NA),
      preference_and_burden_visit_timestamp = replace(
      preference_and_burden_visit_timestamp,
      preference_and_burden_visit_timestamp == "[not completed]", NA)) %>%
    mutate(preference_and_burden_bl_timestamp = lubridate::as_datetime(preference_and_burden_bl_timestamp),
           preference_and_burden_visit_timestamp = lubridate::as_datetime(preference_and_burden_visit_timestamp))
  col_order <- names(rochester)
  if (any(purrr::map(rochester[date_cols], class) == "character")) {
    error_message <- paste("There are some date columns in the rochester clinical",
               "which were cast as character and may prevent them from",
               "being deidentified")
    synSendMessage(list("3342492"), "Error in deidentifying AHPD data",
                   error_message, contentType = "text")
    stop(error_message)
  }
  rochester_perturbed <- perturb_dates(
    df = rochester,
    users = users,
    source_name = "ROCHESTER",
    guid = "guid",
    date_cols = date_cols)
  return(rochester_perturbed[col_order])
}

perturb_at_home_pd_two_dates <- function(users) {
  at_home_pd_two_redcap <- read_syn_csv(AT_HOME_PD_2_REDCAP)
  date_fields <- c(
    "start_dt", "stop_dt", "wddt", "pd_irb_dt", "dd_sd", "dd_ed",
    "disp_date", "device_return_dt", "visit_dt", "reconsent_date",
    "irbapprovedt", "subjsigdt", "onsetdt", "compliance_dttm",
    "ffup_dttm", "ffup_dtime", "ffup_resched_dtime", "ddl_dttm",
    "phoneorientdttm", "visit_dttm", "act_fitbit_dttm", "concldttm",
    "econsent_study_user_timestamp", "icf_dttm", "icl_dttm", "invsigdttm",
    "enrollment_sch_dttm", "device_acc_dttm", "pgis_dttm", "fall_questionnaire_dt",
    "pro_mdsupdrs_dttm", "moca_dttm", "nms_dttm", "mdsupdrs_dttm", "initials_2_v1",
    "nms_summ_date", "sus_fitbit_dttm", "sus_phone_dttm", "sus_actigraph_dttm",
    "system_usability_scale_fitbit_timestamp", "system_usability_scale_mpower_timestamp",
    "system_usability_scale_actigraph_timestamp", "prb_dttm", "cgis_dttm", "revent_dttm",
    "pd_record_date", "pd_date", "determinefall_visitdt", "reslvdt", "comm_dttm",
    "cl_dttm", "datecontact", "inex_dttm", "date_diagnosis", "demo_dttm",
    "initial_dttm", "randomization_dt", "mseadl", "fall_questionnaire_timestamp",
    "patient_global_impression_scale_timestamp", "perceived_research_burden_assessment_timestamp"
  )

  # Join the user day_offset by matching record_id to guid
  ahpd2_joined <- at_home_pd_two_redcap %>%
    left_join(users %>% rename(record_id = guid) %>% filter(source=="ROCHESTER"), by = "record_id")

  # Apply the perturbation for each date column safely
  df_dates_perturbed <- ahpd2_joined %>%
    select(any_of(c(date_fields, "day_offset")))
  df_no_date_cols <- ahpd2_joined %>%
    select(!any_of(c(date_fields, "day_offset")))
  for (col in names(df_dates_perturbed)) {
    if (col != "day_offset" && typeof(df_dates_perturbed[[col]]) != "logical") {
      df_dates_perturbed <- df_dates_perturbed %>%
        mutate(!!col := !!sym(col) + lubridate::days(day_offset))
    }
  }
  df_perturbed <- df_no_date_cols %>%
    bind_cols(df_dates_perturbed) %>%
    select(-day_offset, -source)
  return(df_perturbed)
}

perturb_bridge_dates <- function(users, table_mapping = NULL) {
  bridge <- purrr::map(names(BRIDGE_MAPPING), read_syn_table)
  names(bridge) <- names(BRIDGE_MAPPING)
  if (!is.null(table_mapping)) {
    deidentified_bridge <- purrr::map(table_mapping, read_syn_table)
    bridge_diff <- purrr::map(names(bridge), function(source_id) {
      source_table <- bridge[[source_id]] %>%
        mutate(recordId = as.character(recordId))
      target_table <- deidentified_bridge[[source_id]] %>%
        mutate(recordId = as.character(recordId))
      new_records <- anti_join(source_table, target_table, by = "recordId")
      return(new_records)
    })
  names(bridge_diff) <- names(bridge)
  } else {
    bridge_diff <- bridge
  }
  bridge_perturbed <- purrr::map(bridge_diff, function(df) {
    col_order <- names(df)
    df <- df %>%
      mutate(uploadDate = lubridate::as_date(uploadDate),
             createdOn = lubridate::as_datetime(createdOn))
    if (has_name(df, "metadata.startDate")) {
      df <- df %>%
        mutate(metadata.startDate = lubridate::as_datetime(metadata.startDate))
    }
    if (has_name(df, "metadata.endDate")) {
      df <- df %>%
        mutate(metadata.endDate = lubridate::as_datetime(metadata.endDate))
    }
    bridge_perturbed <- perturb_dates(
      df = df,
      users = users,
      source_name = "MPOWER",
      guid = "externalId",
      date_cols = c("uploadDate", "createdOn", "metadata.startDate",
                    "metadata.endDate", "displacement.timestampDate"))
    bridge_perturbed <- bridge_perturbed %>%
      distinct(recordId, .keep_all = TRUE)
    return(bridge_perturbed[col_order])
    })
  bridge_perturbed <- purrr::map(bridge_perturbed, function(df) { # drop identifiable files
    file_cols <- c("rawData", "left_tapping.samples", "right_tapping.samples",
               "left_motion.json", "right_motion.json", "trackedItems.items",
               "balance_motion.json", "walk_motion.json", "rawMetadata")
    df_with_file_cols_removed <- df %>%
      select_if(!(names(.) %in% file_cols))
    return(df_with_file_cols_removed)
  })
  return(bridge_perturbed)
}

perturb_dates <- function(df, users, source_name, guid, date_cols=NULL) {
  if (nrow(df) == 0) {
    return(df)
  }
  df[[guid]] <- as.character(df[[guid]])
  df_with_offsets <- users %>%
    filter(source == source_name) %>%
    select(guid, day_offset)
  guid_flag <- TRUE
  if (!has_name(df, "guid")) {
    df$guid <- df[[guid]]
    guid_flag <- FALSE
  }
  df_with_offsets <- df %>%
    left_join(df_with_offsets, by = "guid")
  if (!guid_flag) {
    df <- select(df, -guid)
  }
  df_offsets <- lubridate::days(df_with_offsets$day_offset)
  if (is.null(date_cols)) {
    df_dates_perturbed <- df_with_offsets %>%
      select_if(lubridate::is.timepoint)
    df_perturbed <- df_with_offsets %>%
      select_if(~ !lubridate::is.timepoint(.))
  } else {
    df_dates_perturbed <- df_with_offsets %>%
      select_if(names(.) %in% date_cols)
    df_perturbed <- df_with_offsets %>%
      select_if(!(names(.) %in% date_cols))
  }
  df_dates_perturbed <- purrr::map_dfc(df_dates_perturbed, function(col) {
    if (typeof(col) == "logical") { # empty col
      return(col)
    }
    return(col + df_offsets)
  })
  df_perturbed <- df_perturbed %>%
    bind_cols(df_dates_perturbed) %>%
    select(-day_offset)
  return(df_perturbed)
}

store_at_home_pd_perturbed <- function(rochester_dataset) {
  fname <- "deidentified_exported_records.csv"
  write_csv(rochester_dataset, fname, na="")
  f <- synapser::File(fname, parent = ROCHESTER_PARENT)
  synStore(f, used = list(AT_HOME_PD_REDCAP))
  unlink(fname)
}

store_at_home_pd_two_perturbed <- function(rochester_dataset) {
  fname <- "deidentified_exported_records.csv"
  write_csv(rochester_dataset, fname, na="")
  f <- synapser::File(fname, parent = AT_HOME_PD_2_PARENT)
  synStore(f, used = list(AT_HOME_PD_2_REDCAP))
  unlink(fname)
}

store_mjff_perturbed <- function(mjff_dataset, table_mapping=NULL) {
  purrr::map2(names(mjff_dataset), mjff_dataset, function(.x, .y) {
    file_info <- synGet(.x, downloadFile = FALSE)
    fname <- paste0("deidentified_", file_info$properties$name)
    write_csv(.y, fname)
    f <- synapser::File(fname, parent = MJFF_PARENT)
    synStore(f, used = list(.x))
    unlink(fname)
  })
}

store_bridge_perturbed <- function(bridge_dataset, table_mapping=NULL) {
  if (!is.null(table_mapping)) {
    purrr::map2(bridge_dataset, names(bridge_dataset), function(df, source) {
      if (nrow(df) > 0) {
        target <- table_mapping[[source]]
        synStore(Table(target, df), used = list(source))
      }
    })
  }
  else {
    purrr::map2(names(bridge_dataset), bridge_dataset, function(.x, .y) {
      table_info <- synGet(.x)
      tname <- paste("Deidentified",
                      table_info$properties$name)
      table_cols <- synGetTableColumns(.x)$asList()
      schema <- Schema(name = tname, columns = table_cols, parent = BRIDGE_PARENT)
      synStore(Table(schema, .y), used = list(.x))
    })
  }
}

main <- function() {
  # set env variables synapseUsername and synapsePassword before running
  synLogin(Sys.getenv("synapseUsername"), Sys.getenv("synapsePassword"))
  users <- curate_user_list() %>%
    update_user_list()
  at_home_pd_dataset <- perturb_at_home_pd_dates(users)
  at_home_pd_two_dataset <- perturb_at_home_pd_two_dates(users)
  mjff_dataset <- perturb_mjff_dates(users)
  bridge_dataset <- perturb_bridge_dates(users, table_mapping = BRIDGE_MAPPING)
  store_at_home_pd_perturbed(at_home_pd_dataset)
  store_at_home_pd_two_perturbed(at_home_pd_two_dataset)
  store_mjff_perturbed(mjff_dataset)
  store_bridge_perturbed(bridge_dataset, table_mapping = BRIDGE_MAPPING)
}

main()
