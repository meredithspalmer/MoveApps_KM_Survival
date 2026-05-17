library(move2)
library(survival)
library(survminer)
library(ggplot2)   
library(dplyr)
library(lubridate)
library(stringr)
library(sf)
library(forcats)
library(tidyr)
library(purrr) 
library(viridis)
library(ggpubr)

# logger.fatal(), logger.error(), logger.warn(), logger.info(), logger.debug(), logger.trace()

# Survival Function 
rFunction = function(data,  
                     time_period_start, 
                     time_period_end, 
                     censor_capture_mortality,
                     fix_na_start_times, 
                     fix_na_end_times,  
                     subset_condition_1,
                     subset_condition_define_1, 
                     subset_condition_2,
                     subset_condition_define_2, 
                     group_comparison_individual,
                     survival_yr_start,
                     animal_birth_hatch_year_table, 
                     life_table_days, 
                     calc_month_mort, 
                     calc_tracking_history, 
                     add_cis, 
                     zoom_to_plot) {
  
  ## Load auxiliary data ------------------------------------------------------ 
  
  if(!is.null(animal_birth_hatch_year_table)){
    animal_birth_hatch_year_table <- read.csv(getAuxiliaryFilePath("animal_birth_hatch_year_table"))
    logger.info("Auxiliary birth/hatch year data loaded.")
  }
  
  
  ## Cleaning and cropping ----------------------------------------------------
  
  data <- dplyr::filter(data, !sf::st_is_empty(data))       # Exclude empty locations
  data <- mt_filter_unique(data)                            # Exclude marked outliers 
  data <- data %>% filter_track_data(is_test == FALSE)      # Exclude data marked "test"
  
  
  ## Aggregate across multiple deployments (where present) ---
  
  # Extract event-level data 
  events <- data |>
    as_tibble() |>
    dplyr::select(any_of(c(
      "deployment_id",
      "individual_local_identifier",   
      "timestamp")))
  
  # Columns to drop from track data (administrative/redundant metadata)
  exclude_cols <- c("acknowledgements",
                    "citation",
                    "contact_person_name",
                    "deployment_local_identifier",
                    "group_id",
                    "has_quota",
                    "i_am_collaborator",
                    "i_am_owner",
                    "i_can_see_data",
                    "i_have_download_access",
                    "individual_comments",
                    "license_type",
                    "main_location",
                    "principal_investigator_name",
                    "serial_no",
                    "study_number_of_deployments",
                    "study_objective",
                    "study_permission",
                    "study_site",
                    "study_type",
                    "suspend_license_terms",
                    "tag_comments",
                    "tag_failure_comments",
                    "tag_local_identifier",
                    "tag_number_of_deployments",
                    "taxon_canonical_name",
                    "taxon_detail",
                    "taxon_ids",
                    "there_are_data_which_i_cannot_see", 
                    "is_test")
  
  tracks <- mt_track_data(data) |>
    mutate(mortality_location_filled = if_else(
      is.na(mortality_location) | st_is_empty(mortality_location),
      0L, 1L)) |>
    dplyr::select(-any_of(exclude_cols))
  
  # Join track attributes to every event row
  use_deployment_join <- all("deployment_id" %in% names(events),
                             "deployment_id" %in% names(tracks)) &&
    any(!is.na(events$deployment_id)) &&
    any(!is.na(tracks$deployment_id))
  
  if (use_deployment_join) {
    events_with_ind <- events |> left_join(tracks, by = "deployment_id")
    
  } else {
    if (!"individual_local_identifier" %in% names(events)) {
      logger.fatal("Cannot join: neither deployment_id nor individual_local_identifier is available in events")
    }
    logger.info("Joining on individual_local_identifier (deployment_id join not possible)")
    events_with_ind <- events |> left_join(tracks, by = "individual_local_identifier")
  }
  
  events_with_ind <- events_with_ind |>
    relocate(any_of(c("individual_id", "individual_local_identifier",
                      "deployment_id", "timestamp")),
             .before = everything())
  
  # Infer column groups dynamically 
  cols_timestamp_min <- c("timestamp_first_deployed_location", "deploy_on_timestamp")
  cols_timestamp_max <- c("timestamp_last_deployed_location",  "deploy_off_timestamp")
  
  cols_already_handled <- c(
    "individual_id", "individual_local_identifier", "deployment_id",
    "timestamp", "geometry", "mortality_location", "mortality_location_filled",
    cols_timestamp_min, cols_timestamp_max
  )
  
  # Character/factor: collapse unique values
  cols_categorical <- events_with_ind |>
    select(where(~ is.character(.) || is.factor(.))) |>
    select(-any_of(cols_already_handled)) |>
    names()
  
  # Numeric and units: mean
  cols_numeric <- events_with_ind |>
    select(where(~ is.numeric(.) || inherits(., "units"))) |>
    select(-any_of(cols_already_handled)) |>
    names()
  
  # Logical: any() (that is, TRUE if any event flagged TRUE)
  cols_logical <- events_with_ind |>
    select(where(is.logical)) |>
    select(-any_of(cols_already_handled)) |>
    names()
  
  # POSIXct/Date not already in explicit min/max lists: take first non-NA value
  cols_timestamp_first <- events_with_ind |>
    select(where(~ inherits(., "POSIXct") || inherits(., "Date"))) |>
    select(-any_of(c(cols_already_handled, cols_timestamp_min, cols_timestamp_max))) |>
    names()
  
  # sfc geometry columns: take first non-empty geometry per individual
  cols_geometry <- events_with_ind |>
    select(where(~ inherits(., "sfc"))) |>
    select(-any_of(cols_already_handled)) |>
    names()
  
  logger.info(paste("Categorical columns:",    paste(cols_categorical,     collapse = ", ")))
  logger.info(paste("Numeric/units columns:",  paste(cols_numeric,         collapse = ", ")))
  logger.info(paste("Logical columns:",        paste(cols_logical,         collapse = ", ")))
  logger.info(paste("Timestamp (first) cols:", paste(cols_timestamp_first, collapse = ", ")))
  logger.info(paste("Geometry (first) cols:",  paste(cols_geometry,        collapse = ", ")))
  
  
  # Summarise to individual level
  na_like <- function(x) x[NA_integer_]
  
  # Helper function: first non-empty geometry or an empty point if all empty/NA
  first_non_empty_geom <- function(x) {
    non_empty <- x[!is.na(x) & !st_is_empty(x)]
    if (length(non_empty) > 0) non_empty[[1]] else st_point()
  }
  
  summary_table <- events_with_ind |>
    group_by(individual_id, individual_local_identifier) |>
    summarise(
      first_timestamp = min(as.Date(timestamp), na.rm = TRUE),
      last_timestamp  = max(as.Date(timestamp), na.rm = TRUE),
      n_locations     = n(),
      
      n_deployments = if ("deployment_id" %in% names(events_with_ind))
        n_distinct(deployment_id[!is.na(deployment_id)]) else 1L,
      
      mortality_location_filled =
        if ("mortality_location_filled" %in% names(events_with_ind))
          as.integer(any(mortality_location_filled >= 1, na.rm = TRUE)) else NA_integer_,
      
      across(any_of(cols_timestamp_min),
             ~ if (all(is.na(.))) na_like(.) else min(., na.rm = TRUE)),
      
      across(any_of(cols_timestamp_max),
             ~ if (all(is.na(.))) na_like(.) else max(., na.rm = TRUE)),
      
      # Deployment-level timestamp constants — take first non-NA value
      across(any_of(cols_timestamp_first),
             ~ if (all(is.na(.))) na_like(.) else first(na.omit(.))),
      
      across(any_of(cols_categorical),
             ~ if (all(is.na(.))) NA_character_
             else str_c(unique(.[!is.na(.)]), collapse = " | ")),
      
      across(any_of(cols_numeric),
             ~ if (all(is.na(.))) na_like(.) else mean(., na.rm = TRUE)),
      
      # Logical — TRUE if any event is TRUE, NA if all NA
      across(any_of(cols_logical),
             ~ if (all(is.na(.))) NA else any(., na.rm = TRUE)),
      
      # Geometry — first non-empty geometry, or empty point
      across(any_of(cols_geometry),
             ~ st_sfc(first_non_empty_geom(.), crs = st_crs(.))),
      
      .groups = "drop") |>
    
    mutate(
      across(any_of(cols_categorical), ~ if_else(. == "", NA_character_, .)),
      across(any_of(c("deploy_on_timestamp", "deploy_off_timestamp")), as.Date)
    )
  
  
  # Quick checks 
  unhandled_cols <- setdiff(names(events_with_ind),
                            c(cols_already_handled, cols_categorical, cols_numeric,
                              cols_logical, cols_timestamp_first, cols_geometry,
                              "first_timestamp", "last_timestamp", "n_locations", "n_deployments"))
  
  if (length(unhandled_cols) > 0) {
    logger.warn(paste("Columns not summarised (check types):", paste(unhandled_cols, collapse = ", ")))
  }
  
  
  ## Inject NA columns for optional mortality fields absent from this dataset --- 
  if (!"mortality_date" %in% names(summary_table)) summary_table$mortality_date <- NA_Date_
  if (!"mortality_type" %in% names(summary_table)) summary_table$mortality_type <- NA_character_
  if (!"death_comments" %in% names(summary_table)) summary_table$death_comments <- NA_character_
  if (!"deployment_end_comments" %in% names(summary_table)) summary_table$deployment_end_comments <- NA_character_
  if (!"deployment_end_type" %in% names(summary_table)) summary_table$deployment_end_type <- NA_character_
  
  
  ## Clean dates ---
  
  # Start times  
  if(fix_na_start_times == "timestamp"){
    summary_table <- summary_table %>% 
      mutate(missing_timestamp_start = is.na(deploy_on_timestamp))
    n_missing <- sum(is.na(summary_table$deploy_on_timestamp), na.rm = TRUE)
    
    summary_table <- summary_table %>%
      mutate(deploy_on_timestamp = if_else(
        is.na(deploy_on_timestamp),
        as.Date(first_timestamp),
        deploy_on_timestamp)) %>% 
      dplyr::select(-missing_timestamp_start)
    
    if (n_missing > 0) {
      logger.info(
        sprintf("Warning: Replaced %d missing deploy_on_timestamp value%s with first_timestamp.",
                n_missing, if (n_missing == 1) "" else "s"))
    }
  } 
  
  if (fix_na_start_times == "remove"){
    n_missing <- sum(is.na(summary_table$deploy_on_timestamp))
    summary_table <- summary_table %>% filter(!is.na(deploy_on_timestamp))
    
    if (n_missing > 0) {
      logger.info(sprintf("Warning: Removed %d deploy_on_timestamp value%s that were NA.", n_missing,
                          if (n_missing == 1) "" else "s"), call. = FALSE, immediate. = TRUE)
    }
  }
  
  # End times 
  if(fix_na_end_times == "timestamp"){
    summary_table <- summary_table %>%
      mutate(missing_timestamp_end = is.na(deploy_off_timestamp))
    n_missing <- sum(is.na(summary_table$deploy_off_timestamp), na.rm = TRUE)
    
    summary_table <- summary_table %>%
      mutate(deploy_off_timestamp = if_else(is.na(deploy_off_timestamp),
                                            as.Date(last_timestamp),
                                            deploy_off_timestamp)) %>% 
      dplyr::select(-missing_timestamp_end)
    
    if (n_missing > 0) {
      logger.info(
        sprintf("Warning: Replaced %d missing deploy_off_timestamp value%s with last_timestamp.", 
                n_missing, if (n_missing == 1) "" else "s"), call. = FALSE, immediate. = TRUE)
    }
  } 
  
  if (fix_na_end_times == "systime"){
    summary_table <- summary_table %>%
      mutate(missing_timestamp_end = is.na(deploy_off_timestamp))
    n_missing <- sum(is.na(summary_table$deploy_off_timestamp), na.rm = TRUE)
    
    summary_table <- summary_table %>%
      mutate(deploy_off_timestamp = if_else(is.na(deploy_off_timestamp),
                                            Sys.Date(), 
                                            deploy_off_timestamp))%>% 
      dplyr::select(-missing_timestamp_end)
    
    if (n_missing > 0) {
      logger.info(sprintf("Warning: Replaced %d missing deploy_off_timestamp value%s with current date.",
                          n_missing, if (n_missing == 1) "" else "s"), call. = FALSE, immediate. = TRUE)
    }
  } 
  
  if (fix_na_end_times == "remove"){
    n_missing <- sum(is.na(is.na(summary_table$deploy_off_timestamp)))
    summary_table <- summary_table %>% filter(!is.na(deploy_off_timestamp))
    
    if (n_missing > 0) {
      logger.info(
        sprintf("Warning: Removed %d deploy_off_timestamp and/or deploy_on_timestamp value%s that were NA.", 
                n_missing, if (n_missing == 1) "" else "s"), call. = FALSE, immediate. = TRUE)
    }
  } 
  
  # Remove data for individuals where "deploy_off_timestamp" occurs before "deploy_on_timestamp" 
  n_original <- nrow(summary_table) 
  summary_table <- summary_table %>%
    filter(deploy_off_timestamp >= deploy_on_timestamp)
  n_removed <- n_original - nrow(summary_table)  
  
  if (n_removed > 0) {
    logger.info(
      sprintf("Warning: Removed %d individual%s where deploy_off_timestamp < deploy_on_timestamp.",
              n_removed, if (n_removed == 1) "" else "s"))
  }
  
  
  ## Crop data to user-defined temporal windows ---
  
  # Removed censored data (mortalities within set period of capture)
  if(censor_capture_mortality > 0){
    n_before <- nrow(summary_table)
    
    summary_table <- summary_table %>%
      mutate(raw_deploy_on_timestamp = deploy_on_timestamp) %>%  
      mutate(censor_cutoff = deploy_on_timestamp + lubridate::days(censor_capture_mortality)) %>%
      mutate(remove_due_to_early_end = !is.na(deploy_off_timestamp) & deploy_off_timestamp <= censor_cutoff) %>%
      filter(!remove_due_to_early_end) %>%
      mutate(deploy_on_timestamp = censor_cutoff) %>%
      select(-censor_cutoff, -remove_due_to_early_end)
    
    n_after  <- nrow(summary_table)
    n_removed <- n_before - n_after
    
    if (n_removed > 0) {
      logger.info(paste0("Warning: Removed ", n_removed, " individual(s) because deploy_off_timestamp occurred within ", censor_capture_mortality, " day(s) after deploy_on_timestamp"))
    } 
  }
  
  
  # Crop to study period of interest 
  
  # Save original deploy_off_time 
  summary_table <- summary_table %>% mutate(raw_deploy_off_timestamp = deploy_off_timestamp) 
  
  # Define window 
  effective_start <- if (is.null(time_period_start)) {
    min(summary_table$deploy_on_timestamp, na.rm = TRUE)
  } else {
    time_period_start
  }
  
  effective_end <- if (is.null(time_period_end)) {
    max(summary_table$deploy_off_timestamp, na.rm = TRUE)
  } else {
    time_period_end
  }
  
  # Run updates 
  if(!is.null(time_period_start) | !is.null(time_period_end)){
    
    # Crop to window 
    n_original <- nrow(summary_table) 
    summary_table <- summary_table %>%
      
      # Determine if the deployment overlaps study window 
      mutate(overlaps_study = deploy_on_timestamp <= effective_end & 
               deploy_off_timestamp  >= effective_start) %>%
      filter(overlaps_study | is.na(overlaps_study)) %>%    
      
      # Crop to window 
      mutate(first_timestamp = pmax(deploy_on_timestamp, effective_start, na.rm = TRUE),
             last_timestamp  = pmin(deploy_off_timestamp, effective_end,   na.rm = TRUE)) %>%
      
      # Clean 
      select(-overlaps_study) 
    
    n_removed <- n_original - nrow(summary_table)
    if (n_removed > 0) {
      logger.info(
        sprintf("Warning: %d record%s did not overlap the user-defined study window and were removed.",
                n_removed, if (n_removed == 1) "" else "s"))
    }
  } 
  
  
  ## Calculate entry time and exit time (for staggered entry) ---
  origin_date <- if (is.null(time_period_start) || is.na(time_period_start)) {
    min(summary_table$deploy_on_timestamp, na.rm = TRUE)
  } else {
    time_period_start
  }
  
  summary_table <- summary_table %>%
    mutate(analysis_entry_date = pmax(deploy_on_timestamp, effective_start, na.rm = TRUE),
           analysis_exit_date  = pmin(deploy_off_timestamp, effective_end, na.rm = TRUE),
           entry_time_days = as.numeric(difftime(analysis_entry_date, origin_date, units = "days")),
           exit_time_days  = as.numeric(difftime(analysis_exit_date,  origin_date, units = "days"))) 
  
  
  ## Calculate mortality indicator --------------------------------------------
  # Event = 1 if observed death, 0 if censored or survived 
  
  # Death comments to flag
  positive_pattern <- "dead|death|died|cod|predation|predator|vehicle|collision|killed|poach|poached|shot|hunt|harvest|harvested|mortality"
  
  # Check which columns are present 
  has_death_comments          <- "death_comments"          %in% names(summary_table)
  has_deployment_end_comments <- "deployment_end_comments" %in% names(summary_table)
  has_deployment_end_type     <- "deployment_end_type"     %in% names(summary_table)
  has_mortality_type          <- "mortality_type"          %in% names(summary_table)
  has_mortality_location      <- "mortality_location_filled" %in% names(summary_table)
  
  # Construct mortality event
  summary_table <- summary_table %>%
    
    # Initialise mortality event
    mutate(mortality_event = NA_real_) %>%
    
    # Identify survivors (individuals who last beyond study)
    mutate(survived_beyond_study = !is.na(raw_deploy_off_timestamp) &
             raw_deploy_off_timestamp > as.Date(effective_end),
           mortality_event = if_else(survived_beyond_study, 0L, mortality_event)) %>%
    
    # Update optional columns to remove ambiguity for survivors
    { if (has_death_comments)
      mutate(., death_comments = if_else(survived_beyond_study, "survived beyond study", death_comments))
      else . } %>%
    { if (has_deployment_end_comments)
      mutate(., deployment_end_comments = if_else(survived_beyond_study, "survived beyond study", deployment_end_comments))
      else . } %>%
    { if (has_deployment_end_type)
      mutate(., deployment_end_type = if_else(survived_beyond_study, "survived beyond study", deployment_end_type))
      else . } %>%
    { if (has_mortality_location)
      mutate(., mortality_location_filled = if_else(survived_beyond_study, 0L, mortality_location_filled))
      else . } %>%
    
    # A. death_comments keywords
    { if (has_death_comments)
      mutate(., mortality_event = case_when(
        str_detect(tolower(death_comments), "\\bnot\\b")      ~ 0L,
        str_detect(tolower(death_comments), positive_pattern) ~ 1L,
        mortality_event == 1L                                  ~ 1L,
        TRUE                                                   ~ mortality_event))
      else . } %>%
    
    # B. deployment_end_comments keywords
    { if (has_deployment_end_comments)
      mutate(., mortality_event = case_when(
        mortality_event == 1L ~ 1L,
        str_detect(tolower(deployment_end_comments), positive_pattern) ~ 1L,
        TRUE ~ mortality_event))
      else . } %>%
    
    # C. mortality_type is filled (non-NA)
    { if (has_mortality_type)
      mutate(., mortality_event = case_when(
        mortality_event == 1L  ~ 1L,
        !is.na(mortality_type) ~ 1L,
        TRUE                   ~ mortality_event))
      else . } %>%
    
    # D. mortality_location_filled
    { if (has_mortality_location)
      mutate(., mortality_event = case_when(
        mortality_location_filled >= 1 ~ 1L,
        mortality_event == 1L          ~ 1L,
        TRUE                           ~ mortality_event))
      else . } %>%
    
    # E. deployment_end_type
    { if (has_deployment_end_type)
      mutate(., mortality_event = case_when(
        mortality_event == 1L ~ 1L,
        str_detect(tolower(deployment_end_type), "\\bdead\\b|\\bdeath\\b")           ~ 1L,
        tolower(deployment_end_type) %in%
          c("removal", "other", "unknown", "survived beyond study")                  ~ 0L,
        is.na(deployment_end_type) & is.na(mortality_event)                          ~ 0L,
        TRUE                                                                         ~ mortality_event))
      else
        mutate(., mortality_event = if_else(is.na(mortality_event), 0L, mortality_event))
    } %>%
    
    # Final censoring: remaining NA → 0 if deploy_off_timestamp is present
    mutate(mortality_event = if_else(
      is.na(mortality_event) & !is.na(deploy_off_timestamp),
      0L,
      mortality_event)) %>%
    
    select(-survived_beyond_study) |>
    relocate(mortality_event, .after = any_of("deployment_end_type"))
  
  # Error out: no deaths
  n_mort_events <- sum(summary_table$mortality_event == 1L, na.rm = TRUE)
  if (n_mort_events == 0) {
    logger.fatal("Cannot run survival analysis: no mortality events detected.")
  }
  
  # Warning: few deaths
  if (n_mort_events <= 10) {
    logger.warn(sprintf("Few (%d) deaths detected across entire dataset. Particularly if data is further subset, model may have low statistical power. This could potentially result in unreliable estimates and poor predictive capacity", n_mort_events))
  }
  
  
  ## Calculate survival years (if selected) -----------------------------------
  
  if(!is.null(survival_yr_start)){
    
    # Extract survival year start date 
    start_month <- month(survival_yr_start)
    start_day   <- mday(survival_yr_start)     
    
    # Function to handle invalid dates (e.g., Feb 29 in non-leap years)
    safe_make_date <- function(year, month, day) {
      date <- suppressWarnings(make_date(year, month, day))
      if (is.na(date)) {
        ym <- ymd(sprintf("%d-%02d-01", year, month)) %m-% months(1)
        date <- ceiling_date(ym, "month") - days(1)
      }
      date
    }
    
    # Function to assign survival year and period boundaries 
    get_survival_period <- function(date) {
      if (is.na(date)) return(tibble(survival_year = NA_integer_, 
                                     period_start = NA_Date_, 
                                     period_end = NA_Date_))
      y <- year(date)
      period_start_this_year <- safe_make_date(y, start_month, start_day)
      
      if (date >= period_start_this_year) {
        survival_year <- y
        period_start  <- period_start_this_year
        period_end    <- safe_make_date(y + 1, start_month, start_day) - days(1)
      } else {
        survival_year <- y - 1
        period_start  <- safe_make_date(y - 1, start_month, start_day)
        period_end    <- safe_make_date(y, start_month, start_day) - days(1)
      }
      
      tibble(survival_year = survival_year,
             period_start  = period_start,
             period_end    = period_end)
    }
    
    # Vectorized helpers
    get_survival_year  <- function(date) get_survival_period(date)$survival_year
    
    # Determine range of survival years present in the data 
    date_range <- events_with_ind %>%
      summarise(min_ts = min(timestamp, na.rm = TRUE),
                max_ts = max(timestamp, na.rm = TRUE))
    
    min_year       <- get_survival_year(date_range$min_ts)
    max_year       <- get_survival_year(date_range$max_ts)
    possible_years <- seq(min_year, max_year, by = 1)
    
    logger.info(sprintf("Survival years found in data: %d to %d (%d years total)",
                        min_year, max_year, length(possible_years)))
    
    # Create all possible survival periods
    possible_periods <- tibble(survival_year = possible_years) %>%
      mutate(period_info = map(survival_year, ~ {
        start <- safe_make_date(.x, start_month, start_day)
        end   <- safe_make_date(.x + 1, start_month, start_day) - days(1)
        tibble(period_start = start, period_end = end)
      })) %>%
      unnest(period_info)
    
    # Infer which columns to carry into yearly_survival
    has_mortality_date          <- "mortality_date"          %in% names(summary_table)
    has_mortality_type          <- "mortality_type"          %in% names(summary_table)
    has_death_comments          <- "death_comments"          %in% names(summary_table)
    has_deployment_end_comments <- "deployment_end_comments" %in% names(summary_table)
    has_deployment_end_type     <- "deployment_end_type"     %in% names(summary_table)
    
    # Exclude join keys and period-specific columns that are recalculated per year
    cols_exclude_from_carry <- c(
      # Join keys
      "individual_id", "individual_local_identifier",
      # Recalculated per survival year
      "entry_time_days", "exit_time_days",
      "analysis_entry_date", "analysis_exit_date",
      # Raw timestamps kept separately — clipped versions used in yearly output
      "raw_deploy_on_timestamp", "raw_deploy_off_timestamp"
    )
    
    carry_cols <- setdiff(names(summary_table), cols_exclude_from_carry)
    logger.info(paste("Columns carried into yearly_survival:", paste(carry_cols, collapse = ", ")))
    
    # Infer individual-level constants by checking n_distinct per individual
    static_cols <- summary_table %>%
      dplyr::select(individual_id, all_of(carry_cols)) %>%
      group_by(individual_id) %>%
      summarise(across(everything(), n_distinct), .groups = "drop") %>%
      dplyr::select(-individual_id) %>%
      summarise(across(everything(), max)) %>%
      pivot_longer(everything(), names_to = "col", values_to = "n_distinct") %>%
      filter(n_distinct == 1) %>%
      pull(col)
    
    logger.info(paste("Individual-level constant columns (used for crossing):",
                      paste(static_cols, collapse = ", ")))
    
    # Build individual × year table
    yearly_survival <- summary_table %>%
      
      # Cross individual-level constants with all possible survival periods
      dplyr::select(individual_id, individual_local_identifier,
                    any_of(static_cols)) %>%
      distinct() %>%
      crossing(possible_periods) %>%
      
      # Join back all remaining carry columns
      left_join(summary_table %>%
                  dplyr::select(individual_id, individual_local_identifier,
                                all_of(setdiff(carry_cols, static_cols))),
                by = c("individual_id", "individual_local_identifier")) %>%
      
      # Clip monitoring interval to the survival period
      mutate(monitor_start    = pmax(deploy_on_timestamp, period_start, na.rm = TRUE),
             monitor_end      = pmin(deploy_off_timestamp, period_end,  na.rm = TRUE),
             active_in_period = monitor_start <= monitor_end &
               !is.na(monitor_start) & !is.na(monitor_end)) %>%
      filter(active_in_period) %>%
      
      # Inject NA columns for any optional mortality columns that are absent,
      # so that case_when() can reference them without erroring
      { if (!has_mortality_date)          mutate(., mortality_date          = NA_Date_)     else . } %>%
      { if (!has_mortality_type)          mutate(., mortality_type          = NA_character_) else . } %>%
      { if (!has_death_comments)          mutate(., death_comments          = NA_character_) else . } %>%
      { if (!has_deployment_end_comments) mutate(., deployment_end_comments = NA_character_) else . } %>%
      { if (!has_deployment_end_type)     mutate(., deployment_end_type     = NA_character_) else . } %>%
      
      # Mortality flags per survival year
      mutate(
        entry_date = monitor_start,
        exit_date  = monitor_end,
        
        # Did mortality occur inside this survival year?
        # Priority 1: mortality_date present and inside the period
        # Priority 2: mortality_date absent or NA — fall back to deploy_off_timestamp
        died_this_year = case_when(
          mortality_event == 1L &
            !is.na(mortality_date) &
            mortality_date >= period_start &
            mortality_date <= period_end
          ~ TRUE,
          
          mortality_event == 1L &
            is.na(mortality_date) &
            !is.na(deploy_off_timestamp) &
            deploy_off_timestamp >= period_start &
            deploy_off_timestamp <= period_end
          ~ TRUE,
          
          TRUE ~ FALSE),
        
        mortality_event = as.integer(died_this_year),
        censored        = as.integer(!died_this_year),
        
        # Reported mortality date: prefer mortality_date if present, else deploy_off
        mortality_date_reported = case_when(
          died_this_year & !is.na(mortality_date) ~ as.Date(mortality_date),
          died_this_year                           ~ as.Date(deploy_off_timestamp),
          TRUE                                     ~ NA_Date_),
        
        # Carry mortality metadata only when flagging a death
        mortality_type          = if_else(died_this_year, mortality_type,          NA_character_),
        death_comments          = if_else(died_this_year, death_comments,          NA_character_),
        deployment_end_comments = if_else(died_this_year, deployment_end_comments, NA_character_),
        deployment_end_type     = if_else(died_this_year, deployment_end_type,     NA_character_),
        
        days_at_risk = as.integer(exit_date - entry_date) + 1L) %>%
      
      select(-active_in_period, -monitor_start, -monitor_end, -died_this_year) %>%
      
      # Calculate per-year entry/exit time in days (relative to each period_start)
      mutate(analysis_entry_date = entry_date,
             analysis_exit_date  = exit_date,
             entry_time_days     = as.numeric(difftime(entry_date, period_start, units = "days")),
             exit_time_days      = as.numeric(difftime(exit_date,  period_start, units = "days"))) %>%
      
      arrange(individual_id, survival_year)
    
  } else {
    logger.info("Survival years not calculated.")
  }
  
  
  ## Calculate life stages per year (if selected) -----------------------------
  
  # Note: this needs auxiliary file to be loaded (errors earlier in code upon loading) 
  # Note: this needs "survival_yr_start" to be defined 
  if(!is.null(animal_birth_hatch_year_table) && is.null(survival_yr_start)){
    logger.error("Calculating life-stage requires survival years to be defined.")
  } 
  
  if(!is.null(survival_yr_start) && !is.null(animal_birth_hatch_year_table)){
    
    # Confirm data exists 
    if (!"animal_birth_hatch_year" %in% names(yearly_survival)) {
      logger.error("Column 'animal_birth_hatch_year' does not exist in the data frame. Cannot compute life stage.")
    }
    
    # Calculate age and age_class
    yearly_survival <- yearly_survival %>%
      mutate(animal_birth_hatch_year = as.numeric(animal_birth_hatch_year),
             age = survival_year - animal_birth_hatch_year,
             age = as.integer(pmax(0, age)))
    
    # Prepare thresholds from your existing table  
    thresholds <- animal_birth_hatch_year_table %>%
      filter(!is.na(year_at_start)) %>%           
      arrange(year_at_start) %>%
      distinct(year_at_start, animal_life_stage)
    
    # Create a named vector for fast look-up
    stage_lookup <- setNames(thresholds$animal_life_stage,
                             thresholds$year_at_start)
    
    # Pull stage for animals with unknown birth year from lookup table
    na_stage <- animal_birth_hatch_year_table %>%
      filter(is.na(year_at_start)) %>%
      pull(animal_life_stage) %>%
      first(default = "adult")
    
    # Dynamic assignment 
    yearly_survival <- yearly_survival %>%
      mutate(matched_threshold = findInterval(age, thresholds$year_at_start),
             animal_life_stage_new = case_when(
               is.na(age)             ~ na_stage,           
               matched_threshold == 0 ~ "unknown",            
               TRUE                   ~ stage_lookup[matched_threshold]),
             animal_life_stage_new = coalesce(
               animal_life_stage_new,
               na_stage)) %>%
      select(-matched_threshold)
    
    # Log age class summary
    age_class_summary <- yearly_survival %>%
      count(animal_life_stage_new) %>%
      arrange(desc(n))
    
    logger.info(sprintf("Individuals by age class: %s",
                        paste(sprintf("%s (n=%d)", age_class_summary$animal_life_stage_new,
                                      age_class_summary$n),
                              collapse = ", ")))
    
  } else {
    logger.info("Life stages not calculated.")
  }
  
  
  ## Subset based on condition (if selected) ----------------------------------
  
  # Lookup: condition name -> list(col, label, coerce_fn, summary_ok, yearly_ok)
  subset_spec <- list(
    sex = list(col = "sex", label = "sex", coerce = identity, summary_ok = TRUE,  yearly_ok = TRUE),
    attachment_type = list(col = "attachment_type", label = "attachment type", coerce = identity, summary_ok = TRUE,  yearly_ok = TRUE),
    model = list(col = "model", label = "model", coerce = as.integer, summary_ok = TRUE, yearly_ok = TRUE),
    lifestage = list(col = "animal_life_stage_new", label = "life stage", coerce = identity, summary_ok = FALSE, yearly_ok = TRUE),
    survival_year = list(col = "survival_year", label = "survival year", coerce = as.integer, summary_ok = FALSE, yearly_ok = TRUE)
  )
  
  # Helper function 
  apply_subset <- function(condition, define, summary_table, yearly_survival, survival_yr_start) {
    if (is.null(condition)) return(list(summary_table = summary_table, yearly_survival = yearly_survival))
    
    spec <- subset_spec[[condition]]
    if (is.null(spec)) stop(paste("Unknown subset condition:", condition))
    
    logger.info(paste0("Subsetting by ", spec$label, " (", define, ")"))
    value <- spec$coerce(define)
    
    if (is.null(survival_yr_start)) {
      if (!spec$summary_ok) {
        logger.fatal("This subset only makes sense when data are processed by survival year. Please enter survival year start date.")
      } else {
        summary_table <- summary_table %>% filter(.data[[spec$col]] == value)
      }
    } else {
      if (!spec$yearly_ok) {
        logger.fatal("This subset only makes sense when data are processed by survival year. Please enter survival year start date.")
      } else {
        yearly_survival <- yearly_survival %>% filter(.data[[spec$col]] == value)
      }
    }
    
    list(summary_table = summary_table, yearly_survival = yearly_survival)
  }
  
  # SUBSET CONDITION 1 ---
  result <- apply_subset(subset_condition_1, subset_condition_define_1, summary_table,
                         yearly_survival, survival_yr_start)
  summary_table    <- result$summary_table
  yearly_survival  <- result$yearly_survival
  
  if (!is.null(subset_condition_1)) {
    if (is.null(survival_yr_start) && nrow(summary_table) == 0) {
      logger.fatal("There are no individuals meeting the first subsetting condition.")
    } else if (!is.null(survival_yr_start) && nrow(yearly_survival) == 0) {
      logger.fatal("There are no individuals meeting the first subsetting condition.")
    }
  }
  
  # SUBSET CONDITION 2 ---
  result <- apply_subset(subset_condition_2, subset_condition_define_2, summary_table,
                         yearly_survival, survival_yr_start)
  summary_table   <- result$summary_table
  yearly_survival <- result$yearly_survival
  
  if (!is.null(subset_condition_2)) {
    if (is.null(survival_yr_start) && nrow(summary_table) == 0) {
      logger.fatal("There are no individuals meeting both subsetting conditions.")
    } else if (!is.null(survival_yr_start) && nrow(yearly_survival) == 0) {
      logger.fatal("There are no individuals meeting both subsetting conditions.")
    }
  }
  
  
  ## Clean user-defined grouping attributes (if selected) ---------------------
  
  # Note: This sets group_comparison_individual to NULL if no cleaned values are available 
  
  if(!is.null(group_comparison_individual)) {
    logger.info(paste("Comparing across", group_comparison_individual, "..."))
    
    if(is.null(survival_yr_start)){
      logger.info("... with summary table")
      
      non_na_unique <- unique(na.omit(summary_table[[group_comparison_individual]]))
      
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          logger.fatal(paste0("Warning: The grouping variable, ", group_comparison_individual, ", is entirely NA; no comparison is possible."))
        } else {
          logger.fatal(paste0("Warning: There is only one non-NA comparison covariate in ", group_comparison_individual, "; no comparison is possible."))
        }
        group_comparison_individual <- NULL
      }
      
      if (!is.null(group_comparison_individual)) {
        n_original    <- nrow(summary_table)
        summary_table <- summary_table[!is.na(summary_table[[group_comparison_individual]]) & 
                                         trimws(summary_table[[group_comparison_individual]]) != "", ]
        
        unique_values <- sort(unique(summary_table[[group_comparison_individual]]))
        n_values      <- length(unique_values)
        
        logger.info(sprintf("%d values of comparison covariate detected after cleaning: %s", n_values, paste(unique_values, collapse = ", ")))
        
        n_lost <- n_original - nrow(summary_table)
        if (n_lost > 0) {
          logger.info(sprintf("%d individuals with NA covariate value removed from study.", n_lost))
        }
      }
      
    } else {
      logger.info("... with yearly survival")
      
      non_na_unique <- unique(na.omit(yearly_survival[[group_comparison_individual]]))
      
      if (length(non_na_unique) <= 1) {
        if (length(non_na_unique) == 0) {
          logger.fatal(paste0("Warning: The grouping variable, ", group_comparison_individual, ", is entirely NA; no comparison is possible."))
        } else {
          logger.fatal(paste0("Warning: There is only one non-NA comparison covariate in ", group_comparison_individual, "; no comparison is possible."))
        }
        group_comparison_individual <- NULL
      }
      
      if (!is.null(group_comparison_individual)) {
        n_original      <- nrow(yearly_survival)
        yearly_survival <- yearly_survival[!is.na(yearly_survival[[group_comparison_individual]]) & 
                                             trimws(yearly_survival[[group_comparison_individual]]) != "", ]
        
        unique_values <- sort(unique(yearly_survival[[group_comparison_individual]]))
        n_values      <- length(unique_values)
        
        logger.info(sprintf("%d values of comparison covariate detected after cleaning: %s", n_values, paste(unique_values, collapse = ", ")))
        
        n_lost <- n_original - nrow(yearly_survival)
        if (n_lost > 0) {
          logger.info(sprintf("%d individuals with NA covariate value removed from study.", n_lost))
        }
      }
    }
  }
  
  
  ## Basic summaries of data ------------------------------------------------
  
  # Plot each individual's tracking history --- 
  
  if(calc_tracking_history == TRUE) {
    
    # Helper (base R, no extra packages)
    clamp <- function(x, lo, hi) max(lo, min(hi, x))
    
    if (is.null(survival_yr_start)) {
      logger.info("Plotting tracking history using summary table...")
      
      # Create deployment summary
      deployment_summary <- summary_table |>
        mutate(deploy_on  = as.POSIXct(deploy_on_timestamp),
               deploy_off = as.POSIXct(deploy_off_timestamp),
               duration_days = round(as.numeric(difftime(deploy_off, deploy_on, units = "days")), 1)) |>
        
        filter(deploy_off > deploy_on,
               !is.na(deploy_on),
               !is.na(deploy_off)) |>
        
        group_by(individual_id, individual_local_identifier) |>          
        mutate(first_start = min(deploy_on, na.rm = TRUE) ) |>
        ungroup() |>
        
        mutate(individual_label = fct_reorder(
          paste(individual_id, individual_local_identifier, sep = " – "),
          first_start)) |>
        
        arrange(first_start, deploy_on) |>
        mutate(plot_start = deploy_on,
               plot_end   = deploy_off)
      
      # Total location count
      n_locs_total <- summary_table |>
        summarise(total = sum(n_locations, na.rm = TRUE)) |>
        pull(total)
      
      # Gap detection
      deployment_summary_with_gaps <- deployment_summary |>
        group_by(individual_label) |>
        arrange(plot_start) |>
        mutate(prev_end  = lag(plot_end),
               gap_start = prev_end,
               gap_end   = plot_start,
               gap_days  = as.numeric(difftime(gap_end, gap_start, units = "days"))) |>
        filter(gap_days > 30, !is.na(gap_days)) |>
        ungroup()
      
      # Build the plot  
      tracking_history <- ggplot(deployment_summary) +
        geom_segment(aes(x = plot_start, xend = plot_end,
                         y = individual_label, yend = individual_label),
                     linewidth = 3.2, color = "grey") +
        geom_point(aes(x = plot_start, y = individual_label),
                   color = "#1F77B4", size = 3.5) +
        geom_point(aes(x = plot_end, y = individual_label),
                   color = "#9467BD", size = 3.5) +
        geom_segment(data = deployment_summary_with_gaps,
                     aes(x = gap_start + (gap_end - gap_start)/2,
                         xend = gap_start + (gap_end - gap_start)/2,
                         y = as.numeric(individual_label),
                         yend = as.numeric(individual_label) + 0.45),
                     color = "grey50", linewidth = 1.2,
                     arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
        labs(title = "Individual Collared Periods",
             subtitle = sprintf("%d unique individuals • %d visible deployments • %d locations",
                                n_distinct(deployment_summary$individual_id),
                                nrow(deployment_summary),
                                n_locs_total),
             x = "Time",
             y = "Individual") +
        theme_minimal(base_size = 12) +
        theme(axis.text.y = element_text(size = 8, face = "plain"),
              panel.grid.major.y = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(face = "bold", size = 14),
              plot.subtitle = element_text(size = 11, color = "grey50", margin = margin(b = 10)),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 12)) +
        scale_x_datetime(date_breaks = "1 year",
                         date_labels = "%Y",
                         expand = expansion(mult = c(0.01, 0.03)))
      
      # Dynamic dimensions based on data
      n_individuals <- n_distinct(deployment_summary$individual_id)
      time_range_years <- as.numeric(difftime(max(deployment_summary$plot_end), 
                                              min(deployment_summary$plot_start), 
                                              units = "days")) / 365.25
      
      # Height: base + per-individual allowance, clamped to reasonable bounds
      plot_height <- clamp(2 + n_individuals * 0.18, 5, 40)
      
      # Width: base + per-year allowance, clamped to reasonable bounds
      plot_width  <- clamp(4 + time_range_years * 1.2, 7, 24)
      
      # Save plot
      png(appArtifactPath("tracking_history.png"),
          width  = plot_width,
          height = plot_height,
          units  = "in", res = 300)
      print(tracking_history)
      dev.off()
      
    }
    
    if (!is.null(survival_yr_start)) {
      logger.info("Plotting tracking history using yearly survival")
      
      deployment_summary <- yearly_survival |>
        distinct(individual_id,
                 individual_local_identifier,
                 deploy_on_timestamp,
                 deploy_off_timestamp) |>
        mutate(deploy_on  = as.POSIXct(deploy_on_timestamp),
               deploy_off = as.POSIXct(deploy_off_timestamp),
               duration_days = round(as.numeric(difftime(deploy_off, deploy_on, units = "days")), 1)) |>
        filter(deploy_off > deploy_on,
               !is.na(deploy_on),
               !is.na(deploy_off)) |>
        group_by(individual_id) |>
        mutate(first_start = min(deploy_on, na.rm = TRUE),
               individual_label = fct_reorder(paste(individual_id, individual_local_identifier, sep = " – "),
                                              first_start)) |>
        ungroup() |>
        arrange(first_start, deploy_on) |>
        mutate(plot_start = deploy_on,
               plot_end   = deploy_off)
      
      # Total location count 
      n_locs_total <- yearly_survival |>
        distinct(individual_id, n_locations) |>
        summarise(total = sum(n_locations, na.rm = TRUE)) |>
        pull(total)
      
      # Plot 
      tracking_history <- ggplot(deployment_summary) +
        geom_segment(aes(x = plot_start, xend = plot_end,
                         y = individual_label, yend = individual_label),
                     linewidth = 3.2, color = "grey") +
        geom_point(aes(x = plot_start, y = individual_label),
                   color = "#1F77B4", size = 3.5) +
        geom_point(aes(x = plot_end, y = individual_label),
                   color = "#9467BD", size = 3.5) +
        geom_segment(data = deployment_summary |>
                       group_by(individual_label) |>
                       arrange(plot_start) |>
                       mutate(prev_end  = lag(plot_end),
                              gap_start = prev_end,
                              gap_end   = plot_start,
                              gap_days  = as.numeric(difftime(gap_end, gap_start, units = "days"))) |>
                       filter(gap_days > 30, !is.na(gap_days)),
                     aes(x = gap_start + (gap_end - gap_start)/2,
                         xend = gap_start + (gap_end - gap_start)/2,
                         y = as.numeric(individual_label),
                         yend = as.numeric(individual_label) + 0.45),
                     color = "grey50", linewidth = 1.2,
                     arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
        labs(title    = paste0("Individual Collared Periods"),
             subtitle = sprintf("%d unique individuals • %d visible deployments • %d locations",
                                n_distinct(deployment_summary$individual_id),
                                nrow(deployment_summary),
                                n_locs_total),
             x = "Time",
             y = "Individual") +
        theme_minimal(base_size = 12) +
        theme(axis.text.y        = element_text(size = 8, face = "plain"),
              panel.grid.major.y = element_blank(),
              panel.grid.minor   = element_blank(),
              plot.title         = element_text(face = "bold", size = 14),
              plot.subtitle      = element_text(size = 11, color = "grey50", margin = margin(b = 10)),
              axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title         = element_text(size = 12)) +
        scale_x_datetime(date_breaks = "1 year",
                         date_labels = "%Y",
                         expand      = expansion(mult = c(0.01, 0.03)))
      
      # Dynamic dimensions based on data
      n_individuals <- n_distinct(deployment_summary$individual_id)
      time_range_years <- as.numeric(difftime(max(deployment_summary$plot_end), 
                                              min(deployment_summary$plot_start), 
                                              units = "days")) / 365.25
      
      # Height: base + per-individual allowance, clamped to reasonable bounds
      plot_height <- clamp(2 + n_individuals * 0.18, 5, 40)
      
      # Width: base + per-year allowance, clamped to reasonable bounds
      plot_width  <- clamp(4 + time_range_years * 1.2, 7, 24)
      
      # Save plot
      png(appArtifactPath("tracking_history.png"),
          width  = plot_width,
          height = plot_height,
          units  = "in", res = 300)
      print(tracking_history)
      dev.off()
      
    } 
  }
  
  
  # Calculate monthly mortality ---- 
  if (calc_month_mort == TRUE) {
    
    if (is.null(survival_yr_start)) {
      logger.info("Plotting monthly mortality using summary table...")
      
      min_date <- min(summary_table$deploy_on_timestamp, na.rm = TRUE)
      max_date <- max(summary_table$deploy_off_timestamp, na.rm = TRUE)
      full_years <- seq(year(min_date), year(max_date), by = 1)
      
      mortality_data <- summary_table |>
        filter(mortality_event == 1) |>
        mutate(mort_date   = as.Date(mortality_date),
               deploy_off  = as.Date(deploy_off_timestamp),
               death_date  = coalesce(mort_date, deploy_off),
               death_year  = year(death_date),
               death_month = month(death_date, label = TRUE, abbr = TRUE),
               death_month = factor(death_month, levels = month.abb, ordered = TRUE)) |>
        dplyr::select(death_year, death_month)
      
      monthly_morts <- mortality_data %>%
        count(death_year, death_month, name = "n_mortalities") %>%
        complete(death_year  = full_years,
                 death_month = factor(month.abb, levels = month.abb, ordered = TRUE),
                 fill = list(n_mortalities = 0)) %>%
        mutate(death_month_num = as.integer(death_month),
               death_month     = fct_relevel(death_month, month.abb))
      
      # Plot 
      monthly_mort_plot <- ggplot(monthly_morts, aes(x = death_month, y = factor(death_year), 
                                                     fill = factor(n_mortalities))) +
        geom_tile(color = "white", linewidth = 0.5) +
        scale_fill_brewer(palette  = "PuBu",
                          na.value = "grey92",
                          name     = "Number of\nmortality events",
                          drop     = FALSE) + 
        scale_x_discrete(position = "top") +
        labs(title    = "Monthly Distribution of Confirmed Mortality Events",
             subtitle = paste0(
               "Total events: ", sum(monthly_morts$n_mortalities, na.rm = TRUE),
               " • Time span: ", format(min_date, "%b %Y"),
               " to ", format(max_date, "%b %Y")),
             x = NULL,
             y = "Year") +
        theme_minimal(base_size = 14) +
        theme(panel.grid        = element_blank(),
              axis.ticks        = element_blank(),
              legend.position   = "right",
              legend.title      = element_text(size = 11),
              legend.text       = element_text(size = 10),
              plot.title        = element_text(face = "bold", hjust = 0, size = 16),
              plot.subtitle     = element_text(hjust = 0, size = 12),
              axis.text.x       = element_text(size = 11, face = "bold"),
              axis.text.y       = element_text(size = 11))
      
      # Save plot  
      png(appArtifactPath("monthly_mortality.png"), 
          width = 10, height = 8, 
          units = "in", res = 300)
      print(monthly_mort_plot)
      dev.off()
    }
    
    if (!is.null(survival_yr_start)) {
      
      logger.info("Plotting monthly mortality using yearly survival")
      
      min_date <- min(yearly_survival$deploy_on_timestamp, na.rm = TRUE)
      max_date <- max(yearly_survival$deploy_off_timestamp, na.rm = TRUE)
      full_years <- seq(year(min_date), year(max_date), by = 1)
      
      mortality_data <- yearly_survival |>
        filter(mortality_event == 1) |>
        mutate(mort_date   = as.Date(mortality_date),
               deploy_off  = as.Date(deploy_off_timestamp),
               death_date  = coalesce(mort_date, deploy_off),
               death_year  = year(death_date),
               death_month = month(death_date, label = TRUE, abbr = TRUE),
               death_month = factor(death_month, levels = month.abb, ordered = TRUE)) |>
        dplyr::select(death_year, death_month)
      
      monthly_morts <- mortality_data |>
        count(death_year, death_month, name = "n_mortalities") |>
        complete( death_year  = full_years,                                       
                  death_month = factor(month.abb, levels = month.abb, ordered = TRUE),
                  fill = list(n_mortalities = 0)) |>
        mutate(death_month_num = as.integer(death_month),
               death_month     = fct_relevel(death_month, month.abb))
      
      # Plot 
      monthly_mort_plot <- ggplot(monthly_morts, 
                                  aes(x = death_month, 
                                      y = factor(death_year),
                                      fill = factor(n_mortalities))) +
        geom_tile(color = "white", linewidth = 0.5) +
        scale_fill_brewer(palette  = "PuBu",
                          na.value = "grey92",
                          name     = "Number of\nmortality events",
                          drop     = FALSE) + 
        scale_x_discrete(position = "top") +
        labs(title = "Monthly Distribution of Confirmed Mortality Events",
             subtitle = paste0(
               "Total events: ", sum(monthly_morts$n_mortalities, na.rm = TRUE),
               " • Time span: ", format(min_date, "%b %Y"),
               " to ", format(max_date, "%b %Y")),
             x = NULL,
             y = "Year") +
        theme_minimal(base_size = 14) +
        theme(panel.grid = element_blank(),
              axis.ticks = element_blank(),
              legend.position = "right",
              legend.title = element_text(size = 11),
              legend.text = element_text(size = 10),
              plot.title = element_text(face = "bold", hjust = 0, size = 16),
              plot.subtitle = element_text(hjust = 0, size = 12),
              axis.text.x = element_text(size = 11, face = "bold"),
              axis.text.y = element_text(size = 11))
      
      # Save plot  
      png(appArtifactPath("monthly_mortality.png"), 
          width = 10, height = 8, 
          units = "in", res = 300)
      print(monthly_mort_plot)
      dev.off()
    }
  }
  
  
  ## Survival Analysis: No comparisons ----------------------------------------
  
  if (is.null(survival_yr_start)) {
    
    logger.info("Calculating KM using summary table...")
    
    # Warning for not mortality 
    if(sum(summary_table$mortality_event) == 0){
      logger.warn("There are no mortality events in the chosen subset of data.")
    }
    
    # Fit model 
    fitting_data <- summary_table
    km_fit <- survfit(Surv(entry_time_days, exit_time_days, mortality_event) ~ 1, 
                      data = summary_table)
    
    # Data for life table 
    lt.length.out <- ceiling(max(summary_table$exit_time_days)/life_table_days)
    times <- round(seq(min(summary_table$entry_time_days), max(summary_table$exit_time_days), 
                       length.out = lt.length.out))
    
  } 
  
  if(!is.null(survival_yr_start)) {
    
    logger.info("Calculating KM using survival table...")
    
    # Warning for no mortality 
    if(sum(yearly_survival$mortality_event) == 0){
      logger.warn("There are no mortality events in the chosen subset of data.")
    }
    
    # Fit model 
    fitting_data <- yearly_survival
    km_fit <- survfit(Surv(days_at_risk, mortality_event) ~ 1, data = yearly_survival)
    
    # Data for life table 
    max_days <- max(yearly_survival$days_at_risk, na.rm = TRUE)
    lt.length.out <- ceiling(max_days / life_table_days) + 1   
    times <- round(seq(0, max_days, length.out = lt.length.out))
  }
  
  # Generate life table ---
  s <- summary(km_fit, times = times)
  life_table <- data.frame(time_days        = s$time,
                           n_risk           = s$n.risk,
                           n_event          = s$n.event,
                           survival_prob    = s$surv,
                           std_err          = s$std.err,
                           lower_95         = s$lower,
                           upper_95         = s$upper)
  
  # Save 
  write.csv(life_table, file = appArtifactPath("life_table.csv"), row.names = F)
  
  
  # Annual survival probabilities ---
  if (!is.null(survival_yr_start)) {
    
    logger.info("Calculating annual survival probabilities per survival year...")
    
    # Initialize list to store per-year survfit objects
    survfit_list <- list()
    
    # Loop through each survival year
    for (yr in unique(yearly_survival$survival_year)) {
      data_year <- subset(yearly_survival, survival_year == yr)
      
      fit <- survfit(Surv(days_at_risk, mortality_event) ~ 1, data = data_year)
      survfit_list[[as.character(yr)]] <- fit
    }
    
    # Extract end-of-year survival probability from each survfit object
    annual_surv_df <- data.frame(
      Year      = integer(),
      AnnualSurv = numeric(),
      SE        = numeric(),
      LCI       = numeric(),
      UCI       = numeric(),
      n         = integer()
    )
    
    for (yr_chr in names(survfit_list)) {
      fit <- survfit_list[[yr_chr]]
      s   <- summary(fit)
      idx <- length(s$surv)
      
      if (idx == 0) {
        # No mortality events this year — survival stayed at 1.0
        annual_surv_df <- rbind(annual_surv_df, data.frame(
          Year       = as.integer(yr_chr),
          AnnualSurv = 1.0,
          SE         = NA_real_,
          LCI        = NA_real_,
          UCI        = NA_real_,
          n          = fit$n))
      } else {
        annual_surv_df <- rbind(annual_surv_df, data.frame(
          Year       = as.integer(yr_chr),
          AnnualSurv = s$surv[idx],
          SE         = s$std.err[idx],
          LCI        = s$lower[idx],
          UCI        = s$upper[idx],
          n          = fit$n))
      }
    }
    
    annual_surv_df <- annual_surv_df[order(annual_surv_df$Year), ]
    
    write.csv(annual_surv_df, 
              file = appArtifactPath("annual_survival_by_year.csv"), 
              row.names = FALSE)
  }
  
  
  # Plot KM curve ---
  n.ind <- nrow(fitting_data)
  n.events <- nrow(fitting_data[fitting_data$mortality_event == 1,])
  n.days <- as.numeric(summary(km_fit)$table["median"])
  med <- survminer::surv_median(km_fit)
  
  km_curve <- ggsurvplot(
    km_fit,
    data = fitting_data,
    title = "Kaplan-Meier Survival Curve",
    subtitle = paste0("N = ", n.ind, ", Events = ", n.events, ", Median Survival = ", 
                      med$median, " days"),
    xlab = "Time (days)",
    ylab = "Survival Probability",
    risk.table = FALSE,
    conf.int = add_cis,
    censor.shape = "|",
    censor.size = 3,
    legend = "none",
    pval = FALSE,
    surv.median.line = "hv",        
    palette = c("#E69F00", "#56B4E9"),
    ggtheme = theme_classic(base_size = 12) + 
      theme(plot.title         = element_text(face = "bold", size = 14), 
            plot.subtitle      = element_text(size = 12, color = "gray50"),
            axis.text          = element_text(color = "black"),
            panel.grid.major.y = element_line(color = "gray90"), 
            panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
            line               = element_line(linewidth = 0.1),
            plot.margin        = margin(10, 10, 10, 10)))
  
  # Zoom y-axis to data range
  if (zoom_to_plot) {
    y_min   <- min(km_fit$surv, na.rm = TRUE)
    y_floor <- floor(y_min * 10) / 10
    km_curve$plot <- km_curve$plot + 
      coord_cartesian(ylim = c(y_floor, 1)) +
      annotate("text",
               x     = -Inf,
               y     = y_floor,
               label = paste0("* y-axis zoomed to [", round(y_floor, 2), ", 1]"),
               hjust = -0.05,
               vjust = -0.5,
               size  = 3.2,
               color = "gray50",
               fontface = "italic")
  }
  
  # Save plot  
  png(appArtifactPath("km_survival_curve.png"), width = 10, height = 8, units = "in", res = 300)
  print(km_curve)
  dev.off()
  
  
  # Plot cumulative hazard curve --- 
  cum_hazard <- ggsurvplot(
    km_fit,
    data = fitting_data,
    fun = "cumhaz",
    conf.int = add_cis,
    risk.table = FALSE,
    cumevents = FALSE,                 
    pval = FALSE,                     
    xlab = "Time (days)",
    ylab = "Cumulative Hazard",
    title = "Cumulative Hazard",
    subtitle = paste0("N = ", n.ind, ", Events = ", n.events), 
    palette = c("#E69F00", "#56B4E9"),
    legend = "none",
    ggtheme = theme_classic(base_size = 12) + 
      theme(plot.title   = element_text(face = "bold", size = 14), 
            plot.subtitle = element_text(size = 12, color = "gray50"),
            axis.text    = element_text(color = "black"),
            panel.grid.major.y = element_line(color = "gray90"), 
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            plot.margin  = margin(10, 10, 10, 10)))
  
  if (zoom_to_plot) {
    y_max     <- max(-log(km_fit$surv), na.rm = TRUE)
    y_ceiling <- ceiling(y_max * 10) / 10
    cum_hazard$plot <- cum_hazard$plot +
      coord_cartesian(ylim = c(0, y_ceiling)) +
      annotate("text",
               x        = -Inf,
               y        = y_ceiling,
               label    = paste0("* y-axis zoomed to [0, ", round(y_ceiling, 2), "]"),
               hjust    = -0.05,
               vjust    = 1.5,
               size     = 3.2,
               color    = "gray50",
               fontface = "italic")
  }
  
  # Save plot 
  png(appArtifactPath("cumulative_hazard_plot.png"), width = 10, height = 8, units = "in", res = 300)
  print(cum_hazard)
  dev.off()
  
  
  ## Survival Analysis: Group comparisons -------------------------------------
  
  # For yearly_summary data 
  if (!is.null(group_comparison_individual) && !is.null(survival_yr_start)){
    logger.info(paste0("Comparing ", group_comparison_individual, " using yearly survival data"))
    
    # Fit survival object 
    surv_formula <- as.formula(paste("Surv(days_at_risk, mortality_event) ~", group_comparison_individual))
    km_fit_comp <- survfit(surv_formula, data = yearly_survival) 
    
    # Log-Rank test --- 
    test <- survdiff(surv_formula, data = yearly_survival)
    
    # Extract components
    groups       <- names(test$n)
    n_total      <- sum(test$n)
    events_total <- sum(test$obs)
    chisq_val    <- round(test$chisq, 2)
    df_val       <- length(test$n) - 1
    p_val        <- 1 - pchisq(test$chisq, df_val)
    p_formatted  <- ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val))
    
    # Summary table 
    per_group <- tibble(!!group_comparison_individual := sub(".*=", "", groups),
                        `N`               = test$n,
                        `Events`          = test$obs,
                        `Expected events` = round(test$exp, 2),
                        `O/E ratio`       = round(test$obs / test$exp, 2)) %>%
      mutate(`N (events)` = sprintf("%d(%d)", N, Events), .keep = "unused")
    
    summary_row <- tibble(!!group_comparison_individual := "Overall",
                          `N (events)`       = sprintf("%d(%d)", n_total, events_total),
                          `Chisq (log-rank)` = chisq_val,
                          `df`               = df_val,
                          `p-value`          = p_formatted)
    
    logrank_table <- bind_rows(per_group, summary_row)
    logrank_table <- logrank_table %>%
      dplyr::select(any_of(group_comparison_individual), `N (events)`, `Expected events`, 
                    `O/E ratio`, `Chisq (log-rank)`, df, `p-value`)
    
    # Save
    write.csv(logrank_table, file = appArtifactPath("logrank_table_statistics.csv"), row.names = F)
    
    
    # KM comparison plots ---
    km_fit_comp <- surv_fit(surv_formula, data = yearly_survival)  
    
    # Check if any groups have N == 1 and remove 
    filter_singleton_strata_refit <- function(fit, data) {
      grouping_var <- all.vars(formula(fit))[3]    
      n_per_group  <- fit$n
      group_names  <- sub(".*=", "", names(fit$strata) %||% "Overall")
      
      singletons <- n_per_group == 1
      removed    <- group_names[singletons]   
      
      if (any(singletons)) {
        logger.info(paste("Removed the following singleton group(s) (N=1):\n",
                          paste(" •", removed, collapse = "\n ")), call. = FALSE)
        data_clean <- data[!data[[grouping_var]] %in% removed, , drop = FALSE]
        km_clean   <- surv_fit(formula(fit), data = data_clean)
        return(km_clean)
      }
      return(fit)
    }
    
    km_fit_clean <- filter_singleton_strata_refit(km_fit_comp, yearly_survival)
    
    title_text     <- paste0("Kaplan-Meier Survival Curve: ", group_comparison_individual)
    group_rows     <- logrank_table %>% filter(!!sym(group_comparison_individual) != "Overall") 
    group_rows_clean <- group_rows %>% filter(!str_detect(`N (events)`, "^1\\s*\\("))
    removed_groups <- group_rows %>% filter(str_detect(`N (events)`, "^1\\s*\\(")) %>% pull(1)
    
    if (length(removed_groups) > 0) {
      logger.info(paste0("The following group(s) had N=1 and were removed from the table:\n • ",
                         paste(removed_groups, collapse = "\n • ")), call. = FALSE)
    }
    
    if (nrow(group_rows_clean) == 1){
      logger.fatal("There is only one group left; unable to perform comparisons.")
      
    } else {
      subtitle_text <- paste0(paste(sprintf("N_%s: %s", 
                                            group_rows_clean[[group_comparison_individual]],
                                            group_rows_clean$`N (events)`),
                                    collapse = ", "),
                              "\nP-value: ", 
                              logrank_table$`p-value`[logrank_table[[group_comparison_individual]] == "Overall"])
      
      # Wrap each line independently to preserve the \n break before P-value
      subtitle_text_wrapped <- paste(
        sapply(strsplit(subtitle_text, "\n")[[1]], function(line) {
          paste(strwrap(line, width = 90), collapse = "\n")
        }),
        collapse = "\n"
      )
      
      old_strata_names <- names(km_fit_clean$strata)
      new_strata_names <- sub(".*=", "", old_strata_names)
      names(km_fit_clean$strata) <- new_strata_names
      n_groups <- length(names(km_fit_clean$strata))
      
      # Okabe-Ito colorblind-friendly palette
      okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                     "#0072B2", "#D55E00", "#CC79A7", "#000000",
                     "#999999", "#44AA99")
      my_palette <- okabe_ito[seq_len(n_groups)]
      
      # Plot 
      km_comp_curve <- ggsurvplot(km_fit_clean,
                                  data          = yearly_survival,
                                  title         = title_text,
                                  subtitle      = subtitle_text_wrapped,
                                  conf.int      = add_cis,
                                  risk.table    = FALSE,
                                  palette       = my_palette, 
                                  xlab          = "Days at risk",
                                  ylab          = "Survival probability",
                                  legend.title  = group_comparison_individual,
                                  legend        = "bottom",
                                  legend.labs   = levels(yearly_survival[[group_comparison_individual]]),
                                  censor.shape  = "|",
                                  censor.size   = 4,
                                  font.main     = c(14, "bold", "black"),
                                  font.x        = 12, font.y = 12, font.tickslab = 11, 
                                  ggtheme = theme_classic(base_size = 12) +
                                    theme(plot.title         = element_text(face = "bold", size = 14),
                                          plot.subtitle      = element_text(size = 12, color = "gray50"),
                                          axis.text          = element_text(color = "black"),
                                          panel.grid.major.y = element_line(color = "gray90"),
                                          panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
                                          plot.margin        = margin(20, 10, 10, 10),
                                          legend.direction   = "horizontal",
                                          legend.box         = "horizontal",
                                          legend.spacing.x   = unit(0.3, "cm")))
      
      # Spread legend
      km_comp_curve$plot <- km_comp_curve$plot + 
        guides(color = guide_legend(nrow           = ifelse(n_groups <= 8, 1, 2),
                                    byrow          = TRUE,
                                    title.position = "top"),
               fill  = guide_legend(nrow           = ifelse(n_groups <= 4, 1, 2),
                                    byrow          = TRUE,
                                    title.position = "top"))
      
      # Zoom y-axis to data range
      if (zoom_to_plot) {
        y_min   <- min(km_fit_clean$surv, na.rm = TRUE)
        y_floor <- floor(y_min * 10) / 10
        km_comp_curve$plot <- km_comp_curve$plot +
          coord_cartesian(ylim = c(y_floor, 1)) +
          annotate("text",
                   x        = -Inf,
                   y        = y_floor,
                   label    = paste0("* y-axis zoomed to [", round(y_floor, 2), ", 1]"),
                   hjust    = -0.05,
                   vjust    = -0.5,
                   size     = 3.2,
                   color    = "gray50",
                   fontface = "italic")
      }
      
      # Save plot 
      png(appArtifactPath("km_comparison_curves.png"), width = 10, height = 8, units = "in", res = 300)
      print(km_comp_curve)
      dev.off()
    }
    
    
    ## Cumulative hazard comparison plots ---
    
    # Prepare statistics for subtitle 
    n_per_group      <- km_fit_clean$n
    sum_fit          <- summary(km_fit_clean)
    events_per_group <- tapply(sum_fit$n.event, sum_fit$strata, sum, na.rm = TRUE)
    
    # Clean strata names  
    clean_strata <- gsub("^.*=", "", names(km_fit_clean$strata))
    
    # Create subtitle for groups: "N(Group): N(events)"
    subtitle_parts <- mapply(function(group, n, ev) {
      sprintf("N_%s: %d(%d)", group, n, ev)
    },
    clean_strata, n_per_group, events_per_group, SIMPLIFY = FALSE)
    
    groups_line   <- paste(subtitle_parts, collapse = ", ")
    test          <- surv_pvalue(km_fit_clean, data = yearly_survival)
    pval_text     <- sprintf("P-value: %.3f", test$pval)
    subtitle_text <- paste0(groups_line, "\n", pval_text)
    
    # Wrap each line independently to preserve the \n break before P-value
    subtitle_text_wrapped <- paste(
      sapply(strsplit(subtitle_text, "\n")[[1]], function(line) {
        paste(strwrap(line, width = 90), collapse = "\n")
      }),
      collapse = "\n"
    )
    
    # Plot 
    cum_hazard_comp <- ggsurvplot(km_fit_clean,
                                  data         = yearly_survival,
                                  fun          = "cumhaz",
                                  conf.int     = add_cis,
                                  censor.shape = "|",
                                  censor.size  = 4,
                                  title        = paste0("Cumulative Hazard by ", group_comparison_individual),
                                  subtitle     = subtitle_text_wrapped,
                                  xlab         = "Days at risk",
                                  ylab         = "Cumulative Hazard",
                                  legend       = "bottom",
                                  legend.title = group_comparison_individual,
                                  legend.labs  = levels(yearly_survival[[group_comparison_individual]]),
                                  palette      = my_palette,
                                  risk.table   = FALSE,
                                  cumevents    = FALSE,
                                  font.main    = c(14, "bold", "black"),
                                  font.x       = 12, font.y = 12, font.tickslab = 11,
                                  ggtheme = theme_classic(base_size = 12) +
                                    theme(plot.title         = element_text(face = "bold", size = 14),
                                          plot.subtitle      = element_text(size = 12, color = "gray50"),
                                          axis.text          = element_text(color = "black"),
                                          panel.grid.major.y = element_line(color = "gray90"),
                                          panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
                                          plot.margin        = margin(20, 10, 10, 10),
                                          legend.direction   = "horizontal",
                                          legend.box         = "horizontal",
                                          legend.spacing.x   = unit(0.3, "cm")))
    
    # Spread legend
    cum_hazard_comp$plot <- cum_hazard_comp$plot +
      guides(color = guide_legend(nrow           = ifelse(n_groups <= 8, 1, 2),
                                  byrow          = TRUE,
                                  title.position = "top"),
             fill  = guide_legend(nrow           = ifelse(n_groups <= 4, 1, 2),
                                  byrow          = TRUE,
                                  title.position = "top"))
    
    # Zoom y-axis to data range
    if (zoom_to_plot) {
      y_max     <- max(-log(km_fit_clean$surv), na.rm = TRUE)
      y_ceiling <- ceiling(y_max * 10) / 10
      cum_hazard_comp$plot <- cum_hazard_comp$plot +
        coord_cartesian(ylim = c(0, y_ceiling)) +
        annotate("text",
                 x        = -Inf,
                 y        = y_ceiling,
                 label    = paste0("* y-axis zoomed to [0, ", round(y_ceiling, 2), "]"),
                 hjust    = -0.05,
                 vjust    = 1.5,
                 size     = 3.2,
                 color    = "gray50",
                 fontface = "italic")
    }
    
    # Save plot  
    png(appArtifactPath("cum_hazard_comparison_plot.png"), width = 10, height = 8, units = "in", res = 300)
    print(cum_hazard_comp)
    dev.off()
  }
  
  
  # For summary table data 
  if(!is.null(group_comparison_individual) && is.null(survival_yr_start)) {
    logger.info(paste0("Comparing ", group_comparison_individual, " using summary table data"))
    
    # Fit survival object ---
    summary_table$time_at_risk <- summary_table$exit_time_days - summary_table$entry_time_days
    formula_str <- paste("Surv(time_at_risk, mortality_event) ~", group_comparison_individual)
    surv_formula <- as.formula(formula_str)
    km_fit_comp <- survfit(surv_formula, data = summary_table)
    
    
    # Log-Rank test --- 
    test <- survdiff(surv_formula, data=summary_table)
    
    # Extract components
    groups <- names(test$n)
    n_total <- sum(test$n)
    events_total <- sum(test$obs)
    chisq_val <- round(test$chisq, 2)
    df_val <- length(test$n) - 1
    p_val <- 1 - pchisq(test$chisq, df_val)
    p_formatted <- ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val))
    
    # Summary table 
    per_group <- tibble(!!group_comparison_individual := sub(".*=", "", groups),
                        `N`               = test$n,
                        `Events`          = test$obs,
                        `Expected events` = round(test$exp, 2),
                        `O/E ratio`       = round(test$obs / test$exp, 2)) %>%
      mutate(`N (events)` = sprintf("%d (%d)", N, Events), .keep = "unused")
    
    summary_row <- tibble(!!group_comparison_individual := "Overall",
                          `N (events)`             = sprintf("%d (%d)", n_total, events_total),
                          `Chisq (log-rank)`       = chisq_val,
                          `df`                     = df_val,
                          `p-value`                = p_formatted)
    
    logrank_table <- bind_rows(per_group, summary_row)
    logrank_table <- logrank_table %>%
      dplyr::select(any_of(group_comparison_individual), 
                    `N (events)`,
                    `Expected events`,
                    `O/E ratio`,
                    `Chisq (log-rank)`,
                    df,
                    `p-value`)
    
    # Save
    write.csv(logrank_table, file = appArtifactPath("logrank_table_statistics.csv"), row.names = F)
    
    
    # KM comparison plots ---
    km_fit_comp <- surv_fit(surv_formula, data = summary_table)   
    
    # Check if any groups have N == 1 and remove 
    filter_singleton_strata_refit <- function(fit, data) {
      
      grouping_var <- all.vars(formula(fit))[3]   
      n_per_group <- fit$n
      group_names <- sub(".*=", "", names(fit$strata) %||% "Overall")
      
      # Identify groups where N == 1 
      singletons <- n_per_group == 1
      removed <- group_names[singletons]    
      
      if (any(singletons)) {
        logger.info(paste("Removed the following singleton group(s) (N=1):\n",
                          paste(" •", removed, collapse = "\n ")), call. = FALSE)
        
        # Filter data using clean group names
        data_clean <- data[!data[[grouping_var]] %in% removed, , drop = FALSE]
        
        # Refit 
        km_clean <- surv_fit(formula(fit), data = data_clean)
        
        return(km_clean)
      }
      return(fit)
    }
    
    km_fit_clean <- filter_singleton_strata_refit(km_fit_comp, summary_table)
    
    # Titles and formatting 
    title_text <- paste0("Kaplan-Meier Survival Curve: ", group_comparison_individual)
    
    group_rows <- logrank_table %>% 
      filter(!!sym(group_comparison_individual) != "Overall") 
    
    # Remove any rows where N = 1 
    group_rows_clean <- group_rows %>%
      filter(!str_detect(`N (events)`, "^1\\s*\\("))
    
    removed_groups <- group_rows %>%
      filter(str_detect(`N (events)`, "^1\\s*\\(")) %>%
      pull(1)
    
    if (length(removed_groups) > 0) {
      logger.info(paste0("The following group(s) had N=1 and were removed from the table:\n • ",
                         paste(removed_groups, collapse = "\n • ")), call. = FALSE)
    }
    
    if (nrow(group_rows_clean) == 1){
      logger.fatal("There is only one group left; unable to perform comparisons.")
      
    } else { 
      subtitle_text <- paste0(paste(sprintf("N_%s: %s", 
                                            group_rows_clean[[group_comparison_individual]],
                                            group_rows_clean$`N (events)`),
                                    collapse = ", "),
                              "\nP-value: ", 
                              logrank_table$`p-value`[logrank_table[[group_comparison_individual]] == "Overall"])
      
      # Wrap each line independently to preserve the \n break before P-value
      subtitle_text_wrapped <- paste(
        sapply(strsplit(subtitle_text, "\n")[[1]], function(line) {
          paste(strwrap(line, width = 90), collapse = "\n")
        }),
        collapse = "\n")
      
      old_strata_names <- names(km_fit_clean$strata)
      new_strata_names <- sub(".*=", "", old_strata_names)
      names(km_fit_clean$strata) <- new_strata_names
      n_groups <- length(names(km_fit_clean$strata))
      
      # Plot 
      okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                     "#0072B2", "#D55E00", "#CC79A7", "#000000",
                     "#999999", "#44AA99")
      my_palette <- okabe_ito[seq_len(n_groups)]
      
      km_comp_curve <- ggsurvplot(km_fit_clean,
                                  data = summary_table,
                                  title = title_text,
                                  subtitle = subtitle_text_wrapped,
                                  conf.int = add_cis,
                                  risk.table = FALSE,
                                  palette = my_palette,
                                  xlab = "Days at risk",
                                  ylab = "Survival probability",
                                  legend.title = group_comparison_individual,
                                  legend = "bottom",
                                  legend.labs = levels(summary_table[[group_comparison_individual]]),
                                  censor.shape = "|",
                                  censor.size = 4,
                                  font.main = c(14, "bold", "black"),
                                  font.x = 12, font.y = 12, font.tickslab = 11,
                                  ggtheme = theme_classic(base_size = 12) +
                                    theme(
                                      plot.title         = element_text(face = "bold", size = 14),
                                      plot.subtitle      = element_text(size = 12, color = "gray50"),
                                      axis.text          = element_text(color = "black"),
                                      panel.grid.major.y = element_line(color = "gray90"),
                                      panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
                                      plot.margin        = margin(20, 10, 10, 10),
                                      legend.direction   = "horizontal",
                                      legend.box         = "horizontal",
                                      legend.spacing.x   = unit(0.3, "cm")))
      
      # Spread legend 
      km_comp_curve$plot <- km_comp_curve$plot +
        guides(color = guide_legend(nrow    = ifelse(n_groups <= 8, 1, 2),
                                    byrow   = TRUE,
                                    title.position = "top"),
               fill  = guide_legend(nrow    = ifelse(n_groups <= 4, 1, 2),
                                    byrow   = TRUE,
                                    title.position = "top"))
      
      # Zoom y-axis to data range
      if (zoom_to_plot) {
        y_min   <- min(km_fit_clean$surv, na.rm = TRUE)
        y_floor <- floor(y_min * 10) / 10
        km_comp_curve$plot <- km_comp_curve$plot +
          coord_cartesian(ylim = c(y_floor, 1)) +
          annotate("text",
                   x        = -Inf,
                   y        = y_floor,
                   label    = paste0("* y-axis zoomed to [", round(y_floor, 2), ", 1]"),
                   hjust    = -0.05,
                   vjust    = -0.5,
                   size     = 3.2,
                   color    = "gray50",
                   fontface = "italic")
      }
      
      # Save plot 
      png(appArtifactPath("km_comparison_curves.png"), width = 10, height = 8, units = "in", res = 300)
      print(km_comp_curve)
      dev.off()
    }
    
    
    ## Cumulative hazard comparison plots ---
    
    # Prepare statistics for subtitle 
    n_per_group      <- km_fit_clean$n
    sum_fit          <- summary(km_fit_clean)
    events_per_group <- tapply(sum_fit$n.event, sum_fit$strata, sum, na.rm = TRUE)
    
    # Clean strata names  
    clean_strata <- gsub("^.*=", "", names(km_fit_clean$strata))
    
    # Create one-line subtitle for groups: "N(Group): N (events)"
    subtitle_parts <- mapply(function(group, n, ev) {
      sprintf("N_%s: %d(%d)", group, n, ev)
    },
    clean_strata, n_per_group, events_per_group, SIMPLIFY = FALSE)
    
    groups_line <- paste(subtitle_parts, collapse = ", ")
    test        <- surv_pvalue(km_fit_clean, data = summary_table)
    pval_text   <- sprintf("P-value: %.3f", test$pval)
    subtitle_text <- paste0(groups_line, "\n", pval_text)
    
    # Wrap each line independently to preserve the \n break before P-value
    subtitle_text_wrapped <- paste(
      sapply(strsplit(subtitle_text, "\n")[[1]], function(line) {
        paste(strwrap(line, width = 90), collapse = "\n")
      }),
      collapse = "\n")
    
    # Plot 
    okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
                   "#D55E00", "#CC79A7", "#000000", "#999999", "#44AA99")
    my_palette <- okabe_ito[seq_len(n_groups)]
    
    cum_hazard_comp <- ggsurvplot(km_fit_clean,
                                  data         = summary_table,
                                  fun          = "cumhaz",
                                  conf.int     = add_cis,
                                  censor.shape = "|",
                                  censor.size  = 4,
                                  title        = paste0("Cumulative Hazard by ", group_comparison_individual),
                                  subtitle     = subtitle_text_wrapped,
                                  xlab         = "Days at risk",
                                  ylab         = "Cumulative Hazard",
                                  legend       = "bottom",
                                  legend.title = group_comparison_individual,
                                  legend.labs  = levels(summary_table[[group_comparison_individual]]),
                                  palette      = my_palette,
                                  risk.table   = FALSE,
                                  cumevents    = FALSE,
                                  font.main    = c(14, "bold", "black"),
                                  font.x       = 12, font.y = 12, font.tickslab = 11,
                                  ggtheme = theme_classic(base_size = 12) +
                                    theme(plot.title         = element_text(face = "bold", size = 14),
                                          plot.subtitle      = element_text(size = 12, color = "gray50"),
                                          axis.text          = element_text(color = "black"),
                                          panel.grid.major.y = element_line(color = "gray90"),
                                          panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
                                          plot.margin        = margin(20, 10, 10, 10),
                                          legend.direction   = "horizontal",
                                          legend.box         = "horizontal",
                                          legend.spacing.x   = unit(0.3, "cm")))
    
    # Spread legend 
    cum_hazard_comp$plot <- cum_hazard_comp$plot +
      guides(color = guide_legend(nrow           = ifelse(n_groups <= 8, 1, 2),
                                  byrow          = TRUE,
                                  title.position = "top"),
             fill  = guide_legend(nrow           = ifelse(n_groups <= 4, 1, 2),
                                  byrow          = TRUE,
                                  title.position = "top"))
    
    # Zoom y-axis to data range
    if (zoom_to_plot) {
      y_max     <- max(-log(km_fit_clean$surv), na.rm = TRUE)
      y_ceiling <- ceiling(y_max * 10) / 10
      cum_hazard_comp$plot <- cum_hazard_comp$plot +
        coord_cartesian(ylim = c(0, y_ceiling)) +
        annotate("text",
                 x        = -Inf,
                 y        = y_ceiling,
                 label    = paste0("* y-axis zoomed to [0, ", round(y_ceiling, 2), "]"),
                 hjust    = -0.05,
                 vjust    = 1.5,
                 size     = 3.2,
                 color    = "gray50",
                 fontface = "italic")
    }
    
    # Save plot 
    png(appArtifactPath("cum_hazard_comparison_plot.png"), 
        width = 10, height = 8, units = "in", res = 300)
    print(cum_hazard_comp)
    dev.off()
  }
  
  # Pass original to the next app in the MoveApps workflow
  return(data)
} 
