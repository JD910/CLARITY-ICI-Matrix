library(dplyr)
library(purrr)
library(tidyr)
library(stringr)


# Flatten results: Combine all data.frames in results$*bin$* into one long table
flatten_results2 <- function(results) {
  imap_dfr(results, function(one_bin, bin_name) {
    # one_bin is a list: Development_PFS / Development_OS / Ext1_PFS ...
    dfs <- keep(one_bin, is.data.frame)
    bind_rows(dfs) %>%
      mutate(main_bin_id = bin_name)  # Fallback marker: to guard against main_bin column anomalies
  })
}

get_all3_sig <- function(df, endpoint_vec = c("PFS","OS"), p_cut = 0.05) {
  df %>%
    filter(endpoint %in% endpoint_vec, !is.na(p.value), p.value < p_cut) %>%
    group_by(main_bin, main_level, term, endpoint) %>%
    filter(n_distinct(analysis_set) == 3) %>%   # Significant in all three datasets
    ungroup()
}


get_all2_sig <- function(df, endpoint_vec = c("PFS","OS"), p_cut = 0.05) {
  df %>%
    filter(endpoint %in% endpoint_vec, !is.na(p.value), p.value < p_cut) %>%
    group_by(main_bin, main_level, term, endpoint) %>%
    filter(n_distinct(analysis_set) == 2) %>%   # Significant in both datasets
    ungroup()
}


to_wide3 <- function(df, imp_clean, analysis_set) {
  
  endpoint_to_timevar <- function(endpoint) {
    if (endpoint == "PFS") "PFS_Months" else "OS_Months"
  }
  
  # term -> (var, level)
  # Depends on the existing parse_term_var_level(term, var_candidates)
  make_term <- function(var, level) {
    paste0(var, level)
  }
  
  # -------------------------------
  # mean / median survival calculation
  # -------------------------------
  get_all_stats <- function(main_bin, main_level, term, endpoint, var_candidates, analysis_set_value) {
    
    dat0 <- imp_clean$data %>%
      dplyr::filter(analysis_set == analysis_set_value)  # First filter to the specified dataset
    time_var <- endpoint_to_timevar(endpoint)
    
    pv <- parse_term_var_level(term, var_candidates)
    v  <- pv$var
    lv <- pv$level
    
    # Base filtering: time + main_bin
    dsub0 <- dat0 %>%
      dplyr::filter(!is.na(.data[[time_var]])) %>%
      dplyr::filter(!is.na(.data[[main_bin]])) %>%
      dplyr::filter(as.character(.data[[main_bin]]) == as.character(main_level)) 
    
    # Term constraint (if parsable)
    if (!is.na(v) && !is.na(lv) && v %in% names(dsub0)) {
      dsub <- dsub0 %>%
        dplyr::filter(!is.na(.data[[v]])) %>%
        dplyr::filter(as.character(.data[[v]]) == as.character(lv))
    } else {
      return(c(n = 0, mean = NA_real_, median = NA_real_))
    }
    
    tt <- dsub[[time_var]]
    tt <- tt[!is.na(tt)]
    
    if (length(tt) == 0) {
      return(c(n = 0, mean = NA_real_, median = NA_real_))
    }
    
    c(
      n      = length(tt),
      mean   = mean(tt),
      median = stats::median(tt)
    )
  }
  
  
  # 1) First create table for HR / CI / p
  
  base_wide <- df %>%
    dplyr::select(main_bin, main_level, term, endpoint, analysis_set,
                  estimate, conf.low, conf.high, p.value) %>%
    dplyr::mutate(
      HR_CI = sprintf("%.3f (%.3f–%.3f)", estimate, conf.low, conf.high),
      p_fmt = ifelse(is.na(p.value), NA, formatC(p.value, format="f", digits=4))
    ) %>%
    dplyr::select(main_bin, main_level, term, endpoint, analysis_set,
                  HR_CI, p_fmt) %>%
    tidyr::pivot_wider(
      names_from  = analysis_set,
      values_from = c(HR_CI, p_fmt),
      names_glue  = "{analysis_set}_{.value}"
    )
  
  # 2) Generate all categories for "term belonging variable" and add missing rows
  
  key_df <- df %>%
    dplyr::distinct(main_bin, main_level, term, endpoint)
  
  # Candidate variable names: for parse_term_var_level
  var_candidates <- unique(c(
    key_df$main_bin,
    "Age_bin","Patho_bin","Smoking_bin","NLR_bin","ECOG_bin","FGA_bin",
    "Diameter_bin","Tumor_Purity_bin","Stage_bin","P_Treat_bin",
    "IO_Type_bin","TMB_bin","PD_L1_bin","Sex_bin","Region_bin",
    "Line_bin","Real_World_bin"
  ))
  var_candidates <- as.character(var_candidates[!is.na(var_candidates)])
  
  dat0 <- imp_clean$data %>%
    dplyr::filter(analysis_set == !!analysis_set)  # Use analysis_set passed as parameter
  
  # For each key: parse the term's belonging variable v, and take all levels of v (excluding NA), generate complete term list
  expanded_keys <- do.call(
    rbind,
    lapply(seq_len(nrow(key_df)), function(i) {
      mb <- key_df$main_bin[i]
      ml <- key_df$main_level[i]
      tm <- key_df$term[i]
      ep <- key_df$endpoint[i]
      
      pv <- parse_term_var_level(tm, var_candidates)
      v  <- pv$var
      
      if (is.na(v) || !(v %in% names(dat0))) {
        # Term parsing failed: only keep the original key
        return(key_df[i, , drop = FALSE])
      }
      
      lv_all <- as.character(dat0[[v]])
      lv_all <- lv_all[!is.na(lv_all) & lv_all != ""]
      lv_all <- sort(unique(lv_all))
      
      if (length(lv_all) == 0) {
        return(key_df[i, , drop = FALSE])
      }
      
      out <- data.frame(
        main_bin   = mb,
        main_level = ml,
        term       = make_term(v, lv_all),
        endpoint   = ep,
        stringsAsFactors = FALSE
      )
      out
    })
  )
  expanded_keys <- dplyr::as_tibble(expanded_keys) %>%
    dplyr::distinct(main_bin, main_level, term, endpoint)
  
  new_keys <- expanded_keys %>%
    dplyr::anti_join(
      base_wide %>% dplyr::select(main_bin, main_level, term, endpoint),
      by = c("main_bin","main_level","term","endpoint")
    )
  
  # Fill all columns (HR/CI/p etc.) with NA for the newly added keys
  if (nrow(new_keys) > 0) {
    template_cols <- setdiff(names(base_wide), c("main_bin","main_level","term","endpoint"))
    new_rows <- new_keys
    for (cc in template_cols) new_rows[[cc]] <- NA_character_
    
    # Bind: original table + newly added blank rows
    base_wide2 <- dplyr::bind_rows(base_wide, new_rows %>% dplyr::select(names(base_wide)))
  } else {
    base_wide2 <- base_wide
  }
  
  
  # Calculate ALL_n/mean/median uniformly for "all expanded keys"
  
  all_keys_for_stats <- base_wide2 %>%
    dplyr::select(main_bin, main_level, term, endpoint) %>%
    dplyr::distinct()
  
  stats_mat <- t(mapply(
    FUN = function(main_bin, main_level, term, endpoint) {
      get_all_stats(
        main_bin       = main_bin,
        main_level     = main_level,
        term           = term,
        endpoint       = endpoint,
        var_candidates = var_candidates,
        analysis_set_value = analysis_set  # Pass parameter
      )
    },
    main_bin   = all_keys_for_stats$main_bin,
    main_level = all_keys_for_stats$main_level,
    term       = all_keys_for_stats$term,
    endpoint   = all_keys_for_stats$endpoint,
    SIMPLIFY   = TRUE
  ))
  
  stats_df <- as.data.frame(stats_mat, stringsAsFactors = FALSE)
  for (nm in names(stats_df)) stats_df[[nm]] <- as.numeric(stats_df[[nm]])
  
  stats_df <- dplyr::bind_cols(
    all_keys_for_stats,
    stats_df %>%
      dplyr::rename(
        ALL_n      = n,
        ALL_mean   = mean,
        ALL_median = median
      )
  )
  
  
  # Merge into final output
  
  out <- base_wide2 %>%
    dplyr::left_join(
      stats_df,
      by = c("main_bin","main_level","term","endpoint")
    ) %>%
    dplyr::arrange(main_bin, endpoint, main_level, term)
  
  out
}