# -----------------------------------------
# Parse term -> (var, level)
# Find "term prefix" in candidate variable names, select match
# -----------------------------------------


parse_term_var_level <- function(term, var_candidates) {
  
  # ====== Critical fix: Ensure term is a length-1 character ======
  if (length(term) != 1) term <- term[1]
  term <- as.character(term)
  
  if (is.na(term) || !nzchar(term)) {
    return(list(var = NA_character_, level = NA_character_))
  }
  
  # Do not handle splines, interactions, complex terms
  if (grepl("^ns\\(|^splines::ns\\(|\\(|\\)|:", term)) {
    return(list(var = NA_character_, level = NA_character_))
  }
  
  # Defense: candidate variables may contain NA / non-character
  var_candidates <- as.character(var_candidates)
  var_candidates <- var_candidates[!is.na(var_candidates) & nzchar(var_candidates)]
  if (length(var_candidates) == 0) {
    return(list(var = NA_character_, level = NA_character_))
  }
  
  hits <- var_candidates[
    vapply(var_candidates, function(v) startsWith(term, v), logical(1))
  ]
  
  if (length(hits) == 0) {
    return(list(var = NA_character_, level = NA_character_))
  }
  
  # Select longest prefix
  var <- hits[which.max(nchar(hits))]
  level <- substr(term, nchar(var) + 1, nchar(term))
  
  if (!nzchar(level)) level <- NA_character_
  
  list(var = var, level = level)
}


# -----------------------------------------
# For a single row: Subset in imp_clean$data by set + main_level + term_level,
# Calculate mean/median of time_var (remove NA)
# If term parsing fails/subset is empty/variable has no values -> directly return NA (your requested behavior)
# -----------------------------------------
get_surv_stats_one_row <- function(imp_clean, set_name,
                                   main_bin, main_level,
                                   time_var,
                                   term, var_candidates) {
  
  # Fixed return template: NA if not found
  out_na <- c(n_surv = NA_real_, mean_surv = NA_real_, median_surv = NA_real_)
  
  pv <- parse_term_var_level(term, var_candidates)
  v  <- pv$var
  lv <- pv$level
  
  dat0 <- mice::complete(imp_clean, action = 1)
  
  # Basic filtering: set + strata(main_bin==main_level) + time not NA
  dsub <- dat0 %>%
    dplyr::filter(.data$analysis_set == set_name) %>%
    dplyr::filter(!is.na(.data[[time_var]])) %>%
    dplyr::filter(!is.na(.data[[main_bin]])) %>%
    dplyr::filter(as.character(.data[[main_bin]]) == as.character(main_level))
  
  # term parsing failed: directly NA (no statistics)
  if (is.na(v) || is.na(lv) || !(v %in% names(dsub))) {
    return(out_na)
  }
  
  # term filtering
  dsub <- dsub %>%
    dplyr::filter(!is.na(.data[[v]])) %>%
    dplyr::filter(as.character(.data[[v]]) == as.character(lv))
  
  tt <- dsub[[time_var]]
  tt <- tt[!is.na(tt)]
  
  # No data found: return NA (according to your request)
  if (length(tt) == 0) {
    return(out_na)
  }
  
  c(
    n_surv      = length(tt),
    mean_surv   = mean(tt),
    median_surv = stats::median(tt)
  )
}


run_one_analysis <- function(imp_clean, set_name,
                             main_bin,
                             time_var, status_var,
                             covars_bin,
                             spline_vars = character(0),
                             df_spline = 4) {
  
  dat_imp <- imp_clean %>%
    dplyr::filter(.data$analysis_set == set_name) %>%
    dplyr::filter(stats::complete.cases(dplyr::across(dplyr::all_of(c(time_var, status_var, main_bin)))))
  
  ##Core function for Multivariate Cox analsyis
  res <- run_stratified_mi_cox_nw(
    dat_imp,
    main_bin    = main_bin,
    time_var    = time_var,
    status_var  = status_var,
    covariates  = setdiff(covars_bin, main_bin),
    spline_vars = spline_vars,
    df_spline   = df_spline
  )
  
  
  # =========================================================
  # Minimal defense: If no results (NULL or 0 rows), directly return
  # =========================================================
  if (is.null(res) || !is.data.frame(res) || nrow(res) == 0) {
    return(res)
  }
  
  # Add labels, very convenient for subsequent error checking/summarization
  res$analysis_set <- set_name
  res$endpoint     <- ifelse(time_var == "PFS_Months", "PFS", "OS")
  res$main_bin_id  <- main_bin
  
  # =========================================================
  # New addition: Calculate average/median survival time for each row's (set_name + main_level + term) subset
  # Returns NA if no data found (ensured by get_surv_stats_one_row)
  # =========================================================
  var_candidates <- unique(c(setdiff(covars_bin, main_bin), main_bin))
  
  stats_list <- lapply(seq_len(nrow(res)), function(i) {
    get_surv_stats_one_row(
      imp_clean       = imp_clean,
      set_name        = set_name,
      main_bin        = main_bin,
      main_level      = res$main_level[i],
      time_var        = time_var,
      term            = res$term[i],
      var_candidates  = var_candidates
    )
  })
  
  stats_mat <- do.call(rbind, stats_list)
  stats_df  <- as.data.frame(stats_mat, stringsAsFactors = FALSE)
  
  stats_df$n_surv      <- as.numeric(stats_df$n_surv)
  stats_df$mean_surv   <- as.numeric(stats_df$mean_surv)
  stats_df$median_surv <- as.numeric(stats_df$median_surv)
  
  res <- dplyr::bind_cols(res, stats_df)
  res
}

run_module <- function(imp_clean, main_bin, covars_bin,
                       spline_vars_PFS = character(0),  ##Parameter reserved for spline-restricted variables
                       spline_vars_OS  = character(0),  ##Parameter reserved for spline-restricted variables
                       sets = c("Development","External_Test1","External_Test2")) {
  
  out <- list()
  
  for (set_name in sets) {
    out[[paste0(set_name, "_PFS")]] <- run_one_analysis(
      imp_clean, set_name,
      main_bin = main_bin,
      time_var = "PFS_Months", status_var = "PFS_Event",
      covars_bin = covars_bin,
      spline_vars = spline_vars_PFS
    )
    
    out[[paste0(set_name, "_OS")]] <- run_one_analysis(
      imp_clean, set_name,
      main_bin = main_bin,
      time_var = "OS_Months", status_var = "OS_Event",
      covars_bin = covars_bin,
      spline_vars = spline_vars_OS
    )
  }
  
  out
}

run_stratified_mi_cox_nw <- function(imp, 
                                     main_bin,      # e.g., "NLR_bin" / "Age_bin" / "PDL1_bin"
                                     time_var,      # "PFS_Months" / "OS_Months"
                                     status_var,    # "PFS_Event" / "OS_Event"
                                     covariates,    # Other *_bin variables entering Cox
                                     spline_vars = NULL, ##Parameter reserved for spline-restricted variables
                                     df_spline = 4) {
  
  
  # Levels of the main variable
  levs <- levels(imp$data[[main_bin]])
  if (is.null(levs)) levs <- unique(imp$data[[main_bin]])
  
  # Exclude NA
  levs <- levs[!is.na(levs)]
  
  m <- imp$m
  res_all <- list()
  n_patient <- 0
  
  for (lev in levs) {
    cat("Analyzing:", main_bin, "=", lev, "\n")
    
    # Construct Cox right-hand side: covariates + spline terms
    rhs <- covariates
    
    # Check and add spline terms
    if (!is.null(spline_vars) && length(spline_vars) > 0) {
      # Filter out non-existent variables
      existing_spline_vars <- spline_vars[spline_vars %in% names(imp$data)]
      if (length(existing_spline_vars) > 0) {
        # Check if they are numeric variables
        numeric_spline_vars <- character(0)
        for (var in existing_spline_vars) {
          if (is.numeric(imp$data[[var]])) {
            numeric_spline_vars <- c(numeric_spline_vars, var)
          } else {
            warning("Variable ", var, " is not numeric, skipping spline term")
          }
        }
        
        if (length(numeric_spline_vars) > 0) {
          rhs <- c(rhs, paste0("ns(", numeric_spline_vars, ", df = ", df_spline, ")"))
        }
      }
    }
    
    form_str <- paste0(
      "Surv(", time_var, ", ", status_var, ") ~ ",
      paste(rhs, collapse = " + ")
    )
    form_cox <- as.formula(form_str)
    
    # Fit a Cox model on each imputed dataset
    fits <- vector("list", m)
    
    for (k in seq_len(m)) {
      dat_k <- complete(imp, action = k)
      
      dat_sub <- dat_k %>%
        dplyr::filter(
          #analysis_set == "Development",  ##Commented out, already specified earlier
          .data[[main_bin]] == lev,
          !is.na(.data[[time_var]]),
          !is.na(.data[[status_var]])
        )
      
      n_patient<-nrow(dat_sub)
      # If too few samples for this level, can skip
      n_events <- sum(dat_sub[[status_var]] == 1, na.rm = TRUE)
      if (nrow(dat_sub) < 30 || n_events < 5) {
        cat("  Imputation set", k, ": Insufficient samples (n =", nrow(dat_sub), ", events =", n_events, ")\n")
        fits[[k]] <- NULL
      } else {
        tryCatch({
          fits[[k]] <- coxph(form_cox, data = dat_sub)
          cat("  Imputation set", k, ": Success (n =", nrow(dat_sub), ", events =", n_events, ")\n")
        }, error = function(e) {
          cat("  Imputation set", k, ": Failed -", e$message, "\n")
          fits[[k]] <- NULL
        })
      }
    }
    
    # Remove NULL models
    fits <- Filter(Negate(is.null), fits)
    
    if (length(fits) < 1) {
      cat("Level", lev, "has insufficient valid models, skipping\n")
      next
    }
    
    # Construct "mira" object independent of with()
    mira_obj <- list(call = match.call(), analyses = fits)
    class(mira_obj) <- "mira"
    
    pool_fit <- mice::pool(mira_obj)
    tab <- summary(pool_fit, conf.int = TRUE, exponentiate = TRUE)
    
    # Fix: If summary.coxph (non-pooled), skip this strata
    if (inherits(tab, "summary.coxph")) {
      cat("Pooled summary returned summary.coxph (possibly due to only 1 model), skipping\n")
      next
    }
    
    tab_df <- as.data.frame(tab)
    
    # If column names are 2.5 % / 97.5 %, change to conf.low / conf.high for convenience in subsequent use
    if ("2.5 %" %in% names(tab_df)) {
      tab_df$conf.low  <- tab_df$`2.5 %`
      tab_df$conf.high <- tab_df$`97.5 %`
    }
    
    # Add metadata
    tab_df$main_bin   <- main_bin
    tab_df$main_level <- lev
    tab_df$time_var   <- time_var
    tab_df$term       <- tab_df$term  # Ensure term column is retained
    tab_df$n_patients <- n_patient
    
    rownames(tab_df) <- NULL
    
    # Select and reorder columns (here estimate is HR because exponentiate = TRUE)
    final_cols <- c("main_bin", "main_level", "time_var", "term",
                    "estimate", "std.error", "statistic", "p.value","n_patients")
    if ("conf.low" %in% names(tab_df)) {
      final_cols <- c(final_cols, "conf.low", "conf.high")
    }
    
    tab_df <- tab_df[, intersect(final_cols, names(tab_df))]
    
    res_all[[as.character(lev)]] <- tab_df
  }
  
  if (length(res_all) == 0) {
    warning("No results generated")
    return(NULL)
  }
  
  # Combine all results
  dplyr::bind_rows(res_all)
}