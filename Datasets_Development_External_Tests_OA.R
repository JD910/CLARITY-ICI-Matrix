# Function to pick IDs to meet target sample size
pick_by_target_n <- function(id_vec, target_n, sizes_df, prefer = c("large","small")) {
  prefer <- match.arg(prefer)
  ids <- intersect(unique(id_vec), sizes_df$PubMedID)
  if (length(ids) == 0 || target_n <= 0) return(character(0))
  
  df <- sizes_df %>% dplyr::filter(PubMedID %in% ids)
  # Sort by sample size: large first or small first
  df <- if (prefer == "large") df %>% dplyr::arrange(dplyr::desc(n_pat), PubMedID) else df %>% dplyr::arrange(n_pat, PubMedID)
  
  picked <- character(0); acc <- 0
  for (i in seq_len(nrow(df))) {
    if (acc >= target_n) break
    pid <- df$PubMedID[i]
    n_i <- df$n_pat[i]
    picked <- c(picked, pid)
    acc <- acc + n_i
  }
  unique(picked)
}

# Force external test 2 to contain only East Asian cohorts
force_ext2_eas_only <- function(external2_ids, east_ids) {
  external2_ids <- unique(external2_ids)
  unique(intersect(external2_ids, unique(east_ids)))
}

# Get sample counts for each dataset partition
get_set_counts <- function(Sta1, ext1_ids, ext2_ids) {
  all_ids <- unique(Sta1$PubMedID)
  dev_ids <- setdiff(all_ids, c(ext1_ids, ext2_ids))
  
  c(
    Development    = sum(Sta1$PubMedID %in% dev_ids),
    External_Test1 = sum(Sta1$PubMedID %in% ext1_ids),
    External_Test2 = sum(Sta1$PubMedID %in% ext2_ids)
  )
}

# Calculate partition score based on multiple criteria
calc_partition_score <- function(
    Sta1,
    external1_ids,
    external2_ids,
    balance_vars,
    min_n = 100,
    ratio_target = c(Development=0.35, External_Test1=0.33, External_Test2=0.32),
    w_zero = 500,        # Strong penalty for "level count = 0 in a set" (encourages all three sets to have all levels)
    w_min  = 80,         # min_n penalty
    w_balance = 1,       # Subgroup proportion difference penalty
    w_size = 30,         # 2:1:1 size penalty
    verbose = FALSE
){
  sets <- c("Development", "External_Test1", "External_Test2")
  
  # =========================================================
  # Small modification: Hard constraint - external2 must only contain EastAsia (EAS)
  # Depends on externally defined east_ids (PubMedID vector)
  # =========================================================
  if (!exists("east_ids", inherits = TRUE)) {
    stop("east_ids (EastAsia PubMedID vector) must be defined in the global environment first", call. = FALSE)
  }
  external2_ids <- intersect(unique(external2_ids), unique(get("east_ids", inherits = TRUE)))
  
  # Defensive deduplication: ext1 / ext2 should not overlap
  external1_ids <- setdiff(unique(external1_ids), external2_ids)
  
  all_ids <- unique(Sta1$PubMedID)
  dev_ids <- setdiff(all_ids, c(external1_ids, external2_ids))
  
  # Assign analysis_set labels based on partition
  dat <- Sta1 %>%
    dplyr::mutate(
      analysis_set_tmp = dplyr::case_when(
        PubMedID %in% dev_ids         ~ "Development",
        PubMedID %in% external1_ids   ~ "External_Test1",
        PubMedID %in% external2_ids   ~ "External_Test2",
        TRUE                           ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(analysis_set_tmp))
  
  # ---------- size penalty ----------
  n_set <- dat %>%
    dplyr::count(analysis_set_tmp, name = "n") %>%
    tidyr::complete(analysis_set_tmp = sets, fill = list(n = 0)) %>%
    dplyr::arrange(match(analysis_set_tmp, sets))
  
  N <- sum(n_set$n)
  target_n <- round(ratio_target * N)
  size_pen <- sum(((n_set$n - target_n)^2) / (max(1, N)^2))
  
  # ---------- balance/min_n/zero penalty ----------
  score_balance <- 0
  score_min <- 0
  score_zero <- 0
  violations <- list()
  
  for (v in balance_vars) {
    
    if (!v %in% colnames(dat)) next
    
    # Remove NA (do not treat NA as a level)
    dv <- dat %>%
      dplyr::select(analysis_set_tmp, dplyr::all_of(v)) %>%
      dplyr::mutate(value = as.character(.data[[v]])) %>%
      dplyr::filter(!is.na(value) & value != "")
    
    # Skip if too few valid samples for this variable (avoid noise)
    if (nrow(dv) < 10) next
    
    # Counts/proportions per set (only based on non-NA)
    tab <- dv %>%
      dplyr::count(analysis_set_tmp, value, name = "n") %>%
      dplyr::group_by(analysis_set_tmp) %>%
      dplyr::mutate(p = n / sum(n)) %>%
      dplyr::ungroup()
    
    # Ensure all three sets and all levels are expanded
    all_levels <- sort(unique(dv$value))
    tab <- tidyr::complete(
      tab,
      analysis_set_tmp = sets,
      value = all_levels,
      fill = list(n = 0, p = 0)
    )
    
    # --- 1) zero penalty: level count = 0 in a set (encourages all sets to have all levels)
    tab_zero <- tab %>% dplyr::filter(n == 0)
    if (nrow(tab_zero) > 0) {
      score_zero <- score_zero + nrow(tab_zero)
      violations[[paste0(v, "_zero")]] <- tab_zero
    }
    
    # --- 2) min_n penalty: level count in a set < min_n
    tab_min <- tab %>% dplyr::filter(n > 0 & n < min_n)
    if (nrow(tab_min) > 0) {
      score_min <- score_min + sum(min_n - tab_min$n)
      violations[[paste0(v, "_min")]] <- tab_min
    }
    
    # --- 3) balance penalty: proportion difference across three sets for same level (smaller range_p is better)
    tab_bal <- tab %>%
      dplyr::group_by(value) %>%
      dplyr::summarise(range_p = max(p) - min(p), .groups = "drop")
    
    score_balance <- score_balance + sum(tab_bal$range_p^2)
  }
  
  total_score <- w_zero * score_zero + w_min * score_min + w_balance * score_balance + w_size * size_pen
  
  if (isTRUE(verbose)) {
    message(
      "score_zero=", score_zero,
      " score_min=", score_min,
      " score_balance=", signif(score_balance, 4),
      " size_pen=", signif(size_pen, 4),
      " total=", signif(total_score, 6)
    )
  }
  
  list(
    total_score = total_score,
    score_zero = score_zero,
    score_min = score_min,
    score_balance = score_balance,
    size_pen = size_pen,
    n_set = setNames(n_set$n, n_set$analysis_set_tmp),
    target_n = target_n,
    violations = violations
  )
}

# Main optimization function using swap operations
optimize_partition_by_swaps <- function(
    Sta1,
    external1_ids,
    external2_ids,
    balance_vars,
    ext1_pool_ids,
    ext2_pool_ids,
    min_n = 100,
    ratio_target = c(Development=0.5, External_Test1=0.25, External_Test2=0.25),
    max_iter = 600,
    n_try_each_iter = 120,
    seed = 2026,
    w_min = 80,
    w_balance = 1,
    w_size = 30,
    verbose = TRUE
){
  set.seed(seed)
  
  all_ids <- unique(Sta1$PubMedID)
  
  # Initial cleaning and deduplication
  external1_ids <- intersect(unique(external1_ids), unique(ext1_pool_ids))
  external2_ids <- intersect(unique(external2_ids), unique(ext2_pool_ids))
  external1_ids <- setdiff(external1_ids, external2_ids)
  external2_ids <- setdiff(external2_ids, external1_ids)
  
  best_ext1 <- external1_ids
  best_ext2 <- external2_ids
  
  # Calculate initial score
  best_detail <- calc_partition_score(
    Sta1, best_ext1, best_ext2,
    balance_vars = balance_vars,
    min_n = min_n,
    ratio_target = ratio_target,
    w_min = w_min,
    w_balance = w_balance,
    w_size = w_size
  )
  best_score <- best_detail$total_score
  
  if (verbose) {
    message("Initial score: ", best_score)
    message("Initial counts: ", paste(names(best_detail$n_set), best_detail$n_set, collapse = ", "))
    message("Target counts : ", paste(names(best_detail$target_n), best_detail$target_n, collapse = ", "))
  }
  
  # Helper: Determine which set to swap based on current size deviation
  choose_swap_type <- function(n_set, target_n) {
    # Which external set is more "deficient" -> prioritize adding cohorts (from dev)
    # Which external set has more -> prioritize swapping out (back to dev)
    d1 <- n_set["External_Test1"] - target_n["External_Test1"]
    d2 <- n_set["External_Test2"] - target_n["External_Test2"]
    
    # If ext2 is more deficient, higher probability to increase ext2; vice versa
    p_ext2 <- plogis(-3 * d2 / max(1, target_n["External_Test2"])) # deficiency -> p increases
    p_ext1 <- plogis(-3 * d1 / max(1, target_n["External_Test1"]))
    
    # Normalize
    p <- c(ext1 = p_ext1, ext2 = p_ext2)
    p <- p / sum(p)
    sample(c("ext1","ext2"), 1, prob = p)
  }
  
  # Main optimization loop
  for (it in seq_len(max_iter)) {
    improved <- FALSE
    dev_ids <- setdiff(all_ids, c(best_ext1, best_ext2))
    
    # Make swap direction closer to size target
    swap_type <- choose_swap_type(best_detail$n_set, best_detail$target_n)
    
    for (k in seq_len(n_try_each_iter)) {
      
      if (swap_type == "ext1") {
        if (length(best_ext1) == 0) next
        out_id <- sample(best_ext1, 1)
        
        in_cand <- intersect(dev_ids, ext1_pool_ids)
        if (length(in_cand) == 0) next
        in_id <- sample(in_cand, 1)
        
        new_ext1 <- union(setdiff(best_ext1, out_id), in_id)
        new_ext2 <- best_ext2
        
      } else {
        if (length(best_ext2) == 0) next
        out_id <- sample(best_ext2, 1)
        
        in_cand <- intersect(dev_ids, ext2_pool_ids)
        if (length(in_cand) == 0) next
        in_id <- sample(in_cand, 1)
        
        new_ext1 <- best_ext1
        new_ext2 <- union(setdiff(best_ext2, out_id), in_id)
      }
      
      # Pool/deduplication defense
      new_ext1 <- intersect(unique(new_ext1), ext1_pool_ids)
      new_ext2 <- intersect(unique(new_ext2), ext2_pool_ids)
      new_ext1 <- setdiff(new_ext1, new_ext2)
      new_ext2 <- setdiff(new_ext2, new_ext1)
      
      # Calculate score for new partition
      sc_detail <- calc_partition_score(
        Sta1, new_ext1, new_ext2,
        balance_vars = balance_vars,
        min_n = min_n,
        ratio_target = ratio_target,
        w_min = w_min,
        w_balance = w_balance,
        w_size = w_size
      )
      sc <- sc_detail$total_score
      
      # Update if improved
      if (sc < best_score) {
        best_score <- sc
        best_ext1 <- new_ext1
        best_ext2 <- new_ext2
        best_detail <- sc_detail
        improved <- TRUE
        if (verbose) {
          message("Iter ", it, " improved: ", signif(best_score, 6), " (swap ", swap_type, ")")
          message("Counts: ", paste(names(best_detail$n_set), best_detail$n_set, collapse = ", "))
        }
        break
      }
    }
    
    if (!improved) {
      if (verbose) message("No improvement at iter ", it, ", stop.")
      break
    }
  }
  
  list(
    external1_ids = best_ext1,
    external2_ids = best_ext2,
    best_score = best_score,
    score_detail = best_detail
  )
}