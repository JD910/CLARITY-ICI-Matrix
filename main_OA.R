rm(list = ls())
gc()

library("readxl")
library("dplyr")
library("tidyverse")
library("mice")
library("survival")
library("openxlsx")
library("writexl")

setwd("Your personal path")
workbook1<-"Your personal path\\NSCLC_CLARITY_ICI_Matrix_OA.xlsx"
Sta1 = read_excel(workbook1,1,guess_max = 5000)

getwd()

## ===============================
## 1. Preprocess variables for imputation and CLARITY-ICI Matrix modeling
## ===============================

Sta1 %>%
  count(region) %>%
  mutate(prop = n / sum(n))

cohort_meta <- Sta1 %>%
  group_by(PubMedID, region) %>%
  summarise(n_patients = n(), .groups = "drop")


Sta1 <- Sta1 %>%
  mutate(
    # 1=IO_mono, 2=IO_chemo, 3=IO_antiangiogenic,
    # 4=IO_double, 5=IO_other_combo, 6=Unknown
    IO_Type = case_when(
      IOType == 1 ~ "Mono",
      IOType == 2 ~ "Chemo",
      #Drug == 3 ~ "IO_antiangiogenic",
      #Drug == 4 ~ "IO_double",
      #Drug == 5 ~ "IO_other_combo",
      IOType == 3 | IOType == 4 | IOType == 5 ~ "Other",  # IO+Targeted/Combination Therapy
      #Drug == 6 ~ "Unknown",
      TRUE       ~ NA_character_
    ),
    
    # 1 = Retrospective real-world, 2 = Prospective clinical trial
    Real_World = case_when(
      DataSource == 1 ~ "Yes",
      DataSource == 2 ~ "No",
      TRUE              ~ NA_character_
    ),
    
    ## Region: EastAsia / Western (Europe+NorthAmerica) / Other
    Region_EW = case_when(
      region == "EastAsia"                     ~ "EastAsia",
      region %in% c("Europe", "NorthAmerica")  ~ "Western",
      region == "Other"                        ~ "Other",
      TRUE                                     ~ "Unknown"
    ),
    
  )

###Data transformation

vars_for_mice <- c(
  "OS_Months", "OS_Event",
  "PFS_Months", "PFS_Event",
  "Age", "NLR", "PDL1", "TMB", "Diameter", "TumorPurity",
  "Sex", "Smoking", "ECOG", "Stage", "Histo", "FGA",  "P_Treat","IO_Type", "region", "Region_EW", "Line",
  "Real_World",
  "PFS_Months", "PFS_Event", "OS_Months", "OS_Event","PubMedID"
)

dat_mice_prepare <- Sta1 %>%
  select(all_of(vars_for_mice))


## Set variable types (factor / ordered)

dat_mice_prepare <- dat_mice_prepare %>%
  mutate(
    ## Binary/Multiclass to factor
    Sex          = factor(Sex),          
    Smoking      = factor(Smoking),      
    Histo    = factor(Histo),   
    IO_Type = factor(IO_Type), 
    Region_EW = factor(Region_EW), 
    region = factor(region), 
    Real_World = factor(Real_World),
    P_Treat = factor(P_Treat),
    Line_B = factor(Line),
    
    
    ## Ordinal variables → ordered
    ECOG  = ordered(ECOG),               
    Stage = ordered(Stage),             
    
  )

## Categorical variables
dat_mice_bin <- dat_mice_prepare %>%
  mutate(
    # Age grouping
    Age_bin = factor(ifelse(Age > 65, ">65", "≤65"), 
                     levels = c(">65", "≤65")),
    
    # Pathology type (copy)
    Patho_bin = factor(Histo,
                       levels = c("LUAD", "LUSC", "Other"),
                       labels = c("LUAD", "LUSC", "Other")),
    
    
    # Smoking status grouping
    Smoking_bin = factor(ifelse(Smoking == 0, "No", "Yes"),
                         levels = c("No", "Yes")),
    
    # ECOG grouping
    ECOG_bin = factor(ifelse(as.numeric(as.character(ECOG)) >= 1, "≥1", "<1"),
                      levels = c("≥1", "<1")),
    
    # FGA grouping
    FGA_bin = factor(ifelse(FGA < 0.2, "<0.2", "≥0.2"),
                     levels = c("<0.2", "≥0.2")),
    
    # Tumor diameter grouping
    Diameter_bin = factor(ifelse(Diameter > 50, ">50", "≤50"),
                          levels = c(">50", "≤50")),
    
    # Tumor purity grouping
    Tumor_Purity_bin = factor(ifelse(TumorPurity < 50, "<50", "≥50"),
                              levels = c("<50", "≥50")),
    
    # Stage grouping
    Stage_bin = factor(ifelse(Stage == 4, "IV", "I~III"),
                       levels = c("IV", "I~III")),
    
    # Prior treatment grouping
    P_Treat_bin = factor(ifelse(P_Treat == 1, "Yes", "No"),
                         levels = c("Yes", "No")),
    
    # Treatment type grouping
    IO_Type_bin = IO_Type,
    
    # Ethnic group/Race grouping
    Region_bin = factor(region,
                        levels = c("EastAsia", "Europe", "NorthAmerica", "Other"),
                        labels = c("EAS", "EUR", "NAM", "Other")),
    
    # NLR median grouping
    NLR_bin = factor(ifelse(as.numeric(as.character(NLR)) > 
                              median(as.numeric(as.character(NLR)), na.rm = TRUE), 
                            "High", "Low"),
                     levels = c("High", "Low")),
    
    # TMB grouping (High if >10)
    TMB_bin = factor(ifelse(as.numeric(as.character(TMB)) > 10, "High", "Low"),
                     levels = c("High", "Low")),
    
    ####For TMB three class (in a sensitivity analysis)
    
    #TMB_bin = factor(
    #  case_when(
    #    as.numeric(as.character(TMB)) < 5 ~ "Low",
    #    as.numeric(as.character(TMB)) > 15 ~ "High",
    #    TRUE ~ "Intermediate"  # Between 5-15
    #  ),
    #  levels = c("Low", "Intermediate", "High")
    #)
    
    # PD-L1 grouping (High if ≥50%)
    PD_L1_bin = factor(ifelse(as.numeric(as.character(PDL1)) >= 0.5, "High", "Low"),
                       levels = c("High", "Low")),
    
    # Sex grouping
    Sex_bin = factor(ifelse(Sex == 1, "Male", "Female"),
                     levels = c("Male", "Female")),
    
    # Line of therapy grouping
    Line_bin = factor(ifelse(Line == "First", "1L", "≥2L"),
                      levels = c("1L", "≥2L")),
    
    Real_World_bin = Real_World,
  )



set.seed(2026)

balance_vars <- c(
  "Age_bin", "Patho_bin", "Smoking_bin", "ECOG_bin", "FGA_bin", "Diameter_bin","Real_World_bin",
  "Tumor_Purity_bin", "Stage_bin", "P_Treat_bin", "IO_Type_bin", "Region_bin", "NLR_bin","TMB_bin","PD_L1_bin","Sex_bin","Line_bin"
)


## ===============================
## 2. Dataset partitioning
## ensuring the training set integrates data from all regions, 
## external test 1 contains only Western population data, 
## and external test 2 contains only East Asian population data.
## Ensuring the dataset proportions follow the designed `ratio_target`.
## ===============================

source("Datasets_Development_External_Tests_OA.R", encoding="utf-8")

min_n <- 20

# Set target proportions (Development, External1, External2)
ratio_target <- c(Development = 0.35, External_Test1 = 0.33, External_Test2 = 0.32)

N_total <- nrow(dat_mice_bin)
target_n <- round(ratio_target * N_total)
target_ext1 <- target_n["External_Test1"]
target_ext2 <- target_n["External_Test2"]
target_dev  <- target_n["Development"]


# Extract region cohort pools

east_ids    <- cohort_meta %>% filter(region == "EastAsia")      %>% pull(PubMedID)
eur_ids     <- cohort_meta %>% filter(region == "Europe")        %>% pull(PubMedID)
northam_ids <- cohort_meta %>% filter(region == "NorthAmerica")  %>% pull(PubMedID)
other_ids   <- cohort_meta %>% filter(region == "Other")         %>% pull(PubMedID)

ext1_pool <- unique(c(eur_ids, northam_ids, other_ids))        # external1 pool
ext2_pool <- unique(c(east_ids))         # external2 pool

cohort_sizes <- dat_mice_bin %>% dplyr::count(PubMedID, name = "n_pat")

external1_ids_init <- pick_by_target_n(ext1_pool, target_ext1, cohort_sizes, prefer = "large")

east_total_n <- dat_mice_bin %>% filter(PubMedID %in% east_ids) %>% nrow()
target_east_in_ext2 <- floor(east_total_n / 1.5)

external2_east_ids_init <- pick_by_target_n(east_ids, target_east_in_ext2, cohort_sizes, prefer = "large")
external2_ids_init <- external2_east_ids_init

# Remove duplicates
external1_ids_init <- setdiff(external1_ids_init, external2_ids_init)
external2_ids_init <- setdiff(external2_ids_init, external1_ids_init)

# Optimization (including 1:1:1 size penalty)

opt_res <- optimize_partition_by_swaps(
  Sta1 = Sta1,
  external1_ids = external1_ids_init,
  external2_ids = external2_ids_init,
  balance_vars = balance_vars,
  ext1_pool_ids = ext1_pool,
  ext2_pool_ids = ext2_pool,
  min_n = min_n,
  ratio_target = ratio_target,
  max_iter = 800,
  n_try_each_iter = 150,
  seed = 2026,
  w_min = 80,
  w_balance = 1,
  w_size = 30,
  verbose = TRUE
)

external1_ids_final <- opt_res$external1_ids
external2_ids_final <- opt_res$external2_ids
#
# Apply partition labels (analysis_set)
#
all_study_ids <- unique(dat_mice_bin$PubMedID)
development_ids <- setdiff(all_study_ids, c(external1_ids_final, external2_ids_final))

dat_mice_bin <- dat_mice_bin %>%
  mutate(
    analysis_set = case_when(
      PubMedID %in% development_ids      ~ "Development",
      PubMedID %in% external1_ids_final  ~ "External_Test1",
      PubMedID %in% external2_ids_final  ~ "External_Test2",
      TRUE ~ "Unassigned"
    )
  )
table(dat_mice_bin$analysis_set)

# Uncomment if metadata needs saving
complete_partition <- list(
  external1_queues = external1_ids_final,
  external2_queues = external2_ids_final,
  total_patients = c(
    external1 = dat_mice_bin %>% filter(analysis_set == "External_Test1") %>% nrow(),
    external2 = dat_mice_bin %>% filter(analysis_set == "External_Test2") %>% nrow(),
    development = dat_mice_bin %>% filter(analysis_set == "Development") %>% nrow()
  ),
  targets_211 = opt_res$score_detail$target_n,
  best_score = opt_res$best_score,
  note = paste0(
    "Optimized by swaps under pool constraints with 2:1:1 size penalty; min subgroup target=",
    min_n, ". External2 prioritizes EastAsia balance with Development."
  )
)

# save(
#   dat_mice_bin,
#   external1_ids_final,
#   external2_ids_final,
#   development_ids,
#   complete_partition,
#   file = "partition_results_balanceAll_1_1_1_OA_3000_Examples.RData"
# )


# Load RData file
load("partition_results_balanceAll_1_1_1_OA_3000_Examples.RData")



## ===============================
## 3. MICE imputation
## ===============================

## Keep only columns with "_bin" suffix and other necessary columns
dat_mice <- dat_mice_bin %>%
  select(
    # Continuous raw variables (for splines)
    Age, NLR, TMB, PDL1, Diameter, TumorPurity, FGA,
    
    # All columns with "_bin"
    contains("_bin"),
    
    # Other columns to keep (if they don't have "_bin" suffix)
    PFS_Months, PFS_Event, OS_Months, OS_Event,analysis_set
  )

## Force correction of continuous variable types to avoid being treated as ordered factor
dat_mice <- dat_mice %>%
  mutate(
    Age         = as.numeric(Age),
    NLR         = as.numeric(NLR),
    TMB         = as.numeric(TMB),
    PDL1        = as.numeric(PDL1),
    Diameter    = as.numeric(Diameter),
    TumorPurity = as.numeric(TumorPurity),
    FGA         = as.numeric(FGA)
  )


## The following comments outline the MICE imputation process, which is time-consuming. 
## Readers can choose to comment it out and directly load the post-imputation file to continue.

meth <- rep("", ncol(dat_mice))
names(meth) <- names(dat_mice)

for (v in names(dat_mice)) {
  x <- dat_mice[[v]]
  
  # Do not impute outcomes & analysis_set
  if (v %in% c("PFS_Months","PFS_Event","OS_Months","OS_Event","analysis_set")) {
    meth[v] <- ""
    next
  }
  
  # Continuous variables → pmm
  if (is.numeric(x)) {
    meth[v] <- "pmm"
    next
  }
  
  # Factors: Binary → logreg, Multiclass → polyreg
  if (is.factor(x)) {
    if (nlevels(x) == 2) {
      meth[v] <- "logreg"
    } else {
      meth[v] <- "polyreg"
    }
  }
}

meth

pred <- matrix(1, ncol(dat_mice), ncol(dat_mice))
colnames(pred) <- names(dat_mice)
rownames(pred) <- names(dat_mice)

# Do not use outcomes and analysis_set to predict others
pred[c("PFS_Months","PFS_Event","OS_Months","OS_Event","analysis_set"), ] <- 0

# Cannot predict itself
diag(pred) <- 0


set.seed(2025)

imp <- mice(
  dat_mice,
  m              = 10,      # >10 imputed datasets
  method         = meth,
  predictorMatrix = pred,
  maxit          = 20,      # Readers can give more iterations to avoid non-convergence
  printFlag      = TRUE
)


###Load pre-MICE data for drawing Kaplan-Meier survival curves stratified by each main index
load("mice_imputationbalanceAll_1_1_1_OA_3000_Examples.RData") 
##Save data after MICE imputation
load("mice_imputation_cleanbalanceAll_1_1_1_OA_3000_Examples.RData")  


## ===============================
## 4. Calculate the specific two-step construction process of the CLARITY-ICI Matrix
## ===============================

source("Multivariate_Cox_Regression_analysis_OA.R", encoding="utf-8")

##Training set, test set, validation set

covars_common <- c(
  "Age_bin","Patho_bin","Smoking_bin","NLR_bin","ECOG_bin","FGA_bin",
  "Diameter_bin","Tumor_Purity_bin","Stage_bin","P_Treat_bin",
  "IO_Type_bin","TMB_bin","PD_L1_bin","Sex_bin","Region_bin",
  "Line_bin","Real_World_bin"
)



specs <- list(
  "ECOG_bin", "NLR_bin", "FGA_bin", "Age_bin", "Smoking_bin", "Patho_bin","Diameter_bin",
  "Tumor_Purity_bin", "Stage_bin", "P_Treat_bin", "IO_Type_bin", "TMB_bin", "PD_L1_bin",
  "Sex_bin", "Region_bin", "Line_bin", "Real_World_bin"
)

results <- lapply(specs, function(mb) {
  run_module(
    imp_clean,
    main_bin   = mb,
    covars_bin = covars_common
  )
})
names(results) <- specs
#save(results, file = "results_run_modulebalanceAll_1_1_1_OA_3000_Examples.RData")
load("results_run_modulebalanceAll_1_1_1_OA_3000_Examples.RData")




## ===============================
## 5. Output the results of the two-step construction of the CLARITY-ICI Matrix
## ===============================

source("CLARITY-ICI Matrix_Output_For_plotting_OA.R", encoding="utf-8")

#####All CLARITY-ICI Matrix: Obtain multivariate Cox proportional hazards regression analysis results for PFS and OS across the three datasets.
all_res <- flatten_results2(results)

# Standardize fields (ensure keys are consistent for subsequent joins)
all_res <- all_res %>%
  mutate(
    analysis_set = as.character(analysis_set),
    endpoint     = as.character(endpoint),
    main_bin     = as.character(main_bin),
    main_level   = as.character(main_level),
    term         = as.character(term)
  ) %>%
  filter(
    analysis_set %in% c("Development","External_Test1","External_Test2"),
    endpoint %in% c("PFS","OS")
  )

#  Critical duplicate check: the same key should not appear multiple times
dup_keys <- all_res %>%
  count(main_bin, analysis_set, endpoint, main_level, term, name = "n") %>%
  filter(n > 1)

if (nrow(dup_keys) > 0) stop("Duplicated rows detected for same (main_bin, analysis_set, endpoint, main_level, term). Please check copying/overwriting.")


All_CLARITY_ICI_PFS <- get_all3_sig(all_res, endpoint_vec = "PFS", p_cut = 0.05) %>%
  arrange(main_bin, main_level, term, analysis_set)

All_CLARITY_ICI_OS  <- get_all3_sig(all_res, endpoint_vec = "OS", p_cut = 0.05) %>%
  arrange(main_bin, main_level, term, analysis_set)

All_CLARITY_ICI <- bind_rows(
  All_CLARITY_ICI_PFS %>% mutate(sig_type = "All3_PFS"),
  All_CLARITY_ICI_OS  %>% mutate(sig_type = "All3_OS")
) %>% arrange(main_bin, sig_type, main_level, term, analysis_set)

All_CLARITY_ICI_OUTPUT <- to_wide3(All_CLARITY_ICI, imp_clean = imp_clean, analysis_set = c("Development","External_Test1","External_Test2"))
All_CLARITY_ICI_OUTPUT
write.xlsx(All_CLARITY_ICI_OUTPUT, file = "Your personal path\\All_CLARITY_ICI_OUTPUT.xlsx")



#####Western CLARITY-ICI Matrix: Obtain multivariate Cox proportional hazards regression analysis results for PFS and OS across the Development and External_Test1 datasets.
all_res <- flatten_results2(results)
# Standardize fields (ensure keys are consistent for subsequent joins)
all_res <- all_res %>%
  mutate(
    analysis_set = as.character(analysis_set),
    endpoint     = as.character(endpoint),
    main_bin     = as.character(main_bin),
    main_level   = as.character(main_level),
    term         = as.character(term)
  ) %>%
  filter(
    analysis_set %in% c("Development","External_Test1"),
    endpoint %in% c("PFS","OS")
  )

# Critical duplicate check: the same key should not appear multiple times
dup_keys <- all_res %>%
  count(main_bin, analysis_set, endpoint, main_level, term, name = "n") %>%
  filter(n > 1)

if (nrow(dup_keys) > 0) stop("Duplicated rows detected for same (main_bin, analysis_set, endpoint, main_level, term). Please check copying/overwriting.")

Western_CLARITY_ICI_PFS <- get_all2_sig(all_res, endpoint_vec = "PFS", p_cut = 0.05) %>%
  arrange(main_bin, main_level, term, analysis_set)

Western_CLARITY_ICI_OS  <- get_all2_sig(all_res, endpoint_vec = "OS", p_cut = 0.05) %>%
  arrange(main_bin, main_level, term, analysis_set)

Western_CLARITY_ICI <- bind_rows(
  Western_CLARITY_ICI_PFS %>% mutate(sig_type = "All3_PFS"),
  Western_CLARITY_ICI_OS  %>% mutate(sig_type = "All3_OS")
) %>% arrange(main_bin, sig_type, main_level, term, analysis_set)

Western_CLARITY_ICI_OUTPUT <- to_wide3(Western_CLARITY_ICI, imp_clean = imp_clean, analysis_set = c("Development","External_Test1"))
#write.xlsx(Western_CLARITY_ICI_OUTPUT, file = "Your personal path\\Western_CLARITY_ICI_OUTPUT.xlsx")



#####East Asian CLARITY-ICI Matrix: Obtain multivariate Cox proportional hazards regression analysis results for PFS and OS across the Development and External_Test2 datasets.
all_res <- flatten_results2(results)
# Standardize fields (ensure keys are consistent for subsequent joins)
all_res <- all_res %>%
  mutate(
    analysis_set = as.character(analysis_set),
    endpoint     = as.character(endpoint),
    main_bin     = as.character(main_bin),
    main_level   = as.character(main_level),
    term         = as.character(term)
  ) %>%
  filter(
    analysis_set %in% c("Development","External_Test2"),
    endpoint %in% c("PFS","OS")
  )

# Critical duplicate check: the same key should not appear multiple times
dup_keys <- all_res %>%
  count(main_bin, analysis_set, endpoint, main_level, term, name = "n") %>%
  filter(n > 1)

if (nrow(dup_keys) > 0) stop("Duplicated rows detected for same (main_bin, analysis_set, endpoint, main_level, term). Please check copying/overwriting.")

library(dplyr)
library(tidyr)

EastAsian_CLARITY_ICI_PFS <- get_all2_sig(all_res, endpoint_vec = "PFS", p_cut = 0.05) %>%
  arrange(main_bin, main_level, term, analysis_set)

EastAsian_CLARITY_ICI_OS  <- get_all2_sig(all_res, endpoint_vec = "OS", p_cut = 0.05) %>%
  arrange(main_bin, main_level, term, analysis_set)

EastAsian_CLARITY_ICI <- bind_rows(
  EastAsian_CLARITY_ICI_PFS %>% mutate(sig_type = "All3_PFS"),
  EastAsian_CLARITY_ICI_OS  %>% mutate(sig_type = "All3_OS")
) %>% arrange(main_bin, sig_type, main_level, term, analysis_set)

EastAsian_CLARITY_ICI_OUTPUT <- to_wide3(EastAsian_CLARITY_ICI, imp_clean = imp_clean, analysis_set = c("Development","External_Test2"))
#write.xlsx(EastAsian_CLARITY_ICI_OUTPUT, file = "Your personal path\\EastAsian_CLARITY_ICI_OUTPUT.xlsx")


## ===============================
## 6. Plot the CLARITY-ICI Matrix, here only the Western CLARITY-ICI Matrix is shown, the other two are similar.
## ===============================


# Read Western CLARITY-ICI Matrix data, i.e., the data saved in the previous step.

workbook <- "Download the provided excel to your personal path\\Western_CLARITY_ICI_OUTPUT.xlsx"
df <- read_excel(workbook, guess_max = 5000)

need_cols <- c("main_bin","main_level","endpoint","term","ALL_median")
miss <- setdiff(need_cols, names(df))
if (length(miss) > 0) stop("Excel file missing columns: ", paste(miss, collapse = ", "))

# 1) Rows (y-axis)

df2 <- df %>%
  mutate(
    main_level = as.character(main_level),
    endpoint   = as.character(endpoint),
    
    # X-coordinate: term without "_bin"
    term = str_replace_all(as.character(term), "_bin", " "),
    
    ALL_Median = suppressWarnings(as.numeric(ALL_median)),
    
    # Y-coordinate: main_bin without "_bin"
    main_bin_clean = str_replace_all(as.character(main_bin), "_bin", ""),
    
    # Generate Y using the cleaned main_bin
    Y = paste(main_bin_clean, main_level, endpoint, sep = " ")
  ) %>%
  filter(!is.na(term), term != "", !is.na(Y), Y != "")


# Row ordering: by group_prefix + descending row mean
row_stat <- df2 %>%
  group_by(Y) %>%
  summarise(
    row_mean = ifelse(all(is.na(ALL_Median)), NA_real_, mean(ALL_Median, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    group_prefix = str_extract(Y, "^[^ ]+ [^ ]+")
  ) %>%
  arrange(group_prefix, desc(row_mean)) %>%
  mutate(row_num = row_number())

df2 <- df2 %>%
  left_join(row_stat %>% select(Y, group_prefix, row_num), by = "Y")


# Column ordering: by ascending term mean

term_order <- df2 %>%
  group_by(term) %>%
  summarise(term_mean = ifelse(all(is.na(ALL_Median)), NA_real_, mean(ALL_Median, na.rm = TRUE)),
            .groups = "drop") %>%
  arrange(term_mean) %>%
  pull(term)

df2 <- df2 %>%
  mutate(term = factor(term, levels = term_order))

# Calculate divider lines (by group_prefix)

divider_lines <- row_stat %>%
  group_by(group_prefix) %>%
  summarise(y_pos = max(row_num) + 0.5, .groups = "drop") %>%
  filter(!is.na(group_prefix)) %>%
  add_row(y_pos = 0.5) %>%
  add_row(y_pos = max(row_stat$row_num) + 0.5) %>%
  distinct(y_pos) %>%
  arrange(y_pos)

# Plot coloring

color_low  <- "#00468B"
color_high <- "#42B540"

vmin <- min(df2$ALL_Median, na.rm = TRUE)
vmed <- median(df2$ALL_Median, na.rm = TRUE)
vmax <- max(df2$ALL_Median, na.rm = TRUE)


p <- ggplot(df2, aes(x = term, y = reorder(Y, -row_num), fill = ALL_Median)) +
  geom_tile(color = "white", linewidth = 0.7) +
  scale_fill_gradientn(
    colours = c(color_low, "white", color_high),
    values  = scales::rescale(c(vmin, vmed, vmax)),
    na.value = "gray90",
    name = "Months"
  ) +
  geom_text(
    aes(label = ifelse(is.na(ALL_Median), "", round(ALL_Median, 1))),
    color = "black",
    size = 5.5,
    fontface = "bold"
  ) +
  geom_hline(
    data = divider_lines,
    aes(yintercept = y_pos),
    color = "black",
    linewidth = 0.8
  ) +
  labs(
    x = "Step 2: Median survival outcomes",
    y = "Step 1: Identify your interested variables and the target events",
    title = "Clinical Attribute–Based Risk Matrix (Western)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    text = element_text(family = "sans"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.title  = element_text(face = "bold", size = 18),
    plot.title  = element_text(hjust = 0.5, face = "bold", size = 20),
    panel.grid  = element_blank(),
    plot.margin = unit(c(1,1,1,1), "cm")
  ) +
  guides(fill = guide_colorbar(barwidth = 1.5, barheight = 15))

print(p) ##Western CLARITY-ICI Matrix plotting completed
