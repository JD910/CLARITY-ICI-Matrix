# Clinical Attribute–Based Risk Matrix (CLARITY-ICI Matrix)

## Prognostication in Immunotherapy-Treated Non-Small-Cell Lung Cancer Individual-Patient Data Analysis of 17,051 Patients

This repository provides the full, executable R workflow used to construct the CLARITY-ICI Matrix for NSCLC patients receiving immune checkpoint inhibitor (ICI) therapy. 

The code reproduces the two-step modeling and matrix-generation process described in the study and allows users to generate a demonstrative matrix based on 3,000 patient samples.

All required .xlsx and .RData files are included in this repository. After downloading the repository, users can run the workflow end-to-end without modification except for setting a local working path.

---

## 📂 Repository Contents

Main_CLARITY_ICI_Matrix.R — Main pipeline script  
Datasets_Development_External_Tests_OA.R  
Multivariate_Cox_Regression_analysis_OA.R  
CLARITY-ICI Matrix_Output_For_plotting_OA.R  

NSCLC_CLARITY_ICI_Matrix_OA.xlsx  
partition_results_balanceAll_1_1_1_OA_3000_Examples.RData  
mice_imputationbalanceAll_1_1_1_OA_3000_Examples.RData  
mice_imputation_cleanbalanceAll_1_1_1_OA_3000_Examples.RData  
results_run_modulebalanceAll_1_1_1_OA_3000_Examples.RData  

Western_CLARITY_ICI_OUTPUT.xlsx (example matrix for plotting)

---

## 🧠 What the Pipeline Does

This repository demonstrates the full CLARITY-ICI Matrix construction workflow, including:

### 1. Data preprocessing and clinical variable transformation

### 2. Region-stratified dataset partitioning

Development  
External Test 1 (Western)  
External Test 2 (East Asian)

### 3. Multiple imputation using MICE

### 4. Multivariable Cox proportional hazard regression modeling

### 5. Two-step CLARITY-ICI Matrix construction

### 6. Output of matrix tables (PFS / OS)

### 7. Plotting the CLARITY-ICI Matrices

The included example reproduces the matrix based on 3,000 patients for demonstration purposes.

---

## 📦 Software Requirements

Please install the following R packages:

```r
install.packages(c(
  "readxl","dplyr","tidyverse","mice","survival",
  "openxlsx","writexl","stringr","ggplot2"
))
```

---

## ▶️ Step-by-Step Usage

### Step 1 — Download the Repository

Clone or download as ZIP and unzip locally


### Step 2 — Open R / RStudio and Set Working Directory

```r
setwd("Your personal path")
```

### Step 3 — Load Required Excel and RData Files

All required files are already included in this repository.

The main script automatically loads:

NSCLC_CLARITY_ICI_Matrix_OA.xlsx  
Partition metadata (partition_results_*.RData)  
Pre-imputed datasets  
Post-imputation cleaned data  
Cox model outputs  

### Step 4 — Run the Main Pipeline Script

Run:

```bash
Main_CLARITY_ICI_Matrix.R
```

This executes:

preprocessing  
dataset partitioning  
(optional) MICE imputation  
multivariable Cox regression for each stratum 
CLARITY-ICI matrix construction  
export of matrix tables  

Users may skip full imputation and directly load the prepared .RData files (as indicated in the script comments) to reduce runtime.

### Step 5 — Generate CLARITY-ICI Matrix Outputs

The script produces:

All-population CLARITY-ICI Matrix  
Western CLARITY-ICI Matrix  
East Asian CLARITY-ICI Matrix  

Matrix tables are exported as .xlsx files and can be directly used for visualization or interpretation.

### Step 6 — Plot the Western CLARITY-ICI Matrix as example

Load the provided output:

Western_CLARITY_ICI_OUTPUT.xlsx

Then run the plotting section in the script to generate the heatmap illustrating:

Step 1 — stratifying key clinical attributes  
Step 2 — summarizing median survival outcomes  

This reproduces the Western-population CLARITY-ICI Matrix example.

---

## 📝 Notes for Users

The script reproduces the two-step CLARITY-ICI Matrix framework exactly as implemented in the study.

All uploaded datasets are anonymized and provided solely for educational and methodological demonstration.

Users may adapt the workflow to:

apply the matrix framework to new cohorts  
conduct subgroup modeling  
extend to additional endpoints or populations  

---

## 📖 Citation

If you use this code or workflow, please cite:

Clinical Attribute–Based Risk Matrix for Prognostication in Immunotherapy-Treated Non-Small-Cell Lung Cancer:  
An Individual-Patient Data Analysis of 17,051 Patients

---

## 📬 Contact

For questions or collaboration inquiries, please feel free to open an issue or contact the study authors.
