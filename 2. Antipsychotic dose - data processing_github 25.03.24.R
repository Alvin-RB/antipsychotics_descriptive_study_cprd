# ----------------
# "Prescribing of antipsychotics for people diagnosed with severe mental illness in UK primary care"

# 2. Antipsychotic dose - data processing_github 25.03.24 ####
# This script covers the data processing required to analyse antipsychotic dose over the first year of prescribing
# Subsequent analyses are carried out using the "Antipsychotic dose - analysis" script
# ----------------
# Last run: 25/03/2024

# Clear memory
# rm(list = ls())

# Packages
library(dplyr)
library(readr)
library(purrr)
library(gtsummary)
library(gt)
library(forcats)
library(lubridate)
library(tidyr)
library(doseminer)
library(drugprepr)
library(readxl)
library(chlorpromazineR)
library(stringr)
library(tidylog)
library(data.table)

# Set file path
path <- # removed
  
# Set working directory
# removed

# Load files
load("Cohort.Rdata") # file containing individual level data on patients diagnosed with SMI (including time-invariant variables)
load("Antipsychotics.Rdata") # file containing prescription level data on antipsychotics for the cohort

# Define parameters for saving plots
output_dir <- # removed
file_date <- format(Sys.Date(), "%Y-%m-%d")

# SELECT COHORT ####
# Note that these paramaters are as per the "Sample derivation and trends" script

overall_total <- n_distinct(Cohort$patid)

# Remove ineligible patients
Cohort <- Cohort %>%
  mutate(studyend = as.Date("2019-12-31"),
         studystart = as.Date("2000-01-01")) %>%
  rowwise() %>%
  mutate(overallend = pmin(deathdate, regenddate, lcd, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(regstartdate < overallend) %>% # remove those with regstart after death, regend or lcd
  filter(overallend > studystart) %>% # remove those where overallend is before the study start date
  filter(regstartdate <= studyend) %>% # remove those where regstart is after the study end date
  filter(first_diagnosis_date < overallend) # remove those where first diagnosis is after overall end (no potential for follow-up from diagnosis)

potentiallyeligible_total <- n_distinct(Cohort$patid)

# Actively registered (at least 30 days in study period)
Cohort <- Cohort %>%
  mutate(active = ifelse(year(regstartdate) <= 2019 & year(overallend) >= 2000, 1, 0),
         days_covered = ifelse(active == 1, as.numeric(pmin(as.Date(overallend), as.Date("2019-12-31")) - pmax(as.Date(regstartdate), as.Date("2000-01-01")) + 1), 0)) %>% # count the number of days the patient was actively registered in the study period
  filter(days_covered > 30) # remove those actively registered in study period for less than 30 days

activelyregistered_total <- n_distinct(Cohort$patid)

# SMI cohort (first diagnosed 2000-2017)
Cohort <- Cohort %>%
  filter(first_diagnosis_date >= '2000-01-01' & first_diagnosis_date < '2018-01-01')

firstdiag_2000to2017 <- n_distinct(Cohort$patid)

apdose <- Antipsychotics %>%
  select(patid, deliverymethod, issuedate, AP, productname, prodcode = prodcodeid, formulation, dosage_text, strength, duration, quantity) %>%
  filter(AP != "Prochlorperazine") %>% # # excl. Prochlorperazine as very commonly used but not for SMI
  distinct() %>%
  group_by(patid) %>%
  mutate(min = min(issuedate)) %>%
  ungroup() %>%
  mutate(year = year(issuedate)) %>%
  filter(year >= 2000 & year < 2020) # keep only 2000-2019

apdose <- Cohort %>% 
  select(patid, pracid, last_smi_diag) %>%
  inner_join(apdose, by = "patid")

# Select AP prescriptions in people first receiving an oral AP from year 2000-2019
ap_count <- n_distinct(apdose$patid)
consort1 <- data.frame(level = 1, Included = ap_count, Excluded = 0, Reason = "First diagnosed with SMI 2000-2017 and ever prescribed an AP")

apdose <- apdose %>%
  filter(deliverymethod == "Oral") %>%
  select(-deliverymethod)
  
oralap_count <- n_distinct(apdose$patid)
consort2 <- data.frame(level = 2, Included = oralap_count, Excluded = ap_count - oralap_count, Reason = "Received an oral AP")

apdose <- apdose %>%
  filter(min >= '2000-01-01' & min < '2020-01-01') %>% # keep only those first prescribed between 2000-2019
  arrange(., issuedate)

studyperiod_count <- n_distinct(apdose$patid)
consort3 <- data.frame(level = 3, Included = studyperiod_count, Excluded = oralap_count - studyperiod_count, Reason = "First received an oral AP 2000-2019")

# Keep only patients with > 1 prescription date
# Keep only up to first 12 prescriptions per patient
apdose <- apdose %>%
  group_by(patid) %>%
  mutate(date_rank = dense_rank(issuedate)) %>%
  ungroup() %>%
  group_by(patid) %>%
  mutate(sum_row = n_distinct(date_rank)) %>%
  ungroup() %>%
  filter(sum_row > 1) %>%
  filter(date_rank < 13) %>%
  select(-sum_row, -date_rank)

morethan1presc_count <- n_distinct(apdose$patid)
consort4 <- data.frame(level = 4, Included = morethan1presc_count, Excluded = studyperiod_count - morethan1presc_count, Reason = "Received more than prescription")

# CALCULATE NUMBER OF DAILY DOSES #### 

# Clean dosage texts with lookup tables
lookup <- read_excel("lookup.xlsx")
conditional_lookup <- read_excel("conditional_lookup.xlsx")

apdose <- apdose %>%
  mutate(dosage_text_raw = dosage_text) %>%
  mutate(dosage_text = gsub("\\*", "", dosage_text)) %>%
  left_join(conditional_lookup, by = c("dosage_text", "productname")) %>%
  mutate(success = ifelse(is.na(code), 0, 1),
         dosage_text = ifelse(success == 1, code, dosage_text)) %>%
  select(-success, -code) %>%
  left_join(lookup, by = "dosage_text") %>%
  mutate(success = ifelse(is.na(code), 0, 1),
         dosage_text = ifelse(success == 1, code, dosage_text)) %>%
  select(-success, -code, -`Dupe check`)

# Extract free text prescription instructions and convert to numeric values, using doseminer
free_text <- with(apdose, dosage_text[!duplicated(dosage_text) & nchar(dosage_text) > 0])
extracted <- doseminer::extract_from_prescription(free_text)
apdose <- merge(extracted, apdose, by.x = 'raw', by.y = 'dosage_text', all.x = TRUE)

# Create new variables as per doseminer instructions
apdose <- apdose %>%
  separate(dose, c('min_dose', 'max_dose'), sep = '-',
           convert = TRUE, fill = 'right') %>%
  separate(itvl, c('min_itvl', 'max_itvl'), sep = '-',
           convert = TRUE, fill = 'right') %>%
  separate(freq, c('min_freq', 'max_freq'), sep = '-',
           convert = TRUE, fill = 'right') %>%
  mutate(dose = coalesce((min_dose + max_dose) / 2, min_dose),
         itvl = coalesce((min_itvl + max_itvl) / 2, min_itvl, 1),
         freq = coalesce((min_freq + max_freq) / 2, min_freq),
         ndd = freq * dose / itvl,
         dose_num = readr::parse_number(strength), # create numeric value of dose
         dose_num = case_when(grepl("micro", productname) | grepl("micro", strength) ~ dose_num/1000, #Convert mcg to mg
                              TRUE ~ dose_num)) %>%
  select(-min_dose, -max_dose, -min_itvl, -max_itvl, -min_freq, -max_freq, -unit) %>%
  select(patid, pracid, last_smi_diag, issuedate, AP, productname, formulation, dosage_text_raw, dosage_text_processed = raw, strength, dose_num, ndd, everything())

# Manual cleaning
apdose <- apdose %>%
  mutate(quantity = case_when(patid == "removed" & issuedate == "removed" & quantity == "420" ~ 42,
                              TRUE ~ quantity))

apdose <- apdose %>%
  mutate(ndd_calc = case_when(!grepl("ml", strength) & quantity > 1 & quantity == duration & is.na(ndd) ~ round(quantity / duration, 1), # if qty/dur are the same & > 1 & AP is not prescribed in ml, calculate NDD
                              !grepl("ml", strength) & duration > 6 & duration <= 122 & quantity > 1 & quantity < 500 & is.na(ndd) ~ round(quantity / duration, 1), # if qty/dur are not the same, apply some limits
                              TRUE ~ NA_real_), # if criteria not met, set to NA
         ndd_calc = round(ndd_calc/.25)*.25, # round to nearest 0.25
         ndd_calc = case_when(ndd_calc == 0.0 ~ NA_real_, # if rounded value is 0, set to missing
                              TRUE ~ ndd_calc),
         ndd = coalesce(ndd, ndd_calc)) %>% # use newly calculated ndd if original ndd is missing
  select(-ndd_calc) %>%
  mutate(ndd_raw = ndd)

# Impute missing NDDs with various median values, starting more precise
apdose <- impute_ndd(apdose, 'median', group = c("patid", "prodcode"))
apdose <- impute_ndd(apdose, 'median', group = c("patid", "productname"))
apdose <- impute_ndd(apdose, 'median', group = c("patid", "AP", "strength"))
apdose <- impute_ndd(apdose, 'median', group = c("patid", "AP", "dose_num"))
apdose <- impute_ndd(apdose, 'median', group = c("pracid", "last_smi_diag", "prodcode"))
apdose <- impute_ndd(apdose, 'median', group = c("pracid", "last_smi_diag", "productname"))
apdose <- impute_ndd(apdose, 'median', group = c("pracid", "last_smi_diag", "AP", "strength"))
apdose <- impute_ndd(apdose, 'median', group = c("pracid", "last_smi_diag", "AP", "dose_num"))
apdose <- impute_ndd(apdose, 'median', group = c("last_smi_diag", "prodcode"))
apdose <- impute_ndd(apdose, 'median', group = c("last_smi_diag", "productname"))
apdose <- impute_ndd(apdose, 'median', group = c("last_smi_diag", "AP", "strength"))
apdose <- impute_ndd(apdose, 'median', group = c("last_smi_diag", "AP", "dose_num"))
apdose <- impute_ndd(apdose, 'median', group = c("AP", "strength"))
apdose <- impute_ndd(apdose, 'median', group = c("AP", "dose_num"))

apdose <- apdose %>%
  mutate(total_daily_dose = dose_num * ndd) # times dose in mg by number of tablets taken per day

# Convert to chlorpromazine and olanzapine equivalent doses 
# (use DDD method (leucht2016) in the first instance, otherwise use the consensus method i.e for droperidol (gardner2010) & cariprazine (Leuct2020)
DDD_consensus <- chlorpromazineR::add_key(base = leucht2016, add = gardner2010, trim = TRUE)
DDD_consensus <- chlorpromazineR::add_key(base = DDD_consensus, add = leucht2020, trim = TRUE)

apdose <- apdose %>%
  mutate(AP = case_when(AP == "Pericyazine" ~ "Periciazine", # to match spelling in the Leucht DDD key
                        TRUE ~ AP))

apdose <- as.data.frame(apdose) # base file must be a DF

apdose <- chlorpromazineR::to_ap(apdose, convert_to_ap = "olanzapine", convert_to_route = "oral", 
                                      ap_label = "AP", dose_label = "total_daily_dose", key = DDD_consensus)

# Limit to up to 4 prescriptions on a given medication per day, sum total dose, filter to one row per date
apdose <- apdose %>%
  group_by(patid, issuedate, AP) %>%
  mutate(num = row_number()) %>% # number of prescriptions of an individual AP per day per patient
  ungroup()

overthree <- apdose %>%
  filter(num >= 4)

apdose <- apdose %>%
  filter(num < 4) %>% # where multiple prescriptions on the same day, consider up to three prescriptions of each AP per patient, 
  group_by(patid, issuedate) %>%
  mutate(doseknownforall = !any(is.na(ndd)),
         dose_total = ifelse(doseknownforall != FALSE, sum(ap_eq), NA), # sum the total dose per prescription date, if there are none missing
           num = row_number()) %>%
  ungroup() %>%
  filter(num == 1) # filter to one row per prescription date per patient

apdose <- apdose %>%
  arrange(., issuedate) %>%
  group_by(patid) %>%
  mutate(num = row_number()) %>%
  ungroup()

# Create a column for each dose and each dose date (doses 1-12)
dosevars <- apdose %>%
  mutate(num = sprintf("%02d", num)) %>% # add leading zeros to column names
  pivot_wider(
    id_cols = patid,
    names_from = num,
    names_prefix = "dose",
    values_from = c("dose_total", "issuedate"),
    values_fill = NA)

# Format column names
names(dosevars) <- gsub("dose_total_", "", names(dosevars))
names(dosevars)[grepl("issuedate_", names(dosevars))] <- paste0(sub("issuedate_", "", grep("issuedate_", names(dosevars), value = TRUE)), "_date")

# Ungroup dfs
apdose <- apdose %>%
  ungroup()

# Remove dosevars that do no contribute any actual doses
dosevars <- dosevars %>%
  ungroup() %>%
  mutate(not_na_count = rowSums(!is.na(.[, 2:13]))) %>%
  filter(not_na_count > 1) %>%
  select(-not_na_count)

nodoses_count <- n_distinct(dosevars$patid)
consort5 <- data.frame(level = 5, Included = nodoses_count, Excluded = morethan1presc_count - nodoses_count, Reason = "No known doses in the first 12 (even after imputation)")

#Merge dose data with Cohort data
Antipsychotic_dose <- dosevars %>%
  left_join(Cohort, by = "patid")

# If the dose before and after a missing dose are the same, impute it as the same value (likely continued between the time points)
Antipsychotic_dose <- Antipsychotic_dose %>%
  mutate(dose02_calc = case_when(is.na(dose02) & !is.na(dose01) & !is.na(dose03) & dose01 == dose03 ~ dose01,
                                 TRUE ~ NA),
         dose03_calc = case_when(is.na(dose03) & !is.na(dose02) & !is.na(dose04) & dose02 == dose04 ~ dose02,
                                 TRUE ~ NA),
         dose04_calc = case_when(is.na(dose04) & !is.na(dose03) & !is.na(dose05) & dose03 == dose05 ~ dose03,
                                 TRUE ~ NA),
         dose05_calc = case_when(is.na(dose05) & !is.na(dose04) & !is.na(dose06) & dose04 == dose06 ~ dose04,
                                 TRUE ~ NA),
         dose06_calc = case_when(is.na(dose06) & !is.na(dose05) & !is.na(dose07) & dose05 == dose07 ~ dose05,
                                 TRUE ~ NA),
         dose07_calc = case_when(is.na(dose07) & !is.na(dose06) & !is.na(dose08) & dose06 == dose08 ~ dose06,
                                 TRUE ~ NA),
         dose08_calc = case_when(is.na(dose08) & !is.na(dose07) & !is.na(dose09) & dose07 == dose09 ~ dose07,
                                 TRUE ~ NA),
         dose09_calc = case_when(is.na(dose09) & !is.na(dose08) & !is.na(dose10) & dose08 == dose10 ~ dose08,
                                 TRUE ~ NA),
         dose10_calc = case_when(is.na(dose10) & !is.na(dose09) & !is.na(dose11) & dose09 == dose11 ~ dose09,
                                 TRUE ~ NA),
         dose11_calc = case_when(is.na(dose11) & !is.na(dose10) & !is.na(dose12) & dose10 == dose12 ~ dose10,
                                 TRUE ~ NA),
         dose02 = coalesce(dose02, dose02_calc),
         dose03 = coalesce(dose03, dose03_calc),
         dose04 = coalesce(dose04, dose04_calc),
         dose05 = coalesce(dose05, dose05_calc),
         dose06 = coalesce(dose06, dose06_calc),
         dose07 = coalesce(dose07, dose07_calc),
         dose08 = coalesce(dose08, dose08_calc),
         dose09 = coalesce(dose09, dose09_calc),
         dose10 = coalesce(dose10, dose10_calc),
         dose11 = coalesce(dose11, dose11_calc)) %>%
  select(-dose02_calc, dose03_calc, dose04_calc, dose05_calc, dose06_calc, dose07_calc, dose08_calc, dose09_calc, dose10_calc, dose11_calc)

#Save to file
n_distinct(Antipsychotic_dose$patid)
save(Antipsychotic_dose, file = "Antipsychotic_dose.Rdata")

# Misc ##
# Count number of prescription dates with doses available across dataset
consort6 <- Antipsychotic_dose %>%
  select(2:13) %>%
  summarise(Included = sum(!is.na(.))) %>%
  mutate(Reason = "Total number of included prescription dates",
         level = 6)

# Consort
consort <- rbind(consort1, consort2, consort3, consort4, consort5) %>%
  bind_rows(consort6) %>%
  select(-level) %>%
  gt() %>%
  fmt_number(columns = c("Included", "Excluded"), use_seps = TRUE, decimals = 0) %>%
  sub_missing(missing_text = "-") %>%
  gtsave(paste0(output_dir, "doseconsort_", file_date, ".docx"))

# How many NDDs required imputation?
imputednumbers <- apdose %>%
  filter(is.na(ndd_raw) & !is.na(ndd)) %>% # filter to rows where NDD was imputed (where raw ndd is not available)
  select(ndd) %>%
  tbl_summary(label = ndd ~ "Imputed NDD")

n_distinct(imputednumbers$patid)
nonimputednumbers <- apdose %>%
  filter(!is.na(ndd_raw)) %>%
  select(ndd_raw) %>%
  tbl_summary(label = ndd_raw ~ "Not imputed NDD")

allndds_imputedandnot <- apdose %>%
  select(ndd) %>%
  tbl_summary(label = ndd ~ "All NDDs")

no_patients_imputed <- length(unique(apdose$patid[is.na(apdose$ndd_raw) & !is.na(apdose$ndd)]))
total_ndds_imputedandnot <- sum(!is.na(apdose$ndd))
total_imputed_ndds <- sum(!is.na(apdose$ndd) & is.na(apdose$ndd_raw))
imputed_table_footnote <- paste0("In total, ", format(total_imputed_ndds, big.mark = ","), " of ", format(total_ndds_imputedandnot, big.mark = ","), " (",  sprintf("%.1f%%", (total_imputed_ndds / total_ndds_imputedandnot) * 100), ") were imputed. The imputed NDDs were across ", format(no_patients_imputed, big.mark = ","), " patients.")

imputedndds <- tbl_merge(tbls = list(allndds_imputedandnot, imputednumbers, nonimputednumbers), tab_spanner = c("All", "Imputed", "Not imputed")) %>%
  as_gt() %>%
  tab_footnote(imputed_table_footnote) %>%
  gtsave(paste0(output_dir, "imputed_ndds_numbers_", file_date, ".docx"))

# Count the number of prescrptions excluded due to being possible duplicates, in patients that end up in the final cohort
overthree_tab <- Antipsychotic_dose %>%
  select(patid) %>%
  inner_join(overthree, by = "patid") %>%
  summarise(distinct_pts = n_distinct(patid),
            total_prescriptions = n()) %>%
  gt() %>%
  cols_label(distinct_pts = "Number of patients with prescriptions considered erroneous (>3 duplicate)",
             total_prescriptions = "Total number of prescriptions excluded due to being considered erronenous (>3 duplicate)") %>%
  fmt_number(columns = c("distinct_pts", "total_prescriptions"), use_seps = TRUE, decimals = 0) %>%
  gtsave(paste0(output_dir, "erroneous_prescriptions_", file_date, ".docx"))

# Package versions
packageVersion("dplyr") # 1.1.4
packageVersion("readr") # 2.1.5
packageVersion("purrr") # 1.0.2
packageVersion("gtsummary") # 1.7.2
packageVersion("gt") # 0.10.1
packageVersion("forcats") # 1.0.0
packageVersion("lubridate") # 1.9.2
packageVersion("tidyr") # 1.3.0
packageVersion("doseminer") # 0.1.2
packageVersion("drugprepr") # 0.0.4
packageVersion("readxl") # 1.4.3
packageVersion("chlorpromazineR") # 0.2.0
packageVersion("stringr") # 1.5.0
packageVersion("tidylog") # 1.0.2
packageVersion("data.table") # 1.14.10

