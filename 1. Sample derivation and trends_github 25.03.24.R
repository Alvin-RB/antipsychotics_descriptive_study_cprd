# ----------------
# "Prescribing of antipsychotics for people diagnosed with severe mental illness in UK primary care"

# 1. Sample derivation and trends_github 25.03.24.R ####
# This script covers creating the study cohort and analyses of antipsychotic prescribing trends
# ----------------
# Last run: 25/03/2024

# Clear memory
# rm(list = ls())

# Packages
library(dplyr)
library(gtsummary)
library(gt)
library(lubridate)
library(ggplot2)
library(scales)
library(readr)
library(reshape2)
library(gridExtra)
library(tibble)
library(tidylog)
library(tidyr)
library(stringr)
library(forcats)
library(patchwork)
library(viridis)
library(ggrepel)

# Set correct file path
path <- # removed

# Load files
load("Cohort.Rdata") # file containing individual level data on patients diagnosed with SMI (including time-invariant variables)
load("Antipsychotics.Rdata") # file containing prescription level data on antipsychotics for the cohort

# Define parameters for saving plots
output_dir <- # removed
file_date <- format(Sys.Date(), "%Y-%m-%d")
file_type <- ".png"

# DEFINE SAMPLE ####

overall_total <- n_distinct(Cohort$patid)
consort1 <- data.frame(level = 1, Included = overall_total, Excluded = 0, Reason = "Overall sample")

# Remove ineligible patient records
Cohort <- Cohort %>%
  mutate(studyend = as.Date("2019-12-31"), # study end date
       studystart = as.Date("2000-01-01")) %>% # study start date
  rowwise() %>%
  mutate(overallend = pmin(deathdate, regenddate, lcd, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(regstartdate < overallend) %>% # remove those with regstart after death, regend or lcd
  filter(overallend > studystart) %>% # remove those where death, regend or lcd is before the study start date
  filter(regstartdate <= studyend) %>% # remove those where regstart is after the study end date
  filter(first_diagnosis_date < overallend) # remove those where first diagnosis is after overall end (no potential for follow-up from diagnosis)

potentiallyeligible_total <- n_distinct(Cohort$patid)
consort2 <- data.frame(level = 2, Included = potentiallyeligible_total, Excluded = overall_total - potentiallyeligible_total, Reason = "Potentially eligible")

# Actively registered (at least 30 days in study period)
Cohort <- Cohort %>%
  mutate(active = ifelse(year(regstartdate) <= 2019 & year(overallend) >= 2000, 1, 0),
         days_covered = ifelse(active == 1, as.numeric(pmin(as.Date(overallend), as.Date("2019-12-31")) - pmax(as.Date(regstartdate), as.Date("2000-01-01")) + 1), 0)) %>% # count the number of days the patient was actively registered in the study period
  filter(days_covered > 30) # remove those actively registered in study period for less than 30 days

activelyregistered_total <- n_distinct(Cohort$patid)
consort3 <- data.frame(level = 3, Included = activelyregistered_total, Excluded = potentiallyeligible_total - activelyregistered_total, Reason = "Actively registered in study period")

# SMI cohort (first diagnosed 2000-2017)
Cohort <- Cohort %>%
  filter(first_diagnosis_date >= '2000-01-01' & first_diagnosis_date < '2018-01-01') %>%
  mutate(diagnosis_year = year(first_diagnosis_date))

firstdiag_2000to2017 <- n_distinct(Cohort$patid) 
consort4 <- data.frame(level = 4, Included = firstdiag_2000to2017, Excluded = activelyregistered_total - firstdiag_2000to2017, Reason = "First diagnosed 2000-2017")

# Antipsychotics sample (prescriptions 2000-2019)
define_sample_parameters <- function(df) { 
  df %>%
    distinct() %>% # duplicate prescriptions not important here as reporting at individual AP level
    filter(AP != "Prochlorperazine") %>% # exclude as not used as an AP
    mutate(year = year(issuedate)) %>% # mutate year of prescription
    filter(year >= 2000 & year < 2020) %>% # keep only 2000-2019 prescriptions
    inner_join(select(Cohort, patid, first_diagnosis_date, diagnosis_year, last_smi_diag, first_ap_date, ethnicity_cat_cprdhes, gender), by = "patid")} # merge with diagnosis and other characteristics

antipsychotics_sample <- as_tibble(Antipsychotics) %>%
  select(patid, issuedate, deliverymethod, AP) %>% # minimal AP details required, as reporting at the individual AP level
  define_sample_parameters()

total_prescriptions <- nrow(antipsychotics_sample)
consort5 <- data.frame(level = 5, Included = total_prescriptions, Excluded = "-", Reason = "Total number of AP prescriptions")

# LAI sample (prescriptions 2000-2019)
injections_sample <- as_tibble(Antipsychotics) %>%
  select(patid, issuedate, AP, deliverymethod) %>% # minimal AP details required, as reporting at the individual AP level
  filter(deliverymethod == "Injection") %>%
  select(-deliverymethod) %>%
  define_sample_parameters()

total_LAI_prescriptions <- nrow(injections_sample)
consort6 <- data.frame(level = 5, Included = total_LAI_prescriptions, Excluded = "-", Reason = "Total number of AP LAI prescriptions")

# Create/save consort flow
consort <- rbind(consort1, consort2, consort3, consort4, consort5, consort6) %>%
  gt() %>%
  gtsave(paste0(output_dir, "trendsconsort_", today(), ".docx"))

# COMPARISON OF PATIENTS PRESCRIBED VERSUS NOT PRESCRIBED ANTIPSYCHOTICS ####

# Derive additional fields - concomitant medications
process_comeds <- function(df) { 
  df %>%
    select(patid, issuedate) %>%
    distinct() %>%
    mutate(year = year(issuedate)) %>%
    filter(year >= 2000 & year < 2020) %>% # keep only 2000-2019
    select(patid) %>%
    distinct()}
    
load("MoodStabilisers.Rdata") # file containing prescription level data on mood stabilisers

moodstabilisers <- MoodStabilisers %>%
  process_comeds() %>%
  mutate(moodstab_bin = 1)

load("Antidepressants.Rdata") # file containing prescription level data on antidepressants

antidepressants <- Antidepressants %>%
  process_comeds() %>%
  mutate(antid_bin = 1)

# Antipsychotic injection in study period
injection <- injections_sample %>%
  select(patid, issuedate, AP) %>%
  group_by(patid) %>%
  mutate(injectionap_rec_bin_cohort = 1, # received antipsychotic injection prescription during study period
         injection_firstdate_cohort = min(issuedate)) %>% # first date of antipsychotic injection prescription during study period
  select(patid, injectionap_rec_bin_cohort, injection_firstdate_cohort) %>%
  ungroup() %>%
  distinct()

# Oral AP in study period
oral_aps <- antipsychotics_sample %>%
  select(patid, issuedate, AP, deliverymethod) %>%
  filter(deliverymethod == "Oral") %>%
  group_by(patid) %>%
  mutate(oralap_rec_bin_cohort = 1, # received oral antipsychotic during study period
         first_oralap_date_cohort = min(issuedate), # first date of oral antipsychotic prescription during study period
         last_oralap_date_cohort = max(issuedate)) %>% # last date of oral antipsychotic prescription during study period
  ungroup() %>%
  select(patid, oralap_rec_bin_cohort, first_oralap_date_cohort, last_oralap_date_cohort) %>%
  distinct()

# Select and derive variables for inclusion in table
characteristics_table_sample <- Cohort %>%
  left_join(antipsychotics_sample %>%
              mutate(ap_in_study = 1) %>% # received an antipsychotic during study period
              select(patid, ap_in_study) %>%
              distinct(),
            by = "patid") %>%
  left_join(moodstabilisers, by = "patid") %>%
  left_join(antidepressants, by = "patid") %>%
  left_join(injection, by = "patid") %>%
  left_join(oral_aps, by = "patid") %>%
  rowwise() %>%
  mutate(overallstart = pmax(regstartdate, first_diagnosis_date, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(ap_in_study = if_else(is.na(ap_in_study), "Not prescribed antipsychotic", "Prescribed antipsychotic"),
         gender = case_when(gender == "Indeterminate" | gender == "Unknown" ~ NA,
                            TRUE ~ gender),
         ethnicity_cat_cprdhes = str_to_title(ethnicity_cat_cprdhes),
         last_smi_diag = case_when(last_smi_diag == "bipolar" ~ "Bipolar disorder",
                                   last_smi_diag == "other psychosis" ~ "Other non-organic psychoses",
                                   last_smi_diag == "schizophrenia" ~ "Schizophrenia"),
         oralap_rec_bin_cohort = case_when(oralap_rec_bin_cohort == 1 ~ "Yes",
                                           TRUE ~ "No"),
         injectionap_rec_bin_cohort = case_when(injectionap_rec_bin_cohort == 1 ~ "Yes",
                                                TRUE ~ "No"),
         moodstab_bin = case_when(moodstab_bin == 1 ~ "Yes",
                                  TRUE ~ "No"),
         antid_bin = case_when(antid_bin == 1 ~ "Yes",
                               TRUE ~ "No"),
         noantidormoodstab_bin = case_when(moodstab_bin == "No" & antid_bin == "No" ~ "No antidepressant/mood stabiliser",
                                           TRUE ~ "At least one antidepressant or mood stabiliser"),
         diagnosis_to_endfup = as.numeric((overallend - first_diagnosis_date)/365.25), # years
         timeactivelyregistered = days_covered / 365.25, # years
         time_diagtooral_cohort = first_oralap_date_cohort - first_diagnosis_date, # days
         time_diagtolai_cohort = (injection_firstdate_cohort - first_diagnosis_date)/365.25, # years
         age_atfirstoralap_cohort = year(first_oralap_date_cohort) - yob,
         age_atfirstoralap_cohort_cat = case_when(age_atfirstoralap_cohort < 30 ~ "<30",
                                                  age_atfirstoralap_cohort >= 30 & age_atfirstoralap_cohort < 40 ~ "30-49",
                                                  age_atfirstoralap_cohort >= 40 & age_atfirstoralap_cohort < 65 ~ "40-64",
                                                  age_atfirstoralap_cohort >= 65 ~ "65+",
                                                  TRUE ~ "Unknown"),
         age_atfirstinjectionap_cohort = year(injection_firstdate_cohort) - yob,
         ap_priorto_diag = case_when(first_ap_date < first_diagnosis_date ~ "Yes",
                                     TRUE ~ "No"),
         firstap_timeperiod = case_when(first_ap_date < '2000-01-01' ~ "<2000", # first antipsychotic in or out of study period (note that first_ap_date is a different field from Cohort, not the one created above for in the study period)
                                        first_ap_date >= '2000-01-01' & first_ap_date < '2010-01-01' ~ "2000-2009",
                                        first_ap_date >= '2010-01-01' & first_ap_date < '2020-01-01' ~ "2010-2019",
                                        TRUE ~ NA)) %>%
  select(ap_in_study,
         patid, gender, ethnicity_cat_cprdhes, region, patprac_2019imd_quintile, timeactivelyregistered,
         last_smi_diag, age_atfirstdiag, age_atfirstdiag_cat, diagnosis_year, moodstab_bin, antid_bin, noantidormoodstab_bin, diagnosis_to_endfup,
         firstap_timeperiod, ap_priorto_diag, oralap_rec_bin_cohort, time_diagtooral_cohort, age_atfirstoralap_cohort, age_atfirstoralap_cohort_cat,
         injectionap_rec_bin_cohort, time_diagtolai_cohort, age_atfirstinjectionap_cohort, firsttolastap) %>%
  select(-patid)

# Specify table parameters
# Labels
label_list <- list(
  ethnicity_cat_cprdhes = "Ethnicity, n (%)",
  gender = "Sex, n (%)",
  region = "Geographic region, n (%)",
  patprac_2019imd_quintile = "IMD quintile, n (%)",
  last_smi_diag	= "SMI diagnosis, n (%)",
  age_atfirstdiag = "Age at first SMI diagnosis, median (IQR)",
  age_atfirstdiag_cat = "Age at first SMI diagnosis (category), n (%)",
  diagnosis_to_endfup = "Time from SMI diagnosis to end of follow-up (years), median (IQR)",
  timeactivelyregistered = "Time actively registered in study period (years), median (IQR)",
  time_diagtooral_cohort = "Time from SMI diagnosis to first oral antipsychotic (days), median (IQR)",
  oralap_rec_bin_cohort = "Ever prescribed oral antipsychotic, n (%)",
  injectionap_rec_bin_cohort = "Ever prescribed LAI antipsychotic, n (%)",
  time_diagtolai_cohort = "Time from SMI diagnosis to first LAI antipsychotic (years), median (IQR)",
  age_atfirstoralap_cohort = "Age at first oral antipsychotic, median (IQR)",
  age_atfirstoralap_cohort_cat = "Age at first oral antipsychotic category, n (%)",
  age_atfirstinjectionap_cohort = "Age at first LAI, median (IQR)",
  firstap_timeperiod = "Antipsychotic initiation time-period, n (%)",
  firsttolastap = "Time from first to last antipsychotic (years), median (IQR)",
  moodstab_bin = "Prescribed a mood stabiliser, n (%)",
  antid_bin = "Prescribed an antidepressant, n (%)",
  noantidormoodstab_bin = "No mood stabiliser or antidepressant, n (%)",
  diagnosis_year = "Year of SMI diagnosis, median (IQR)",
  ap_priorto_diag = "Prescribed antipsychotic prior to SMI diagnosis date, n (%)")

# Digits
digits_list <- list(
  diagnosis_year ~ 0,
  all_categorical() ~ c(0, 1))

# Footnotes
noap <- sum(characteristics_table_sample$ap_in_study == "Not prescribed antipsychotic")
prior <- sum(characteristics_table_sample$ap_in_study == "Not prescribed antipsychotic" & characteristics_table_sample$firstap_timeperiod == "<2000", na.rm = TRUE)

footnote_text <- paste0(
  "Amongst the ", format(noap, big.mark = ","), " patients not prescribed an antipsychotic during the study period, ", format(prior, big.mark = ","), 
  " (", trimws(sprintf("%.2f%%", prior / noap * 100)), ") were prescribed an antipsychotic prior to the year 2000.")

# Build table using gtsummary > gt
table <- as.data.frame(characteristics_table_sample %>%
                         tbl_summary(by = ap_in_study,
                                     label = label_list,
                                     sort = everything() ~ "alphanumeric",
                                     digits = digits_list)) %>%
                         rename_all(~gsub("\\*", "", .)) %>%
                         mutate(row_num = row_number()) %>%
                         filter(!row_num %in% c(4, 5, 56, 60, 62, 68, 71, 73, 75)) %>% # remove 'unknown' rows that are not relevant
                         mutate_at(vars(-row_num), ~ifelse(row_num == 45, gsub("2,", "2", .), .)) %>% # remove commas in calendar years
                         add_row(row_num = 0, Characteristic = "DEMOGRAPHICS") %>% # add table sub-headings
                         add_row(row_num = 34.5, Characteristic = "MENTAL HEALTH") %>%
                         add_row(row_num = 51.5, Characteristic = "ANTIPSYCHOTICS") %>%
                         arrange(row_num) %>%
                         mutate(across(2, ~ifelse(row_num >= 52, "-", .))) %>% # set all antipsychotic variables to '-' as NA for people not prescribed
                         mutate(row_num = row_number()) %>%
                         mutate(across(2, ~ifelse(row_num %in% c(53, 61), NA, .))) %>% # categorical variable labels do not need the '-'
                         gt() %>% # convert to gt
                         cols_hide(row_num) %>%
                         sub_missing(missing_text = "") %>% # use blank for empty cells
                         cols_align(
                           align = "center", # center allign columns
                           columns = 2:3) %>%
                         tab_style(style = cell_text(weight = "bold"), # format spanners bold
                                   locations = cells_column_labels(everything())) %>%
                         tab_style(style = cell_text(weight = "bold"), # format spanners bold
                                   locations = cells_body(columns = 1, rows = c(1, 2, 6, 13, 26, 33, 34, 35, 39, 40, 45, 46, 47, 48, 51, 52, 53, 57, 58, 59, 60, 61, 66, 67, 68, 69))) %>%
                         tab_footnote("IMD, index of multiple deprivation, severe mental illness; LAI, long-acting injectable.") %>%
                         tab_footnote(footnote_text, locations = cells_body(columns = 1, rows = 53)) %>%
                         tab_footnote("Amongst patients registered at primary care practices in England only.", locations = cells_body(columns = 1, rows = 26)) %>%
                         tab_footnote("During the study period, 2000-2019.", locations = cells_body(columns = 1, rows = c(46, 47, 48)))
                       
# Generate a characteristics table for each SMI diagnosis (requires minor changes to the above overall code given that SMI diagnosis was a variable in the above)
generate_diagnosis_table <- function(diagnosis) {
  table <- as.data.frame(characteristics_table_sample %>%
                            filter(last_smi_diag == diagnosis) %>%
                            tbl_summary(by = ap_in_study,
                                        label = label_list,
                                        sort = everything() ~ "alphanumeric",
                                        digits = digits_list)) %>%
    rename_all(~gsub("\\*", "", .)) %>%
    mutate(row_num = row_number()) %>%
    filter(!row_num %in% c(4, 5, 35, 36, 54, 58, 60, 66, 69, 71, 73)) %>% # remove 'unknown' rows that are not relevant
    mutate_at(vars(-row_num), ~ifelse(row_num == 43, gsub("2,", "2", .), .)) %>% # remove commas in calendar years
    add_row(row_num = 0, Characteristic = "DEMOGRAPHICS") %>% # add table sub headings
    add_row(row_num = 34.5, Characteristic = "MENTAL HEALTH") %>%
    add_row(row_num = 49.5, Characteristic = "ANTIPSYCHOTICS") %>%
    arrange(row_num) %>%
    mutate(across(2, ~ifelse(row_num >= 50, "-", .))) %>% # set all antipsychotic variables to '-' as NA for people not prescribed
    mutate(row_num = row_number()) %>%
    mutate(across(2, ~ifelse(row_num %in% c(49, 57), NA, .))) %>% # categorical variable labels do not need the '-'
    gt() %>%
   cols_hide(row_num) %>%
    sub_missing(missing_text = "") %>% # use - for empty cells
    cols_align(
      align = "center", # center allign columns
      columns = 2:3) %>%
    tab_style(style = cell_text(weight = "bold"), # format spanners bold
              locations = cells_column_labels(everything())) %>%
    tab_style(style = cell_text(weight = "bold"), # format spanners bold
              locations = cells_body(columns = 1, rows = c(1, 2, 6, 13, 26, 33, 34, 35, 36, 41, 42, 43, 44, 47, 48, 49, 53, 54, 55, 56, 57, 62, 63, 64, 65))) %>%
    tab_footnote("IMD, index of multiple deprivation, severe mental illness; LAI, long-acting injectable. Row percentages are shown for comparisons of between-group characteristics. Column percentages are shown for characteristics relating to just the patients that were prescribed antipsychotics.") %>%
    tab_footnote("Amongst patients registered at primary care practices in England only.", locations = cells_body(columns = 1, rows = 26)) %>%
    tab_footnote("During the study period, 2000-2019.", locations = cells_body(columns = 1, rows = c(42, 43, 44)))}

schz_table <- generate_diagnosis_table("Schizophrenia")
bipolar_table <- generate_diagnosis_table("Bipolar disorder")
otherpsych_table <- generate_diagnosis_table("Other non-organic psychoses")

# LAI information stratified by ethnicity
ethnicity_lais <- characteristics_table_sample %>%
  select(ethnicity_cat_cprdhes, ap_in_study, injectionap_rec_bin_cohort, age_atfirstinjectionap_cohort, time_diagtolai_cohort) %>%
  tbl_summary(by = ethnicity_cat_cprdhes,
              digits = all_categorical() ~ c(0, 1)) %>%
  bold_labels()

# Time from antipsychotic to diagnosis
timefromaptodiag <- antipsychotics_sample %>%
  select(patid, first_ap_date, first_diagnosis_date) %>%
  distinct() %>%
  filter(first_ap_date < first_diagnosis_date) %>%
  mutate(diff_in_months = as.numeric(difftime(first_diagnosis_date, first_ap_date, units = "days")) / 30.44) %>%
  select(diff_in_months) %>%
  tbl_summary(label = diff_in_months ~ "Number of months between first AP and first SMI diagnosis (amongst those receiving an AP prior to diagnosis), median (IQR)")

#Save to file
table %>%
  gtsave(paste0(output_dir, "patientscharacteristics_ap_vs_noap_", today(), ".docx"))

schz_table %>%
  gtsave(paste0(output_dir, "schz_patientscharacteristics_ap_vs_noap_", today(), ".docx"))

bipolar_table %>%
  gtsave(paste0(output_dir, "bipolar_patientscharacteristics_ap_vs_noap_", today(), ".docx"))
  
otherpsych_table %>%
  gtsave(paste0(output_dir, "otherpsych_patientscharacteristics_ap_vs_noap_", today(), ".docx"))
  
ethnicity_lais %>%
  as_gt() %>%
  gtsave(paste0(output_dir, "ethniciy_lais_", today(), ".docx"))

timefromaptodiag %>%
  as_gt() %>%
  gtsave(paste0(output_dir, "timefromaptodiag_", today(), ".docx"))

# Standardised mean differences 
smd_table <- as.data.frame(characteristics_table_sample %>%
                             select(ap_in_study, gender, ethnicity_cat_cprdhes, region, patprac_2019imd_quintile,
                                    last_smi_diag, age_atfirstdiag, diagnosis_year, moodstab_bin, antid_bin, timeactivelyregistered) %>%
                             tbl_summary(by = ap_in_study,
                                         label = list(
                                           gender ~ "Sex", 
                                           ethnicity_cat_cprdhes ~ "Ethnicity",
                                           region ~ "Geographic region", 
                                           patprac_2019imd_quintile ~ "Deprivation quintile",
                                           last_smi_diag ~ "SMI diagnosis",
                                           age_atfirstdiag ~ "Age at first SMI diagnosis",
                                           diagnosis_year ~ "Year of first SMI diagnosis",
                                           moodstab_bin ~ "Prescribed a mood stabiliser", 
                                           antid_bin ~ "Prescribed an antidepressant", 
                                           timeactivelyregistered ~ "Time actively registered in primary care during study period"),
                                         sort = everything() ~ "alphanumeric") %>%
                             add_difference(test = list(
                               all_continuous() ~ "cohens_d",
                               all_categorical() ~ "smd", 
                               group = "ap_in_study"))) %>%
  select(1, 4, 5) %>%
  rename_all(~gsub("\\*", "", .)) %>%
  filter(!is.na(Difference)) %>%
  separate(3, into = c("lower_ci", "upper_ci"), sep = ", ") %>%
  mutate_at(vars(2:4), as.numeric)

# Create plot
smd_table$Characteristic <- factor(smd_table$Characteristic, levels = unique(sort(smd_table$Characteristic)))

smd_plot <- ggplot(smd_table, aes(y = factor(Characteristic, levels = rev(levels(Characteristic))), x = Difference)) +
  geom_point(color = "darkred", size = 2) +   # Add points for the estimates with a different color
  geom_errorbarh(aes(xmax = upper_ci, xmin = lower_ci), size = 0.5, height = 0.2, color = "steelblue") +   # Add error bars for confidence intervals
  geom_vline(xintercept = c(-0.1, 0, 0.1), linetype = "dashed", color = "grey50", alpha = 0.5) +   # Add vertical lines at -0.1, 0, and 0.1
  labs(x = "Standardized nean difference", y = NULL) + # Add axis labels
  scale_x_continuous(breaks = c(-0.30, -0.10, 0, 0.10, 0.30), limits = c(-0.5, 0.5)) +  # Customize x-axis
  ggtitle("SMD (95% CI)") + # Add a heading for estimates column
  theme_minimal() + 
  theme(panel.grid.major.y = element_blank(),  # Remove grid lines
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(size = 12, hjust = 1),  # Adjust title
        plot.margin = margin(1, 3, 1, 1, "cm"),  # Adjust margins
        panel.border = element_rect(color = "black", fill = NA, size = 1)) +  # Add border
  geom_text(aes(label = sprintf("%.2f (%.2f, %.2f)", Difference, lower_ci, upper_ci), x = 0.40), 
            hjust = 0, vjust = 0.5, size = 3) +
  coord_cartesian(xlim = c(-0.35, 0.60))   # Set x-axis limits

ggsave(smd_plot,
       path = output_dir,
       filename = paste0("smd_forest_plot_", file_date, file_type),
       dpi = 300, width = 10, height = 12, bg = "white")

# ABSOLUTE NUMBER RECEIVING EACH ANTIPSYCHOTIC 2000-2019 ####

# Function to create plot
create_number_receiving_plot <- function(data) {
  
  absolute_num <- data %>%
    select(patid, AP) %>%
    group_by(AP) %>%
    summarize(num_patients = n_distinct(patid)) %>%
    arrange(num_patients) %>%
    filter(num_patients >= 50) # only include if at least 50 participants
  
  ggplot(absolute_num, aes(x = reorder(AP, desc(num_patients), FUN = sum), y = num_patients)) +
  geom_bar(stat = "identity", position = "stack", fill = "#32CD32") +
  xlab("") +
  ylab("Number of patients") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  coord_flip() +
  scale_y_continuous(labels = comma)}

number_receiving_each_ap <- create_number_receiving_plot(antipsychotics_sample)
number_receiving_each_lai_ap <- create_number_receiving_plot(injections_sample)

ggsave(number_receiving_each_ap,
       path = output_dir,
       filename = paste0("number_receiving_each_ap_", file_date, file_type),
       dpi = 300, width = 8, height = 9, bg = "white")

ggsave(number_receiving_each_lai_ap,
       path = output_dir,
       filename = paste0("number_receiving_each_lai_ap_", file_date, file_type),
       dpi = 300, width = 8, height = 9, bg = "white")

# Percentages of patients receiving each antipsychotic during study period
total_pts <- sum(characteristics_table_sample$ap_in_study == "Prescribed antipsychotic", na.rm = TRUE)

pcts <- antipsychotics_sample %>%
  select(patid, AP) %>%
  group_by(AP) %>%
  summarize(num_patients = n_distinct(patid)) %>%
  arrange(num_patients) %>%
  filter(num_patients >= 50) %>% # only include if at least 50 participants
  mutate(pct = round(num_patients / total_pts * 100, 1)) %>%
  gt() %>%
  gtsave(paste0(output_dir, "ap_pcts_", today(), ".docx"))

# What proportion of prescriptions are accounted for by olanzapine/risperidone/aripiprazole/quetiapine?
proportion_top4 <- antipsychotics_sample %>%
  select(patid, issuedate, AP) %>%
  mutate(main_four = case_when(AP %in% c("Aripiprazole", "Olanzapine", "Risperidone", "Quetiapine") ~ 1,
                               TRUE ~ 0)) %>%
  select(main_four) %>%
  tbl_summary(digits = all_categorical() ~ c(0, 1)) %>%
  as_gt() %>%
  gtsave(paste0(output_dir, "top4_proportionofaps_", today(), ".docx"))

# TRENDS OVER TIME ####

# Calculate the number of active SMI patients in each year 2000:2019 
# (active = diagnosed at the year (allowing up to 2 calendars years for recording delay), still registered, alive and not after LCD)
active_by_year <- Cohort %>%
  select(patid, last_smi_diag, first_diagnosis_date, regstartdate, overallend) %>%
  mutate(diagnosis_start_year = case_when(first_diagnosis_date > regstartdate ~ year(first_diagnosis_date) - 2, # if diagnosis after regstart, then use diagnosis year (-2, allowing for the diagnosis to be up to 2 years prior)
                                          TRUE ~ NA),
         reg_start_year = year(regstartdate),
         start_year = case_when(reg_start_year > diagnosis_start_year & !is.na(diagnosis_start_year) ~ reg_start_year,
                                is.na(diagnosis_start_year) ~ reg_start_year,
                                TRUE ~ diagnosis_start_year),
         end_year = year(overallend)) %>% # if regstart is after diagnosis start, then use reg start as this will be when f-up can start
  select(patid, last_smi_diag, start_year, end_year)

active_by_year_schizophrenia <- active_by_year %>%
  filter(last_smi_diag == "schizophrenia")

active_by_year_bipolar <- active_by_year %>%
  filter(last_smi_diag == "bipolar")

active_by_year_otherpsychosis <- active_by_year %>%
  filter(last_smi_diag == "other psychosis")

active_by_year <- active_by_year %>%
  select(patid, start_year, end_year)

# Define a function to create a binary indicator for if patient was active in a given study year
add_year_indicators <- function(df, start_year, end_year) {
  years <- 2000:2019 # Create a vector of study years
  for (year in years) {
    df <- df %>%
      mutate(!!paste0("year_", year) := as.integer(year >= start_year & year <= end_year))} # Add column for each year with active status
  return(df)}

# Call function for each data frame
active_by_year <- add_year_indicators(active_by_year, start_year, end_year)
active_by_year_schizophrenia <- add_year_indicators(active_by_year_schizophrenia, start_year, end_year)
active_by_year_bipolar <- add_year_indicators(active_by_year_bipolar, start_year, end_year)
active_by_year_otherpsychosis <- add_year_indicators(active_by_year_otherpsychosis, start_year, end_year)

# Convert from wide to long
active_by_year <- tidyr::gather(active_by_year, key = "year", value = "active", starts_with("year_"))
active_by_year_schizophrenia <- tidyr::gather(active_by_year_schizophrenia, key = "year", value = "active", starts_with("year_"))
active_by_year_bipolar <- tidyr::gather(active_by_year_bipolar, key = "year", value = "active", starts_with("year_"))
active_by_year_otherpsychosis <- tidyr::gather(active_by_year_otherpsychosis, key = "year", value = "active", starts_with("year_"))

# Convert year to numeric
active_by_year$year <- as.numeric(sub("year_", "", active_by_year$year))
active_by_year_schizophrenia$year <- as.numeric(sub("year_", "", active_by_year_schizophrenia$year))
active_by_year_bipolar$year <- as.numeric(sub("year_", "", active_by_year_bipolar$year))
active_by_year_otherpsychosis$year <- as.numeric(sub("year_", "", active_by_year_otherpsychosis$year))

# Aggregate count of active patients per year
aggregate_count <- function(df) { 
  df %>%
    group_by(year) %>%
    summarize(count_active = sum(as.numeric(active), na.rm = TRUE))}

active_counts <- active_by_year %>%
  aggregate_count()

active_counts_schizophrenia <- active_by_year_schizophrenia %>%
  aggregate_count()

active_counts_bipolar <- active_by_year_bipolar %>%
  aggregate_count()

active_counts_otherpsychosis <- active_by_year_otherpsychosis %>%
  aggregate_count()

# Prescribing rates per 1,000 active SMI patients ####

# Overall prescribing rate ####

# Count the number of people prescribed an AP each year
overall_presc_rate_counts <- antipsychotics_sample %>% 
  select(patid, diagnosis_year, issuedate, year) %>%
  mutate(count_year = diagnosis_year - 2) %>% 
  filter(year >= count_year) %>% # remove prescriptions that are more than 2 years prior to diagnosis
  select(-count_year) %>%
  group_by(year) %>% 
  summarise(n_exp = n_distinct(patid)) %>% 
  ungroup() %>%
  left_join(active_counts, by = "year") %>%
  mutate(last_smi_diag = "Overall")

# Standardised to 1,000 patients and add CIs
calculate_CI <- function(exposed, total, confidence_level = 0.95) {
  
  rate <- exposed / total * 1000 # Calculate rate per 1000
  SE <- sqrt(rate * (1000 - rate) / total) # Calculate standard error
  z <- qnorm((1 + confidence_level) / 2) # Calculate z value based on the confidence level
  MOE <- z * SE # Calculate margin of error
  lower_bound <- rate - MOE # Calculate CI
  upper_bound <- rate + MOE # Calculate CI
  
  return(paste(round(lower_bound, 2), ",", round(upper_bound, 2))) # Return CI}

# Add CI to data frame
overall_presc_rate_counts <- overall_presc_rate_counts %>%
  mutate(rate_obs_1000 = n_exp / count_active * 1000,
         CI = calculate_CI(n_exp, count_active)) %>%
  separate(CI, into = c("lower_ci", "upper_ci"), sep = " , ") %>%
  mutate_at(vars(5:7), as.numeric) %>%
  mutate_at(vars(5:7), ~ round(., digits = 0))

overall_presc_rate_counts %>%
  gt() %>%
  cols_hide(last_smi_diag) %>%
  cols_label( # tidy column labels
    n_exp = "Number prescribed AP",
    count_active = "Number active",
    year = "Year",
    rate_obs_1000 = "Rate/1000 (95% CI)") %>%
  cols_merge(columns = c("rate_obs_1000", "lower_ci", "upper_ci"), pattern = "{1} ({2}, {3})") %>%
  gtsave(paste0(output_dir, "overallAPprescribingrates_", today(), ".docx"))
  
proportion_overtime_overall <- ggplot(overall_presc_rate_counts, aes(x = factor(year))) +
  geom_bar(aes(y = count_active, fill = "Active patients"), stat = "identity", alpha = 0.5) +  # Total number
  geom_bar(aes(y = n_exp, fill = "Prescribed antipsychotics"), stat = "identity", alpha = 0.5) +  # Number exposed
  labs(x = "Year",
       y = "Number of patients") +
  scale_y_continuous(labels = scales::comma) +  # Add thousand separators for y-axis labels
  scale_fill_manual(values = c("Active patients" = "orange", "Prescribed antipsychotics" = "#32CD32"),
                    name = NULL) +  # Remove legend title
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")  # Position legend under the plot

ggsave(proportion_overtime_overall,
       path = output_dir,
       filename = paste0("proportion_overtime_overall_", file_date, file_type),
       dpi = 300, width = 8, height = 9, bg = "white")

# Individual antipsychotics ##### 
# prescribing rates per 1000 (with active SMI patients as denominator) graphs

# Remove very infrequently prescribed antipsychotics for individual trends
top_APs <- antipsychotics_sample %>%
  select(patid, AP) %>%
  group_by(AP) %>%
  summarize(num_patients = n_distinct(patid)) %>%
  arrange(desc(num_patients)) %>%
  slice(1:15)

top_APs <- unique(top_APs$AP)

trends_sample <- antipsychotics_sample %>%
  mutate(count_year = diagnosis_year - 2) %>%
  filter(year >= count_year) %>% 
  filter(AP %in% top_APs)

# Create DF for each SMI diagnosis
descriptivesample_schiz <- trends_sample %>%
  filter(last_smi_diag == "schizophrenia")

descriptivesample_bipolar <- trends_sample %>%
  filter(last_smi_diag == "bipolar")

descriptivesample_other <- trends_sample %>%
  filter(last_smi_diag == "other psychosis")

df_list <- list(trends_sample, descriptivesample_schiz, descriptivesample_bipolar, descriptivesample_other)
active_dfs <- list(active_counts, active_counts_schizophrenia, active_counts_bipolar, active_counts_otherpsychosis)

# Iterate over each data frame
individual_ap_rate_plot_list <- list()

colours <- c("black", "blue", "red", "#B0BF1A", "#915C83", "purple", "#FF9966", "#9E1B32","darkorange", "darkgreen", "#E68FAC", "magenta", "darkgrey", "#A67B5B", "#1DACD6")

for (i in seq_along(df_list)) {
  
  # Count the number of people for each medication and year
  df_count <- df_list[[i]] %>% 
    group_by(AP, year) %>% 
    summarise(n_exp = n_distinct(patid)) %>% 
    filter(year >= 2000 & year < 2020) %>% 
    ungroup()
  
  # Merge with numbers active in each year, calculate rate
  df_total <- df_count %>%
    left_join(active_dfs[[i]], by = "year") %>%
    mutate(rate = (n_exp / count_active) * 1000)
  
  # Create plot title
  plot_title <- paste("Rate of prescribing per 1,000 active SMI patients:")
  
  # Add the subtitle to the plot title
  plot_title <- paste(plot_title, subtitles[i])
  
  plot <- ggplot(df_total, aes(x = year, y = rate, group = AP, color = AP, label = AP)) +
    geom_line(aes(linetype = ifelse(AP == "Total", "solid", "dashed")), stat = "identity", size = 1.5, alpha = 0.9) +
    expand_limits(x = c(2000, NA)) +
    expand_limits(x = c(NA, 2020)) +
    coord_cartesian(clip = "off") +
    geom_label_repel(
      data = subset(df_total, year == max(year)),
      nudge_x = 1,
      nudge_y = 1,
      segment.curvature = -0.2,
      segment.ncp = 2,
      box.padding = 0.25,  # Adjust box padding to control label distance
      segment.size = 0.75,
      size = 4,
      aes(label = AP),
      min.segment.length = 0.6, 
      max.overlaps = Inf,
      arrow = arrow(length = unit(0.01, "npc")),
      point.padding = 0.9, 
    ) +
    labs(
      #    title = plot_title,
      x = "Year",
      y = "Prescribing rate per 1,000",
      color = "Medication"
    ) +
    scale_x_continuous(breaks = seq(2000, 2020, 2), expand = c(0, 0.4)) +
    scale_color_manual(values = setNames(colours, unique(df_total$AP))) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_y_continuous(labels = scales::comma_format(scale = 1/1))
  
  # Add plot to plot list
  individual_ap_rate_plot_list[[i]] <- plot}

# Save each plot as PNG
rates_active_file_names <- paste0(seq_along(individual_ap_rate_plot_list), "_individualaprates_", file_date, "_", gsub("[^A-Za-z0-9]", "", subtitles), file_type)
for (i in seq_along(individual_ap_rate_plot_list)) {
  file_path <- file.path(output_dir, rates_active_file_names[i])
  ggsave(file_path, individual_ap_rate_plot_list[[i]], width = 15, height = 8, dpi = 300)}


# FIRST- AND SECOND-LINE ANTIPSYCHOTICS ####

# Identify antipsychotics prescribed on the first prescription date (could be more than one antipsychotic per person)

# Order by issue date (oldest to newest)
incident_ap <- antipsychotics_sample %>%
  filter(first_ap_date >= '2000-01-01') %>%
  filter(first_ap_date == issuedate)

table_overall_incidentap <- incident_ap %>%
  select(AP) %>%
  tbl_summary(include = AP,
              digits = all_categorical() ~ c(0, 1))

table_2000to2009_incidentap <- incident_ap %>%
  filter(issuedate >= '2000-01-01' & issuedate < '2010-01-01') %>%
  select(AP) %>%
  tbl_summary(include = AP,
              digits = all_categorical() ~ c(0, 1))

table_2010plus_incidentap <- incident_ap %>%
  filter(issuedate >= '2010-01-01') %>%
  select(AP) %>%
  tbl_summary(include = AP,
              digits = all_categorical() ~ c(0, 1))

mergedtable_incidentap <- tbl_merge(
  tbls = list(table_overall_incidentap, table_2000to2009_incidentap, table_2010plus_incidentap),
  tab_spanner = c("**Overall**", "**2000-2009**", "**2010+**"))

# Second-line antipsychotics

# Order by issue date (oldest to newest)
secondline_ap <- antipsychotics_sample %>%
  filter(first_ap_date >= '2000-01-01') %>%
  select(patid, AP, issuedate) %>%
  group_by(patid, AP) %>%
  mutate(min = min(issuedate)) %>%
  ungroup() %>%
  select(patid, AP, min) %>%
  distinct() %>%
  group_by(patid) %>%
  mutate(overallmin = min(min)) %>%
  ungroup() %>%
  filter(overallmin != min) %>%
  group_by(patid) %>%
  mutate(secondmin = min(min)) %>%
  ungroup() %>%
  filter(secondmin == min) %>%
  select(patid, AP, issuedate = secondmin)

table_overall_secondlineap <- secondline_ap %>%
  select(AP) %>%
  tbl_summary(include = AP,
              digits = all_categorical() ~ c(0, 1))

table_2000to2009_secondlineap <- secondline_ap %>%
  filter(issuedate >= '2000-01-01' & issuedate < '2010-01-01') %>%
  select(AP) %>%
  tbl_summary(include = AP,
              digits = all_categorical() ~ c(0, 1))

table_2010plus_secondlineap <- secondline_ap %>%
  filter(issuedate >= '2010-01-01') %>%
  select(AP) %>%
  tbl_summary(include = AP,
              digits = all_categorical() ~ c(0, 1))

mergedtable_secondlineap <- tbl_merge(
  tbls = list(table_overall_secondlineap, table_2000to2009_secondlineap, table_2010plus_secondlineap),
  tab_spanner = c("**Overall**", "**2000-2009**", "**2010+**"))

#First- and second-line (note that time periods are not shown in the output)
incidentandsecondlineap <- tbl_merge(tbls = list(mergedtable_incidentap, mergedtable_secondlineap),
                                     tab_spanner = c("**First-line**", "**Second-line**")) %>%
  as_gt() %>%
  gtsave(paste0(output_dir, "incidentandsecondlineap_", today(), ".docx"))

# Package versions
packageVersion("dplyr") # 1.1.4
packageVersion("gtsummary") # 1.7.2
packageVersion("gt") # 0.10.1
packageVersion("lubridate") # 1.9.2
packageVersion("ggplot2") # 3.4.4
packageVersion("scales") # 1.2.1
packageVersion("readr") # 2.1.5
packageVersion("reshape2") # 1.4.4
packageVersion("gridExtra") # 2.3
packageVersion("tibble") # 3.2.1
packageVersion("tidylog") # 1.0.2
packageVersion("tidyr") # 1.3.0
packageVersion("stringr") # 1.5.0
packageVersion("forcats") # 1.0.0
packageVersion("patchwork") # 1.1.3
packageVersion("viridis") # 0.6.4
packageVersion("ggrepel") # 0.9.5
