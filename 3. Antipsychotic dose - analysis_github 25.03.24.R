# ------------------------
# "Prescribing of antipsychotics for people diagnosed with severe mental illness in UK primary care"

# 3. Antipsychotic dose - analysis_github 25.03.24.R ####
# This script covers analyses of antipsychotic dose over the first year of prescribing
# The dataset used in these analyses is previously created using the "Antipsychotic dose - data processing" script
# ------------------------
# Last run: 25/03/2024

# Clear memory
# rm(list = ls())

# Packages
library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(gtsummary)
library(forcats)
library(gridExtra)
library(lubridate)
library(tidyr)
library(readxl)
library(stringr)
library(ggpubr)
library(grid)
library(gt)
library(gtable)
library(patchwork)
library(tidylog)

# Set file path
path <- # removed
  
# Set working directory
# removed
  
# Load data
load("Antipsychotic_dose.Rdata") # as created by the antipsychotic dose data processing script

# FUNCTIONS

# Function to calculate summary statistics for graphing (can be used on grouped and ungrouped data)
calculate_dose_data <- function(dose_data, group_var = NULL) {
  
  if (!is.null(group_var)) {
    # if a grouping variable is provided, split the data by group and calculate the summary statistics for each group
    dose_data <- split(dose_data, dose_data[[group_var]])
    summary <- lapply(dose_data, function(sub_data) {
      calculate_dose_data(sub_data, group_var = NULL)})
    summary <- do.call(rbind, summary)
    summary$Group <- rownames(summary)
    rownames(summary) <- NULL
    
    # Remove characters after and including the first '.' character in the group column
    summary$Group <- gsub("\\..*", "", summary$Group)
    
  } else {
    # create an empty data frame to store results
    summary <- data.frame(dose = factor(character(), levels = paste0("dose", sprintf("%02d", 1:12))),
                          mean = numeric(),
                          sd = numeric(),
                          n_obs = numeric(),
                          lower_ci = numeric(),
                          upper_ci = numeric())
    
    # loop over each dose column
    for (i in 1:12) {
      
      # get column name and values
      dose_col <- paste0("dose", sprintf("%02d", i))
      dose_values <- dose_data[[dose_col]]
      
      # check that dose_values is a numeric/integer
      if (!is.numeric(dose_values) && !is.integer(dose_values)) {
        stop(paste0("Column ", dose_col, " is not numeric or integer"))}
      
      # convert dose_values to numeric if it's a character column
      if (is.character(dose_values)) {
        dose_values <- as.numeric(dose_values)}
      
      # calculate summary statistics
      n_obs <- sum(!is.na(dose_values))
      mean_val <- mean(dose_values, na.rm = TRUE)
      sd_val <- sd(dose_values, na.rm = TRUE)
      se_val <- sd_val / sqrt(n_obs)
      z_val <- qnorm(0.975)
      ci_val <- z_val * se_val
      lower_ci <- mean_val - ci_val
      upper_ci <- mean_val + ci_val
      
      # add results to summary data frame
      summary[i, ] <- list(dose_col, mean_val, sd_val, n_obs, lower_ci, upper_ci)}}
  
  return(summary)}

# Create a function to plot each summary data frame (ungrouped)
plot_summary_ungrouped <- function(summary, title) {

  # Create the plot
  u <- ggplot(summary, aes(x = dose, y = mean)) +
    geom_line(color = "blue", linewidth = 1, group = 1) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.3, fill = "lightblue", group = 1) +
    labs(x = "Prescription number", y = "Dose (mg)") +
    ggtitle("Mean by Dose") +
    coord_cartesian(ylim=olanz_yaxis) +
    scale_x_discrete(labels=xlabels) +
    ggtitle(paste("Olanzapine equivalent dose", title)) +
    theme_classic() +
    theme(axis.line = element_line(color = "grey", linewidth = 1), # format axes
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom", # move legend below the plot
        legend.title = element_blank()) # remove legend title

  # Return the plot
  return(u)}

# Create a function to plot each summary data frame (grouped)
plot_summary <- function(summary, title) {
  
  # Create the plot
  p <- ggplot(summary, aes(x = dose, y = mean, group = Group, color = Group)) +
    geom_line(linetype = "solid", linewidth = 1.5) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.3, fill = "lightblue", color = NA) +
    labs(x = "Prescription number", y = "Dose (mg)") +
    theme_classic() +
    coord_cartesian(ylim=olanz_yaxis) +
    scale_x_discrete(labels=xlabels) +
    ggtitle(paste("Olanzapine equivalent dose", title)) +
    theme(axis.line = element_line(color = "grey", linewidth = 1), # format axes
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom", # move legend below the plot
          legend.title = element_blank()) # remove legend title

  # Return the plot
  return(p)}

# Create a function to create the dose denominator table
create_dose_table <- function(summary, title) {
  # select relevant columns
  if ("Group" %in% colnames(summary)) {
    dose_table <- summary %>% 
      select(dose, n_obs, mean, sd, Group)
  } else {
    dose_table <- summary %>% 
      select(dose, n_obs, sd, mean)}
  
  dose_table <- dose_table %>%
    filter(grepl("01|06|12", dose))
  
  # format n column
  dose_table$n <- paste(sprintf("%.1f", dose_table$mean), " (", sprintf("%.1f", dose_table$sd), ") [n=", format(dose_table$n_obs, big.mark = ","), "]", sep = "")
  dose_table$sd <- NULL
  dose_table$mean <- NULL
  dose_table$n_obs <- NULL
  
  # format dose column
  dose_table$dose <- str_replace(dose_table$dose, "^dose0", "Dose ")
  dose_table$dose <- str_replace(dose_table$dose, "^dose", "Dose ")
  
  # pivot wider
  dose_table_wide <- tidyr::pivot_wider(dose_table, names_from = dose, values_from = "n") %>%
    tibble::rownames_to_column(var = "n") %>% # rename row name from 1 to n
    dplyr::mutate(n = if_else(row_number() == 1, "n", n)) %>%
    tibble::column_to_rownames(var = "n")
  
  # return table
  return(dose_table_wide)}

# Data management
Antipsychotic_dose <- Antipsychotic_dose %>% 
  mutate(last_smi_diag = case_when(last_smi_diag == "bipolar" ~ "Bipolar disorder",
                                   last_smi_diag == "other psychosis" ~ "Other non-organic psychoses",
                                            TRUE ~ last_smi_diag),
         last_smi_diag = str_to_title(last_smi_diag), # capitalisation
         ethnicity_cat_cprdhes = fct_na_value_to_level(ethnicity_cat_cprdhes, level = "Unknown"), #Set missing ethnicities to own category
         ethnicity_cat_cprdhes = str_to_title(ethnicity_cat_cprdhes),
         apinitiation_timeperiod = case_when(apinitiation_timeperiod == "<2005" ~ "2000-2004",
                                             apinitiation_timeperiod == "2015+" ~ "2015-2019",
                                             TRUE ~ apinitiation_timeperiod)) # capitalisation after "-"
  
# CALCULATE SUMMARY DATA TO BE GRAPHED
# Run the calculate_dose_data function on the data frames (ungrouped and by grouping variable)

# Overall
summary <- calculate_dose_data(Antipsychotic_dose)
summary_diagnosis <- calculate_dose_data(Antipsychotic_dose, group_var = "last_smi_diag")
summary_ethnicity <- calculate_dose_data(Antipsychotic_dose, group_var = "ethnicity_cat_cprdhes")
summary_gender <- calculate_dose_data(Antipsychotic_dose, group_var = "gender") %>% filter(Group == "Female" | Group == "Male")
summary_age <- calculate_dose_data(Antipsychotic_dose, group_var = "age_atfirstap_cat")
summary_timeperiod <- calculate_dose_data(Antipsychotic_dose, group_var = "apinitiation_timeperiod")
summary_region <- calculate_dose_data(Antipsychotic_dose, group_var = "region")
summary_deprivation <- calculate_dose_data(Antipsychotic_dose, group_var = "patprac_2019imd_quintile")

# PLOTS

# Settings

#Define y axes 
olanz_yaxis <- c(4.0, 11.5) # for dose mg graphs

#Define x axis labels
xlabels <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12')

#Ungrouped

#Define grid names, for naming of tables and plots
ungrouped_grid_names <- c("Overall")

# Create a list of ungrouped summary data frames
summary_ungrouped_list <- list(summary)

# Use purrr to iterate over each summary data frame and plot it with specified title
ungroupedplots_list <- map2(summary_ungrouped_list,  c(""), plot_summary_ungrouped)

#Grouped

#Define grid names, for naming of tables and plots
grid_names <- c("diagnosis", "ethnicity", "gender", "age", "timeperiod", "region", "imd")

# Create list of grouped summary data frames
summary_list <- list(summary_diagnosis, summary_ethnicity, summary_gender, summary_age, summary_timeperiod, summary_region, summary_deprivation)

# Use purrr to iterate over each summary data frame and plot it with specified title
plots_list <- map2(summary_list, c("by severe mental illness diagnosis", "by ethnicity", "by sex", "by age", "by time period", "by region", "by Index of Multiple Deprivation"), plot_summary)

# loop through plots and dynamically name each plot
for (i in seq_along(plots_list)) {
  assign(paste0("plot_by_", grid_names[i]), plots_list[[i]])}

# TABLES OF DOSES/DENOMINATORS

#Ungrouped

# Use purrr to iterate over each ungrouped summary data frame (from summary_ungrouped_list) and create a dose table with specified title
ungrouped_dose_table_list <- map2(summary_ungrouped_list, ungrouped_grid_names, create_dose_table)

#Grouped

# Use purrr to iterate over each summary data frame (from summary_list) and create a dose table with specified title
dose_table_list <- map2(summary_list, grid_names, create_dose_table)

# PLOT GRAPH AND TABLE

# Combine all lists
all_plots <- c(plots_list, ungroupedplots_list) # plots
all_dose_tables <- c(dose_table_list, ungrouped_dose_table_list) # tables
all_grid_names <- c(grid_names, ungrouped_grid_names) # grid names

for (i in seq_along(all_plots)) {
  
  # Convert plot and table to grob objects
  p_grob <- ggplotGrob(all_plots[[i]])
  table_grob <- tableGrob(all_dose_tables[[i]], theme = ttheme_minimal(), rows = NULL)
  
  # Remove cell background color
  for (j in seq_len(nrow(all_dose_tables[[i]]))) {
    for (k in seq_len(ncol(all_dose_tables[[i]]))) {
      table_grob <- editGrob(table_grob,
                             gPath(sprintf("cell-%d-%d", j, k)),
                             gp = gpar(fill = NA, col = "black", lwd = 1))}}
  
  # Combine plot and table using grid.arrange
  grid <- grid.arrange(p_grob, table_grob, layout_matrix = rbind(c(1,1), c(2,2)), 
                       heights = unit.c(unit(0.8, "npc"), unit(0.2, "npc")), 
                       widths = unit(c(0.8, 0.2), c("npc", "npc")))
  
  # Save each grid as a separate object with a unique name
  assign(paste0(all_grid_names[i]), grid)}

# SAVE TO FILE

# Specify output directory
output_dir <- # removed

# Iterate over each named grid
for (i in seq_along(all_grid_names)) {
  
  # Create the file path and name
  file_name <- paste0(all_grid_names[[i]], "_", format(Sys.Date(), "%Y-%m-%d"), ".png")
  file_path <- file.path(output_dir, file_name)
  
  # Get the grid object
  my_grid <- eval(parse(text = all_grid_names[[i]]))
  
  # Save the grid as a PNG file with high dpi
  ggsave(file_path, my_grid, dpi = 300, width = 12, height = 10)}

# Combined SMI diagnosis and ethnicity plot

# Convert plot and table to grob objects
p_grob <- ggplotGrob(all_plots[[1]])
table_grob <- tableGrob(all_dose_tables[[1]], theme = ttheme_minimal(), rows = NULL)
p_grob2 <- ggplotGrob(all_plots[[2]])
table_grob2 <- tableGrob(all_dose_tables[[2]], theme = ttheme_minimal(), rows = NULL)

# Remove cell background color
for (j in seq_len(nrow(all_dose_tables[[1]]))) {
  for (k in seq_len(ncol(all_dose_tables[[1]]))) {
    table_grob <- editGrob(table_grob,
                           gPath(sprintf("cell-%d-%d", j, k)),
                           gp = gpar(fill = NA, col = "black", lwd = 1))}}

for (j in seq_len(nrow(all_dose_tables[[2]]))) {
  for (k in seq_len(ncol(all_dose_tables[[2]]))) {
    table_grob2 <- editGrob(table_grob2,
                            gPath(sprintf("cell-%d-%d", j, k)),
                            gp = gpar(fill = NA, col = "black", lwd = 1))}}

# Combine plot and table using grid.arrange
combined_diagnosis_ethnicity_grid <- grid.arrange(p_grob, p_grob2, table_grob, table_grob2, nrow = 2, heights = c(3.5, 1.5))
text <- "combined_ethnicity_diagnosis"
file_name <- paste0(text, "_", format(Sys.Date(), "%Y-%m-%d"), ".png")
file_path <- file.path(output_dir, file_name)

ggsave(file_path, combined_diagnosis_ethnicity_grid, dpi = 300, width = 18, height = 10)

# TABLES

#Table of differences between dose time points
timediff <- Antipsychotic_dose %>%
  ungroup() %>%
  select(last_smi_diag, dose01_date, dose02_date, dose03_date, dose04_date, dose05_date, dose06_date, 
         dose07_date, dose08_date, dose09_date, dose10_date, dose11_date, dose12_date) %>%
  mutate(dose1to2_datediff = dose02_date - dose01_date, 
         dose2to3_datediff = dose03_date - dose02_date,
         dose3to4_datediff = dose04_date - dose03_date,
         dose4to5_datediff = dose05_date - dose04_date,
         dose5to6_datediff = dose06_date - dose05_date,
         dose6to7_datediff = dose07_date - dose06_date,
         dose7to8_datediff = dose08_date - dose07_date,
         dose8to9_datediff = dose09_date - dose08_date,
         dose9to10_datediff = dose10_date - dose09_date,
         dose10to11_datediff = dose11_date - dose10_date,
         dose11to12_datediff = dose12_date - dose11_date) %>%
  select(dose1to2_datediff, dose2to3_datediff, dose3to4_datediff, dose4to5_datediff, dose5to6_datediff, dose6to7_datediff, dose7to8_datediff, 
         dose8to9_datediff, dose9to10_datediff, dose10to11_datediff, dose11to12_datediff) %>%
  tbl_summary(
    type = all_continuous() ~ "continuous2",
    statistic = all_continuous() ~ c("{median} ({p25}, {p75})","{mean} ({sd})", "{min} - {max}")) %>%
  as_gt() %>%
  gtsave(paste0(output_dir, "timedifferencebetweendates_indoseanalysis_", today(), ".docx"))

# Package versions
packageVersion("dplyr") # 1.1.4
packageVersion("readr") # 2.1.5
packageVersion("purrr") # 1.0.2
packageVersion("ggplot2") # 3.4.4
packageVersion("gtsummary") # 1.7.2
packageVersion("forcats") # 1.0.0
packageVersion("gridExtra") # 2.3
packageVersion("lubridate") # 1.9.2
packageVersion("tidyr") # 1.3.0
packageVersion("readxl") # 1.4.3
packageVersion("stringr") # 1.5.0
packageVersion("ggpubr") # 0.6.0
packageVersion("grid") # 4.3.1
packageVersion("gt") # 0.10.1
packageVersion("gtable") # 0.3.4
packageVersion("patchwork") # 1.1.3
packageVersion("tidylog") # 1.0.2
