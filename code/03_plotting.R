pacman::p_load(tidyverse, ggthemes, ggpubr)

source("code/functions.R")

# list of scenario names: case 1
scenarios <- c("26", "65")
for (s in scenarios) {
  input_file  <- file.path("output", "scenarios", paste0("scenario_", s, ".csv"))
  output_file <- file.path("output", "figures", paste0("scenario_", s, "_plots.png"))
  
  results <- read.csv(input_file)
  plot_and_save_scenario(results, output_file)
}

scenarios <- c("26", "65")
for (s in scenarios) {
  input_file  <- file.path("output", "scenarios", paste0("scenario_", s, ".csv"))
  output_file <- file.path("output", "figures", paste0("ag_", s, "_plots.png"))
  
  results <- read.csv(input_file)
  
  ag_data <- transform_data(results %>% rename(value = ag_prev))
  ag_plot <- plot_antigen_prevalence(ag_data, NULL)
  ggsave(filename = output_file, 
    plot = ag_plot, width = 6, height = 3.5, bg = "white" )
}