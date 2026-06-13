pacman::p_load(Rcpp, parallel, tidyverse, dplyr,
               readxl, ggplot2, ggpubr, ggthemes)

start.time <- Sys.time()

## MDA offset
convert_year_to_time <- function(year, start_year = 2001) {
  # Convert year to time in months since the start of the simulation
  return((year - start_year) * 12)
}

sourceCpp(file.path("code", "Transfil.cpp")) # the TRANSFIL model

posterior_df <- read.csv(file.path("output", "posterior_samples.csv"))
posterior_k <- posterior_df$k
posterior_vth <- posterior_df$vth

n_runs <- 100  # number of simulation runs

int_df <- read.csv(file.path("data", "cleaned", "interventions.csv")) %>%
  mutate(time = convert_year_to_time(year))

set.seed(123)  # For reproducibility
posterior_ks <- sample(posterior_k, n_runs, replace = FALSE)
posterior_vths <- sample(posterior_vth, n_runs, replace = FALSE)

output_path <- file.path("output", "scenarios")
if(!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

simulate_scenario <- function(k, vth, sys_comp, int_df) {
  mda_settings <- list("time" = int_df$time + 36, ## lf prevalence run for 3 years
                       "coverage" = int_df$MDA_coverage )

  bednet_settings <- list("time" = int_df$time + 36,
                          "coverage" = int_df$bednet_coverage)

  x <- new(Data, 10000)
  x$set_parameter("k", k)
  x$set_parameter("vth", vth)  
  x$set_parameter("systematic_compliance", sys_comp)
  x$set_parameter("importation_rate", 0.0003)

  # Burn-in period
  x$burn_in()
  
  # Set MDA and bednet coverage
  x$set_MDA(mda_settings)
  x$set_bednet(bednet_settings)
  
  # Simulate through to 2040
  x$timestep(convert_year_to_time(2040) + 36) 

  out <- data.frame(
    sim_time = 0:(length(x$mf_prev)-1), # 2001 to 2040 in months
    mf_prev = x$mf_prev,
    ag_prev = x$worm_prev
  )
  
  return(out)
}


# interventions data
forward_interventions <- data.frame(year = 2026:2040, 
                                    MDA_coverage = c(rep(0,15)),
                                    bednet_coverage = rep(0.64, 15)) %>% 
  mutate(time = convert_year_to_time(year))

all_int_df <- int_df %>%
  bind_rows(forward_interventions)

# ##############################################################################
# Modeliing with 65% compliance rate
# ##############################################################################


results_65 <- lapply(1:n_runs, function(i) {
  out <- simulate_scenario(k = posterior_ks[i], vth = posterior_vths[i], 
                           sys_comp = 0.35, all_int_df)
  out$sim_id <- i
  return(out)}) %>% bind_rows()

write.csv(results_65, file.path(output_path, "scenario_65.csv"), row.names = FALSE)

# ##############################################################################
# Modeliing with 26% compliance rate
# ##############################################################################

results_26 <- lapply(1:n_runs, function(i) {
  out <- simulate_scenario(k = posterior_ks[i], vth = posterior_vths[i],
                           sys_comp = 0.74, all_int_df)
  out$sim_id <- i
  return(out)}) %>% bind_rows()

write.csv(results_26, file.path(output_path, "scenario_26.csv"), row.names = FALSE)

end.time <- Sys.time()
cat("Simulation completed in:", 
    round(difftime(end.time, start.time, units = "hours"), 2), "hours\n")