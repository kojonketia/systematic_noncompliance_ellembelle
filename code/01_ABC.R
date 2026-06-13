rm(list = ls())

pacman::p_load(Rcpp, parallel, tidyverse, abc, dplyr, Hmisc)

sourceCpp(file.path("code", "Transfil.cpp"))   # The TRANSFIL model
set.seed(123)

start.time <- Sys.time()

prev_df <- read.csv(file.path("data", "cleaned", "prevalence.csv"))
prevalence_obs <- prev_df$prevalence[prev_df$year == 2001]


simulate_baseline_prevalence <- function(k, vth) { 
	pop_size <- 1000
	x <- new(Data, pop_size)
	x$set_parameter("k", k)
	x$set_parameter("vth", vth)  
	x$set_parameter("importation_rate", 0.0003)
	x$burn_in()
	return(x$worm_prev)
}

# Simulation function
simulate_results <- function(i) {
	k_prior <- runif(1, 0.001, 0.4)
  vth_prior <- runif(1, 30, 150)

	sim_data <- simulate_baseline_prevalence(k = k_prior, vth = vth_prior)   

	# combine generated k parameter and simulated prevalence 
	out <- data.frame(
		k = k_prior, 
		vth = vth_prior,
		prev = sim_data
	)
	return(out)
}


# --- ABC Setup ---
N <- 200000 # number of runs
n_cores <- 4

output_file <- file.path("output", "simulated_baseline_prevalence.csv")

if (file.exists(output_file)) {
  print("Output file already exists. Reading existing results...")
  results <- read.csv(output_file)
} else {
  print(sprintf("Running ABC with %d simulations using %d cores...", N, n_cores))
  
  cl <- makeCluster(n_cores)
  
  clusterEvalQ(cl, {
    library(Rcpp)
    Rcpp::sourceCpp(file.path("code", "Transfil.cpp"))
  })
  
  # Export needed objects/functions to workers
  clusterExport(cl, varlist = c("simulate_results", 
                                "simulate_baseline_prevalence",
                                "Data"), 
                envir = environment())
  
  results_list <- parLapply(cl, 1:N, simulate_results)
  stopCluster(cl)
  
  results <- dplyr::bind_rows(results_list)
  
  write.csv(results, output_file, row.names = FALSE)
}

tol <- 0
posterior_df <- results %>% 
	filter(abs(prev - prevalence_obs) <= tol)


posterior_path <- file.path("output", "posterior_samples.csv")
write.csv(posterior_df, posterior_path, row.names = FALSE)


postk_vs_prev <- posterior_df %>% 
	ggplot(aes(x = k, y = prev)) + 
	geom_point() +
	geom_hline(yintercept = prevalence_obs, linetype = "dashed", color = "blue") +
	labs(title = "Posterior Samples of k Parameter",
	     x = "k Parameter", y = "Simulated Baseline Prevalence") +
	theme_minimal()

postvth_vs_prev <- posterior_df %>% 
	ggplot(aes(x = vth, y = prev)) + 
	geom_point() +
	geom_hline(yintercept = prevalence_obs, linetype = "dashed", color = "blue") +
	labs(title = "Posterior Samples of k Parameter",
	     x = "k Parameter", y = "Simulated Baseline Prevalence") +
	theme_minimal()



## Simulate from the posterior predictive
pp_df <- Map(function(k, vth) {
  simulated_prev <- simulate_baseline_prevalence(k, vth)
  
  data.frame(
    k = k, 
    vth = vth,
    prevalence = simulated_prev
  )
}, posterior_df$k, posterior_df$vth) %>% 
  dplyr::bind_rows()


ppd_plot <- pp_df  %>% 
	ggplot(aes(x = prevalence)) + 
	geom_histogram() +
	geom_vline(xintercept = 18, linetype = "dashed", color = "blue") +
	labs(title = "Posterior Predictive distribution",
		x = "Prevalence") +
	theme_minimal()

figures_path <- file.path("output", "figures")
if (!dir.exists(figures_path)) {
	dir.create(figures_path, recursive = TRUE)
}

ggsave(
	filename = file.path(figures_path, "posterior_predictive_plot.png"), 
	plot = ppd_plot, width = 7, height = 5, dpi = 300
)


# Extract prior and posterior samples to plot
prior_density <- density(results$k)
posterior_density <- density(posterior_df$k)
posterior_mean <- mean(posterior_df$k)
posterior_ci <- smean.cl.normal(posterior_df$k) 


prior_df <- data.frame(
	k = seq(0.001, 0.3, length.out = 1000),
	prior_density = dunif(seq(0.001, 0.3, length.out = 1000), min = 0.001, max = 0.3)
)

posterior_and_prior_plot <- ggplot() +
	geom_density(data = posterior_df, aes(x = k, y = ..density..), fill = "red", alpha = 0.4, color = "red", size = 1) +
	geom_line(data = prior_df, aes(x = k, y = prior_density), color = "blue", size = 1) +
	geom_area(data = prior_df, aes(x = k, y = prior_density), fill = "blue", alpha = 0.3) +
	geom_vline(xintercept = posterior_mean, color = "red", linetype = "dashed", size = 1) +
	annotate("text", x = posterior_mean, y = max(density(posterior_df$k)$y) * 0.9,
			 label = sprintf("Mean: %.3f, 95%% CI: (%.4f, %.4f)", 
							 posterior_mean, quantile(posterior_df$k, 0.025), quantile(posterior_df$k, 0.975)),
			 color = "red", hjust = -1, vjust = 1, size = 4) +
	labs(title = "Posterior and Prior Density of k",
			 x = "k",
			 y = "Density") +
	theme_minimal() +
	scale_x_continuous(limits = c(0, 0.3))


ggsave(
	filename = file.path(figures_path, "prior_and_posterior.png"),
	plot = posterior_and_prior_plot, width = 7, height = 5, dpi = 300
)


end.time <- Sys.time()
cat("Simulation completed in:", 
    round(difftime(end.time, start.time, units = "hours"), 2), "hours\n")



