
transform_data <- function(sim_data){
  summary_data <- sim_data %>% group_by(sim_time) %>%
    summarise(median = median(value),
              lower = quantile(value, 0.025), # 0.025
              upper = quantile(value, 0.975)) # 0.0975
  
  return (summary_data)
}

# function to plot antigen prevalence
plot_antigen_prevalence <- function(data, plot_title="Antigen Prevalence"){
  #x-axis year labels
  years <- c(
    "1998", rep("", 2),
    "2001", rep("", 3), "2005", rep("", 4),
    "2010", rep("", 4), "2015", rep("", 4),
    "2020", rep("", 4), "2025", rep("", 4),
    "2030", rep("", 4), "2035", rep("", 4),
    "2040"
  )
  
  ggplot(data, aes(x = sim_time)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, color = "95% CI"), 
                fill = "lightblue", alpha = 0.8) +
    geom_line(aes(y = median, color = "Median prevalence"),
              linewidth = 0.5) +
    geom_hline(aes(yintercept = c(2), color = "2% threshold"), 
               lty="dashed") +
    labs(title = plot_title, x = "Year", y = "Prevalence (%)") +
    theme_tufte(base_size=16) + 
    geom_rangeframe(data=data.frame(x=c(1,505), y=c(0,50)), aes(x=x, y=y)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_text(face = "bold"), 
          axis.title.y = element_text(face = "bold"),
          axis.line = element_line(linewidth = 0.2),
          axis.ticks.length = unit(0.2, "cm"),
          plot.title = element_text(face = "bold", hjust = 0.5), 
          legend.position = c(0.8,0.9) ) +
    scale_x_continuous(breaks = seq(1, 505, 12), 
                       labels = years, 
                       expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0,30,5), 
                       labels = seq(0,30,5), 
                       expand = c(0,0), 
                       limits = c(0,30)) +
    scale_color_manual(name = "", 
                       breaks = c("Median prevalence","95% CI","2% threshold"),
                       values = c("Median prevalence" = "blue",
                                  "95% CI" = 'lightblue',
                                  "2% threshold" = "red"))
}

# function to plot mf prevalence
plot_mf_prevalence <- function(data, plot_title="Mf Prevalence"){
  #x-axis year labels
  years <- c(
    "1998", rep("", 2),
    "2001", rep("", 3), "2005", rep("", 4),
    "2010", rep("", 4), "2015", rep("", 4),
    "2020", rep("", 4), "2025", rep("", 4),
    "2030", rep("", 4), "2035", rep("", 4),
    "2040"
  )
  
  ggplot(data, aes(x = sim_time)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, color = "95% CI"), 
                fill = "lightblue", alpha=0.8) +
    geom_line(aes(y = median, color="Median prevalence"), 
              linewidth = 0.5) +
    geom_hline(aes(yintercept=c(1), 
                   color="1% threshold"), lty="dashed") +
    labs(title=plot_title, x="Year", y="Prevalence (%)") +
    theme_tufte(base_size = 16) +
    geom_rangeframe(data=data.frame(x=c(1,505), y=c(0,45)), aes(x=x, y=y)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          axis.line = element_line(linewidth = 0.2),
          axis.ticks.length=unit(0.2, "cm"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = c(0.8,0.9) ) +
    scale_x_continuous(breaks = seq(1, 506, 12), 
                       labels = years, 
                       expand = c(0,0) ) +
    scale_y_continuous(breaks = seq(0,35,5), 
                       labels = seq(0,35,5), 
                       expand = c(0,0), 
                       limits = c(0,35)) +
    scale_color_manual(name = "",
                       breaks = c("Median prevalence","95% CI","1% threshold"),
                       values = c("Median prevalence" = "blue",
                                  "95% CI" = 'lightblue',
                                  "1% threshold" = "red"))
}

plot_and_save_scenario <- function(scenario_data, output_file) {
  ag_data <- transform_data(scenario_data %>% rename(value = ag_prev))
  ag_plot <- plot_antigen_prevalence(ag_data)
  
  mf_data <- transform_data(scenario_data %>% rename(value = mf_prev))
  mf_plot <- plot_mf_prevalence(mf_data)
  
  combined_plot <- ggarrange(
    ag_plot, mf_plot, 
    ncol = 1, font.label = list(size = 16)
  )
  
  ggsave(
    filename = output_file, 
    plot = combined_plot, width = 6, height = 7, bg = "white"
  )
}


plot_both_prevalence <- function(data1, data2, 
                                    label1 = "Dataset 1", 
                                    label2 = "Dataset 2",
                                    plot_title="Antigen Prevalence"){
  
  # Add dataset labels
  data1$Dataset <- label1
  data2$Dataset <- label2
  
  # Combine datasets
  data_all <- rbind(data1, data2)
  
  # x-axis year labels
  years <- c(
    "1998", rep("", 2),
    "2001", rep("", 3), "2005", rep("", 4),
    "2010", rep("", 4), "2015", rep("", 4),
    "2020", rep("", 4), "2025", rep("", 4),
    "2030", rep("", 4), "2035", rep("", 4),
    "2040"
  )
  
  plot <- ggplot(data_all, aes(x = sim_time, group = Dataset)) +
    
    # Confidence intervals
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = Dataset), 
                alpha = 0.3) +
    
    # Median lines
    geom_line(aes(y = median, color = Dataset),
              linewidth = 0.7) +
    
    # Threshold line
    geom_hline(yintercept = 2, color = "red", linetype="dashed") +
    
    labs(title = plot_title, x = "Year", y = "Prevalence (%)") +
    
    theme_tufte(base_size=16) + 
    geom_rangeframe(data=data.frame(x=c(1,505), y=c(0,50)), aes(x=x, y=y)) +
    
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_text(face = "bold"), 
          axis.title.y = element_text(face = "bold"),
          axis.line = element_line(linewidth = 0.2),
          axis.ticks.length = unit(0.2, "cm"),
          plot.title = element_text(face = "bold", hjust = 0.5), 
          legend.position = c(0.8,0.9)) +
    
    scale_x_continuous(breaks = seq(1, 505, 12), 
                       labels = years, 
                       expand = c(0,0)) +
    
    scale_y_continuous(breaks = seq(0,35,5), 
                       labels = seq(0,35,5), 
                       expand = c(0,0), 
                       limits = c(0,35)) +
    
    scale_color_manual(values = c(label1 = "blue", label2 = "darkgreen")) +
    scale_fill_manual(values = c(label1 = "lightblue", label2 = "lightgreen"))
  
  print(plot)
}



