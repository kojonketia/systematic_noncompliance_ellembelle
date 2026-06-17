rm(list = ls())

pacman::p_load(dplyr, tidyr, ggh4x, stats, ggpubr, readxl, writexl, ggplot2,
               officer, ggthemes, ggridges, reshape2, graphics, gtsummary,
               flextable, hrbrthemes, purrr, binom, gt)

theme_gtsummary_compact()
set_flextable_defaults(table.layout = "autofit")
set_flextable_defaults(background.color = "white")

# set up -----------------------------------------------------------------------
study_data <- read_excel("study_MDA_data.xlsx")
df <- filter(study_data, Age > 4)

df$`Age group` <- cut(df$Age,
                    breaks = c(4, 10, 15, 20, 30, 40, 50, 60, Inf),
                    labels = c('5-10', '11-15', '16-20', '21-30', '31-40', 
                               '41-50', '51-60', '>60'),
                    right = TRUE)
df$`Level of compliance` <- factor(df$`Level of compliance`, 
                                   levels=c("Systematic compliant",
                                            "Systematic non-compliant",
                                            "Random compliant"))

# summary tables ---------------------------------------------------------------
summary_counts <- function(columns, label_list = NULL){
  s_tab <- tbl_summary(df |> select(columns),
                       include = columns,
                       missing = "no",
                       percent = "column",
                       type = all_dichotomous() ~ "categorical",
                       statistic = list(all_categorical() ~ "{n}"),
                       digits = list(all_categorical() ~ c(0, 1)),
                       label = label_list) |>
    modify_header(label ~ "**Variable**") |> bold_labels()
  return (s_tab)
}

summary_percent <- function(columns, label_list = NULL){
  p_tab <- tbl_summary(df |> select(columns),
                       include = columns,
                       missing = "no",
                       percent = "column",
                       type = all_dichotomous() ~ "categorical",
                       digits = list(all_categorical() ~ 1),
                       statistic = list(all_categorical() ~ "{p}"),
                       label = label_list) |> 
    bold_labels()
  return(p_tab)
}

summary_table <- tbl_merge(
  tbls = list(summary_percent(c("2019 A"), list(`2019 A` = "Variable")), 
              summary_percent(c("2021 A"), list(`2021 A` = "Variable")),
              summary_percent(c("2022 A"), list(`2022 A` = "Variable"))),
  tab_spanner = c("**2019**", "**2021**",  "**2022**"))
summary_table

characteristics_counts <- summary_counts(c("Age group", "Sex"))
characteristics_percentage <- summary_percent(c("Age group", "Sex"))

characteristics_table <- tbl_merge(
  tbls = list(characteristics_counts, characteristics_percentage),
  tab_spanner = c("**Frequency**", "**Percentage**"))
characteristics_table

# frequency tables -------------------------------------------------------------
freq_table <- function(columns, by_var, include_var){
  f_tab <- tbl_summary(df |> select(columns),
                       by = by_var,
                       include = include_var,
                       missing = "no",
                       percent = "row",
                       type = all_dichotomous() ~ "categorical",
                       digits = list(all_categorical() ~ c(0, 1))) |>
    add_p(simulate.p.value=TRUE) |> bold_labels() |> bold_p() |>
    #add_overall(col_label = "**Total**, N = {style_number(N)}") |>
    modify_header(label ~ "**Variable**") |>
    modify_footnote(everything() ~ NA)
  return(f_tab)
}

# demographics by different levels of compliance
dc_table <- freq_table(c("Age group", "Sex", "Level of compliance"),
                       c("Level of compliance"),
                       c("Age group", "Sex") )
dc_table

# participation in MDA ---------------------------------------------------------
participation_data <- df %>%
  select(`Only MDA 1`, `Only MDA 2`, `Only MDA 3`,
         `Only MDA 1 and 2`, `Only MDA 1 and 3`,
         `Only MDA 2 and 3`, `MDA 1, 2 and 3`, Never)

results <- map_dfr(names(participation_data), function(col) {
  col_values <- na.omit(participation_data[[col]])
  total <- length(col_values)
  
  freq_counts <- table(col_values)
  
  map_dfr(names(freq_counts), function(value) {
    count <- as.numeric(freq_counts[[value]])
    ci <- binom.confint(count, total, conf.level = 0.95, methods = "wilson")
    
    tibble(
      `Participated in MDA` = col,
      Value = value,
      `Frequency (n = 1957)` = count,
      `Percentage % (95% CI)` = sprintf("%.1f (%.1f-%.1f)",
                                        100 * count / total,
                                        100 * ci$lower,
                                        100 * ci$upper)) })
})

participation_frequency_data <- results %>%
  filter(Value == "Yes") %>%
  select(-Value)

participation_frequency_data

participation_frequency_data |>
  gt() |>
  tab_header(
    title = "MDA participation") 

pft <- flextable(participation_frequency_data)

# Correlation between MDA coverage and rounds ----------------------------------
df$mda1 <- ifelse(df$`Only MDA 1` == "Yes", 1, 0)
df$mda2 <- ifelse(df$`Only MDA 2` == "Yes", 1, 0)
df$mda3 <- ifelse(df$`Only MDA 3` == "Yes", 1, 0)

print(paste0("Correlation between MDA 1 and 2: ", cor(df$mda1, df$mda2)))
print(paste0("Correlation between MDA 1 and 3: ", cor(df$mda1, df$mda3)))
print(paste0("Correlation between MDA 2 and 3: ", cor(df$mda2, df$mda3)))

attendance <- df %>% select(mda1, mda2, mda3)

estimate_rho_transitions <- function(attendance) {
  attendance <- as.matrix(attendance)
  
  # overall attendance probability
  p <- mean(attendance)
  
  # transitions
  prev <- c(attendance[,1], attendance[,2])
  next_round <- c(attendance[,2], attendance[,3])
  
  P11 <- mean(next_round[prev == 1])  # P(1 | 1)
  P10 <- mean(next_round[prev == 0])  # P(1 | 0)
  
  rho1 <- (P11 - p) / (1 - p)
  rho2 <- (p - P10) / p
  
  rho <- mean(c(rho1, rho2), na.rm = TRUE)
  return(rho)
}

estimate_rho_transitions(attendance)

# plots ------------------------------------------------------------------------
drug_taken_data <- df %>%
  pivot_longer(
    cols = c(`Only MDA 1`, `Only MDA 2`, `Only MDA 3`),
    names_to = "MDA",
    values_to = "Status"
  ) %>%
  mutate(
    MDA = recode(MDA,
                 "Only MDA 1" = "2019",
                 "Only MDA 2" = "2021",
                 "Only MDA 3" = "2022"),
    Status = recode(Status,
                    "Yes" = "Taken",
                    "No"  = "Not Taken")
  ) %>%
  group_by(MDA, Status) %>%
  summarise(Value = n(), .groups = "drop") %>%
  group_by(MDA) %>%
  mutate(
    Percentage = (Value / sum(Value)) * 100
  ) %>%
  ungroup() %>%
  select(MDA, Value, Percentage, Status)

drug_taken_data

drug_taken_data$Status <- factor(drug_taken_data
                                 $Status, levels = c("Taken", "Not Taken"))
drug_taken_data$MDA <- as.numeric(drug_taken_data$MDA)
status_colors <- c("Taken" = "gray30", "Not Taken" = "lightgray")

dtp <- ggplot(drug_taken_data, aes(x = MDA, y = Percentage, fill = Status)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.6), 
           color = "black",
           width = 0.5) +
  scale_fill_manual(values = status_colors) +
  labs(x = "Year", y = "Percentage (%)") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.title = element_blank(),
    axis.ticks.length = unit(0.25, "cm"),  
    axis.ticks = element_line(color = "black")) +
  geom_label(data = data.frame(
    Year = c(2020),
    Mean = c(50),
    label = c("No MDA")),
    aes(x = Year, y = Mean, label = label),
    inherit.aes = FALSE,
    size = 3,
    hjust = 0.5) +
  guides(y = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(breaks = seq(0, 100, 20),           
                     minor_breaks = seq(0, 100, 10),
                     labels = c("0","20","40","60", "80","100"),
                     limits = c(0,100), expand = c(0,0)) 
dtp
ggsave("dt.png", plot = dtp, width = 5, height = 3)


# saving to word ---------------------------------------------------------------
doc <- officer::read_docx()
doc <- body_add_par(doc, value = "")
doc <- body_add_flextable(doc, as_flex_table(summary_table))
doc <- body_add_par(doc, value = "")
doc <- body_add_flextable(doc, as_flex_table(characteristics_table))
doc <- body_add_par(doc, value = "")
doc <- body_add_flextable(doc, as_flex_table(dc_table))
doc <- body_add_par(doc, value = "")
doc <- body_add_flextable(doc, pft)
doc <- body_add_par(doc, value = "")
print(doc, "Results_updated.docx")
