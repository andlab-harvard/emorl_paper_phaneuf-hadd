####################################################
### Analyze EmoRL Data with Computational Models ###
####################################################

# Written by: Camille Phaneuf-Hadd (cphaneuf@g.harvard.edu)
# Last updated: 4/23/26

# Inputs: processed learning task data
# Computes: visualizations and statistical assessments
# Outputs: png plots and text files into results/learn

### Set up Script -----

# Load needed libraries
require(pacman) # for p_load()
p_load(tidyverse, # for df manipulation
       dplyr, # for %>% and other operators
       plotrix, # for std.error()
       ggplot2, # for plotting
       sjPlot, # for plot_model()
       performance, # for check_predictions()
       GGally, # for ggcorr()
       ggeffects, # for ggpredict()
       glmmTMB, # for glmmTMB()
       DHARMa) # for simulateResiduals()

# Load shared EmoRL functions
source("utilities.R")

# Set paths
in_path <- '../data/'
out_path <- '../results/learn/'

# Read in model comparison (only files requiring further analysis)
aic_comp_long <- read.csv(paste0(in_path, "model_comp/aic_comp_long.csv"))
bic_comp_long <- read.csv(paste0(in_path, "model_comp/bic_comp_long.csv"))
pseudor_comp_long <- read.csv(paste0(in_path, "model_comp/pseudor_comp_long.csv"))

# Read in model recovery
summary_aic <- read.csv(paste0(in_path, "model_recov/summary_aic.csv"))
summary_bic <- read.csv(paste0(in_path, "model_recov/summary_bic.csv"))
mdn_aic_long <- read.csv(paste0(in_path, "model_recov/mdn_aic_long.csv"))
mdn_bic_long <- read.csv(paste0(in_path, "model_recov/mdn_bic_long.csv"))
best_aic <- read.csv(paste0(in_path, "model_recov/bestAICmodelfits.csv"))
best_bic <- read.csv(paste0(in_path, "model_recov/bestBICmodelfits.csv"))

# Read in parameter fits
vanilla_fitParams <- read.csv(paste0(in_path, "param_fits/vanilla_fitParams.csv"))
vanilla_allSubFiles <- read.csv(paste0(in_path, "param_fits/vanilla_allSubFiles.csv"))
emoBias_fitParams <- read.csv(paste0(in_path, "param_fits/emoBias_fitParams.csv"))
emoBias_allSubFiles <- read.csv(paste0(in_path, "param_fits/emoBias_allSubFiles.csv"))
cond_fitParams <- read.csv(paste0(in_path, "param_fits/cond_fitParams.csv"))
cond_allSubFiles <- read.csv(paste0(in_path, "param_fits/cond_allSubFiles.csv"))
condBeta_fitParams <- read.csv(paste0(in_path, "param_fits/condBeta_fitParams.csv"))
condBeta_allSubFiles <- read.csv(paste0(in_path, "param_fits/condBeta_allSubFiles.csv"))

# Read in parameter recovery
vanilla_fullParams <- read.csv(paste0(in_path, "param_recov/vanilla_fullParams.csv"))
emoBias_fullParams <- read.csv(paste0(in_path, "param_recov/emoBias_fullParams.csv"))
cond_fullParams <- read.csv(paste0(in_path, "param_recov/cond_fullParams.csv"))
condBeta_fullParams <- read.csv(paste0(in_path, "param_recov/condBeta_fullParams.csv"))

### STEP 1: MODEL COMPARISON -----

### Final data processing -----

# Change aic variable types
str(aic_comp_long)
aic_comp_long$Model <- factor(aic_comp_long$Model, levels = mods)
aic_comp_long$AgeGroup <- factor(aic_comp_long$AgeGroup, levels = grps)
str(aic_comp_long)

# Create summary aic data frames
summary_aic_comp <- aic_comp_long %>%
  group_by(Model) %>%
  summarise(m = mean(AIC_Value), se = std.error(AIC_Value))
summary_aic_comp_age <- aic_comp_long %>%
  group_by(Model, AgeGroup) %>%
  summarise(m = mean(AIC_Value), se = std.error(AIC_Value))
summary_aic_comp_age_mdn <- aic_comp_long %>%
  group_by(Model, AgeGroup) %>%
  summarise(med = median(AIC_Value))

# Change bic variable types
str(bic_comp_long)
bic_comp_long$Model <- factor(bic_comp_long$Model, levels = mods)
bic_comp_long$AgeGroup <- factor(bic_comp_long$AgeGroup, levels = grps)
str(bic_comp_long)

# Create long and summary bic data frames
summary_bic_comp <- bic_comp_long %>%
  group_by(Model) %>%
  summarise(m = mean(BIC_Value), se = std.error(BIC_Value))
summary_bic_comp_age <- bic_comp_long %>%
  group_by(Model, AgeGroup) %>%
  summarise(m = mean(BIC_Value), se = std.error(BIC_Value))
summary_bic_comp_age_mdn <- bic_comp_long %>%
  group_by(Model, AgeGroup) %>%
  summarise(med = median(BIC_Value))

# Change pseudor variable types
str(pseudor_comp_long)
pseudor_comp_long$Model <- factor(pseudor_comp_long$Model, levels = c("Random", mods))
pseudor_comp_long$AgeGroup <- factor(pseudor_comp_long$AgeGroup, levels = grps)
str(pseudor_comp_long)

### AIC values across age -----

# Raw AIC value scatterplot
ggplot(data = aic_comp_long, aes(x = ExactAge, y = AIC_Value, group = Model, color = Model)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = seq(50, 250, by = 50), colour = 'grey80') +
  scale_y_continuous(breaks = seq(50, 250, by = 50)) +  
  geom_point(alpha = 1, size = 2.75) +
  labs(x = "Age (Years)", y = "AIC Values", colour = "Model") +
  scale_color_manual(values = nine_color_scheme_mods) +
  emorl_theme + theme(legend.position = "right")
ggsave(paste0(out_path, "model_comp/aic_age.png"), plot = last_plot(), width = 9, height = 5)

# Raw AIC value histogram
ggplot(data = aic_comp_long, aes(x = AIC_Value, fill = Model, color = Model)) +
  scale_fill_manual(values = nine_color_scheme_mods) + 
  scale_color_manual(values = nine_color_scheme_mods) + 
  scale_x_continuous(breaks = seq(50, 250, 100)) +
  geom_hline(yintercept = seq(0, 240, by = 60), colour = 'grey80') +
  scale_y_continuous(breaks = seq(0, 240, by = 60)) +  
  geom_histogram(binwidth = 25, alpha = .75) +
  facet_grid(cols = vars(AgeGroup)) +
  labs(x = "AIC Values", y = "Number of Participants", colour = "Model") +
  emorl_theme + theme(legend.position = "right")
ggsave(paste0(out_path, "model_comp/aic_dist.png"), plot = last_plot(), width = 9, height = 5)

### BIC values across age -----

# Raw BIC value scatterplot
ggplot(data = bic_comp_long, aes(x = ExactAge, y = BIC_Value, group = Model, color = Model)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = seq(75, 275, by = 50), colour = 'grey80') +
  scale_y_continuous(breaks = seq(75, 275, by = 50)) +  
  geom_point(alpha = 1, size = 2.75) +
  labs(x = "Age (Years)", y = "BIC Values", colour = "Model") +
  scale_color_manual(values = nine_color_scheme_mods) +
  emorl_theme + theme(legend.position = "right")
ggsave(paste0(out_path, "model_comp/bic_age.png"), plot = last_plot(), width = 9, height = 5)

# Raw BIC value histogram
ggplot(data = bic_comp_long, aes(x = BIC_Value, fill = Model, color = Model)) +
  scale_fill_manual(values = nine_color_scheme_mods) + 
  scale_color_manual(values = nine_color_scheme_mods) + 
  scale_x_continuous(breaks = seq(75, 275, 100)) +
  geom_hline(yintercept = seq(0, 120, by = 30), colour = 'grey80') +
  scale_y_continuous(breaks = seq(0, 120, by = 30)) +  
  geom_histogram(binwidth = 25, alpha = .75) +
  facet_grid(cols = vars(AgeGroup)) +
  labs(x = "BIC Values", y = "Number of Participants", colour = "Model") +
  emorl_theme + theme(legend.position = "right")
ggsave(paste0(out_path, "model_comp/bic_dist.png"), plot = last_plot(), width = 9, height = 5)

### AIC summary -----

# Summary AIC values across the sample
ggplot(data = summary_aic_comp, aes(x = Model, y = m, color = Model)) +
  geom_hline(yintercept = seq(195, 235, by = 10), colour = 'grey80') +
  scale_y_continuous(breaks = seq(195, 235, by = 10)) +  
  geom_point(alpha = 1, size = 2.75) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = .1, size = 1) +
  annotate("text", x = "Vanilla", y = median(aic_comp_long$AIC_Value[aic_comp_long$Model == "Vanilla"]), label = "x", color = vanilla, size = 4) +
  annotate("text", x = "EmoQ", y = median(aic_comp_long$AIC_Value[aic_comp_long$Model == "EmoQ"]), label = "x", color = emo_q, size = 4) +
  annotate("text", x = "EmoBias", y = median(aic_comp_long$AIC_Value[aic_comp_long$Model == "EmoBias"]), label = "x", color = emo_bias, size = 4) +
  annotate("text", x = "EmoHier", y = median(aic_comp_long$AIC_Value[aic_comp_long$Model == "EmoHier"]), label = "x", color = emo_heir, size = 4) +
  annotate("text", x = "EmoAvg", y = median(aic_comp_long$AIC_Value[aic_comp_long$Model == "EmoAvg"]), label = "x", color = emo_avg, size = 4) +
  annotate("text", x = "Cong", y = median(aic_comp_long$AIC_Value[aic_comp_long$Model == "Cong"]), label = "x", color = cong, size = 4) +
  annotate("text", x = "Cond", y = median(aic_comp_long$AIC_Value[aic_comp_long$Model == "Cond"]), label = "x", color = cond, size = 4) +
  annotate("text", x = "CongBeta", y = median(aic_comp_long$AIC_Value[aic_comp_long$Model == "CongBeta"]), label = "x", color = cong_beta, size = 4) +
  annotate("text", x = "CondBeta", y = median(aic_comp_long$AIC_Value[aic_comp_long$Model == "CondBeta"]), label = "x", color = cond_beta, size = 4) +
  labs(x = NULL, y = "Summary AIC Values") +
  scale_color_manual(values = nine_color_scheme_mods) +
  emorl_theme + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(out_path, "model_comp/figS2.png"), plot = last_plot(), width = 9, height = 7)

# Mean AIC values by age group
ggplot(data = summary_aic_comp_age, aes(x = Model, y = m, color = Model)) +
  geom_hline(yintercept = seq(165, 245, by = 20), colour = 'grey80') +
  scale_y_continuous(breaks = seq(165, 245, by = 20)) +  
  geom_point(alpha = 1, size = 2.75) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = .1, size = 1) +
  facet_grid(cols = vars(AgeGroup)) +
  labs(x = "Model", y = "Mean AIC Values") +
  scale_color_manual(values = nine_color_scheme_mods) +
  emorl_theme + theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(out_path, "model_comp/mean_aic_age.png"), plot = last_plot(), width = 9, height = 7)

# Median AIC values by age group
summary_aic_comp_age_mdn

### BIC summary -----

# Summary BIC values across the sample
ggplot(data = summary_bic_comp, aes(x = Model, y = m, color = Model)) +
  geom_hline(yintercept = seq(205, 245, by = 10), colour = 'grey80') +
  scale_y_continuous(breaks = seq(205, 245, by = 10)) +  
  geom_point(alpha = 1, size = 2.75) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = .1, size = 1) +
  annotate("text", x = "Vanilla", y = median(bic_comp_long$BIC_Value[bic_comp_long$Model == "Vanilla"]), label = "x", color = vanilla, size = 4) +
  annotate("text", x = "EmoQ", y = median(bic_comp_long$BIC_Value[bic_comp_long$Model == "EmoQ"]), label = "x", color = emo_q, size = 4) +
  annotate("text", x = "EmoBias", y = median(bic_comp_long$BIC_Value[bic_comp_long$Model == "EmoBias"]), label = "x", color = emo_bias, size = 4) +
  annotate("text", x = "EmoHier", y = median(bic_comp_long$BIC_Value[bic_comp_long$Model == "EmoHier"]), label = "x", color = emo_heir, size = 4) +
  annotate("text", x = "EmoAvg", y = median(bic_comp_long$BIC_Value[bic_comp_long$Model == "EmoAvg"]), label = "x", color = emo_avg, size = 4) +
  annotate("text", x = "Cong", y = median(bic_comp_long$BIC_Value[bic_comp_long$Model == "Cong"]), label = "x", color = cong, size = 4) +
  annotate("text", x = "Cond", y = median(bic_comp_long$BIC_Value[bic_comp_long$Model == "Cond"]), label = "x", color = cond, size = 4) +
  annotate("text", x = "CongBeta", y = median(bic_comp_long$BIC_Value[bic_comp_long$Model == "CongBeta"]), label = "x", color = cong_beta, size = 4) +
  annotate("text", x = "CondBeta", y = median(bic_comp_long$BIC_Value[bic_comp_long$Model == "CondBeta"]), label = "x", color = cond_beta, size = 4) +
  labs(x = NULL, y = "Summary BIC Values") +
  scale_color_manual(values = nine_color_scheme_mods) +
  emorl_theme + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(out_path, "model_comp/summ_bic.png"), plot = last_plot(), width = 9, height = 7)

# Mean BIC values by age group
ggplot(data = summary_bic_comp_age, aes(x = Model, y = m, color = Model)) +
  geom_hline(yintercept = seq(175, 255, by = 20), colour = 'grey80') +
  scale_y_continuous(breaks = seq(175, 255, by = 20)) +  
  geom_point(alpha = 1, size = 2.75) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = .1, size = 1) +
  facet_grid(cols = vars(AgeGroup)) +
  labs(x = "Model", y = "Mean BIC Values") +
  scale_color_manual(values = nine_color_scheme_mods) +
  emorl_theme + theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(out_path, "model_comp/mean_bic_age.png"), plot = last_plot(), width = 9, height = 7)

# Median BIC values by age group
summary_bic_comp_age_mdn

### Best-fit counts across age -----

# Prepare for AIC counts by age group
aic_comp_long_simp <- aic_comp_long[aic_comp_long$Model %in% best_mods, ]
aic_comp_wide_simp <- pivot_wider(aic_comp_long_simp, names_from = "Model", values_from = "AIC_Value")
aic_comp_wide_simp$minAIC <- apply(aic_comp_wide_simp[, 4:7], 1, function(row) {
  colnames(aic_comp_wide_simp[, 4:7])[which.min(row)]
})
aic_comp_wide_simp$minAIC <- factor(aic_comp_wide_simp$minAIC, levels = best_mods)

# Display AIC counts in plot
ggplot(aic_comp_wide_simp, aes(x = minAIC, fill = minAIC, color = minAIC)) +
  scale_fill_manual(values = c(vanilla, emo_bias, cond, cond_beta)) + 
  scale_color_manual(values = c(vanilla, emo_bias, cond, cond_beta)) + 
  geom_hline(yintercept = seq(0, 18, by = 6), colour = 'grey80') +
  scale_y_continuous(breaks = seq(0, 18, by = 6)) +  
  geom_bar(alpha = .75, show.legend = FALSE) + 
  facet_grid(cols = vars(AgeGroup)) +
  emorl_theme +
  labs(x = NULL, y = "Number of Participants") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(out_path, "model_comp/figS5.png"), plot = last_plot(), width = 9, height = 7)

# Display AIC counts in table
table(aic_comp_wide_simp$AgeGroup, aic_comp_wide_simp$minAIC)

# Prepare for BIC counts by age group
bic_comp_long_simp <- bic_comp_long[bic_comp_long$Model %in% best_mods, ]
bic_comp_wide_simp <- pivot_wider(bic_comp_long_simp, names_from = "Model", values_from = "BIC_Value")
bic_comp_wide_simp$minBIC <- apply(bic_comp_wide_simp[, 4:7], 1, function(row) {
  colnames(bic_comp_wide_simp[, 4:7])[which.min(row)]
})
bic_comp_wide_simp$minBIC <- factor(bic_comp_wide_simp$minBIC, levels = best_mods)

# Display BIC counts in plot
ggplot(bic_comp_wide_simp, aes(x = minBIC, fill = minBIC, color = minBIC)) +
  scale_fill_manual(values = c(vanilla, emo_bias, cond, cond_beta)) + 
  scale_color_manual(values = c(vanilla, emo_bias, cond, cond_beta)) + 
  geom_hline(yintercept = seq(0, 24, by = 8), colour = 'grey80') +
  scale_y_continuous(breaks = seq(0, 24, by = 8)) +  
  geom_bar(alpha = .75, show.legend = FALSE) + 
  facet_grid(cols = vars(AgeGroup)) +
  emorl_theme +
  labs(x = NULL, y = "Number of Participants") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(out_path, "model_comp/bic_counts.png"), plot = last_plot(), width = 9, height = 7)

# Display BIC counts in table
table(bic_comp_wide_simp$AgeGroup, bic_comp_wide_simp$minBIC)

### Likelihood probability values -----

# Pseudor values
pseudor_comp_long <- pseudor_comp_long[pseudor_comp_long$Model != "Random", ]
ggplot(data = pseudor_comp_long, aes(x = ExactAge, y = Pseudor_Value, group = Model, color = Model)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = c(0,1), colour = 'black', linetype = "longdash") +
  geom_hline(yintercept = c(.25,.5,.75), colour = 'grey80') +
  scale_y_continuous(breaks = c(0,.25,.5,.75,1)) +  
  geom_point(alpha = 1, size = 2.75) +
  annotate("text", x = 6.6, y = median(pseudor_comp_long$Pseudor_Value[pseudor_comp_long$Model == "Vanilla"]), label = "x", color = vanilla, size = 4) +
  annotate("text", x = 6.7, y = median(pseudor_comp_long$Pseudor_Value[pseudor_comp_long$Model == "EmoQ"]), label = "x", color = emo_q, size = 4) +
  annotate("text", x = 6.8, y = median(pseudor_comp_long$Pseudor_Value[pseudor_comp_long$Model == "EmoBias"]), label = "x", color = emo_bias, size = 4) +
  annotate("text", x = 6.9, y = median(pseudor_comp_long$Pseudor_Value[pseudor_comp_long$Model == "EmoHier"]), label = "x", color = emo_heir, size = 4) +
  annotate("text", x = 7, y = median(pseudor_comp_long$Pseudor_Value[pseudor_comp_long$Model == "EmoAvg"]), label = "x", color = emo_avg, size = 4) +
  annotate("text", x = 7.1, y = median(pseudor_comp_long$Pseudor_Value[pseudor_comp_long$Model == "Cong"]), label = "x", color = cong, size = 4) +
  annotate("text", x = 7.2, y = median(pseudor_comp_long$Pseudor_Value[pseudor_comp_long$Model == "Cond"]), label = "x", color = cond, size = 4) +
  annotate("text", x = 7.3, y = median(pseudor_comp_long$Pseudor_Value[pseudor_comp_long$Model == "CongBeta"]), label = "x", color = cong_beta, size = 4) +
  annotate("text", x = 7.4, y = median(pseudor_comp_long$Pseudor_Value[pseudor_comp_long$Model == "CondBeta"]), label = "x", color = cond_beta, size = 4) +
  labs(x = "Age (Years)", y = "Pseudor Values") +
  scale_color_manual(values = nine_color_scheme_mods) +
  emorl_theme + theme(legend.position = "right")
ggsave(paste0(out_path, "model_comp/summ_pseudor_age.png"), plot = last_plot(), width = 9, height = 5)

### STEP 2: MODEL RECOVERY -----

### Vanilla model, AIC/BIC -----

# Summary AIC values
vanilla_aic_summary <- summary_aic[grep("*vanillaM", summary_aic$Data_Model),]
vanilla_aic_summary <- vanilla_aic_summary %>%
  mutate(Data_Model = case_when(
    Data_Model == "aic_vanillaM_vanillaD" ~ "Vanilla",
    Data_Model == "aic_vanillaM_emoBiasD" ~ "EmoBias",
    Data_Model == "aic_vanillaM_condD" ~ "Cond",
    Data_Model == "aic_vanillaM_condBetaD" ~ "CondBeta"
  ))
vanilla_aic_summary$Data_Model <- factor(vanilla_aic_summary$Data_Model, levels = best_mods)
ggplot(data = vanilla_aic_summary, aes(x = Data_Model, y = m, color = Data_Model)) +
  geom_hline(yintercept = summ_aicbic_range, colour = 'grey80') +
  scale_y_continuous(breaks = summ_aicbic_range) +  
  geom_point(alpha = 1, size = 2.75) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = .1, size = 1) +
  annotate("text", x = "Vanilla", y = mdn_aic_long[mdn_aic_long$Data_Model == "mdnAIC_vanillaM_vanillaD", "AIC_Value"], label = "x", color = vanilla, size = 4) +
  annotate("text", x = "EmoBias", y = mdn_aic_long[mdn_aic_long$Data_Model == "mdnAIC_vanillaM_emoBiasD", "AIC_Value"], label = "x", color = emo_bias, size = 4) +
  annotate("text", x = "Cond", y = mdn_aic_long[mdn_aic_long$Data_Model == "mdnAIC_vanillaM_condD", "AIC_Value"], label = "x", color = cond, size = 4) +
  annotate("text", x = "CondBeta", y = mdn_aic_long[mdn_aic_long$Data_Model == "mdnAIC_vanillaM_condBetaD", "AIC_Value"], label = "x", color = cond_beta, size = 4) +
  labs(x = "Input Data", y = "Summary AIC Values", title = "Vanilla Model") +
  scale_color_manual(values = four_color_scheme_mods) +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "model_recov/summ_aic_vanilla.png"), plot = last_plot(), width = 9, height = 7)

# Summary BIC values
vanilla_bic_summary <- summary_bic[grep("*vanillaM", summary_bic$Data_Model),]
vanilla_bic_summary <- vanilla_bic_summary %>%
  mutate(Data_Model = case_when(
    Data_Model == "bic_vanillaM_vanillaD" ~ "Vanilla",
    Data_Model == "bic_vanillaM_emoBiasD" ~ "EmoBias",
    Data_Model == "bic_vanillaM_condD" ~ "Cond",
    Data_Model == "bic_vanillaM_condBetaD" ~ "CondBeta"
  ))
vanilla_bic_summary$Data_Model <- factor(vanilla_bic_summary$Data_Model, levels = best_mods)
ggplot(data = vanilla_bic_summary, aes(x = Data_Model, y = m, color = Data_Model)) +
  geom_hline(yintercept = summ_aicbic_range, colour = 'grey80') +
  scale_y_continuous(breaks = summ_aicbic_range) +  
  geom_point(alpha = 1, size = 2.75) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = .1, size = 1) +
  annotate("text", x = "Vanilla", y = mdn_bic_long[mdn_bic_long$Data_Model == "mdnBIC_vanillaM_vanillaD", "BIC_Value"], label = "x", color = vanilla, size = 4) +
  annotate("text", x = "EmoBias", y = mdn_bic_long[mdn_bic_long$Data_Model == "mdnBIC_vanillaM_emoBiasD", "BIC_Value"], label = "x", color = emo_bias, size = 4) +
  annotate("text", x = "Cond", y = mdn_bic_long[mdn_bic_long$Data_Model == "mdnBIC_vanillaM_condD", "BIC_Value"], label = "x", color = cond, size = 4) +
  annotate("text", x = "CondBeta", y = mdn_bic_long[mdn_bic_long$Data_Model == "mdnBIC_vanillaM_condBetaD", "BIC_Value"], label = "x", color = cond_beta, size = 4) +
  labs(x = "Input Data", y = "Summary BIC Values", title = "Vanilla Model") +
  scale_color_manual(values = four_color_scheme_mods) +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "model_recov/summ_bic_vanilla.png"), plot = last_plot(), width = 9, height = 7)

### EmoBias model, AIC/BIC -----

# Summary AIC values
emoBias_aic_summary <- summary_aic[grep("*emoBiasM", summary_aic$Data_Model),]
emoBias_aic_summary <- emoBias_aic_summary %>%
  mutate(Data_Model = case_when(
    Data_Model == "aic_emoBiasM_vanillaD" ~ "Vanilla",
    Data_Model == "aic_emoBiasM_emoBiasD" ~ "EmoBias",
    Data_Model == "aic_emoBiasM_condD" ~ "Cond",
    Data_Model == "aic_emoBiasM_condBetaD" ~ "CondBeta"
  ))
emoBias_aic_summary$Data_Model <- factor(emoBias_aic_summary$Data_Model, levels = best_mods)
ggplot(data = emoBias_aic_summary, aes(x = Data_Model, y = m, color = Data_Model)) +
  geom_hline(yintercept = summ_aicbic_range, colour = 'grey80') +
  scale_y_continuous(breaks = summ_aicbic_range) +  
  geom_point(alpha = 1, size = 2.75) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = .1, size = 1) +
  annotate("text", x = "Vanilla", y = mdn_aic_long[mdn_aic_long$Data_Model == "mdnAIC_emoBiasM_vanillaD", "AIC_Value"], label = "x", color = vanilla, size = 4) +
  annotate("text", x = "EmoBias", y = mdn_aic_long[mdn_aic_long$Data_Model == "mdnAIC_emoBiasM_emoBiasD", "AIC_Value"], label = "x", color = emo_bias, size = 4) +
  annotate("text", x = "Cond", y = mdn_aic_long[mdn_aic_long$Data_Model == "mdnAIC_emoBiasM_condD", "AIC_Value"], label = "x", color = cond, size = 4) +
  annotate("text", x = "CondBeta", y = mdn_aic_long[mdn_aic_long$Data_Model == "mdnAIC_emoBiasM_condBetaD", "AIC_Value"], label = "x", color = cond_beta, size = 4) +
  labs(x = "Input Data", y = "Summary AIC Values", title = "EmoBias Model") +
  scale_color_manual(values = four_color_scheme_mods) +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "model_recov/summ_aic_emoBias.png"), plot = last_plot(), width = 9, height = 7)

# Summary BIC values
emoBias_bic_summary <- summary_bic[grep("*emoBiasM", summary_bic$Data_Model),]
emoBias_bic_summary <- emoBias_bic_summary %>%
  mutate(Data_Model = case_when(
    Data_Model == "bic_emoBiasM_vanillaD" ~ "Vanilla",
    Data_Model == "bic_emoBiasM_emoBiasD" ~ "EmoBias",
    Data_Model == "bic_emoBiasM_condD" ~ "Cond",
    Data_Model == "bic_emoBiasM_condBetaD" ~ "CondBeta"
  ))
emoBias_bic_summary$Data_Model <- factor(emoBias_bic_summary$Data_Model, levels = best_mods)
ggplot(data = emoBias_bic_summary, aes(x = Data_Model, y = m, color = Data_Model)) +
  geom_hline(yintercept = summ_aicbic_range, colour = 'grey80') +
  scale_y_continuous(breaks = summ_aicbic_range) +  
  geom_point(alpha = 1, size = 2.75) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = .1, size = 1) +
  annotate("text", x = "Vanilla", y = mdn_bic_long[mdn_bic_long$Data_Model == "mdnBIC_emoBiasM_vanillaD", "BIC_Value"], label = "x", color = vanilla, size = 4) +
  annotate("text", x = "EmoBias", y = mdn_bic_long[mdn_bic_long$Data_Model == "mdnBIC_emoBiasM_emoBiasD", "BIC_Value"], label = "x", color = emo_bias, size = 4) +
  annotate("text", x = "Cond", y = mdn_bic_long[mdn_bic_long$Data_Model == "mdnBIC_emoBiasM_condD", "BIC_Value"], label = "x", color = cond, size = 4) +
  annotate("text", x = "CondBeta", y = mdn_bic_long[mdn_bic_long$Data_Model == "mdnBIC_emoBiasM_condBetaD", "BIC_Value"], label = "x", color = cond_beta, size = 4) +
  labs(x = "Input Data", y = "Summary BIC Values", title = "EmoBias Model") +
  scale_color_manual(values = four_color_scheme_mods) +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "model_recov/summ_bic_emoBias.png"), plot = last_plot(), width = 9, height = 7)

### Cond model, AIC/BIC -----

# Summary AIC values
cond_aic_summary <- summary_aic[grep("*condM", summary_aic$Data_Model),]
cond_aic_summary <- cond_aic_summary %>%
  mutate(Data_Model = case_when(
    Data_Model == "aic_condM_vanillaD" ~ "Vanilla",
    Data_Model == "aic_condM_emoBiasD" ~ "EmoBias",
    Data_Model == "aic_condM_condD" ~ "Cond",
    Data_Model == "aic_condM_condBetaD" ~ "CondBeta"
  ))
cond_aic_summary$Data_Model <- factor(cond_aic_summary$Data_Model, levels = best_mods)
ggplot(data = cond_aic_summary, aes(x = Data_Model, y = m, color = Data_Model)) +
  geom_hline(yintercept = summ_aicbic_range, colour = 'grey80') +
  scale_y_continuous(breaks = summ_aicbic_range) +  
  geom_point(alpha = 1, size = 2.75) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = .1, size = 1) +
  annotate("text", x = "Vanilla", y = mdn_aic_long[mdn_aic_long$Data_Model == "mdnAIC_condM_vanillaD", "AIC_Value"], label = "x", color = vanilla, size = 4) +
  annotate("text", x = "EmoBias", y = mdn_aic_long[mdn_aic_long$Data_Model == "mdnAIC_condM_emoBiasD", "AIC_Value"], label = "x", color = emo_bias, size = 4) +
  annotate("text", x = "Cond", y = mdn_aic_long[mdn_aic_long$Data_Model == "mdnAIC_condM_condD", "AIC_Value"], label = "x", color = cond, size = 4) +
  annotate("text", x = "CondBeta", y = mdn_aic_long[mdn_aic_long$Data_Model == "mdnAIC_condM_condBetaD", "AIC_Value"], label = "x", color = cond_beta, size = 4) +
  labs(x = "Input Data", y = "Summary AIC Values", title = "Cond Model") +
  scale_color_manual(values = four_color_scheme_mods) +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "model_recov/summ_aic_cond.png"), plot = last_plot(), width = 9, height = 7)

# Summary BIC values
cond_bic_summary <- summary_bic[grep("*condM", summary_bic$Data_Model),]
cond_bic_summary <- cond_bic_summary %>%
  mutate(Data_Model = case_when(
    Data_Model == "bic_condM_vanillaD" ~ "Vanilla",
    Data_Model == "bic_condM_emoBiasD" ~ "EmoBias",
    Data_Model == "bic_condM_condD" ~ "Cond",
    Data_Model == "bic_condM_condBetaD" ~ "CondBeta"
  ))
cond_bic_summary$Data_Model <- factor(cond_bic_summary$Data_Model, levels = best_mods)
ggplot(data = cond_bic_summary, aes(x = Data_Model, y = m, color = Data_Model)) +
  geom_hline(yintercept = summ_aicbic_range, colour = 'grey80') +
  scale_y_continuous(breaks = summ_aicbic_range) +  
  geom_point(alpha = 1, size = 2.75) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = .1, size = 1) +
  annotate("text", x = "Vanilla", y = mdn_bic_long[mdn_bic_long$Data_Model == "mdnBIC_condM_vanillaD", "BIC_Value"], label = "x", color = vanilla, size = 4) +
  annotate("text", x = "EmoBias", y = mdn_bic_long[mdn_bic_long$Data_Model == "mdnBIC_condM_emoBiasD", "BIC_Value"], label = "x", color = emo_bias, size = 4) +
  annotate("text", x = "Cond", y = mdn_bic_long[mdn_bic_long$Data_Model == "mdnBIC_condM_condD", "BIC_Value"], label = "x", color = cond, size = 4) +
  annotate("text", x = "CondBeta", y = mdn_bic_long[mdn_bic_long$Data_Model == "mdnBIC_condM_condBetaD", "BIC_Value"], label = "x", color = cond_beta, size = 4) +
  labs(x = "Input Data", y = "Summary BIC Values", title = "Cond Model") +
  scale_color_manual(values = four_color_scheme_mods) +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "model_recov/summ_bic_cond.png"), plot = last_plot(), width = 9, height = 7)

### CondBeta model, AIC/BIC -----

# Summary AIC values
condBeta_aic_summary <- summary_aic[grep("*condBetaM", summary_aic$Data_Model),]
condBeta_aic_summary <- condBeta_aic_summary %>%
  mutate(Data_Model = case_when(
    Data_Model == "aic_condBetaM_vanillaD" ~ "Vanilla",
    Data_Model == "aic_condBetaM_emoBiasD" ~ "EmoBias",
    Data_Model == "aic_condBetaM_condD" ~ "Cond",
    Data_Model == "aic_condBetaM_condBetaD" ~ "CondBeta"
  ))
condBeta_aic_summary$Data_Model <- factor(condBeta_aic_summary$Data_Model, levels = best_mods)
ggplot(data = condBeta_aic_summary, aes(x = Data_Model, y = m, color = Data_Model)) +
  geom_hline(yintercept = summ_aicbic_range, colour = 'grey80') +
  scale_y_continuous(breaks = summ_aicbic_range) +  
  geom_point(alpha = 1, size = 2.75) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = .1, size = 1) +
  annotate("text", x = "Vanilla", y = mdn_aic_long[mdn_aic_long$Data_Model == "mdnAIC_condBetaM_vanillaD", "AIC_Value"], label = "x", color = vanilla, size = 4) +
  annotate("text", x = "EmoBias", y = mdn_aic_long[mdn_aic_long$Data_Model == "mdnAIC_condBetaM_emoBiasD", "AIC_Value"], label = "x", color = emo_bias, size = 4) +
  annotate("text", x = "Cond", y = mdn_aic_long[mdn_aic_long$Data_Model == "mdnAIC_condBetaM_condD", "AIC_Value"], label = "x", color = cond, size = 4) +
  annotate("text", x = "CondBeta", y = mdn_aic_long[mdn_aic_long$Data_Model == "mdnAIC_condBetaM_condBetaD", "AIC_Value"], label = "x", color = cond_beta, size = 4) +
  labs(x = "Input Data", y = "Summary AIC Values", title = "CondBeta Model") +
  scale_color_manual(values = four_color_scheme_mods) +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "model_recov/summ_aic_condBeta.png"), plot = last_plot(), width = 9, height = 7)

# Summary BIC values
condBeta_bic_summary <- summary_bic[grep("*condBetaM", summary_bic$Data_Model),]
condBeta_bic_summary <- condBeta_bic_summary %>%
  mutate(Data_Model = case_when(
    Data_Model == "bic_condBetaM_vanillaD" ~ "Vanilla",
    Data_Model == "bic_condBetaM_emoBiasD" ~ "EmoBias",
    Data_Model == "bic_condBetaM_condD" ~ "Cond",
    Data_Model == "bic_condBetaM_condBetaD" ~ "CondBeta"
  ))
condBeta_bic_summary$Data_Model <- factor(condBeta_bic_summary$Data_Model, levels = best_mods)
ggplot(data = condBeta_bic_summary, aes(x = Data_Model, y = m, color = Data_Model)) +
  geom_hline(yintercept = summ_aicbic_range, colour = 'grey80') +
  scale_y_continuous(breaks = summ_aicbic_range) +  
  geom_point(alpha = 1, size = 2.75) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = .1, size = 1) +
  annotate("text", x = "Vanilla", y = mdn_bic_long[mdn_bic_long$Data_Model == "mdnBIC_condBetaM_vanillaD", "BIC_Value"], label = "x", color = vanilla, size = 4) +
  annotate("text", x = "EmoBias", y = mdn_bic_long[mdn_bic_long$Data_Model == "mdnBIC_condBetaM_emoBiasD", "BIC_Value"], label = "x", color = emo_bias, size = 4) +
  annotate("text", x = "Cond", y = mdn_bic_long[mdn_bic_long$Data_Model == "mdnBIC_condBetaM_condD", "BIC_Value"], label = "x", color = cond, size = 4) +
  annotate("text", x = "CondBeta", y = mdn_bic_long[mdn_bic_long$Data_Model == "mdnBIC_condBetaM_condBetaD", "BIC_Value"], label = "x", color = cond_beta, size = 4) +
  labs(x = "Input Data", y = "Summary BIC Values", title = "CondBeta Model") +
  scale_color_manual(values = four_color_scheme_mods) +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "model_recov/summ_bic_condBeta.png"), plot = last_plot(), width = 9, height = 7)

### Best-fit counts -----

# Save best-fit counts
fit_counts <- paste0(out_path, 'model_recov/fit_counts.txt')
sink(fit_counts)
cat("Vanilla Model AIC Counts\n\n")
best_aic[, c(grep("*vanillaM", names(best_aic)))]
cat("\nVanilla Model BIC Counts\n\n")
best_bic[, c(grep("*vanillaM", names(best_bic)))]
cat("\nEmoBias Model AIC Counts\n\n")
best_aic[, c(grep("*emoBiasM", names(best_aic)))]
cat("\nEmoBias Model BIC Counts\n\n")
best_bic[, c(grep("*emoBiasM", names(best_bic)))]
cat("\nCond Model AIC Counts\n\n")
best_aic[, c(grep("*condM", names(best_aic)))]
cat("\nCond Model BIC Counts\n\n")
best_bic[, c(grep("*condM", names(best_bic)))]
cat("\nCondBeta Model AIC Counts\n\n")
best_aic[, c(grep("*condBetaM", names(best_aic)))]
cat("\nCondBeta Model BIC Counts\n\n")
best_bic[, c(grep("*condBetaM", names(best_bic)))]
sink()

### STEP 3: PARAMETER RECOVERY AND FITS -----

### Define q-value estimate variables -----

# Column to legend label mapping
deck_map <- c(CONcon_HR_1 = "H+_H+F-",
              conCON_FZ_2 = "F-_H+F-",
              CONincon_HR_3 = "H+_H+H-",
              conINCON_HZ_4 = "H-_H+H-",
              INCONcon_FR_5 = "F+_F+F-",
              inconCON_FZ_6 = "F-_F+F-",
              INCONincon_FR_7 = "F+_F+H-",
              inconINCON_HZ_8 = "H-_F+H-")

# Legend levels
deck_levels <- unname(deck_map)

### Vanilla q-value estimates -----

# Get average q-value per trial
vanilla_allQMeans <- vanilla_allSubFiles %>% select(trial_labels, all_of(names(deck_map))) %>%
  pivot_longer(cols = -trial_labels,
               names_to = "deck",
               values_to = "q") %>%
  mutate(CardDeckType = factor(deck_map[deck], levels = deck_levels)) %>%
  group_by(trial_labels, CardDeckType) %>%
  summarise(m = mean(q, na.rm = TRUE), se = std.error(q), .groups = "drop")

# Plot average q-value per trial across trials
ggplot(data = vanilla_allQMeans, aes(x = trial_labels, y = m, group = CardDeckType, color = CardDeckType)) +
  geom_point(alpha = 1, size = 2) +
  geom_line(size = 1) +
  geom_hline(yintercept = c(.25, .75), colour = 'grey80') +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed") +
  scale_y_continuous(breaks = c(.25, .5, .75)) +
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  labs(x = "Trial Number", y = "Mean Q-Value Estimate", colour = "Card Deck", title = "Vanilla") +
  emorl_theme + theme(legend.position = "bottom")
ggsave(paste0(out_path, "param_fits/vanilla_qval.png"), plot = last_plot(), width = 9, height = 7)

### EmoBias q-value estimates -----

# Get average q-value per trial
emoBias_allQMeans <- emoBias_allSubFiles %>% select(trial_labels, all_of(names(deck_map))) %>%
  pivot_longer(cols = -trial_labels,
               names_to = "deck",
               values_to = "q") %>%
  mutate(CardDeckType = factor(deck_map[deck], levels = deck_levels)) %>%
  group_by(trial_labels, CardDeckType) %>%
  summarise(m = mean(q, na.rm = TRUE), se = std.error(q), .groups = "drop")

# Plot average q-value per trial across trials
ggplot(data = emoBias_allQMeans, aes(x = trial_labels, y = m, group = CardDeckType, color = CardDeckType)) +
  geom_point(alpha = 1, size = 2) +
  geom_line(size = 1) +
  geom_hline(yintercept = c(.25, .75), colour = 'grey80') +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed") +
  scale_y_continuous(breaks = c(.25, .5, .75)) +
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  labs(x = "Trial Number", y = "Mean Q-Value Estimate", colour = "Card Deck", title = "EmoBias") +
  emorl_theme + theme(legend.position = "bottom")
ggsave(paste0(out_path, "param_fits/emoBias_qval.png"), plot = last_plot(), width = 9, height = 7)

# EmoBias q-differences -----

# Prepare to plot q-differences
diff_df <- emoBias_allSubFiles
diff_df$`H+F-` <- diff_df$CONcon_HR_1 - diff_df$conCON_FZ_2
diff_df$`H+H-` <- diff_df$CONincon_HR_3 - diff_df$conINCON_HZ_4
diff_df$`F+F-` <- diff_df$INCONcon_FR_5 - diff_df$inconCON_FZ_6
diff_df$`F+H-` <- diff_df$INCONincon_FR_7 - diff_df$inconINCON_HZ_8
diff_df_simp <- diff_df[, c("StudyID", "AgeGroup", "trial_labels", "H+F-", "H+H-", "F+F-", "F+H-")]
diff_df_simp_long <- pivot_longer(diff_df_simp, cols = c("H+F-", "H+H-", "F+F-", "F+H-"), names_to = "Condition", values_to = "Q-Value Difference")
diff_df_simp_long$Condition <- factor(diff_df_simp_long$Condition, levels = c("H+F-", "H+H-", "F+F-", "F+H-"))
diff_df_simp_long$AgeGroup <- factor(diff_df_simp_long$AgeGroup, levels = grps)

# Plot q-differences at the sample level
ggplot(data = diff_df_simp_long, aes(x = trial_labels, y = `Q-Value Difference`, color = Condition, fill = Condition)) +
  geom_hline(yintercept = c(-.25, 0, .25), colour = 'grey80') +
  geom_hline(yintercept = .5, linetype = "dashed", colour = gold) +
  scale_y_continuous(breaks = c(-.25, 0, .25, .5)) +
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  scale_fill_manual(values = four_color_scheme_cond) + 
  scale_color_manual(values = four_color_scheme_cond) + 
  geom_smooth(method = "gam", size = 1) +
  labs(x = "Trial Number") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/qdiff.png"), plot = last_plot(), width = 9, height = 7)

# Plot q-differences by age group
ggplot(data = diff_df_simp_long, aes(x = trial_labels, y = `Q-Value Difference`, color = Condition, fill = Condition)) +
  geom_hline(yintercept = c(-.25, 0, .25), colour = 'grey80') +
  geom_hline(yintercept = .5, linetype = "dashed", colour = gold) +
  scale_y_continuous(breaks = c(-.25, 0, .25, .5)) +
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  scale_fill_manual(values = four_color_scheme_cond) + 
  scale_color_manual(values = four_color_scheme_cond) + 
  geom_smooth(method = "gam", size = 1) +
  facet_grid(cols = vars(AgeGroup)) +
  labs(x = "Trial Number") +
  emorl_theme + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(out_path, "param_fits/fig3.png"), plot = last_plot(), width = 9, height = 7)

### Cond q-value estimates -----

# Get average q-value per trial
cond_allQMeans <- cond_allSubFiles %>% select(trial_labels, all_of(names(deck_map))) %>%
  pivot_longer(cols = -trial_labels,
               names_to = "deck",
               values_to = "q") %>%
  mutate(CardDeckType = factor(deck_map[deck], levels = deck_levels)) %>%
  group_by(trial_labels, CardDeckType) %>%
  summarise(m = mean(q, na.rm = TRUE), se = std.error(q), .groups = "drop")

# Plot average q-value per trial across trials
ggplot(data = cond_allQMeans, aes(x = trial_labels, y = m, group = CardDeckType, color = CardDeckType)) +
  geom_point(alpha = 1, size = 2) +
  geom_line(size = 1) +
  geom_hline(yintercept = c(.25, .75), colour = 'grey80') +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed") +
  scale_y_continuous(breaks = c(.25, .5, .75)) +
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  labs(x = "Trial Number", y = "Mean Q-Value Estimate", colour = "Card Deck", title = "Cond") +
  emorl_theme + theme(legend.position = "bottom")
ggsave(paste0(out_path, "param_fits/cond_qval.png"), plot = last_plot(), width = 9, height = 7)

### CondBeta q-value estimates -----

# Get average q-value per trial
condBeta_allQMeans <- condBeta_allSubFiles %>% select(trial_labels, all_of(names(deck_map))) %>%
  pivot_longer(cols = -trial_labels,
               names_to = "deck",
               values_to = "q") %>%
  mutate(CardDeckType = factor(deck_map[deck], levels = deck_levels)) %>%
  group_by(trial_labels, CardDeckType) %>%
  summarise(m = mean(q, na.rm = TRUE), se = std.error(q), .groups = "drop")

# Plot average q-value per trial across trials
ggplot(data = condBeta_allQMeans, aes(x = trial_labels, y = m, group = CardDeckType, color = CardDeckType)) +
  geom_point(alpha = 1, size = 2) +
  geom_line(size = 1) +
  geom_hline(yintercept = c(.25, .75), colour = 'grey80') +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed") +
  scale_y_continuous(breaks = c(.25, .5, .75)) +
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  labs(x = "Trial Number", y = "Mean Q-Value Estimate", colour = "Card Deck", title = "CondBeta") +
  emorl_theme + theme(legend.position = "bottom")
ggsave(paste0(out_path, "param_fits/condBeta_qval.png"), plot = last_plot(), width = 9, height = 7)

### Vanilla RPE estimates -----

# Get means and standard errors for each trial's RPEs
vanilla_allSubFiles$RPE[vanilla_allSubFiles$RPE == -999999] <- NA # Recode missing values
range(vanilla_allSubFiles$RPE, na.rm = TRUE)
vanilla_mean_rpe <- vanilla_allSubFiles %>% 
  group_by(trial_labels) %>% 
  summarise(m = mean(RPE, na.rm = TRUE))

# Plot average RPE per trial across trials
ggplot(data = vanilla_mean_rpe, aes(x = trial_labels, y = m)) +
  geom_point(alpha = 1, size = 2) +
  geom_line(size = 1) +
  geom_hline(yintercept = c(-.2, .2), colour = 'grey80') +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  scale_y_continuous(breaks = c(-.2, 0, .2)) +
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  labs(x = "Trial Number", y = "Mean RPE", title = "Vanilla") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/vanilla_rpe.png"), plot = last_plot(), width = 9, height = 7)

### EmoBias RPE estimates -----

# Get means and standard errors for each trial's RPEs
emoBias_allSubFiles$RPE[emoBias_allSubFiles$RPE == -999999] <- NA # Recode missing values
range(emoBias_allSubFiles$RPE, na.rm = TRUE)
emoBias_mean_rpe <- emoBias_allSubFiles %>% 
  group_by(trial_labels) %>% 
  summarise(m = mean(RPE, na.rm = TRUE))

# Plot average RPE per trial across trials
ggplot(data = emoBias_mean_rpe, aes(x = trial_labels, y = m)) +
  geom_point(alpha = 1, size = 2) +
  geom_line(size = 1) +
  geom_hline(yintercept = c(-.2, .2), colour = 'grey80') +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  scale_y_continuous(breaks = c(-.2, 0, .2)) +
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  labs(x = "Trial Number", y = "Mean RPE", title = "EmoBias") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/emoBias_rpe.png"), plot = last_plot(), width = 9, height = 7)

### Cond RPE estimates -----

# Get means and standard errors for each trial's RPEs
cond_allSubFiles$RPE[cond_allSubFiles$RPE == -999999] <- NA # Recode missing values
range(cond_allSubFiles$RPE, na.rm = TRUE)
cond_mean_rpe <- cond_allSubFiles %>% 
  group_by(trial_labels) %>% 
  summarise(m = mean(RPE, na.rm = TRUE))

# Plot average RPE per trial across trials
ggplot(data = cond_mean_rpe, aes(x = trial_labels, y = m)) +
  geom_point(alpha = 1, size = 2) +
  geom_line(size = 1) +
  geom_hline(yintercept = c(-.2, .2), colour = 'grey80') +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  scale_y_continuous(breaks = c(-.2, 0, .2)) +
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  labs(x = "Trial Number", y = "Mean RPE", title = "Cond") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/cond_rpe.png"), plot = last_plot(), width = 9, height = 7)

### CondBeta RPE estimates -----

# Get means and standard errors for each trial's RPEs
condBeta_allSubFiles$RPE[condBeta_allSubFiles$RPE == -999999] <- NA # Recode missing values
range(condBeta_allSubFiles$RPE, na.rm = TRUE)
condBeta_mean_rpe <- condBeta_allSubFiles %>% 
  group_by(trial_labels) %>% 
  summarise(m = mean(RPE, na.rm = TRUE))

# Plot average RPE per trial across trials
ggplot(data = condBeta_mean_rpe, aes(x = trial_labels, y = m)) +
  geom_point(alpha = 1, size = 2) +
  geom_line(size = 1) +
  geom_hline(yintercept = c(-.2, .2), colour = 'grey80') +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  scale_y_continuous(breaks = c(-.2, 0, .2)) +
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  labs(x = "Trial Number", y = "Mean RPE", title = "CondBeta") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/condBeta_rpe.png"), plot = last_plot(), width = 9, height = 7)

### EmoBias weight estimates

# Get means and standard errors for each trial's weights
mean_wgt <- emoBias_allSubFiles %>% 
  group_by(trial_labels) %>% 
  summarise(m = mean(Weight, na.rm = TRUE))

# Plot average weight per trial across trials
ggplot(data = mean_wgt, aes(x = trial_labels, y = m)) +
  geom_point(alpha = 1, size = 2) +
  geom_line(size = 1) +
  geom_hline(yintercept = c(0, .25), colour = 'grey80') +
  geom_hline(aes(yintercept = .5), linetype = "dashed") +
  scale_y_continuous(breaks = c(0, .25, .5)) +
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  labs(x = "Trial Number", y = "Mean Weight", title = "EmoBias") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/vanilla_wgt.png"), plot = last_plot(), width = 9, height = 7)

### Vanilla, graphing beta parameter estimates -----

# Plot beta distribution
ggplot(data = vanilla_fitParams, aes(x = beta)) +
  scale_x_continuous(breaks = beta_axis_labels) +
  geom_hline(yintercept = seq(0, 36, by = 9), colour = 'grey90') +
  scale_y_continuous(breaks = seq(0, 36, by = 9)) + 
  geom_histogram(binwidth = 1, alpha = .75, color = vanilla, fill = vanilla) +
  labs(x = "Beta", y = "Number of Participants") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/vanilla_beta_dist.png"), plot = last_plot(), width = 7, height = 7)

# Plot beta across age
ggplot(data = vanilla_fitParams, aes(x = ExactAge, y = beta)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = beta_axis_labels, colour = 'grey90') +
  scale_y_continuous(breaks = beta_axis_labels) +
  geom_point(alpha = 1, size = 2.75, colour = vanilla) +
  geom_smooth(alpha = .25, method = "lm", colour = vanilla, fill = vanilla) +
  labs(x = "Age (Years)", y = "Beta") +
  emorl_theme 
ggsave(paste0(out_path, "param_fits/vanilla_beta_age.png"), plot = last_plot(), width = 7, height = 7)

### Vanilla, graphing alpha parameter estimates -----

# Plot alpha distribution
ggplot(data = vanilla_fitParams, aes(x = alpha)) +
  scale_x_continuous(breaks = alpha_axis_labels) +
  geom_hline(yintercept = seq(0, 60, by = 15), colour = 'grey90') +
  scale_y_continuous(breaks = seq(0, 60, by = 15)) + 
  geom_histogram(binwidth = .05, alpha = .75, color = vanilla, fill = vanilla) +
  labs(x = "Alpha", y = "Number of Participants") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/vanilla_alpha_dist.png"), plot = last_plot(), width = 7, height = 7)

# Plot alpha across age
ggplot(data = vanilla_fitParams, aes(x = ExactAge, y = alpha)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = alpha_axis_labels, colour = 'grey90') +
  scale_y_continuous(breaks = alpha_axis_labels) +
  geom_point(alpha = 1, size = 2.75, colour = vanilla) +
  geom_smooth(alpha = .25, method = "lm", colour = vanilla, fill = vanilla) +
  labs(x = "Age (Years)", y = "Alpha") +
  emorl_theme 
ggsave(paste0(out_path, "param_fits/vanilla_alpha_age.png"), plot = last_plot(), width = 7, height = 7)

### Vanilla, beta parameter comparisons -----

# Compare input and output betas
cor.test(vanilla_fullParams$beta_in, vanilla_fullParams$beta_out, method = "pearson")
ggplot(data = vanilla_fullParams, aes(x = beta_in, y = beta_out, color = factor(StudyID))) + 
  geom_point(size = 2.75) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = 'grey80') +
  geom_hline(yintercept = beta_axis_labels, colour = 'grey80') +
  geom_vline(xintercept = beta_axis_labels, colour = 'grey80') +
  scale_x_continuous(breaks = beta_axis_labels) +
  labs(x = "Randomly Sampled Input", y = "Fitted Output Estimate") +
  ggtitle("Beta") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "param_recov/vanilla_beta_cor.png"), plot = last_plot(), width = 7, height = 7)

### Vanilla, alpha parameter comparisons -----

# Compare input and output alphas
cor.test(vanilla_fullParams$alpha_in, vanilla_fullParams$alpha_out, method = "pearson")
ggplot(data = vanilla_fullParams, aes(x = alpha_in, y = alpha_out, color = factor(StudyID))) + 
  geom_point(size = 2.75) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = 'grey80') +
  geom_hline(yintercept = alpha_axis_labels, colour = 'grey80') +
  geom_vline(xintercept = alpha_axis_labels, colour = 'grey80') +
  labs(x = "Randomly Sampled Input", y = "Fitted Output Estimate") +
  ggtitle("Alpha") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "param_recov/vanilla_alpha_cor.png"), plot = last_plot(), width = 7, height = 7)

### Vanilla, all parameter comparisons -----

# Correlations between all parameters
vanilla_fullParams_Simp <- subset(vanilla_fullParams, select = -c(StudyID, lik_in, lik_out))
names(vanilla_fullParams_Simp) <- c("B_in", "a_in", "B_out", "a_out")
vanilla_p <- ggcorr(vanilla_fullParams_Simp, label = TRUE, label_size = 6)
vanilla_p <- vanilla_p + theme(axis.text.x = element_text(size = 12), 
                               axis.text.y = element_text(size = 12),
                               legend.text = element_text(size = 12))
vanilla_p
ggsave(paste0(out_path, "param_recov/vanilla.png"), plot = vanilla_p, width = 5, height = 4)

### EmoBias, graphing beta parameter estimates -----

# Plot beta distribution
ggplot(data = emoBias_fitParams, aes(x = beta)) +
  scale_x_continuous(breaks = beta_axis_labels) +
  geom_hline(yintercept = seq(0, 36, by = 9), colour = 'grey90') +
  scale_y_continuous(breaks = seq(0, 36, by = 9)) + 
  geom_histogram(binwidth = 1, alpha = .75, color = emo_bias, fill = emo_bias) +
  labs(x = "Beta", y = "Number of Participants") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/emoBias_beta_dist.png"), plot = last_plot(), width = 7, height = 7)

# Plot beta across age (note lm vs. glm model)
ggplot(data = emoBias_fitParams, aes(x = ExactAge, y = beta)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = beta_axis_labels, colour = 'grey90') +
  scale_y_continuous(breaks = beta_axis_labels) +
  geom_point(alpha = 1, size = 2.75, colour = emo_bias) +
  geom_smooth(alpha = .25, method = "lm", colour = emo_bias, fill = emo_bias) +
  labs(x = "Age (Years)", y = "Beta") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/emoBias_beta_age.png"), plot = last_plot(), width = 7, height = 7)

### EmoBias, graphing alpha parameter estimates -----

# Plot alpha distribution
ggplot(data = emoBias_fitParams, aes(x = alpha)) +
  scale_x_continuous(breaks = alpha_axis_labels) +
  geom_hline(yintercept = seq(0, 60, by = 15), colour = 'grey90') +
  scale_y_continuous(breaks = seq(0, 60, by = 15)) + 
  geom_histogram(binwidth = .05, alpha = .75, color = emo_bias, fill = emo_bias) +
  labs(x = "Alpha", y = "Number of Participants") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/emoBias_alpha_dist.png"), plot = last_plot(), width = 7, height = 7)

# Plot alpha across age (note lm vs. glm model)
ggplot(data = emoBias_fitParams, aes(x = ExactAge, y = alpha)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = alpha_axis_labels, colour = 'grey90') +
  scale_y_continuous(breaks = alpha_axis_labels) +
  geom_point(alpha = 1, size = 2.75, colour = emo_bias) +
  geom_smooth(alpha = .25, method = "lm", colour = emo_bias, fill = emo_bias) +
  labs(x = "Age (Years)", y = "Alpha") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/emoBias_alpha_age.png"), plot = last_plot(), width = 7, height = 7)

### EmoBias, graphing omega parameter estimates -----

# Plot omega distribution
ggplot(data = emoBias_fitParams, aes(x = omega)) +
  scale_x_continuous(breaks = alpha_axis_labels) +
  geom_hline(yintercept = seq(0, 60, by = 15), colour = 'grey90') +
  scale_y_continuous(breaks = seq(0, 60, by = 15)) + 
  geom_histogram(binwidth = .05, alpha = .75, color = emo_bias, fill = emo_bias) +
  labs(x = "Omega", y = "Number of Participants") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/emoBias_omega_dist.png"), plot = last_plot(), width = 7, height = 7)

# Plot omega across age (note lm vs. glm model)
ggplot(data = emoBias_fitParams, aes(x = ExactAge, y = omega)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = alpha_axis_labels, colour = 'grey90') +
  scale_y_continuous(breaks = alpha_axis_labels) +
  geom_point(alpha = 1, size = 2.75, colour = emo_bias) +
  geom_smooth(alpha = .25, method = "lm", colour = emo_bias, fill = emo_bias) +
  labs(x = "Age (Years)", y = "Omega") +
  emorl_theme 
ggsave(paste0(out_path, "param_fits/emoBias_omega_age.png"), plot = last_plot(), width = 7, height = 7)

### EmoBias, analyzing parameter estimates -----

# Run and evaluate beta regression
beta_mod <- lm(beta ~ ExactAge, data = emoBias_fitParams)
set.seed(123)
check_predictions(beta_mod) # Poor
qqnorm(residuals(beta_mod), pch = 1, frame = FALSE)
qqline(residuals(beta_mod), col = "red", lwd = 2) # Decent, but see if Beta distribution better fulfills model assumptions
emoBias_fitParams$beta_fract <- emoBias_fitParams$beta / 30 # Normalize (values need to fall between 0 and 1)
table(emoBias_fitParams$beta_fract) # No values at 0 or 1 bounds
beta_mod_beta <- glmmTMB(beta_fract ~ ExactAge, data = emoBias_fitParams, family = beta_family)
set.seed(123)
check_predictions(beta_mod_beta) # Good
set.seed(123)
simres <- simulateResiduals(fittedModel = beta_mod_beta, n = 250)
hist(simres$scaledResiduals, xlab = "scaled residuals", main = "Histogram Residuals") # Somewhat uniform distribution
plotQQunif(simres, testUniformity = FALSE, testOutliers = FALSE) # Noticeably better PPC --> proceed with Beta distribution

# Effects plot for beta regression
plot_model(beta_mod_beta, type = "pred", terms = "ExactAge [all]", show.data = TRUE)
beta_preds <- ggpredict(beta_mod_beta, terms = "ExactAge [all]")
ggplot(data = beta_preds, aes(x = x, y = predicted)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  scale_y_continuous(breaks = seq(.05, .25, .05)) +
  geom_line(color = emo_bias, size = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.25, fill = emo_bias) +
  labs(x = "Age (Years)", y = "Model Predictions of\nBeta Estimates") +
  emorl_theme 
ggsave(paste0(out_path, "param_fits/figS4b.png"), plot = last_plot(), width = 7, height = 7)

# Run and evaluate alpha regression
alpha_mod <- lm(alpha ~ ExactAge, data = emoBias_fitParams)
set.seed(123)
check_predictions(alpha_mod) # Poor
qqnorm(residuals(alpha_mod), pch = 1, frame = FALSE)
qqline(residuals(alpha_mod), col = "red", lwd = 2) # Poor, see if Beta distribution better fulfills model assumptions
table(emoBias_fitParams$alpha) # Lots of values at 0 and 1 bounds
emoBias_fitParams$alpha_squeeze <- betaSqueeze(emoBias_fitParams$alpha)
table(emoBias_fitParams$alpha_squeeze) # No values at 0 or 1 bounds
alpha_mod_beta <- glmmTMB(alpha_squeeze ~ ExactAge, data = emoBias_fitParams, family = beta_family)
set.seed(123)
check_predictions(alpha_mod_beta) # Decent
set.seed(123)
simres <- simulateResiduals(fittedModel = alpha_mod_beta, n = 250)
hist(simres$scaledResiduals, xlab = "scaled residuals", main = "Histogram Residuals") # Still not uniform distribution
plotQQunif(simres, testUniformity = FALSE, testOutliers = FALSE) # Noticeably better PPC --> proceed with Beta distribution

# Effects plot for alpha regression
plot_model(alpha_mod_beta, type = "pred", terms = "ExactAge [all]", show.data = TRUE)
alpha_preds <- ggpredict(alpha_mod_beta, terms = "ExactAge [all]")
ggplot(data = alpha_preds, aes(x = x, y = predicted)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  scale_y_continuous(breaks = seq(.15, .55, .1), limits = c(.15, .55)) +
  geom_line(color = emo_bias, size = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.25, fill = emo_bias) +
  labs(x = "Age (Years)", y = "Model Predictions of\nAlpha Estimates") +
  emorl_theme 
ggsave(paste0(out_path, "param_fits/figS4a.png"), plot = last_plot(), width = 7, height = 7)

# Run and evaluate omega regression
omega_mod <- lm(omega ~ ExactAge, data = emoBias_fitParams)
set.seed(123)
check_predictions(omega_mod) # Poor
qqnorm(residuals(omega_mod), pch = 1, frame = FALSE)
qqline(residuals(omega_mod), col = "red", lwd = 2) # Poor, see if Beta distribution better fulfills model assumptions
table(emoBias_fitParams$omega) # Lots of values at 0 bound and some values at 1 bound
emoBias_fitParams$omega_squeeze <- betaSqueeze(emoBias_fitParams$omega)
table(emoBias_fitParams$omega_squeeze) # No values at 0 or 1 bounds
omega_mod_beta <- glmmTMB(omega_squeeze ~ ExactAge, data = emoBias_fitParams, family = beta_family)
set.seed(123)
check_predictions(omega_mod_beta) # Decent
set.seed(123)
simres <- simulateResiduals(fittedModel = omega_mod_beta, n = 250)
hist(simres$scaledResiduals, xlab = "scaled residuals", main = "Histogram Residuals") # Still not uniform distribution
plotQQunif(simres, testUniformity = FALSE, testOutliers = FALSE) # Noticeably better PPC --> proceed with Beta distribution

# Effects plot for omega regression
plot_model(omega_mod_beta, type = "pred", terms = "ExactAge [all]", show.data = TRUE)
omega_preds <- ggpredict(omega_mod_beta, terms = "ExactAge [all]")
ggplot(data = omega_preds, aes(x = x, y = predicted)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  scale_y_continuous(breaks = seq(.25, .65, .1), limits = c(.25, .65)) +
  geom_line(color = emo_bias, size = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.25, fill = emo_bias) +
  labs(x = "Age (Years)", y = "Model Predictions of\nOmega Estimates") +
  emorl_theme 
ggsave(paste0(out_path, "param_fits/figS4c.png"), plot = last_plot(), width = 7, height = 7)

# Save regression model outputs
param_mods <- paste0(out_path, 'param_fits/param_glms.txt')
sink(param_mods)
summary(beta_mod_beta)
cat("\n\n")
summary(alpha_mod_beta)
cat("\n\n")
summary(omega_mod_beta)
sink()

### EmoBias, beta parameter comparisons -----

# Compare input and output betas
cor.test(emoBias_fullParams$beta_in, emoBias_fullParams$beta_out, method = "pearson")
ggplot(data = emoBias_fullParams, aes(x = beta_in, y = beta_out, color = factor(StudyID))) + 
  geom_point(size = 2.75) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = 'grey80') +
  geom_hline(yintercept = beta_axis_labels, colour = 'grey80') +
  geom_vline(xintercept = beta_axis_labels, colour = 'grey80') +
  scale_x_continuous(breaks = beta_axis_labels) +
  labs(x = "Randomly Sampled Input", y = "Fitted Output Estimate") +
  ggtitle("Beta") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "param_recov/emoBias_beta_cor.png"), plot = last_plot(), width = 7, height = 7)

### EmoBias, alpha parameter comparisons -----

# Compare input and output alphas
cor.test(emoBias_fullParams$alpha_in, emoBias_fullParams$alpha_out, method = "pearson")
ggplot(data = emoBias_fullParams, aes(x = alpha_in, y = alpha_out, color = factor(StudyID))) + 
  geom_point(size = 2.75) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = 'grey80') +
  geom_hline(yintercept = alpha_axis_labels, colour = 'grey80') +
  geom_vline(xintercept = alpha_axis_labels, colour = 'grey80') +
  labs(x = "Randomly Sampled Input", y = "Fitted Output Estimate") +
  ggtitle("Alpha") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "param_recov/emoBias_alpha_cor.png"), plot = last_plot(), width = 7, height = 7)

### EmoBias, omega parameter comparisons -----

# Compare input and output omegas
cor.test(emoBias_fullParams$omega_in, emoBias_fullParams$omega_out, method = "pearson")
ggplot(data = emoBias_fullParams, aes(x = omega_in, y = omega_out, color = factor(StudyID))) + 
  geom_point(size = 2.75) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = 'grey80') +
  geom_hline(yintercept = alpha_axis_labels, colour = 'grey80') +
  geom_vline(xintercept = alpha_axis_labels, colour = 'grey80') +
  labs(x = "Randomly Sampled Input", y = "Fitted Output Estimate") +
  ggtitle("Omega") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "param_recov/emoBias_omega_cor.png"), plot = last_plot(), width = 7, height = 7)

### EmoBias, all parameter comparisons -----

# Correlations between all parameters
emoBias_fullParams_Simp <- subset(emoBias_fullParams, select = -c(StudyID, lik_in, lik_out))
names(emoBias_fullParams_Simp) <- c("B_in", "a_in", "w_in", "B_out", "a_out", "w_out")
emoBias_p <- ggcorr(emoBias_fullParams_Simp, label = TRUE, label_size = 6)
emoBias_p <- emoBias_p + theme(axis.text.x = element_text(size = 12), 
                               axis.text.y = element_text(size = 12),
                               legend.text = element_text(size = 12))
emoBias_p
ggsave(paste0(out_path, "param_recov/emoBias.png"), plot = emoBias_p, width = 5, height = 4)

### Cond, graphing beta parameter estimates -----

# Plot beta distribution
ggplot(data = cond_fitParams, aes(x = beta)) +
  scale_x_continuous(breaks = beta_axis_labels) +
  geom_hline(yintercept = seq(0, 36, by = 9), colour = 'grey90') +
  scale_y_continuous(breaks = seq(0, 36, by = 9)) + 
  geom_histogram(binwidth = 1, alpha = .75, color = cond, fill = cond) +
  labs(x = "Beta", y = "Number of Participants") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/cond_beta_dist.png"), plot = last_plot(), width = 7, height = 7)

# Plot beta across age
ggplot(data = cond_fitParams, aes(x = ExactAge, y = beta)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = beta_axis_labels, colour = 'grey90') +
  scale_y_continuous(breaks = beta_axis_labels) +
  geom_point(alpha = 1, size = 2.75, colour = cond) +
  geom_smooth(alpha = .25, method = "lm", colour = cond, fill = cond) +
  labs(x = "Age (Years)", y = "Beta") +
  emorl_theme 
ggsave(paste0(out_path, "param_fits/cond_beta_age.png"), plot = last_plot(), width = 7, height = 7)

### Cond, graphing alpha H+F- parameter estimates -----

# Plot alpha H+F- distribution
ggplot(data = cond_fitParams, aes(x = alphaCC)) +
  scale_x_continuous(breaks = alpha_axis_labels) +
  geom_hline(yintercept = seq(0, 60, by = 15), colour = 'grey90') +
  scale_y_continuous(breaks = seq(0, 60, by = 15)) + 
  geom_histogram(binwidth = .05, alpha = .75, color = cond, fill = cond) +
  labs(x = "Alpha H+F-", y = "Number of Participants") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/cond_alpha_h+f-_dist.png"), plot = last_plot(), width = 7, height = 7)

# Plot alpha H+F- across age
ggplot(data = cond_fitParams, aes(x = ExactAge, y = alphaCC)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = alpha_axis_labels, colour = 'grey90') +
  scale_y_continuous(breaks = alpha_axis_labels) +
  geom_point(alpha = 1, size = 2.75, colour = cond) +
  geom_smooth(alpha = .25, method = "lm", colour = cond, fill = cond) +
  labs(x = "Age (Years)", y = "Alpha H+F-") +
  emorl_theme 
ggsave(paste0(out_path, "param_fits/cond_alpha_h+f-_age.png"), plot = last_plot(), width = 7, height = 7)

### Cond, graphing alpha H+H- parameter estimates -----

# Plot alpha H+H- distribution
ggplot(data = cond_fitParams, aes(x = alphaCI)) +
  scale_x_continuous(breaks = alpha_axis_labels) +
  geom_hline(yintercept = seq(0, 60, by = 15), colour = 'grey90') +
  scale_y_continuous(breaks = seq(0, 60, by = 15)) + 
  geom_histogram(binwidth = .05, alpha = .75, color = cond, fill = cond) +
  labs(x = "Alpha H+H-", y = "Number of Participants") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/cond_alpha_h+h-_dist.png"), plot = last_plot(), width = 7, height = 7)

# Plot alpha H+H- across age
ggplot(data = cond_fitParams, aes(x = ExactAge, y = alphaCI)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = alpha_axis_labels, colour = 'grey90') +
  scale_y_continuous(breaks = alpha_axis_labels) +
  geom_point(alpha = 1, size = 2.75, colour = cond) +
  geom_smooth(alpha = .25, method = "lm", colour = cond, fill = cond) +
  labs(x = "Age (Years)", y = "Alpha H+H-") +
  emorl_theme 
ggsave(paste0(out_path, "param_fits/cond_alpha_h+h-_age.png"), plot = last_plot(), width = 7, height = 7)

### Cond, graphing alpha F+F- parameter estimates -----

# Plot alpha F+F- distribution
ggplot(data = cond_fitParams, aes(x = alphaIC)) +
  scale_x_continuous(breaks = alpha_axis_labels) +
  geom_hline(yintercept = seq(0, 60, by = 15), colour = 'grey90') +
  scale_y_continuous(breaks = seq(0, 60, by = 15)) + 
  geom_histogram(binwidth = .05, alpha = .75, color = cond, fill = cond) +
  labs(x = "Alpha F+F-", y = "Number of Participants") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/cond_alpha_f+f-_dist.png"), plot = last_plot(), width = 7, height = 7)

# Plot alpha F+F- across age
ggplot(data = cond_fitParams, aes(x = ExactAge, y = alphaIC)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = alpha_axis_labels, colour = 'grey90') +
  scale_y_continuous(breaks = alpha_axis_labels) +
  geom_point(alpha = 1, size = 2.75, colour = cond) +
  geom_smooth(alpha = .25, method = "lm", colour = cond, fill = cond) +
  labs(x = "Age (Years)", y = "Alpha F+F-") +
  emorl_theme 
ggsave(paste0(out_path, "param_fits/cond_alpha_f+f-_age.png"), plot = last_plot(), width = 7, height = 7)

### Cond, graphing alpha F+H- parameter estimates -----

# Plot alpha F+H- distribution
ggplot(data = cond_fitParams, aes(x = alphaII)) +
  scale_x_continuous(breaks = alpha_axis_labels) +
  geom_hline(yintercept = seq(0, 60, by = 15), colour = 'grey90') +
  scale_y_continuous(breaks = seq(0, 60, by = 15)) + 
  geom_histogram(binwidth = .05, alpha = .75, color = cond, fill = cond) +
  labs(x = "Alpha F+H-", y = "Number of Participants") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/cond_alpha_f+h-_dist.png"), plot = last_plot(), width = 7, height = 7)

# Plot alpha F+H- across age
ggplot(data = cond_fitParams, aes(x = ExactAge, y = alphaII)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = alpha_axis_labels, colour = 'grey90') +
  scale_y_continuous(breaks = alpha_axis_labels) +
  geom_point(alpha = 1, size = 2.75, colour = cond) +
  geom_smooth(alpha = .25, method = "lm", colour = cond, fill = cond) +
  labs(x = "Age (Years)", y = "Alpha F+H-") +
  emorl_theme 
ggsave(paste0(out_path, "param_fits/cond_alpha_f+h-_age.png"), plot = last_plot(), width = 7, height = 7)

### Cond, beta parameter comparisons -----

# Compare input and output betas
cor.test(cond_fullParams$beta_in, cond_fullParams$beta_out, method = "pearson")
ggplot(data = cond_fullParams, aes(x = beta_in, y = beta_out, color = factor(StudyID))) + 
  geom_point(size = 2.75) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = 'grey80') +
  geom_hline(yintercept = beta_axis_labels, colour = 'grey80') +
  geom_vline(xintercept = beta_axis_labels, colour = 'grey80') +
  scale_x_continuous(breaks = beta_axis_labels) +
  labs(x = "Randomly Sampled Input", y = "Fitted Output Estimate") +
  ggtitle("Beta") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "param_recov/cond_beta_cor.png"), plot = last_plot(), width = 7, height = 7)

### Cond, alpha parameter comparisons -----

# Compare input and output alphas
cor.test(cond_fullParams$alpha_cc_in, cond_fullParams$alpha_cc_out, method = "pearson")
cor.test(cond_fullParams$alpha_ci_in, cond_fullParams$alpha_ci_out, method = "pearson")
cor.test(cond_fullParams$alpha_ic_in, cond_fullParams$alpha_ic_out, method = "pearson")
cor.test(cond_fullParams$alpha_ii_in, cond_fullParams$alpha_ii_out, method = "pearson")
ggplot(data = cond_fullParams, aes(x = alpha_cc_in, y = alpha_cc_out, color = factor(StudyID))) + 
  geom_point(size = 2.75) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = 'grey80') +
  geom_hline(yintercept = alpha_axis_labels, colour = 'grey80') +
  geom_vline(xintercept = alpha_axis_labels, colour = 'grey80') +
  labs(x = "Randomly Sampled Input", y = "Fitted Output Estimate") +
  ggtitle("Alpha H+F-") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "param_recov/cond_alpha_h+f-_cor.png"), plot = last_plot(), width = 7, height = 7)
ggplot(data = cond_fullParams, aes(x = alpha_ci_in, y = alpha_ci_out, color = factor(StudyID))) + 
  geom_point(size = 2.75) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = 'grey80') +
  geom_hline(yintercept = alpha_axis_labels, colour = 'grey80') +
  geom_vline(xintercept = alpha_axis_labels, colour = 'grey80') +
  labs(x = "Randomly Sampled Input", y = "Fitted Output Estimate") +
  ggtitle("Alpha H+H-") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "param_recov/cond_alpha_h+h-_cor.png"), plot = last_plot(), width = 7, height = 7)
ggplot(data = cond_fullParams, aes(x = alpha_ic_in, y = alpha_ic_out, color = factor(StudyID))) + 
  geom_point(size = 2.75) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = 'grey80') +
  geom_hline(yintercept = alpha_axis_labels, colour = 'grey80') +
  geom_vline(xintercept = alpha_axis_labels, colour = 'grey80') +
  labs(x = "Randomly Sampled Input", y = "Fitted Output Estimate") +
  ggtitle("Alpha F+F-") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "param_recov/cond_alpha_f+f-_cor.png"), plot = last_plot(), width = 7, height = 7)
ggplot(data = cond_fullParams, aes(x = alpha_ii_in, y = alpha_ii_out, color = factor(StudyID))) + 
  geom_point(size = 2.75) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = 'grey80') +
  geom_hline(yintercept = alpha_axis_labels, colour = 'grey80') +
  geom_vline(xintercept = alpha_axis_labels, colour = 'grey80') +
  labs(x = "Randomly Sampled Input", y = "Fitted Output Estimate") +
  ggtitle("Alpha F+H-") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "param_recov/cond_alpha_f+h-_cor.png"), plot = last_plot(), width = 7, height = 7)

### Cond, all parameter comparisons -----

# Correlations between all parameters
cond_fullParams_Simp <- subset(cond_fullParams, select = -c(StudyID, lik_in, lik_out))
names(cond_fullParams_Simp) <- c("B_in", "a_1_in", "a_2_in", "a_3_in", "a_4_in", 
                                 "B_out", "a_1_out", "a_2_out", "a_3_out", "a_4_out")
cond_p <- ggcorr(cond_fullParams_Simp, label = TRUE, label_size = 6)
cond_p <- cond_p + theme(axis.text.x = element_text(size = 12), 
                         axis.text.y = element_text(size = 12),
                         legend.text = element_text(size = 12))
cond_p
ggsave(paste0(out_path, "param_recov/cond.png"), plot = cond_p, width = 10, height = 8)

### CondBeta, graphing beta H+F- parameter estimates -----

# Plot beta H+F- distribution
ggplot(data = condBeta_fitParams, aes(x = betaCC)) +
  scale_x_continuous(breaks = beta_axis_labels) +
  geom_hline(yintercept = seq(0, 36, by = 9), colour = 'grey90') +
  scale_y_continuous(breaks = seq(0, 36, by = 9)) + 
  geom_histogram(binwidth = 1, alpha = .75, color = cond_beta, fill = cond_beta) +
  labs(x = "Beta H+F-", y = "Number of Participants") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/condBeta_beta_h+f-_dist.png"), plot = last_plot(), width = 7, height = 7)

# Plot beta H+F- across age
ggplot(data = condBeta_fitParams, aes(x = ExactAge, y = betaCC)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = beta_axis_labels, colour = 'grey90') +
  scale_y_continuous(breaks = beta_axis_labels) +
  geom_point(alpha = 1, size = 2.75, colour = cond_beta) +
  geom_smooth(alpha = .25, method = "lm", colour = cond_beta, fill = cond_beta) +
  labs(x = "Age (Years)", y = "Beta H+F-") +
  emorl_theme 
ggsave(paste0(out_path, "param_fits/condBeta_beta_h+f-_age.png"), plot = last_plot(), width = 7, height = 7)

### CondBeta, graphing beta H+H- parameter estimates -----

# Plot beta H+H- distribution
ggplot(data = condBeta_fitParams, aes(x = betaCI)) +
  scale_x_continuous(breaks = beta_axis_labels) +
  geom_hline(yintercept = seq(0, 36, by = 9), colour = 'grey90') +
  scale_y_continuous(breaks = seq(0, 36, by = 9)) + 
  geom_histogram(binwidth = 1, alpha = .75, color = cond_beta, fill = cond_beta) +
  labs(x = "Beta H+H-", y = "Number of Participants") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/condBeta_beta_h+h-_dist.png"), plot = last_plot(), width = 7, height = 7)

# Plot beta H+H- across age
ggplot(data = condBeta_fitParams, aes(x = ExactAge, y = betaCI)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = beta_axis_labels, colour = 'grey90') +
  scale_y_continuous(breaks = beta_axis_labels) +
  geom_point(alpha = 1, size = 2.75, colour = cond_beta) +
  geom_smooth(alpha = .25, method = "lm", colour = cond_beta, fill = cond_beta) +
  labs(x = "Age (Years)", y = "Beta H+H-") +
  emorl_theme 
ggsave(paste0(out_path, "param_fits/condBeta_beta_h+h-_age.png"), plot = last_plot(), width = 7, height = 7)

### CondBeta, graphing beta F+F- parameter estimates -----

# Plot beta F+F- distribution
ggplot(data = condBeta_fitParams, aes(x = betaIC)) +
  scale_x_continuous(breaks = beta_axis_labels) +
  geom_hline(yintercept = seq(0, 36, by = 9), colour = 'grey90') +
  scale_y_continuous(breaks = seq(0, 36, by = 9)) + 
  geom_histogram(binwidth = 1, alpha = .75, color = cond_beta, fill = cond_beta) +
  labs(x = "Beta F+F-", y = "Number of Participants") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/condBeta_beta_f+f-_dist.png"), plot = last_plot(), width = 7, height = 7)

# Plot beta F+F- across age
ggplot(data = condBeta_fitParams, aes(x = ExactAge, y = betaIC)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = beta_axis_labels, colour = 'grey90') +
  scale_y_continuous(breaks = beta_axis_labels) +
  geom_point(alpha = 1, size = 2.75, colour = cond_beta) +
  geom_smooth(alpha = .25, method = "lm", colour = cond_beta, fill = cond_beta) +
  labs(x = "Age (Years)", y = "Beta F+F-") +
  emorl_theme 
ggsave(paste0(out_path, "param_fits/condBeta_beta_f+f-_age.png"), plot = last_plot(), width = 7, height = 7)

### CondBeta, graphing beta F+H- parameter estimates -----

# Plot beta F+H- distribution
ggplot(data = condBeta_fitParams, aes(x = betaII)) +
  scale_x_continuous(breaks = beta_axis_labels) +
  geom_hline(yintercept = seq(0, 36, by = 9), colour = 'grey90') +
  scale_y_continuous(breaks = seq(0, 36, by = 9)) + 
  geom_histogram(binwidth = 1, alpha = .75, color = cond_beta, fill = cond_beta) +
  labs(x = "Beta F+H-", y = "Number of Participants") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/condBeta_beta_f+h-_dist.png"), plot = last_plot(), width = 7, height = 7)

# Plot beta F+H- across age
ggplot(data = condBeta_fitParams, aes(x = ExactAge, y = betaII)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = beta_axis_labels, colour = 'grey90') +
  scale_y_continuous(breaks = beta_axis_labels) +
  geom_point(alpha = 1, size = 2.75, colour = cond_beta) +
  geom_smooth(alpha = .25, method = "lm", colour = cond_beta, fill = cond_beta) +
  labs(x = "Age (Years)", y = "Beta F+H-") +
  emorl_theme 
ggsave(paste0(out_path, "param_fits/condBeta_beta_f+h-_age.png"), plot = last_plot(), width = 7, height = 7)

### CondBeta, graphing alpha parameter estimates -----

# Plot alpha distribution
ggplot(data = condBeta_fitParams, aes(x = alpha)) +
  scale_x_continuous(breaks = alpha_axis_labels) +
  geom_hline(yintercept = seq(0, 60, by = 15), colour = 'grey90') +
  scale_y_continuous(breaks = seq(0, 60, by = 15)) + 
  geom_histogram(binwidth = .05, alpha = .75, color = cond_beta, fill = cond_beta) +
  labs(x = "Alpha", y = "Number of Participants") +
  emorl_theme
ggsave(paste0(out_path, "param_fits/condBeta_alpha_dist.png"), plot = last_plot(), width = 7, height = 7)

# Plot alpha across age
ggplot(data = condBeta_fitParams, aes(x = ExactAge, y = alpha)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = alpha_axis_labels, colour = 'grey90') +
  scale_y_continuous(breaks = alpha_axis_labels) +
  geom_point(alpha = 1, size = 2.75, colour = cond_beta) +
  geom_smooth(alpha = .25, method = "lm", colour = cond_beta, fill = cond_beta) +
  labs(x = "Age (Years)", y = "Alpha") +
  emorl_theme 
ggsave(paste0(out_path, "param_fits/condBeta_alpha_age.png"), plot = last_plot(), width = 7, height = 7)

### CondBeta, beta parameter comparisons -----

# Compare input and output betas
cor.test(condBeta_fullParams$beta_cc_in, condBeta_fullParams$beta_cc_out, method = "pearson")
cor.test(condBeta_fullParams$beta_ci_in, condBeta_fullParams$beta_ci_out, method = "pearson")
cor.test(condBeta_fullParams$beta_ic_in, condBeta_fullParams$beta_ic_out, method = "pearson")
cor.test(condBeta_fullParams$beta_ii_in, condBeta_fullParams$beta_ii_out, method = "pearson")
ggplot(data = condBeta_fullParams, aes(x = beta_cc_in, y = beta_cc_out, color = factor(StudyID))) + 
  geom_point(size = 2.75) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = 'grey80') +
  geom_hline(yintercept = beta_axis_labels, colour = 'grey80') +
  geom_vline(xintercept = beta_axis_labels, colour = 'grey80') +
  scale_x_continuous(breaks = beta_axis_labels) +
  labs(x = "Randomly Sampled Input", y = "Fitted Output Estimate") +
  ggtitle("Beta H+F-") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "param_recov/condBeta_beta_h+f-_cor.png"), plot = last_plot(), width = 7, height = 7)
ggplot(data = condBeta_fullParams, aes(x = beta_ci_in, y = beta_ci_out, color = factor(StudyID))) + 
  geom_point(size = 2.75) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = 'grey80') +
  geom_hline(yintercept = beta_axis_labels, colour = 'grey80') +
  geom_vline(xintercept = beta_axis_labels, colour = 'grey80') +
  scale_x_continuous(breaks = beta_axis_labels) +
  labs(x = "Randomly Sampled Input", y = "Fitted Output Estimate") +
  ggtitle("Beta H+H-") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "param_recov/condBeta_beta_h+h-_cor.png"), plot = last_plot(), width = 7, height = 7)
ggplot(data = condBeta_fullParams, aes(x = beta_ic_in, y = beta_ic_out, color = factor(StudyID))) + 
  geom_point(size = 2.75) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = 'grey80') +
  geom_hline(yintercept = beta_axis_labels, colour = 'grey80') +
  geom_vline(xintercept = beta_axis_labels, colour = 'grey80') +
  scale_x_continuous(breaks = beta_axis_labels) +
  labs(x = "Randomly Sampled Input", y = "Fitted Output Estimate") +
  ggtitle("Beta F+F-") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "param_recov/condBeta_beta_f+f-_cor.png"), plot = last_plot(), width = 7, height = 7)
ggplot(data = condBeta_fullParams, aes(x = beta_ii_in, y = beta_ii_out, color = factor(StudyID))) + 
  geom_point(size = 2.75) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = 'grey80') +
  geom_hline(yintercept = beta_axis_labels, colour = 'grey80') +
  geom_vline(xintercept = beta_axis_labels, colour = 'grey80') +
  scale_x_continuous(breaks = beta_axis_labels) +
  labs(x = "Randomly Sampled Input", y = "Fitted Output Estimate") +
  ggtitle("Beta F+H-") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "param_recov/condBeta_beta_f+h-_cor.png"), plot = last_plot(), width = 7, height = 7)

### CondBeta, alpha parameter comparisons -----

# Compare input and output alphas
cor.test(condBeta_fullParams$alpha_in, condBeta_fullParams$alpha_out, method = "pearson")
ggplot(data = condBeta_fullParams, aes(x = alpha_in, y = alpha_out, color = factor(StudyID))) + 
  geom_point(size = 2.75) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = 'grey80') +
  geom_hline(yintercept = alpha_axis_labels, colour = 'grey80') +
  geom_vline(xintercept = alpha_axis_labels, colour = 'grey80') +
  labs(x = "Randomly Sampled Input", y = "Fitted Output Estimate") +
  ggtitle("Alpha") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(out_path, "param_recov/condBeta_alpha_cor.png"), plot = last_plot(), width = 7, height = 7)

### CondBeta, all parameter comparisons -----

# Correlations between all parameters
condBeta_fullParams_Simp <- subset(condBeta_fullParams, select = -c(StudyID, lik_in, lik_out))
names(condBeta_fullParams_Simp) <- c("B_1_in", "B_2_in", "B_3_in", "B_4_in", "a_in",
                                     "B_1_out", "B_2_out", "B_3_out", "B_4_out", "a_out")
condBeta_p <- ggcorr(condBeta_fullParams_Simp, label = TRUE, label_size = 6)
condBeta_p <- condBeta_p + theme(axis.text.x = element_text(size = 12), 
                                 axis.text.y = element_text(size = 12),
                                 legend.text = element_text(size = 12))
condBeta_p
ggsave(paste0(out_path, "param_recov/condBeta.png"), plot = condBeta_p, width = 10, height = 8)
