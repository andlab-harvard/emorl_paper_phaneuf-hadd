##################################################
### Analyze EmoRL Data with Statistical Models ###
##################################################

# Written by: Camille Phaneuf-Hadd (cphaneuf@g.harvard.edu)
# Last updated: 4/30/26

# Inputs: processed learning task, test phase, and scale data
# Computes: visualizations and statistical assessments
# Outputs: png plots and text files into results/learn, test, and scales

### Set up Script -----

# Load needed libraries
require(pacman) # for p_load()
p_load(tidyverse, # for df manipulation
       dplyr, # for %>% and other operators
       plotrix, # for std.error()
       ggplot2, # for plotting
       sjPlot, # for plot_model()
       lmerTest, # for mixed effects models
       performance, # for check_predictions()
       car, # for Anova()
       ggeffects, # for ggpredict()
       optimx, # for optimizers
       emmeans, # for emmeans()
       glmmTMB, # for glmmTMB()
       DHARMa) # for simulateResiduals()

# Load shared EmoRL functions
source("utilities.R")

# Set paths
in_path <- '../data/behav/'
learn_out_path <- '../results/learn/'
test_out_path <- '../results/test/'
scales_out_path <- '../results/scales/'

# Read in task data
learn <- read.csv(paste0(in_path, "learn.csv"))
test <- read.csv(paste0(in_path, "test.csv"))
avg_test <- read.csv(paste0(in_path, "avg_test.csv"))
avg_test_rt <- read.csv(paste0(in_path, "avg_test_rt.csv"))

# Read in scale data
subj_emo_val <- read.csv(paste0(in_path, "subj_emo_val.csv"))
rein_rep <- read.csv(paste0(in_path, "rein_rep.csv"))
subj_val_money <- read.csv(paste0(in_path, "subj_val_money.csv"), check.names = FALSE)

### Final *learning task* data processing -----

# Change variable types
str(learn)
learn$StudyID <- as.factor(learn$StudyID)
learn$condition <- factor(learn$condition, levels = c("H+F-", "H+H-", "F+F-", "F+H-"))
str(learn)

### Final *test phase* frequency data processing -----

# Pivot avg_test longer
avg_test_long <- pivot_longer(avg_test, cols = grep("rate_", colnames(avg_test)), names_to = "deck_type")

# Add reward and emotion labels
reward_labels <- colnames(avg_test)[grepl("Reward", colnames(avg_test))]
zero_labels <- colnames(avg_test)[grepl("Zero", colnames(avg_test))]
avg_test_long <- avg_test_long %>% mutate(Reward = case_when(deck_type %in% reward_labels ~ "75% reinforced",
                                                             deck_type %in% zero_labels ~ "25% reinforced"))
happy_labels <- colnames(avg_test)[grepl("Happy", colnames(avg_test))]
fear_labels <- colnames(avg_test)[grepl("Fear", colnames(avg_test))]
avg_test_long <- avg_test_long %>% mutate(Emotion = case_when(deck_type %in% happy_labels ~ "Happy",
                                                              deck_type %in% fear_labels ~ "Fear"))

# Change variable types
str(avg_test_long)
avg_test_long$StudyID <- as.factor(avg_test_long$StudyID)
avg_test_long$Reward <- factor(avg_test_long$Reward, levels = c("75% reinforced", "25% reinforced"))
avg_test_long$Emotion <- factor(avg_test_long$Emotion, levels = c("Happy", "Fear"))
str(avg_test_long)

### Final *test phase* RT data processing -----

# Pivot avg_test longer
avg_test_rt_long <- pivot_longer(avg_test_rt, cols = grep("rate_", colnames(avg_test_rt)), names_to = "deck_type")

# Add reward and emotion labels
reward_labels_rt <- colnames(avg_test_rt)[grepl("Reward", colnames(avg_test_rt))]
zero_labels_rt <- colnames(avg_test_rt)[grepl("Zero", colnames(avg_test_rt))]
avg_test_rt_long <- avg_test_rt_long %>% mutate(Reward = case_when(deck_type %in% reward_labels_rt ~ "75% reinforced",
                                                                   deck_type %in% zero_labels_rt ~ "25% reinforced"))
happy_labels_rt <- colnames(avg_test_rt)[grepl("Happy", colnames(avg_test_rt))]
fear_labels_rt <- colnames(avg_test_rt)[grepl("Fear", colnames(avg_test_rt))]
avg_test_rt_long <- avg_test_rt_long %>% mutate(Emotion = case_when(deck_type %in% happy_labels_rt ~ "Happy",
                                                                    deck_type %in% fear_labels_rt ~ "Fear"))

# Change variable types
str(avg_test_rt_long)
avg_test_rt_long$StudyID <- as.factor(avg_test_rt_long$StudyID)
avg_test_rt_long$Reward <- factor(avg_test_rt_long$Reward, levels = c("75% reinforced", "25% reinforced"))
avg_test_rt_long$Emotion <- factor(avg_test_rt_long$Emotion, levels = c("Happy", "Fear"))
str(avg_test_rt_long)

# Count how many times a deck was never chosen
no_choice1 <- avg_test_long[(avg_test_long$value == 0), ]
table(no_choice1$deck_type) # Matches below
no_choice2 <- avg_test_rt_long[is.na(avg_test_rt_long$value), ]
table(no_choice2$deck_type) # Matches above

### Final *test phase* (learning pairs) similarity data processing -----

# Subset test into data frames of each pair_class
from_learn <- test[test$pair_class == "from_learn", ]
all_else <- test[test$pair_class != "from_learn", ]

# Change variable types
str(from_learn)
from_learn$StudyID <- as.factor(from_learn$StudyID)
from_learn$pair_type_simp <- factor(from_learn$pair_type_simp, levels = c("H+F-", "H+H-", "F+F-", "F+H-"))
str(from_learn)

# Make summary data frames for from_learn
overall_from_learn_rt <- from_learn %>%
  group_by(StudyID, pair_type_simp) %>%
  dplyr::summarise(ind_avg_rt = mean(rt, na.rm = TRUE), ind_rt_se = std.error(rt, na.rm = TRUE), n = n())

### Final *test phase* (novel pairs) similarity data processing -----

# Change variable types
str(all_else)
all_else$StudyID <- as.factor(all_else$StudyID)
all_else$pair_class <- factor(all_else$pair_class, levels = c("same_rew_same_emo",
                                                              "same_rew_diff_emo",
                                                              "diff_rew_same_emo",
                                                              "diff_rew_diff_emo"))
str(all_else)

### Final *scale* data processing -----

# Add emotion labels
subj_emo_val <- subj_emo_val %>% mutate(emotion = case_when(grepl("Happy", face) ~ "Happy",
                                                            grepl("Fear", face) ~ "Fear"))
rein_rep <- rein_rep %>% mutate(emotion = case_when(grepl("Happy", face) ~ "Happy",
                                                    grepl("Fear", face) ~ "Fear"))

# Change variable types -- subjective ratings
str(subj_emo_val)
subj_emo_val$StudyID <- as.factor(subj_emo_val$StudyID)
subj_emo_val$reward <- factor(subj_emo_val$reward, levels = c("75% reinforced", "25% reinforced"))
subj_emo_val$emotion <- factor(subj_emo_val$emotion, levels = c("Happy", "Fear"))
str(subj_emo_val)

# Change variable types -- reinforcement reports
str(rein_rep)
rein_rep$StudyID <- as.factor(rein_rep$StudyID)
rein_rep$reward <- factor(rein_rep$reward, levels = c("75% reinforced", "25% reinforced"))
rein_rep$emotion <- factor(rein_rep$emotion, levels = c("Happy", "Fear"))
str(rein_rep)

# Pivot subj_val_money longer
amt_cols <- grep("[$]", names(subj_val_money))
subj_val_money_long <- subj_val_money %>% pivot_longer(cols = amt_cols, names_to = "Amount", values_to = "Report")

# Change variable types
str(subj_val_money_long)
subj_val_money_long$StudyID <- as.factor(subj_val_money_long$StudyID)
subj_val_money_long$Amount <- factor(subj_val_money_long$Amount, levels = c("$0.00", "$0.25", "$0.50", "$0.75", "$1.00"))
subj_val_money_long$Report <- as.numeric(subj_val_money_long$Report)
str(subj_val_money_long)

### Learning task analyses ----- 

# Are there age-related changes in the use of reward and emotion information to guide
# learning across time?

##### Accuracy -----

# Check for appropriate model distribution and remove NA values for Bernoulli
hist(learn$acc) # Bernoulli distribution is the correct choice
unique(learn$acc)
learn_noNA <- learn[!is.na(learn$acc), ]
unique(learn_noNA$acc)

# Run and evaluate learning accuracy model
# learn_acc_mod <- glmer(acc ~ condition * instance * ExactAge + (1 | StudyID), data = learn_noNA, family = binomial) # Binomial distribution does not converge; re-fit with new optimizer
# learn_acc_mod <- glmer(acc ~ condition * instance * ExactAge + (1 | StudyID), data = learn_data_noNA, family = binomial, 
#                                                                               control = glmerControl(optimizer = c("bobyqa"))) # Does not do any better, and the outputs are the same, so stick with the simpler model
# learn_acc_mod <- glmer(acc ~ condition * instance * ExactAge + (1 | StudyID), data = learn_data_noNA, family = binomial, 
#                                                                               control = glmerControl(optimizer = c("Nelder_Mead"))) # Does not do any better, and the outputs are the same, so stick with the simpler model
learn_noNA <- learn_noNA %>% mutate(ExactAge_z = scale(ExactAge), ExactAge_z2 = ExactAge_z^2, instance_z = scale(instance)) # Try and fix convergence issues by scaling variables
learn_acc_mod <- glmer(acc ~ condition * instance_z * ExactAge_z + (1 | StudyID), data = learn_noNA, family = binomial) # Success! No convergence issues
set.seed(123)
check_predictions(learn_acc_mod) # Great!

# Explore whether quadratic age model is better fitting
learn_acc_mod2 <- glmer(acc ~ condition * instance_z * (ExactAge_z + ExactAge_z2) + (1 | StudyID), data = learn_noNA, family = binomial)
anova(learn_acc_mod, learn_acc_mod2) # p = 0.1222
# No significant difference --> continue with more parsimonious model

# Inspect learn_acc_mod's significant main effects and interactions
Anova(learn_acc_mod, type = "II")
plot_model(learn_acc_mod, type = "pred", terms = "condition")
plot_model(learn_acc_mod, type = "pred", terms = "instance_z [all]")
plot_model(learn_acc_mod, type = "pred", terms = "ExactAge_z [all]")
plot_model(learn_acc_mod, type = "pred", terms = c("instance_z [all]", "condition"))
plot_model(learn_acc_mod, type = "pred", terms = c("instance_z [all]", "ExactAge_z"))

# Figure out time_pts for post-hoc tests
t1 <- learn_noNA[learn_noNA$instance == 15, ]
unique(t1$instance_z) # -0.6159607
t2 <- learn_noNA[learn_noNA$instance == 30, ]
unique(t2$instance_z) # 0.5389722
t3 <- learn_noNA[learn_noNA$instance == 45, ]
unique(t3$instance_z) # 1.693905

# Conduct post-hoc tests of instance x condition interaction
time_pts <- c(-0.6159607, 0.5389722, 1.693905)
em_inst_pair_acc <- emmeans(learn_acc_mod, specs = c("condition"), by = c("instance_z"), at = list(instance_z = time_pts))
em_inst_pair_acc # Marginal means for each point of interest                       
comp_em_inst_pair_acc <- pairs(em_inst_pair_acc, adjust = "holm")
comp_em_inst_pair_acc 
# After correcting for multiple comparisons:
# --> H+F- acc is significantly greater than acc for other conditions at instance = 15
# --> H+H- and F+F- acc is also significantly greater than F+H- acc at instance = 15
# --> H+F-, H+H-, and F+F- acc is still significantly greater than F+H- acc at instance = 30
# --> H+F- acc is still significantly greater than F+H- acc at instance = 45
# --> H+H- and F+F- acc never differ

# Save learn_acc_mod's outputs
learn_acc_mod_filename <- paste0(learn_out_path, 'learn_acc_mod.txt')
sink(learn_acc_mod_filename)
print(Anova(learn_acc_mod, type = "II"))
cat("\nEMMEANS FOR GLMER ACCURACY MODEL\n\n")
cat("Time points:\n")
print(time_pts)
cat("\nHolm-Corrected P-Values\n\n")
comp_em_inst_pair_acc 
sink()

# Get values for transforming plot axes from z-scored to original scale
inst_mean <- mean(learn_noNA$instance, na.rm = TRUE)
inst_sd <- sd(learn_noNA$instance, na.rm = TRUE)
age_mean <- mean(learn_noNA$ExactAge, na.rm = TRUE)
age_sd <- sd(learn_noNA$ExactAge, na.rm = TRUE)

# Make plots for learn_acc_mod's significant main effects and interactions...

# ...condition main effect
condition_main <- ggpredict(learn_acc_mod, terms = "condition")
ggplot(data = condition_main, aes(x = x, y = predicted, color = x, fill = x)) +
  ylim(.55, .75) +
  geom_point(alpha = 1, size = 5.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = .2, size = 2) +
  scale_color_manual(values = c(cc_tomato, ci_lilac, ic_water, ii_sage)) +
  scale_fill_manual(values = c(cc_tomato, ci_lilac, ic_water, ii_sage)) +
  labs(x = "Condition", y = "Model Predictions of Accuracy") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(learn_out_path, 'condition_main.png'), width = 4.5, height = 7)

# ...instance main effect
instance_main <- ggpredict(learn_acc_mod, terms = c("instance_z [all]"))
ggplot(data = instance_main, aes(x = x, y = predicted)) +
  ylim(.5, .9) +
  scale_x_continuous(labels = function(x) round(x * inst_sd + inst_mean), breaks = time_pts) +
  geom_line(color = cc_tomato, size = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.25, fill = cc_tomato) +
  labs(x = "Trial Number", y = "Model Predictions of Accuracy") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(learn_out_path, 'instance_main.png'), width = 4.5, height = 7)

# ...age main effect
age_main <- ggpredict(learn_acc_mod, terms = c("ExactAge_z [all]"))
ggplot(data = age_main, aes(x = x, y = predicted)) +
  ylim(.5, .9) +
  scale_x_continuous(labels = function(x) round(x * age_sd + age_mean), 
                     breaks = c((8 - age_mean) / age_sd, (12 - age_mean) / age_sd, (16 - age_mean) / age_sd, (20 - age_mean) / age_sd)) +
  geom_line(color = cc_tomato, size = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.25, fill = cc_tomato) +
  labs(x = "Age (Years)", y = "Model Predictions of Accuracy") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(learn_out_path, 'age_main.png'), width = 4.5, height = 7)

# ...instance x condition interaction
instance_condition_interact <- ggpredict(learn_acc_mod, terms = c("instance_z [all]", "condition"))
ggplot(instance_condition_interact, aes(x, predicted)) +
  ylim(.4, .9) +
  scale_x_continuous(labels = function(x) round(x * inst_sd + inst_mean), breaks = time_pts) +
  geom_line(aes(color = group), size = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  geom_vline(xintercept = time_pts, linetype = 'dashed', size = 2) +
  scale_color_manual(values = c(cc_tomato, ci_lilac, ic_water, ii_sage)) +
  scale_fill_manual(values = c(cc_tomato, ci_lilac, ic_water, ii_sage)) +
  labs(x = "Trial Number", y = "Model Predictions of Accuracy", color = "Condition", fill = "Condition") +
  emorl_theme
ggsave(paste0(learn_out_path, 'fig2b.png'), width = 7, height = 7)

# ...instance x age interaction
instance_age_interact <- ggpredict(learn_acc_mod, terms = c("instance_z [all]", "ExactAge_z"))
instance_age_interact$group <- as.numeric(as.character(instance_age_interact$group))
instance_age_interact$group2 <- round(instance_age_interact$group * age_sd + age_mean, 1)
instance_age_interact$group2 <- as.character(instance_age_interact$group2)
instance_age_interact$group2 <- gsub("20", "20.0", instance_age_interact$group2)
instance_age_interact$group2 <- as.factor(instance_age_interact$group2)
ggplot(instance_age_interact, aes(x, predicted)) +
  ylim(.4, .9) +
  scale_x_continuous(labels = function(x) round(x * inst_sd + inst_mean), breaks = time_pts) +
  geom_line(aes(color = group2), size = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group2), alpha = 0.25) +
  scale_color_manual(values = c(fem_tomato, masc_tomato, not_capt_tomato)) +
  scale_fill_manual(values = c(fem_tomato, masc_tomato, not_capt_tomato)) +
  labs(x = "Trial Number", y = "Model Predictions of Accuracy", color = "Age (Years)", fill = "Age (Years)") +
  emorl_theme
ggsave(paste0(learn_out_path, 'fig2a.png'), width = 7, height = 7)

### Test phase analyses ----- 

# How does reward and emotion information persist to bias choice across age?

##### Frequency -----

# Check for appropriate model distribution
hist(avg_test_long$value) # Gaussian distribution is the correct choice

# Run and evaluate test frequency ("accuracy") model
test_acc_mod <- lmer(value ~ Reward * Emotion * ExactAge + (1 | StudyID), data = avg_test_long)
set.seed(123)
check_predictions(test_acc_mod) # Fine
qqnorm(residuals(test_acc_mod), pch = 1, frame = FALSE)
qqline(residuals(test_acc_mod), col = "red", lwd = 2) # Great!

# Inspect test_acc_mod's significant main effects and interactions
Anova(test_acc_mod, type = "II")
plot_model(test_acc_mod, type = "pred", terms = "Reward")
plot_model(test_acc_mod, type = "pred", terms = "Emotion")
plot_model(test_acc_mod, type = "pred", terms = c("ExactAge [all]", "Reward"))

# Save test_acc_mod's outputs
test_acc_mod_filename <- paste0(test_out_path, 'acc/test_acc_mod.txt')
sink(test_acc_mod_filename)
print(Anova(test_acc_mod, type = "II"))
sink()

# Make plots for test_acc_mod's significant main effects and interactions...

# ...reward main effect
rew_main <- ggpredict(test_acc_mod, terms = "Reward")
ggplot(data = rew_main, aes(x = x, y = predicted, color = x, fill = x)) +
  ylim(.4, .7) +
  geom_point(alpha = 1, size = 5.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = .2, size = 2) +
  scale_color_manual(values = c(gold, silver)) +
  scale_fill_manual(values = c(gold, silver)) +
  labs(x = "Reward", y = "Model Predictions of Choice") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(test_out_path, 'acc/rew_acc_main.png'), width = 7, height = 7)

# ...emotion main effect
emo_main <- ggpredict(test_acc_mod, terms = "Emotion")
ggplot(data = emo_main, aes(x = x, y = predicted, color = x, fill = x)) +
  ylim(.4, .7) +
  geom_point(alpha = 1, size = 5.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = .2, size = 2) +
  scale_color_manual(values = c(ci_lilac, ic_water)) +
  scale_fill_manual(values = c(ci_lilac, ic_water)) +
  labs(x = "Emotion", y = "Model Predictions of Choice") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(test_out_path, 'acc/emo_acc_main.png'), width = 7, height = 7)

# ...reward x age interaction
age_rew_interact <- ggpredict(test_acc_mod, terms = c("ExactAge [all]", "Reward"))
ggplot(age_rew_interact, aes(x, predicted)) +
  ylim(.3, .8) +
  geom_line(aes(color = group), size = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  scale_color_manual(values = c(gold, silver)) +
  scale_fill_manual(values = c(gold, silver)) +
  labs(x = "Age (Years)", y = "Model Predictions of Choice", color = "Reward", fill = "Reward") +
  emorl_theme
ggsave(paste0(test_out_path, 'acc/fig4a.png'), width = 7, height = 7)

##### RT -----

# Check for appropriate model distribution
hist(avg_test_rt_long$value) # Gaussian distribution seems to be the correct choice

# Run and evaluate test RT model
test_rt_mod <- lmer(value ~ Reward * Emotion * ExactAge + (1 | StudyID), data = avg_test_rt_long)
set.seed(123)
check_predictions(test_rt_mod) # Fine
qqnorm(residuals(test_rt_mod), pch = 1, frame = FALSE)
qqline(residuals(test_rt_mod), col = "red", lwd = 2) # Decent, but see if Gamma distribution better fulfills model assumptions
avg_test_rt_long_cleanrt <- avg_test_rt_long[!is.na(avg_test_rt_long$value), ] # Gamma distribution cannot take non-positive values
avg_test_rt_long_cleanrt <- avg_test_rt_long_cleanrt[avg_test_rt_long_cleanrt$value != 0, ] # Gamma distribution cannot take non-positive values
test_rt_mod_gamma <- glmer(value ~ Reward * Emotion * ExactAge + (1 | StudyID), family = Gamma(link = log), data = avg_test_rt_long_cleanrt)
set.seed(123)
check_predictions(test_rt_mod_gamma) # Fine
qqnorm(residuals(test_rt_mod_gamma), pch = 1, frame = FALSE)
qqline(residuals(test_rt_mod_gamma), col = "red", lwd = 2) # Convergence issues + not dramatically better --> proceed with Gaussian model

# Inspect test_rt_mod's significant main effects and interactions
Anova(test_rt_mod, type = "II")
plot_model(test_rt_mod, type = "pred", terms = "Reward")
plot_model(test_rt_mod, type = "pred", terms = "Emotion")
plot_model(test_rt_mod, type = "pred", terms = c("ExactAge [all]", "Reward"))

# Save test_rt_mod's outputs
test_rt_mod_filename <- paste0(test_out_path, 'rt/test_rt_mod.txt')
sink(test_rt_mod_filename)
print(Anova(test_rt_mod, type = "II"))
sink()

# Make plots for test_rt_mod's significant main effects and interactions...

# ...reward main effect
rew_rt_main <- ggpredict(test_rt_mod, terms = "Reward")
ggplot(data = rew_rt_main, aes(x = x, y = predicted, color = x, fill = x)) +
  ylim(1100, 1700) +
  geom_point(alpha = 1, size = 5.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = .2, size = 2) +
  scale_color_manual(values = c(gold, silver)) +
  scale_fill_manual(values = c(gold, silver)) +
  labs(x = "Reward", y = "Model Predictions of RT (ms)") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(test_out_path, 'rt/rew_rt_main.png'), width = 7, height = 7)

# ...emotion main effect
emo_rt_main <- ggpredict(test_rt_mod, terms = "Emotion")
ggplot(data = emo_rt_main, aes(x = x, y = predicted, color = x, fill = x)) +
  ylim(1100, 1700) +
  geom_point(alpha = 1, size = 5.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = .2, size = 2) +
  scale_color_manual(values = c(ci_lilac, ic_water)) +
  scale_fill_manual(values = c(ci_lilac, ic_water)) +
  labs(x = "Emotion", y = "Model Predictions of RT (ms)") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(test_out_path, 'rt/emo_rt_main.png'), width = 7, height = 7)

# ...reward x age interaction
age_rew_rt_interact <- ggpredict(test_rt_mod, terms = c("ExactAge [all]", "Reward"))
ggplot(age_rew_rt_interact, aes(x, predicted)) +
  geom_line(aes(color = group), size = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  scale_color_manual(values = c(gold, silver)) +
  scale_fill_manual(values = c(gold, silver)) +
  labs(x = "Age (Years)", y = "Model Predictions of RT (ms)", color = "Reward", fill = "Reward") +
  emorl_theme
ggsave(paste0(test_out_path, 'rt/fig4b.png'), width = 7, height = 7)

##### Similarity (learning pairs) -----

# Make descriptive plot of from_learn choices
ggplot(data = overall_from_learn_rt, aes(x = pair_type_simp, y = ind_avg_rt, color = pair_type_simp, fill = pair_type_simp)) +
  geom_hline(yintercept = proportion_rt_wide_format_y_axis, colour = 'grey90') +
  scale_y_continuous(breaks = proportion_rt_wide_format_y_axis) +
  geom_violin(alpha = 0.25, size = 1) +
  geom_point(position = position_jitter(width = .1), alpha = .25, size = 2.75) +
  geom_boxplot(width = 0.2, alpha = 0, size = 1) +
  scale_color_manual(values = c(cc_tomato, ci_lilac, ic_water, ii_sage)) +
  scale_fill_manual(values = c(cc_tomato, ci_lilac, ic_water, ii_sage)) +
  labs(x = "Condition from Learning Task", y = "Mean RT (ms)") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(test_out_path, 'delib/figS8.png'), width = 9, height = 7)

# Check for appropriate model distribution
hist(from_learn$rt) # Gaussian distribution should be interrogated

# Run and evaluate test deliberation model for learning pairs across age
test_delib_from_learn_mod <- lmer(rt ~ pair_type_simp * ExactAge + (1 | StudyID), data = from_learn)
set.seed(123)
check_predictions(test_delib_from_learn_mod) # Not ideal
qqnorm(residuals(test_delib_from_learn_mod), pch = 1, frame = FALSE)
qqline(residuals(test_delib_from_learn_mod), col = "red", lwd = 2) # Not ideal, see if Gamma distribution better fulfills model assumptions
from_learn_cleanrt <- from_learn[!is.na(from_learn$rt), ] # Gamma distribution cannot take non-positive values
from_learn_cleanrt <- from_learn_cleanrt[from_learn_cleanrt$rt != 0, ] # Gamma distribution cannot take non-positive values
test_delib_from_learn_mod_gamma <- glmer(rt ~ pair_type_simp * ExactAge + (1 | StudyID), family = Gamma(link = log), data = from_learn_cleanrt)
set.seed(123)
check_predictions(test_delib_from_learn_mod_gamma) # Not ideal
qqnorm(residuals(test_delib_from_learn_mod_gamma), pch = 1, frame = FALSE)
qqline(residuals(test_delib_from_learn_mod_gamma), col = "red", lwd = 2) # Convergence issues + not dramatically better --> proceed with Gaussian model

# Inspect test_delib_from_learn_mod's significant main effects and interactions
Anova(test_delib_from_learn_mod, type = "II")

# Save test_delib_from_learn_mod's outputs
test_delib_from_learn_mod_filename <- paste0(test_out_path, 'delib/test_delib_from_learn_mod.txt')
sink(test_delib_from_learn_mod_filename)
print(Anova(test_delib_from_learn_mod, type = "II"))
sink()

##### Similarity (novel pairs) -----

# Check for appropriate model distribution
hist(all_else$rt) # Gaussian distribution should be interrogated

# Run and evaluate test deliberation model for all other pair types across age
test_delib_all_else_mod <- lmer(rt ~ pair_class * ExactAge + (1 | StudyID), data = all_else)
set.seed(123)
check_predictions(test_delib_all_else_mod) # Not ideal
qqnorm(residuals(test_delib_all_else_mod), pch = 1, frame = FALSE)
qqline(residuals(test_delib_all_else_mod), col = "red", lwd = 2) # Not ideal, see if Gamma distribution better fulfills model assumptions
all_else_cleanrt <- all_else[!is.na(all_else$rt), ] # Gamma distribution cannot take non-positive values
all_else_cleanrt <- all_else_cleanrt[all_else_cleanrt$rt != 0, ] # Gamma distribution cannot take non-positive values
test_delib_all_else_mod_gamma <- glmer(rt ~ pair_class * ExactAge + (1 | StudyID), family = Gamma(link = log), data = all_else_cleanrt)
set.seed(123)
check_predictions(test_delib_all_else_mod_gamma) # Not ideal
qqnorm(residuals(test_delib_all_else_mod_gamma), pch = 1, frame = FALSE)
qqline(residuals(test_delib_all_else_mod_gamma), col = "red", lwd = 2) # Convergence issues + not dramatically better --> proceed with Gaussian model

# Inspect test_delib_all_else_mod's significant main effects and interactions
Anova(test_delib_all_else_mod, type = "II")
plot_model(test_delib_all_else_mod, type = "pred", terms = "pair_class")
plot_model(test_delib_all_else_mod, type = "pred", terms = c("ExactAge", "pair_class"))

# Conduct post-hoc tests of ExactAge x pair_class interaction
age_pts <- c(10, 15.5, 21) # Age points of interest
em_age_pair_rt <- emmeans(test_delib_all_else_mod, specs = c("pair_class"), by = c("ExactAge"), at = list(ExactAge = age_pts))
em_age_pair_rt # Marginal means for each point of interest                       
comp_em_age_pair_rt <- contrast(em_age_pair_rt, 
                                method = "trt.vs.ctrl", 
                                ref = "same_rew_same_emo",
                                adjust = "holm")
comp_em_age_pair_rt 
# After correcting for multiple comparisons:
# --> Full Match is significantly slower than Full Mismatch and Reward Match at age = 10
# --> Full Match slowness is driven by emotions at age = 10
# --> Full Match is significantly slower than all other match types at age = 15.5
# --> Full Match slowness is driven by an emotion-reward combination at age = 15.5
# --> Full Match is significantly slower than Full Mismatch and Emotion Match at age = 21
# --> Full Match slowness is driven by rewards at age = 21

# Save test_delib_all_else_mod's outputs
test_delib_all_else_mod_filename <- paste0(test_out_path, 'delib/test_delib_all_else_mod.txt')
sink(test_delib_all_else_mod_filename)
print(Anova(test_delib_all_else_mod, type = "II"))
cat("\nEMMEANS FOR LMER RT MODEL\n\n")
cat("Age points:\n")
print(age_pts)
cat("\nHolm-Corrected P-Values\n\n")
comp_em_age_pair_rt
sink()

# Make plots for test_delib_all_else_mod's significant main effects and interactions...

# ...pair class main effect
pair_class_main <- ggpredict(test_delib_all_else_mod, terms = "pair_class")
pair_class_main$x <- as.character(pair_class_main$x)
pair_class_main$x[pair_class_main$x == "diff_rew_diff_emo"] <- "Full Mismatch"
pair_class_main$x[pair_class_main$x == "diff_rew_same_emo"] <- "Emotion Match"
pair_class_main$x[pair_class_main$x == "same_rew_diff_emo"] <- "Reward Match"
pair_class_main$x[pair_class_main$x == "same_rew_same_emo"] <- "Full Match"
pair_class_main$x <- factor(pair_class_main$x, levels = c("Full Match", "Reward Match", "Emotion Match", "Full Mismatch"))
ggplot(data = pair_class_main, aes(x = x, y = predicted, color = x, fill = x)) +
  geom_point(alpha = 1, size = 5.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = .2, size = 2) +
  scale_color_manual(values = c(cc_tomato, gold, silver, ii_sage)) +
  scale_fill_manual(values = c(cc_tomato, gold, silver, ii_sage)) +
  labs(x = "Degree of Option Match in Test Phase", y = "Model Predictions of RT (ms)") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(test_out_path, 'delib/pair_class_main.png'), width = 9, height = 7)

# ...pair class x age interaction
age_pair_class_interact <- ggpredict(test_delib_all_else_mod, terms = c("ExactAge [all]", "pair_class"))
age_pair_class_interact$group <- as.character(age_pair_class_interact$group)
age_pair_class_interact$group[age_pair_class_interact$group == "diff_rew_diff_emo"] <- "Full Mismatch"
age_pair_class_interact$group[age_pair_class_interact$group == "diff_rew_same_emo"] <- "Emotion Match"
age_pair_class_interact$group[age_pair_class_interact$group == "same_rew_diff_emo"] <- "Reward Match"
age_pair_class_interact$group[age_pair_class_interact$group == "same_rew_same_emo"] <- "Full Match"
age_pair_class_interact$group <- factor(age_pair_class_interact$group, levels = c("Full Match", "Reward Match", "Emotion Match", "Full Mismatch"))
ggplot(age_pair_class_interact, aes(x, predicted)) +
  geom_line(aes(color = group), size = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  scale_color_manual(values = c(cc_tomato, gold, silver, ii_sage)) +
  scale_fill_manual(values = c(cc_tomato, gold, silver, ii_sage)) +
  labs(x = "Age (Years)", y = "Model Predictions of RT (ms)") +
  emorl_theme + theme(legend.title = element_blank()) +
  geom_vline(xintercept = age_pts, linetype = 'dashed', size = 2)
ggsave(paste0(test_out_path, 'delib/fig4c.png'), width = 9, height = 7)

### Scale analyses -----

##### Subjective ratings -----

# How does reward and emotion information persist to bias subjective emotion valence ratings across age?

# Check for appropriate model distribution
hist(subj_emo_val$rating_diff) # Gaussian distribution is the correct choice

# Run and evaluate subjective emotion valence change model
subj_emo_val_mod <- lmer(rating_diff ~ reward * emotion * ExactAge + (1 | StudyID), data = subj_emo_val)
set.seed(123)
check_predictions(subj_emo_val_mod) # Poor
qqnorm(residuals(subj_emo_val_mod), pch = 1, frame = FALSE)
qqline(residuals(subj_emo_val_mod), col = "red", lwd = 2) # Okay
unique(subj_emo_val$rating_diff) # 11 values --> continue with Gaussian model

# Inspect subj_emo_val_mod's significant main effects and interactions
Anova(subj_emo_val_mod, type = "II")
plot_model(subj_emo_val_mod, type = "pred", terms = c("reward"))
plot_model(subj_emo_val_mod, type = "pred", terms = c("emotion"))
plot_model(subj_emo_val_mod, type = "pred", terms = c("ExactAge [all]")) # Marginal
plot_model(subj_emo_val_mod, type = "pred", terms = c("ExactAge", "reward"))
plot_model(subj_emo_val_mod, type = "pred", terms = c("ExactAge", "emotion"))

# Save subj_emo_val_mod's outputs
subj_emo_val_mod_filename <- paste0(scales_out_path, 'subj_emo_val/subj_emo_val_mod.txt')
sink(subj_emo_val_mod_filename)
print(Anova(subj_emo_val_mod, type = "II"))
sink()

# Make plots for subj_emo_val_mod's significant main effects and interactions...

# ...reward main effect
reward_main <- ggpredict(subj_emo_val_mod, terms = "reward")
ggplot(data = reward_main, aes(x = x, y = predicted, color = x, fill = x)) +
  ylim(-.5, .5) +
  geom_hline(yintercept = c(0), colour = 'black', linetype = 'dashed') +
  geom_point(alpha = 1, size = 5.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = .2, size = 2) +
  scale_color_manual(values = c(gold, silver)) +
  scale_fill_manual(values = c(gold, silver)) +
  labs(x = "Reward", y = "Model Predictions of\nValence Change") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(scales_out_path, 'subj_emo_val/rew_main.png'), width = 7, height = 7)

# ...emotion main effect
emotion_main <- ggpredict(subj_emo_val_mod, terms = "emotion")
ggplot(data = emotion_main, aes(x = x, y = predicted, color = x, fill = x)) +
  ylim(0, 1) +
  geom_hline(yintercept = c(0), colour = 'black', linetype = 'dashed') +
  geom_point(alpha = 1, size = 5.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = .2, size = 2) +
  scale_color_manual(values = c(ci_lilac, ic_water)) +
  scale_fill_manual(values = c(ci_lilac, ic_water)) +
  labs(x = "Emotion", y = "Model Predictions of\nValence Change") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(scales_out_path, 'subj_emo_val/emo_main.png'), width = 7, height = 7)

# ...age main effect
age_main <- ggpredict(subj_emo_val_mod, terms = c("ExactAge [all]"))
ggplot(data = age_main, aes(x = x, y = predicted)) +
  ylim(-.5, 1) +
  geom_line(color = cc_tomato, size = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.25, fill = cc_tomato) +
  labs(x = "Age (Years)", y = "Model Predictions of\nValence Change") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(scales_out_path, 'subj_emo_val/age_main.png'), width = 4.5, height = 7)

# ...age x reward interaction
age_reward_interact <- ggpredict(subj_emo_val_mod, terms = c("ExactAge [all]", "reward"))
ggplot(age_reward_interact, aes(x, predicted)) +
  ylim(-1, 1) +
  geom_hline(yintercept = c(0), colour = 'black', linetype = 'dashed') +
  geom_line(aes(color = group), size = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  scale_color_manual(values = c(gold, silver)) +
  scale_fill_manual(values = c(gold, silver)) +
  labs(x = "Age (Years)", y = "Model Predictions of\nValence Change", color = "Reward", fill = "Reward") +
  emorl_theme
ggsave(paste0(scales_out_path, 'subj_emo_val/figS6.png'), width = 7, height = 7)

# ...age x emotion interaction
age_emotion_interact <- ggpredict(subj_emo_val_mod, terms = c("ExactAge [all]", "emotion"))
ggplot(age_emotion_interact, aes(x, predicted)) +
  ylim(-.5, 1.5) +
  geom_hline(yintercept = c(0), colour = 'black', linetype = 'dashed') +
  geom_line(aes(color = group), size = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  scale_color_manual(values = c(ci_lilac, ic_water)) +
  scale_fill_manual(values = c(ci_lilac, ic_water)) +
  labs(x = "Age (Years)", y = "Model Predictions of\nValence Change", color = "Emotion", fill = "Emotion") +
  emorl_theme
ggsave(paste0(scales_out_path, 'subj_emo_val/figS7.png'), width = 7, height = 7)

##### Reinforcement reports -----

# How are reward outcomes tracked across age?

# Check for appropriate model distribution
hist(rein_rep$ReinReport) # Gaussian distribution should be interrogated

# Run and evaluate reinforcement report model
rein_rep_mod <- lmer(ReinReport ~ reward * emotion * ExactAge + (1 | StudyID), data = rein_rep)
check_predictions(rein_rep_mod) # Okay
qqnorm(residuals(rein_rep_mod), pch = 1, frame = FALSE)
qqline(residuals(rein_rep_mod), col = "red", lwd = 2) # Decent, but see if Beta distribution better fulfills model assumptions
rein_rep$report_fract <- rein_rep$ReinReport / 45 # Normalize (values need to fall between 0 and 1)
table(rein_rep$report_fract) # Lots of values at 0 bound and some values at 1 bound
rein_rep$report_squeeze <- betaSqueeze(rein_rep$report_fract)
table(rein_rep$report_squeeze) # No values at 0 or 1 bounds
rein_rep_mod_beta <- glmmTMB(report_squeeze ~ reward * emotion * ExactAge + (1 | StudyID), data = rein_rep, family = beta_family)
set.seed(123)
check_predictions(rein_rep_mod_beta) # Decent
set.seed(123)
simres <- simulateResiduals(fittedModel = rein_rep_mod_beta, n = 250)
hist(simres$scaledResiduals, xlab = "scaled residuals", main = "Histogram Residuals") # Somewhat uniform distribution
plotQQunif(simres, testUniformity = FALSE, testOutliers = FALSE) # Not dramatically better --> proceed with Gaussian model

# Inspect rein_rep_mod's significant main effects and interactions
Anova(rein_rep_mod, type = "II")
plot_model(rein_rep_mod, type = "pred", terms = c("reward"))
plot_model(rein_rep_mod, type = "pred", terms = c("emotion"))
plot_model(rein_rep_mod, type = "pred", terms = c("ExactAge [all]", "reward"))

# Save rein_rep_mod's outputs
rein_rep_filename <- paste0(scales_out_path, 'rein_rep/rein_rep_mod.txt')
sink(rein_rep_filename)
print(Anova(rein_rep_mod, type = "II"))
sink()

# Make plots for rein_rep_mod's significant main effects and interactions...

# ...reward main effect
rew_main <- ggpredict(rein_rep_mod, terms = "reward")
ggplot(data = rew_main, aes(x = x, y = predicted, color = x, fill = x)) +
  ylim(4, 20) +
  geom_point(alpha = 1, size = 5.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = .2, size = 2) +
  scale_color_manual(values = c(gold, silver)) +
  scale_fill_manual(values = c(gold, silver)) +
  labs(x = "Reward", y = "Model Predictions of\nReinforcement Reported") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(scales_out_path, 'rein_rep/rew_main.png'), width = 7, height = 7)

# ...emotion main effect
emo_main <- ggpredict(rein_rep_mod, terms = "emotion")
ggplot(data = emo_main, aes(x = x, y = predicted, color = x, fill = x)) +
  ylim(4, 20) +
  geom_point(alpha = 1, size = 5.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = .2, size = 2) +
  scale_color_manual(values = c(ci_lilac, ic_water)) +
  scale_fill_manual(values = c(ci_lilac, ic_water)) +
  labs(x = "Emotion", y = "Model Predictions of\nReinforcement Reported") +
  emorl_theme + theme(legend.position = "none")
ggsave(paste0(scales_out_path, 'rein_rep/emo_main.png'), width = 7, height = 7)

# ...age x reward interaction
age_rew_interact <- ggpredict(rein_rep_mod, terms = c("ExactAge [all]", "reward"))
ggplot(age_rew_interact, aes(x, predicted)) +
  ylim(4, 20) +
  geom_line(aes(color = group), size = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25) +
  scale_color_manual(values = c(gold, silver)) +
  scale_fill_manual(values = c(gold, silver)) +
  labs(x = "Age (Years)", y = "Model Predictions of\nReinforcement Reported", color = "Reward", fill = "Reward") +
  emorl_theme
ggsave(paste0(scales_out_path, 'rein_rep/figS9.png'), width = 7, height = 7)

##### Subjective value of $ -----

# How are small amounts of money differentiated across age?

# Make descriptive plot of amounts
ggplot(data = subj_val_money_long, aes(x = ExactAge, y = Report, fill = Amount, color = Amount)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = c(0, 100), colour = 'black', linetype = "dashed") +
  geom_hline(yintercept = seq(20, 80, 20), colour = 'grey80') +
  scale_y_continuous(breaks = seq(0, 100, 20)) +  
  geom_point(alpha = 1, size = 1.5) +
  geom_smooth(method = "lm") + 
  scale_color_manual(values = c(silver, not_capt_tomato, masc_tomato, fem_tomato, gold)) +
  scale_fill_manual(values = c(silver, not_capt_tomato, masc_tomato, fem_tomato, gold)) +
  labs(x = "Age (Years)", y = "How much is this amount?\n(Not much money to A lot of money)") +
  emorl_theme
ggsave(paste0(scales_out_path, 'subj_val_money/amts_age.png'), width = 9, height = 7.075)

# Check for appropriate model distribution
hist(subj_val_money_long$Report) # Gaussian distribution is probably incorrect, but try it anyways

# Run and evaluate subjective value of $ model
subj_val_age <- lmer(Report ~ ExactAge * Amount + (1 | StudyID), data = subj_val_money_long)
set.seed(123)
check_predictions(subj_val_age) # Poor
qqnorm(residuals(subj_val_age), pch = 1, frame = FALSE)
qqline(residuals(subj_val_age), col = "red", lwd = 2) # Poor; test Beta distribution to see if model assumptions can be better met
subj_val_money_long$report_fract <- subj_val_money_long$Report / 100 # Normalize (values need to fall between 0 and 1)
table(subj_val_money_long$report_fract) # Lots of values at 0 bound and some values at 1 bound
subj_val_money_long$report_squeeze <- betaSqueeze(subj_val_money_long$report_fract)
table(subj_val_money_long$report_squeeze) # No values at 0 or 1 bounds
subj_val_age_beta <- glmmTMB(report_squeeze ~ ExactAge * Amount + (1 | StudyID), data = subj_val_money_long, family = beta_family)
set.seed(123)
check_predictions(subj_val_age_beta) # Decent
set.seed(123)
simres <- simulateResiduals(fittedModel = subj_val_age_beta, n = 250)
hist(simres$scaledResiduals, xlab = "scaled residuals", main = "Histogram Residuals") # Somewhat uniform distribution
plotQQunif(simres, testUniformity = FALSE, testOutliers = FALSE) # Noticeably better --> proceed with Beta distribution

# Inspect subj_val_age_beta's significant main effects and interactions
Anova(subj_val_age_beta, type = "II")
plot_model(subj_val_age_beta, type = "pred", terms = c("ExactAge")) # Subjective value of $ scores increase with age
plot_model(subj_val_age_beta, type = "pred", terms = c("Amount")) # Subjective value of $ scores increase with value amount
plot_model(subj_val_age_beta, type = "pred", terms = c("ExactAge", "Amount")) # Subjective value of $ scores become less differentiated with age

# Save subjective value of $ model output 
subj_val_model <- paste0(scales_out_path, 'subj_val_money/subj_val_age_lm.txt')
sink(subj_val_model)
Anova(subj_val_age_beta, type = "II")
sink()
