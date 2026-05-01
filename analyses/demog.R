###########################################################
### Calculate EmoRL Demographic and WASI Data Summaries ###
###########################################################

# Written by: Camille Phaneuf-Hadd (cphaneuf@g.harvard.edu)
# Last updated: 3/14/26

# Inputs: demographic and WASI data
# Computes: 
# - visual summary of sample
# - statistical summary of sample
# - visual summary of WASI scores
# - statistical summary of WASI scores
# Outputs: 
# - png plot of age-gender distribution into results/demog
# - txt file of demographic stats into results/demog
# - png plot of WASI scores across age into results/wasi
# - text file of WASI stats into results/wasi

### Set up Script -----

# Load needed libraries
require(pacman) # for p_load()
p_load(tidyverse, # for df manipulation
       dplyr, # for %>% and other operators
       performance, # for check_predictions()
       sjPlot) # for plot_model()

# Load shared EmoRL functions
source("utilities.R")

# Set paths
in_path <- '../data/'
demog_out_path <- '../results/demog/'
wasi_out_path <- '../results/wasi/'

# Read in demographic data
demog <- read.csv(paste0(in_path, "demog.csv"))

# Read in WASI data
wasi <- read.csv(paste0(in_path, "behav/wasi.csv"))

### Save Age-Gender Distribution -----

# With vertical lines at age group boundaries
demog$Gender <- factor(demog$Gender, levels = c("Women", "Men", "Not Captured by Options"))
ggplot(data = demog, aes(x = FlooredAge, fill = Gender, color = Gender)) +
  scale_fill_manual(values = c(fem_tomato, masc_tomato, not_capt_tomato)) + 
  scale_color_manual(values = c(fem_tomato, masc_tomato, not_capt_tomato)) + 
  scale_x_continuous(breaks = c(8:22)) +
  geom_hline(yintercept = seq(0, 10, by = 2), colour = 'grey90') +
  scale_y_continuous(breaks = seq(0, 10, by = 2)) + 
  geom_histogram(binwidth = 1, alpha = .75) +
  geom_vline(xintercept = c(12.5, 17.5), colour = 'black', linetype = 'dashed') +
  labs(x = "Age (Years)", y = "Number of Participants") +
  emorl_theme
ggsave(paste0(demog_out_path, "figS1.png"), plot = last_plot(), width = 9, height = 7)

### Save Demographic Summaries -----

N = length(demog$StudyID)

# Write to file instead of the terminal
sink(paste0(demog_out_path, 'demog.txt'))

cat("All Participants\n")
cat("Total N:", N, "\n")
cat("\nAGE\n")
cat("Total Children:", sum(demog$AgeGroup == "Children"), "\n")
cat("Total Adolescents:", sum(demog$AgeGroup == "Adolescents"), "\n")
cat("Total Adults:", sum(demog$AgeGroup == "Adults"), "\n")
cat("Mean Age:", mean(demog$ExactAge), "\n")
cat("SD Age:", sd(demog$ExactAge), "\n")
cat("\nGENDER\n")
cat("Total Women:", sum(demog$Gender == "Women"), "\n")
cat("Total Men:", sum(demog$Gender == "Men"), "\n")
cat("Total Not Captured by Options:", sum(demog$Gender == "Not Captured by Options"), "\n")
cat("\nRACE\n")
cat("Percentage American Indian or Alaska Native:", round((sum(demog$Race == "American Indian or Alaska Native", na.rm = TRUE) / N) * 100, 2), "\n")
cat("Percentage Asian:", round((sum(demog$Race == "Asian", na.rm = TRUE) / N) * 100, 2), "\n")
cat("Percentage Black or African American:", round((sum(demog$Race == "Black or African American", na.rm = TRUE) / N) * 100, 2), "\n")
cat("Percentage White:", round((sum(demog$Race == "White", na.rm = TRUE) / N) * 100, 2), "\n")
cat("Percentage Not Captured by Options:", round((sum(demog$Race == "Not Captured by Options") / N) * 100, 2), "\n")
cat("\nETHNICITY\n")
cat("Percentage Hispanic:", round((sum(demog$Ethnicity == "Hispanic or Latino") / N) * 100, 2), "\n")
cat("Percentage Not Hispanic:", round((sum(demog$Ethnicity == "Not Hispanic or Latino") / N) * 100, 2), "\n")

# Stop writing to file
sink()

### Save WASI Age Curve -----

# Merge demographic and WASI data
wasi <- merge(wasi, demog[, c("StudyID", "ExactAge")], by = "StudyID")

# Matrix reasoning t-scores across age
ggplot(data = wasi, aes(x = ExactAge, y = t_score)) +
  scale_x_continuous(breaks = seq(8, 22, 2)) +
  geom_hline(yintercept = seq(25, 85, by = 10), colour = 'grey90') +
  scale_y_continuous(breaks = seq(25, 85, by = 10)) +  
  geom_point(fill = not_capt_tomato, color = not_capt_tomato, alpha = 1, size = 1.5) +
  geom_smooth(fill = not_capt_tomato, color = not_capt_tomato, method = "lm") + 
  labs(x = "Age (Years)", y = "T-Score") +
  emorl_theme
ggsave(paste0(wasi_out_path, "wasi_t_age.png"), width = 7, height = 7)

### Save WASI Summary -----

# Check for appropriate model distribution
hist(wasi$t_score) # Gaussian distribution is the correct choice

# Run and evaluate matrix reasoning t-score model
wasi_t_age <- lm(t_score ~ ExactAge, data = wasi)
set.seed(123)
check_predictions(wasi_t_age) # Great!
qqnorm(residuals(wasi_t_age), pch = 1, frame = FALSE)
qqline(residuals(wasi_t_age), col = "red", lwd = 2) # Great!

# Inspect wasi_t_age's significant main effects and interactions
summary(wasi_t_age)
plot_model(wasi_t_age, type = "pred", terms = c("ExactAge"))
# --> Sanity check: matrix reasoning t-scores are stable with age

# Save wasi_t_age's outputs
wasi_model <- paste0(wasi_out_path, 'wasi_age_lm.txt')
sink(wasi_model)
print(summary(wasi_t_age))
sink()
