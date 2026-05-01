#########################################
### Define Shared Functions for EmoRL ###
#########################################

# Written by: Camille Phaneuf-Hadd (cphaneuf@g.harvard.edu)
# Last updated: 4/13/26

#########################
### General Utilities ###
#########################

# Define betaSqueeze() function (credit: Patrick Mair)
betaSqueeze <- function(y) {
  n <- length(y)
  y2 <- (y * (n-1) + 0.5) / n
  return(y2)
}

# Set variable levels
mods <- c("Vanilla", "EmoQ", "EmoBias", "EmoHier", "EmoAvg", "Cong", "Cond", "CongBeta", "CondBeta")
best_mods <- c("Vanilla", "EmoBias", "Cond", "CondBeta")
grps <- c("Children", "Adolescents", "Adults")

##########################
### Plotting Utilities ###
##########################

# Set common plotting variables -- behavioral responses
proportion_acc_format_y_axis <- c(0, .25, .5, .75, 1)
proportion_rt_format_y_axis <- c(0, 600, 1200, 1800, 2400)
proportion_rt_wide_format_y_axis <- c(0, 1100, 2200, 3300, 4400)
proportion_rat_format_y_axis <- c(-1, -.5, 0, .5, 1)

# Set common plotting variables -- computational modeling
beta_axis_labels <- c(0, 6, 12, 18, 24, 30)
alpha_axis_labels <- c(0, .25, .5, .75, 1)
summ_aicbic_range <- c(180, 200, 220, 240, 260)

# Set plotting color scheme -- individual colors
gold <- "#BE8D00"
silver <- "#A6A6A6"
cc_tomato <- "#FE6847"
ci_lilac <- "#D698C1"
ic_water <- "#6494AA"
ii_sage <- "#83B692"
not_capt_tomato <- cc_tomato
masc_tomato <- "#B74F3A"
fem_tomato <- "#742F21"
vanilla <- "#E0D4C1"
emo_q <- ci_lilac
emo_bias <- "#8D6580"
emo_heir <- ic_water
emo_avg <- "#3B5765"
cong <- ii_sage
cond <- cc_tomato
cong_beta <- "#53745F"
cond_beta <- masc_tomato

# Set plotting color scheme -- common combos
four_color_scheme_cond <- c(cc_tomato, ci_lilac, ic_water, ii_sage)
four_color_scheme_mods <- c(vanilla, emo_bias, cond, cond_beta)
nine_color_scheme_mods <- c(vanilla, emo_q, emo_bias, emo_heir, emo_avg, cong, cond, cong_beta, cond_beta)

# Set plotting theme -- RLDM and publication version
emorl_theme <- theme(title = element_text(size = 24, face = "bold", family = "Avenir"),
                     plot.title = element_text(hjust = .5),
                     axis.title.x = element_text(size = 24, family = "Avenir"),
                     axis.title.y = element_text(size = 24, family = "Avenir"),
                     axis.text.x = element_text(size = 18, colour = "black", family = "Avenir"),
                     axis.text.y = element_text(size = 18, colour = "black", family = "Avenir"),
                     legend.text = element_text(size = 18, colour = "black", family = "Avenir"),
                     legend.position = "bottom",
                     legend.key = element_rect(fill = "transparent", color = NA),
                     strip.text.x = element_text(size = 18, colour = "black", family = "Avenir"),
                     strip.text.y = element_text(size = 18, colour = "black", family = "Avenir"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"))
