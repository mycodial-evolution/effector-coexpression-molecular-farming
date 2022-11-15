### Bar charts of relative band density from Western blots
### using the paired t-test for comparison (Figures 8B, S4B, D and F)
### By Freddie King
### 2022/03/31

### Load packages ##############################################################
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(PMCMRplus)
library(MASS)

### Loading and processing data ################################################
# Uses quantified band data from ImageJ in a .csv table with 5 columns:
# 'peak', 'area', 'percentage', 'treatment' and 'bio_rep'
# (see .csv)

# Choose formatted microplate reader data in .csv file
raw <- read.csv(file.choose())
raw$treatment <- as.factor(raw$treatment) # Converts treatments to factors
raw$bio_rep <- as.factor(raw$bio_rep) # Converts biological replicates to factors

# Define the standard error function
se <- function(x){sd(x)/sqrt(length(x))}

# Summarises the individual data points by calculating mean percentage of the 
# control (Ev) within each biological replicate, 
# relative density compared to the control for each treatment and standard error 
# of the relative density values
bands <- raw %>% group_by(bio_rep) %>%
  mutate(standard = mean(subset(percentage, treatment=="Ev"))) %>%
  group_by(treatment, bio_rep) %>%
  mutate(relative_density = percentage / standard) %>%
  ungroup() %>%
  group_by(treatment) %>%
  mutate(percentage_se = se(relative_density))

# Summarises treatment data with percentage mean, relative density mean and 
# standard error of the percentage values
summ_band <- bands %>% group_by(treatment) %>%
  summarise(percentage_mean = mean(percentage),
         relative_density = mean(relative_density),
         percentage_se = percentage_se) %>%
  distinct()

### Statistical testing: paired t-test #########################################

# Select treatment (effector) for statistical comparison
test <- bands %>% subset(treatment == "Ev" | treatment == "E6 P1")

# Paired t-test for testing selected treatment 
model <- t.test(percentage ~ treatment, data = test, 
                alternative="greater", paired = TRUE)

# Converting model factor levels to factors
model$model$treatment <- as.factor(model$model$treatment)

# Extracting factor level order (factor order) from microplate data
model$factor_order <- unique(bands$treatment)

# Ordering factor levels according to row order
model$model$treatment <- factor(model$model$treatment, levels=model$factor_order)

# Extract the p-value
pvalues <- as.data.frame(model$p.value)
names(pvalues) <- "Ev"

# Converting numeric p-values to asterisks
pvalues$Ev <- symnum(pvalues$Ev, corr = FALSE, na = FALSE, 
                      cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                      symbols = c("***", "** ", "*  ", " "))

# Adding p-value asterisks for the selected treatment to the summary dataframe
summ_band$Ev <- c(NA, pvalues$Ev, NA)

### Plotting the relative band density bar chart ###############################
# The grey control bar for comparison (Ev) is set to 1 and relative density
# (calculated within each biological replicate) is plotted with error bars 
# representing the standard errors of the percentage values 

# Re-ordering treatment factor levels according to original data
summ_band$row_order <- unique(bands$treatment)
summ_band$treatment <- factor(summ_band$treatment, 
                              levels=summ_band$row_order)

# Plotting the bar chart
p <- ggplot() +
  # Plots bar chart of relative density for each treatment
  geom_bar(data=summ_band, aes(y=relative_density, x=treatment, fill=treatment),
           stat="identity", width=0.3, alpha=0.75, show.legend = FALSE) +
  # Manually define treatment colours with Ev control set to grey
  scale_fill_manual(values = c("grey50", "#00BFC4", "#00BFC4", "#00BFC4")) +
  # Annotates p-value asterisk from the t-test
  # May need to change position_nudge to get the asterisk in line with the bar
  geom_text(data=summ_band, aes(label=Ev, x=treatment, y = relative_density*1.5), 
            position = position_nudge(x = -0.09), size=7) +
  # Plots error bars representing the percentage standard errors
  geom_errorbar(data=summ_band, aes(y=relative_density, x=treatment, 
                ymin=relative_density-percentage_se, ymax=relative_density+percentage_se),
                width=0.1) +
  # Plots individual data points of relative density for each treatment
  geom_point(data=bands, aes(y=relative_density, x=treatment), colour="grey30", 
             size=1.8, alpha=0.6) +
  # Removes axes labels
  labs(x="", y="") +
  # Flips the chart to make it horizontal 
  coord_flip() +
  # Sets black and white theme then sets font size to 12, hairlines to 0.5 pt
  theme_bw() +
  theme(text = element_text(size = 8, color="black"))
p

# Exports the plot as a .png
png(file="filename.png", units="in", width=8, height=4, res=400)
p
dev.off()
