### Horizontal bar chart of AgroLux, sRFP and GFP values for selected RXLR 
### effectors (Figures 2, 3 and 5B-D)
### By Freddie King
### 2021/12/09

### Load packages ##############################################################
library(ggplot2)
library(dplyr)
library(ggpubr)
library(PMCMRplus)
library(MASS)

### Loading and processing data ################################################
# Uses microplate data in a table with 5 columns: 'treatment', 'rfp', 'lux', 
# 'gfp' and 'experiment' (see effector_screen_data.csv)

# Choose formatted microplate reader data in a .csv file
microplate_dat <- read.csv(file.choose()) # Load csv data of choice

# Subset relevant experiments to plot
microplate_dat <- subset(microplate_dat, 
                         microplate_dat$experiment == " " |
                         microplate_dat$experiment == " ")

# Assigning empty vector (Ev, control) mean values per experiment
microplate_dat <- microplate_dat %>% group_by(experiment) %>%
  mutate(rfp_base = mean(subset(rfp, ?..treatment=="Ev"))) %>%
  mutate(gfp_base = mean(subset(gfp, ?..treatment=="Ev"))) %>%
  mutate(lux_base = mean(subset(lux, ?..treatment=="Ev")))

# Adding standardized values relative to the controls
microplate_dat <- microplate_dat %>% group_by(experiment) %>%
  mutate(rfp2 = rfp - rfp_base) %>%
  mutate(gfp2 = gfp - gfp_base) %>%
  mutate(lux2 = lux - lux_base)

# Define the standard error function and add mean and SE values to the data
se <- function(x){sd(x)/sqrt(length(x))}
summ_dat <- summarise(group_by(microplate_dat, ?..treatment), 
                      rfp_mean = mean(rfp2),
                      rfp_se = se(rfp2),
                      lux_mean = mean(lux2),
                      lux_se = se(lux2),
                      gfp_mean = mean(gfp2),
                      gfp_se = se(gfp2))

# Extracts unique treatments as row order column
summ_dat$row_order <- unique(microplate_dat$?..treatment) 

# Reorders the factor levels according to row order column
summ_dat$?..treatment <- factor(summ_dat$?..treatment, 
                                levels=summ_dat$row_order)

### Statistical testing: blocked two-way ANOVA with post-hoc Dunnett's test ####

# Two-way ANOVAs with blocking variable, the experiment replicate
rfp_model <- aov(rfp ~ ?..treatment*experiment, data = microplate_dat)
lux_model <- aov(lux ~ ?..treatment*experiment, data = microplate_dat)

# Converting model factor levels to factors for both models
rfp_model$model$?..treatment <- as.factor(rfp_model$model$?..treatment)
rfp_model$model$experiment <- as.factor(rfp_model$model$experiment)

lux_model$model$?..treatment <- as.factor(lux_model$model$?..treatment)
lux_model$model$experiment <- as.factor(lux_model$model$experiment)

# Extracting factor level order (factor order) from the microplate data
rfp_model$factor_order <- unique(microplate_dat$?..treatment)
lux_model$factor_order <- unique(microplate_dat$?..treatment)

# Ordering factor levels according to row order
rfp_model$model$?..treatment <- factor(rfp_model$model$?..treatment, 
                                       levels=rfp_model$factor_order)
lux_model$model$?..treatment <- factor(lux_model$model$?..treatment, 
                                       levels=lux_model$factor_order)

# Dunnett's test on the blocked two-way ANOVAs
rfp_pvalues <- as.data.frame(dunnettTest(rfp_model)$p.value)
lux_pvalues <- as.data.frame(dunnettTest(lux_model)$p.value)

# Converting numeric p-values to asterisks
rfp_pvalues$Ev <- symnum(rfp_pvalues$Ev, corr = FALSE, na = FALSE, 
                      cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                      symbols = c(" *** ", " ** ", " *  ", " "))
lux_pvalues$Ev <- symnum(lux_pvalues$Ev, corr = FALSE, na = FALSE, 
                      cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                      symbols = c(" *** ", " ** ", " *  ", " "))

# Adding p-value asterisks to the microplate dataframe
rfp_pvalues <- rbind(c(" "), rfp_pvalues) # Extracts p-values with space
rfp_pvalues$?..treatment <- sort(summ_dat$?..treatment) # Sort by dataframe order

lux_pvalues <- rbind(c(" "), lux_pvalues)
lux_pvalues$?..treatment <- sort(summ_dat$?..treatment)

# Joins the p-values with the summarised microplate dataframe
summ_dat <- left_join(summ_dat, rfp_pvalues, by="?..treatment")
summ_dat <- left_join(summ_dat, lux_pvalues, by="?..treatment")

### Plotting the RFP, luminescence and GFP graph ###############################
# Horizontal bar charts with the control (Ev) set to 0
# Plots SE error bars and individual points in grey
# Bars are ordered in size order

# Subset data according to effector co-expression treatments of interest
EOI<- c("Ev", "AvrPto", " ") # Insert effector names into this list to display
summ_dat <- subset(summ_dat, ?..treatment %in% EOI)
microplate_dat <- subset(microplate_dat, ?..treatment %in% EOI)

### Plotting the sRFP graph ###
p_rfp <- ggplot() +
  # Bar chart displaying RFP values for selected effectors
  geom_bar(data=summ_dat, aes(y=rfp_mean, x=reorder(?..treatment, -rfp_mean), 
           fill=rfp_mean), stat="identity", 
           width=0.5, alpha=0.75, show.legend = FALSE) +
  # Adjusts y-axis values to prevent cut-off axis labels
  scale_y_continuous(expand=expansion(mult=c(0.1,0.1))) +
  # Fills the bars according to mean RFP values
  scale_fill_gradient(low="#ffcccc", high="#ff0000", space="Lab") +
  # Adds asterisk annotations according to p-values
  # Adjust the y-axis multiplier to move the asterisks along the y-axis
  geom_text(data=summ_dat, aes(label=Ev.x, x=reorder(?..treatment, -rfp_mean),
            y = rfp_mean*1.7), vjust = 0.75, 
            position = position_dodge(width = 0.5), size=6) +
  # Plots SE error bars onto the bar charts
  geom_errorbar(data=summ_dat, aes(y=rfp_mean, x=reorder(?..treatment, -rfp_mean), 
            ymin=rfp_mean-rfp_se, ymax=rfp_mean+rfp_se), width=0.25) +
  # Plots the data points with shape according to experiment
  geom_point(data=microplate_dat, aes(y=rfp2, x=?..treatment, shape=experiment), 
             colour="grey30", size=1, alpha=0.6) +
  # Removes axes labels
  labs(shape=" ", x=" ", y=" ") +
  # Sets black and white theme then sets font size to 8, hairlines to 0.5 pt
  theme_bw() +
  theme(text = element_text(size = 8, color="black"), 
        panel.grid.minor = element_line(size = 0.5),
        panel.grid.major = element_line(size = 0.5), legend.position="none") +
  # Flips the chart to make it horizontal 
  coord_flip()
p_rfp

### Plotting the AgroLux graph ###
p_lux <- ggplot() +
  geom_bar(data=summ_dat, aes(y=lux_mean, x=reorder(?..treatment, -lux_mean), 
           fill=lux_mean), stat="identity", 
           width=0.5, alpha=0.75, show.legend = FALSE) +
  scale_y_continuous(expand=expansion(mult=c(0.1,0.1))) +
  scale_fill_gradient(low="#ffe9cc", high="#ff9500", space="Lab") +
  geom_text(data=summ_dat, aes(label=Ev.y, x=reorder(?..treatment, -lux_mean),
            y = lux_mean+(lux_mean*1.4)), vjust = 0.75,
            position = position_dodge(width = 0.5), size=6) +
  geom_errorbar(data=summ_dat, aes(y=lux_mean, x=reorder(?..treatment, -lux_mean),
           ymin=lux_mean-lux_se, ymax=lux_mean+lux_se), width=0.25) +
  geom_point(data=microplate_dat, aes(y=lux2, x=?..treatment, shape=experiment),
             colour="grey30", size=1, alpha=0.6) +
  labs(shape=" ", x=" ", y=" ") +
  theme_bw() +
  theme(text = element_text(size = 8, color="black"), 
        panel.grid.minor = element_line(size = 0.5),
        panel.grid.major = element_line(size = 0.5), legend.position="none") +
  coord_flip()
p_lux

### Plotting the GFP graph ###
p_gfp <- ggplot() +
  geom_bar(data=summ_dat, aes(y=gfp_mean, x=reorder(?..treatment, -gfp_mean), 
           fill=gfp_mean), stat="identity", 
           width=0.5, alpha=0.75, show.legend=FALSE) +
  scale_y_continuous(expand=expansion(mult=c(0.1,0.1))) +
  scale_fill_gradient(low="#e0ffd6", high="#66ff33", space="Lab") +
  geom_errorbar(data=summ_dat, aes(y=gfp_mean, x=reorder(?..treatment, -gfp_mean),
           ymin=gfp_mean-gfp_se, ymax=gfp_mean+gfp_se), width=0.25) +
  geom_point(data=microplate_dat, aes(y=gfp2, x=?..treatment, shape=experiment),
             colour="grey30", size=1, alpha=0.6) +
  labs(shape=" ", x=" ", y=" ") +
  theme_bw() +
  theme(text = element_text(size = 8, color="black"), 
        panel.grid.minor = element_line(size = 0.5),
        panel.grid.major = element_line(size = 0.5), legend.position="none") +
  coord_flip()
p_gfp

# Arranges AgroLux, sRFP and GFP bar charts in a single column on one plot
rfp_lux_gfp_plots <- ggarrange(p_rfp, p_lux, p_gfp, ncol=1, nrow=3)

# Exports the plot as a .png
png(file="filename.png", units="in", width=6, height=12, res=400)
rfp_lux_gfp_plots
dev.off()
