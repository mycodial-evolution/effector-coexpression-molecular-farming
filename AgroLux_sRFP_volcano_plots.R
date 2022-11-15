### Volcano plot of AgroLux and sRFP reporters for all RXLR effectors
### (Figure 4)
### By Freddie King
### 22/01/23

### Load packages ##############################################################
library(tidyverse)
library(PMCMRplus)
library(MASS)
library(ggrepel)
library(ggpubr)

### Loading and processing data ################################################
# Uses microplate data in a .csv table with 5 columns: 'treatment', 'rfp', 'lux', 
# 'gfp' and 'experiment' (see effector_screen_data.csv)

# Choose formatted microplate reader data in .csv file
raw <- read_csv(file.choose())

# Select the required columns from the microplate dataframe
raw <- raw %>%
  dplyr::select(treatment:experiment)

# Extract default row order
row_order <- unique(raw$treatment)

# Set the factors and their levels
raw <- raw %>% mutate(treatment = factor(treatment, levels=row_order),
                    experiment = factor(experiment))

# Assign an ID value to specify linear model per experimental groups
raw <- raw %>%
  mutate(ID = case_when(
    experiment == "FK_009" | experiment == "FK_010" | experiment == "FK_011" ~ "1",
    experiment == "FK_014" | experiment == "FK_015" ~ "2",
    experiment == "FK_016" | experiment == "FK_017" ~ "3",
    experiment == "FK_019" | experiment == "FK_020" ~ "4",
    experiment == "FK_023" | experiment == "FK_024" ~ "5"))

### Calculating reporter fold-change values ####################################

# Set base values within each ID group of Ev RFP and Lux values
df <- raw %>% group_by(ID) %>%
 mutate(rfp_base = mean(subset(rfp, treatment=="Ev"))) %>%
 mutate(lux_base = mean(subset(lux, treatment=="Ev")))

# Calculate fold change by mean over base mean
df <- df %>% group_by(ID, treatment) %>%
  mutate(rfp_fc = log2(mean(rfp) / rfp_base)) %>%
  mutate(lux_fc = log2(mean(lux) / lux_base))

# Remove Ev and AvrBlb2 values and duplicates
df <- df %>%
  filter(treatment!="Ev" & treatment!="AvrBlb2") %>%
  distinct(treatment, .keep_all = TRUE)
  
### Calculate p-values #########################################################

# Define model functions that fit two-way ANOVAs according to ID value
rfp_model <- function(df, id) {
  result <- as.data.frame(dunnettTest(
    aov(rfp ~ treatment*experiment, data=subset(df, ID==id)))$p.value)
  return(result)
}

lux_model <- function(df, id) {
  result <- as.data.frame(dunnettTest(
    aov(lux ~ treatment*experiment, data=subset(df, ID==id)))$p.value)
  return(result)
}

# Calculate p-values using RFP models
rfp_pvalues <- rbind((rfp_model(raw, 1)), (rfp_model(raw, 2)), (rfp_model(raw, 3)),
                     (rfp_model(raw, 4)), (rfp_model(raw, 5)))

# Extract effector treatment names under treatment column
rfp_pvalues <- rfp_pvalues %>% rownames_to_column("treatment")

# Calculate p-values using AgroLux models
lux_pvalues <- rbind((lux_model(raw, 1)), (lux_model(raw, 2)), (lux_model(raw, 3)),
                     (lux_model(raw, 4)), (lux_model(raw, 5)))

# Extract effector treatment names under treatment column
lux_pvalues <- lux_pvalues %>% rownames_to_column("treatment")

# Combine p-values from RFP and AgroLux models 
pvalues <- left_join(rfp_pvalues, lux_pvalues, by = "treatment")

# Convert rownames to column and remove control p-values
# In this case, controls are Ev and AvrBlb2 treatments
pvalues <- pvalues %>%
  rename_at("Ev.x", ~"rfp_pvalues") %>%
  rename_at("Ev.y", ~"lux_pvalues") %>%
  filter(treatment != "Ev") %>%
  filter(!str_detect(treatment, "AvrBlb"))

### Combine p-values and fold-change values ####################################
results <- left_join(df, pvalues)

### Plotting volcano plots #####################################################

# Creating variables to indicate significant changes
# Fold-change and p-value threshold for 'significant' changes can be set
# by changing the values in the inequalities below
results <- results %>%
  mutate(rfp_change = case_when(
    rfp_fc >= 0.3 & rfp_pvalues <= 0.05 ~ "INCREASE",
    rfp_fc <= -0.3 & rfp_pvalues <= 0.05 ~ "DECREASE",
    TRUE ~ "NO"
    )
  ) %>%
  mutate(lux_change = case_when(
    lux_fc >= 0.5 & lux_pvalues <= 0.05 ~ "INCREASE",
    lux_fc <= -0.5 & lux_pvalues <= 0.05 ~ "DECREASE",
    TRUE ~ "NO"
    )
  )

# Creating variables to label significant treatments
results <- results %>%
  mutate(rfp_labels = case_when(
    rfp_change != "NO" ~ treatment,
    TRUE ~ ""
    )
  ) %>%
  mutate(lux_labels = case_when(
    lux_change != "NO" ~ treatment,
    TRUE ~ ""
    )
  )

# Plotting the RFP volcano plot
# Plots fold-changes on the x-axis and -log10(p-values) on the y-axis
# Points are coloured according to significance label and label gives
# treatment name
rfp_vol <- ggplot(results, aes(x=rfp_fc, y=-log10(rfp_pvalues),
  col=rfp_change, label=rfp_labels)) +
  geom_point(size=2) + # Set data point sizes
  # Plots 4 segments as lines indicating thresholds for 'significant' changes
  geom_segment(x=-2, xend=-0.3, y=1.3, yend=1.3, color="grey80", 
               alpha=0.3, size=0.75) +
  geom_segment(x=-0.3, xend=-0.3, y=1.3, yend=4.2, color="grey80", 
               alpha=0.3, size=0.75) +
  geom_segment(x=0.3, xend=1.2, y=1.3, yend=1.3, color="grey80", 
               alpha=0.3, size=0.75) +
  geom_segment(x=0.3, xend=0.3, y=1.3, yend=4.2, color="grey80", 
               alpha=0.3, size=0.75) +
  # Sets black and white theme
  theme_bw() +
  # Sets data point labels to repel each other
  geom_text_repel(na.rm = TRUE) +
  # Removes axes labels
  labs(x=" ",
       y=" ",
       title=" ") + 
  # Sets black and white theme then sets font size to 12, hairlines to 0.5 pt
  theme(text = element_text(size = 12, color="black"), 
        panel.grid.minor = element_line(size = 0.5),
        panel.grid.major = element_line(size = 0.5), legend.position="none")  
rfp_vol

# Plotting the AgroLux volcano plot
lux_vol <- ggplot(results, aes(x=lux_fc, y=-log10(lux_pvalues), 
                                col=lux_change, label=lux_labels)) +
  geom_point(size=2) +
  geom_segment(x=-2, xend=-0.5, y=1.3, yend=1.3, color="grey80", 
               alpha=0.3, size=0.75) +
  geom_segment(x=-0.5, xend=-0.5, y=1.3, yend=15, color="grey80", 
               alpha=0.3, size=0.75) +
  geom_segment(x=0.5, xend=2.1, y=1.3, yend=1.3, color="grey80", 
               alpha=0.3, size=0.75) +
  geom_segment(x=0.5, xend=0.5, y=1.3, yend=15, color="grey80", 
               alpha=0.3, size=0.75) +
  theme_bw() +
  geom_text_repel(na.rm = TRUE) +
  labs(x="",
       y="",
       colour="",
       title="") + 
  theme(text = element_text(size = 12, color="black"), panel.grid.minor = element_line(size = 0.5),
        panel.grid.major = element_line(size = 0.5), legend.position="none")
lux_vol

# Combine the RFP and AgroLux volcano plots in two columns on one plot
vols <- ggarrange(rfp_vol, lux_vol,
          ncol=2, nrow=1)
vols

# Exports the plot as a .png
png(file="filename.png", units="in", width=10, height=8, res=400)
vols
dev.off()