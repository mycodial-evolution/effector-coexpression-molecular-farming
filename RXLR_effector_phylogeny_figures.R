### Phylogeny figures of the screened RXLR effectors and characterised RXLR
### RXLR effectors (Figures 6, S3)
### By Freddie King
### 22/02/07

### Load packages ##############################################################
library(tidyverse)
library(ape)
library(ggplot2)
library(ggtree)
library(ggrepel)
library(treeio)
library(ggpubr)

### Loading and processing data ################################################
# Choose phylogeny with bootstraps in .tre file
tree <- read.tree(file.choose())

# Choose .csv file for converting gene names to treatment names
conversion <- read.csv(file.choose())

# Renames columns and selects relevant columns 
conversion <- conversion %>%
  rename(gene_name = ï..gene_name) %>%
  dplyr::select(gene_name, plate_coordinate, plate_number) %>%
  unite(effector_name, c(plate_coordinate, plate_number),
        sep = " ", remove = FALSE) %>%
  dplyr::select(gene_name, effector_name)

# Lists characterised effectors included in the phylogeny
extras <- c("PiPexRD2", "PiAvr2", "PiPexRD41", "PiAvrblb2", "PiAvrblb1",
            "PiAvr3a", "PsAvr1b", "PsPSR2", "PiPSR2", "PiPexRD1", "PiAvr4")

# Adds gene and effector names to the conversion table
conversion <- conversion %>% add_row(gene_name=extras,
            effector_name=extras)

# Merge the tree file and conversion table to generate temporary labels
temp_labels <- as.data.frame(tree$tip.label)
colnames(temp_labels) <- "labels"
temp_labels <- left_join(x=temp_labels, y=conversion, 
              by=c("labels" = "gene_name"))

# Append new labels to the tree
tree$tip.label <- temp_labels$effector_name

# Define a function to process bootstrap values for plotting using ggtree
bootstrap_function <- function(input_tree) {
  input <- ggtree(input_tree) # Turns tree into ggtree object
  b=input$data
  b=b[!b$isTip,] # Extracts data and removes tip labels
  b=b[!is.na(b$label),] # Removes NAs from labels
  b$label=as.numeric(b$label) # Makes the labels numeric
  return(b)
  }

### Plotting the trees #########################################################

### Plot phylogeny of all 130 screened RXLRs for Figure S3A ###
boots <- bootstrap_function(tree) # Run bootstrap function on total tree
total <- ggtree(tree) +
  geom_tiplab(size=2) + # Adds tip labels at size 2
  # Plots bootstrap values repelling each other 
  geom_text_repel(data=boots, aes(label=label), size=2)
total

# Exports the total tree as a .png
# May need to adjust width and height values to fit the entire tree on
# one page
png(file="File_name.png", units="in", width=12, height=12, res=400)
total
dev.off()

# Define a subset tree plotting function
# Requires 4 inputs: tree file (tree), the effector of interest (target),
# number of levels back from the effector of interest to plot (levels),
# and an x-axis adjustment (x_adjust) to prevent labels being cut off
subset_tree <- function(tree, target, levels, x_adjust) {
  # Subset tree according to levels set
  subset <- tree_subset(tree, target, levels_back=levels) 
  bootstrap <- bootstrap_function(subset) # Process bootstrap values as labels
  plot <- ggtree(subset) + # Plot the tree with tip labels and bootstrap labels
    geom_tiplab(size=2.5) +
    geom_text_repel(data=bootstrap, aes(label=label), size=2.5) +
    xlim(0, max(subset$edge.length)+x_adjust) # Adjust xlim according to input
  return(plot)
}

### Plot the subsetted trees for Figure 6 ###
# Avrblb2 tree
avrblb2 <- subset_tree(tree, "PiAvrblb2", 4, 0.4)
avrblb2

# Avr3a tree
avr3a <- subset_tree(tree, "PiAvr3a", 3, 1.1)
avr3a

# Avr2 and C8 P3 tree
avr2 <- subset_tree(tree, "PiAvr2", 3, 1.2)
avr2

### Plot the subsetted trees of the five candidates for Figure S3B-F ###
# C11 P3 tree
c11p3 <- subset_tree(tree, "C11 P3", 2, 1.5)
c11p3

# E6 P3 tree
e6p3 <- subset_tree(tree, "E6 P3", 2, 0.2)
e6p3

# B4 P2 tree
b4p2 <- subset_tree(tree, "B4 P2", 2, 1.05)
b4p2

# E6 P1 tree
e6p1 <- subset_tree(tree, "E6 P1", 3, 1.5)
e6p1

# F3 P1 tree
f3p1 <- subset_tree(tree, "F3 P1", 3, 0.5)
f3p1

# Combine subsetted trees onto one plot
trees <- ggarrange(avrblb2, avr3a, avr2,
          ncol=2, nrow=2)
trees

# Exports the plot as a .png
png(file="file_name.png", units="in", width=6, height=12, res=400)
trees
dev.off()
