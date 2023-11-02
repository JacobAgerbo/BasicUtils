#!/usr/bin/env Rscript
# R Script for Generating PCA Plot
# Author: Jacob Agerbo Rasmussen
# Load required libraries

# Define a function to generate a boxplot
generate_boxplot <- function(data, 
                             group_var=NULL, 
                             plot_title=NULL, 
                             palette="Dark2", 
                             violin = FALSE,
                             jitter = FALSE,
                             test = FALSE) {
  
  # Make theme
  suppressPackageStartupMessages(library(ggpubr))
  suppressPackageStartupMessages(library(rstatix))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(reshape2))
  suppressPackageStartupMessages(library(RColorBrewer))
  theme_ridges <- source("https://raw.githubusercontent.com/wilkelab/ggridges/master/R/theme.R")
  # Read data
  # Add group variable if provided
  if (!is.null(group_var)) {
    data$group <- group_var
  }
  

  suppressWarnings({
    data.melt <- melt(data)  
  })
  
  if (!is.null(plot_title)) {
    plot_title <- ""
  }
  
  # Create a boxplot using ggplot
  
  if (isTRUE(violin)) {
    boxplot <- ggplot(data.melt, aes(x = group, y = value, fill=group)) +
      geom_violin(alpha = 0.75) +
      ggtitle(plot_title) +
      theme_ridges() +
      scale_fill_brewer(palette = palette) +
      xlab("")
    
    if (isTRUE(test)) {
      # Statistical test
      p <- data.melt %>% levene_test(value ~ group) %>%
        select(p) %>%
        as_vector()
      
      # Choose between Tukey HSD test and Dunn test based on p-value
      if (p > 0.05) {
        test <- paste("Tukey HSD")
        stat.test <- data.melt %>% tukey_hsd(value ~ group) %>%
          adjust_pvalue(method = "bonferroni") %>%
          add_significance("p.adj") %>%
          add_xy_position(x = "group") %>%
          mutate(p.format = p_format(p.adj, accuracy = 0.001, leading.zero = FALSE))
      } else {
        test <- paste("Dunn Test")
        stat.test <- data.melt %>% dunn_test(value ~ group) %>%
          adjust_pvalue(method = "bonferroni") %>%
          add_significance("p.adj") %>%
          add_xy_position(x = "group") %>%
          mutate(p.format = p_format(p.adj, accuracy = 0.001, leading.zero = FALSE))
      }
      
      # Visualization
      boxplot <- ggviolin(data.melt, x = "group", 
                           y = "value", 
                           fill = "group",
                           alpha = 0.75)  +
        theme_ridges() +
        scale_fill_brewer(palette = palette) +
        xlab("") +
        labs(caption = test) +
        stat_pvalue_manual(stat.test, label = "p.format") +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.10)))
    }
  
    }  else {
    boxplot <- ggplot(data.melt, aes(x = group, y = value, fill=group)) +
      geom_boxplot(alpha = 0.75) +
      ggtitle(plot_title) +
      theme_ridges() +
      scale_fill_brewer(palette = palette) +
      xlab("")
    
    if (isTRUE(test)) {
      # Statistical test
      p <- data.melt %>% levene_test(value ~ group) %>%
        select(p) %>%
        as_vector()
      
      # Choose between Tukey HSD test and Dunn test based on p-value
      if (p > 0.05) {
        test <- paste("Tukey HSD")
        stat.test <- data.melt %>% tukey_hsd(value ~ group) %>%
          adjust_pvalue(method = "bonferroni") %>%
          add_significance("p.adj") %>%
          add_xy_position(x = "group") %>%
          mutate(p.format = p_format(p.adj, accuracy = 0.001, leading.zero = FALSE))
      } else {
        test <- paste("Dunn Test")
        stat.test <- data.melt %>% dunn_test(value ~ group) %>%
          adjust_pvalue(method = "bonferroni") %>%
          add_significance("p.adj") %>%
          add_xy_position(x = "group") %>%
          mutate(p.format = p_format(p.adj, accuracy = 0.001, leading.zero = FALSE))
      }
      
      # Visualization
      boxplot <- ggboxplot(data.melt, x = "group", 
                           y = "value", 
                           fill = "group",
                           alpha = 0.75)  +
        theme_ridges() +
        scale_fill_brewer(palette = palette) +
        xlab("") +
        labs(caption = test) +
        stat_pvalue_manual(stat.test, label = "p.format") +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.10)))
    }
  } 
  
  if (isTRUE(jitter)) {
    boxplot <- boxplot +
      geom_jitter()
  }  

  
  # Save the plot as a PNG file
  return(boxplot)
}

# Test function
#set.seed(1234)
#data <- rbind(matrix(rnorm(100, mean = 20, sd = 2), nrow = 100),
#              matrix(rnorm(100, mean = 25, sd = 5), nrow = 100)) %>%
#  as_tibble()
#group_var <- rep(c("Group A", "Group B"), each = 100)

# Test function
#generate_boxplot(data, group_var = group_var, jitter = TRUE,
#                 violin = FALSE, test = TRUE)

