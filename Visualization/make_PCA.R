#!/usr/bin/env Rscript
# R Script for Generating PCA Plot
# Author: Jacob Agerbo Rasmussen
# Load necessary packages
# Function to generate PCA plot
 generate_pca_plot <- function(data, sample_data = NULL, group_var = NULL, method = "euclidean", 
                               palette = "Dark2", alpha = 0.75, lg.position = "bottom", plot.centroids = FALSE,
                               plot_title = "PCA Plot") {
                              
   suppressWarnings({ 
     library(tidyverse, quietly = TRUE,
             warn.conflicts = FALSE)
     theme_ridges <- source("https://raw.githubusercontent.com/wilkelab/ggridges/master/R/theme.R")
   })
  # Perform PCA
  dist_df <- dist(data, method = method)
  pca <- prcomp(data, scale. = TRUE)
  pca_df <- as.data.frame(pca$x)
  
  # Get variance explained by each PC
  per_var_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100

  # Add group variable if provided
  if (!is.null(group_var)) {
    #data$group <- sample_data[,group_var]
    pca_df <- pca_df %>%
    t() %>%
    as_tibble() %>%
    mutate(group = sample_data[,group_var], .before = 1)
  }
  
  # Create PCA plot using ggplot
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
    geom_point(alpha = alpha, size = 3)  +
    xlab(paste("PC 1 (",round(per_var_explained[1],2)," %)", sep = "")) +
    ylab(paste("PC 2 (",round(per_var_explained[2],2)," %)", sep = "")) +
    ggtitle(plot_title) +
    scale_color_brewer(palette = palette) +
    theme_ridges() +
    theme(legend.position = lg.position)

  if (isTRUE(plot.centroids)) {
    # calculate centroids
    centroid <- pca_df %>% 
      dplyr::select("PC1", "PC2", "group") %>%
      group_by(group) %>%
      summarise(mean_PC1 = mean(PC1),
                mean_PC2 = mean(PC2)) 
    
    pca_plot <- pca_plot + geom_point(data = centroid,
                                      aes(x = mean_PC1,
                                          y = mean_PC2,
                                          fill = group),
                                      shape=18,
                                      alpha = 0.85,
                                      size = 7.5) +
      scale_fill_brewer(palette = palette)
  }
  
  
  # Return the plot
  pdf(paste(group_var, '.pdf', sep=""), width=8, height=8)
    pca_plot
    dev.off()
}


# Main script execution
# Your code goes here
#data <- rbind(matrix(rnorm(500, mean = 20, sd = 2), nrow = 100),
#              matrix(rnorm(500, mean = 25, sd = 5), nrow = 100))
#group_var <- rep(c("Group A", "Group B"), each = 100)
# Call the function without group variable

#generate_pca_plot(data, 
#                  group_var = group_var, 
#                  palette ="Dark2", alpha = 0.5, 
#                  plot.centroids = TRUE)