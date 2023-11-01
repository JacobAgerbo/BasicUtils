#!/usr/bin/env Rscript
# R Script for Generating PCA Plot
# Author: Jacob Agerbo Rasmussen
# Load necessary packages
suppressWarnings({ 
  suppressPackageStartupMessages(library(tidyverse, quietly = TRUE,
          warn.conflicts = FALSE))
  theme_ridges <- function(font_size = 14, font_family = "", line_size = .5, grid = TRUE, center_axis_labels = FALSE) {
    half_line <- font_size / 2
    small_rel <- 0.857
    small_size <- small_rel * font_size
    color <- "grey90"
    
    if (grid) {
      panel.grid.major <- element_line(colour = color, linewidth = line_size)
      axis.ticks       <- element_line(colour = color, linewidth = line_size)
      axis.ticks.y     <- axis.ticks
    }
    else {
      panel.grid.major <- element_blank()
      axis.ticks       <- element_line(colour = "black", linewidth = line_size)
      axis.ticks.y     <- element_blank()
    }
    
    if (center_axis_labels) {
      axis_just <- 0.5
    }
    else {
      axis_just <- 1.0
    }
    
    theme_grey(base_size = font_size, base_family = font_family) %+replace%
      theme(
        rect              = element_rect(fill = "transparent", colour = NA, color = NA, linewidth = 0, linetype = 0),
        text              = element_text(family = font_family, face = "plain", colour = "black",
                                         size = font_size, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = .9,
                                         margin = margin(), debug = FALSE),
        axis.text         = element_text(colour = "black", size = small_size),
        #axis.title        = element_text(face = "bold"),
        axis.text.x       = element_text(margin = margin(t = small_size / 4), vjust = 1),
        axis.text.y       = element_text(margin = margin(r = small_size / 4), hjust = 1, vjust = 0),
        axis.title.x      = element_text(
          margin = margin(t = small_size / 2, b = small_size / 4),
          hjust = axis_just
        ),
        axis.title.y      = element_text(
          angle = 90,
          margin = margin(r = small_size / 2, l = small_size / 4),
          hjust = axis_just
        ),
        axis.ticks        = axis.ticks,
        axis.ticks.y      = axis.ticks.y,
        axis.line         = element_blank(),
        legend.key        = element_blank(),
        legend.key.size   = grid::unit(1, "lines"),
        legend.text       = element_text(size = rel(small_rel)),
        legend.justification = c("left", "center"),
        panel.background  = element_blank(),
        panel.border      = element_blank(),
        # make grid lines
        panel.grid.major  = panel.grid.major,
        panel.grid.minor  = element_blank(),
        strip.text        = element_text(size = rel(small_rel)),
        strip.background  = element_rect(fill = "grey80", colour = "grey50", linewidth = 0),
        plot.background   = element_blank(),
        plot.title        = element_text(face = "bold",
                                         size = font_size,
                                         margin = margin(b = half_line), hjust = 0),
        plot.subtitle     = element_text(size = rel(small_rel),
                                         hjust = 0, vjust = 1,
                                         margin = margin(b = half_line * small_rel)),
        plot.caption      = element_text(size = rel(small_rel),
                                         hjust = 1, vjust = 1,
                                         margin = margin(t = half_line * small_rel)),
        plot.margin       = margin(half_line, font_size, half_line, half_line),
        
        complete = TRUE
      )
  }
})

# Function to generate PCA plot
 generate_pca_plot <- function(data, group_var = NULL, method = "euclidean", 
                               palette = "Dark2", alpha = 0.75, lg.position = "bottom", plot.centroids = FALSE,
                               plot_title = "PCA Plot") {
   suppressWarnings({ 
     library(tidyverse, quietly = TRUE,
             warn.conflicts = FALSE)
     theme_ridges <- function(font_size = 14, font_family = "", line_size = .5, grid = TRUE, center_axis_labels = FALSE) {
       half_line <- font_size / 2
       small_rel <- 0.857
       small_size <- small_rel * font_size
       color <- "grey90"
       
       if (grid) {
         panel.grid.major <- element_line(colour = color, linewidth = line_size)
         axis.ticks       <- element_line(colour = color, linewidth = line_size)
         axis.ticks.y     <- axis.ticks
       }
       else {
         panel.grid.major <- element_blank()
         axis.ticks       <- element_line(colour = "black", linewidth = line_size)
         axis.ticks.y     <- element_blank()
       }
       
       if (center_axis_labels) {
         axis_just <- 0.5
       }
       else {
         axis_just <- 1.0
       }
       
       theme_grey(base_size = font_size, base_family = font_family) %+replace%
         theme(
           rect              = element_rect(fill = "transparent", colour = NA, color = NA, linewidth = 0, linetype = 0),
           text              = element_text(family = font_family, face = "plain", colour = "black",
                                            size = font_size, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = .9,
                                            margin = margin(), debug = FALSE),
           axis.text         = element_text(colour = "black", size = small_size),
           #axis.title        = element_text(face = "bold"),
           axis.text.x       = element_text(margin = margin(t = small_size / 4), vjust = 1),
           axis.text.y       = element_text(margin = margin(r = small_size / 4), hjust = 1, vjust = 0),
           axis.title.x      = element_text(
             margin = margin(t = small_size / 2, b = small_size / 4),
             hjust = axis_just
           ),
           axis.title.y      = element_text(
             angle = 90,
             margin = margin(r = small_size / 2, l = small_size / 4),
             hjust = axis_just
           ),
           axis.ticks        = axis.ticks,
           axis.ticks.y      = axis.ticks.y,
           axis.line         = element_blank(),
           legend.key        = element_blank(),
           legend.key.size   = grid::unit(1, "lines"),
           legend.text       = element_text(size = rel(small_rel)),
           legend.justification = c("left", "center"),
           panel.background  = element_blank(),
           panel.border      = element_blank(),
           # make grid lines
           panel.grid.major  = panel.grid.major,
           panel.grid.minor  = element_blank(),
           strip.text        = element_text(size = rel(small_rel)),
           strip.background  = element_rect(fill = "grey80", colour = "grey50", linewidth = 0),
           plot.background   = element_blank(),
           plot.title        = element_text(face = "bold",
                                            size = font_size,
                                            margin = margin(b = half_line), hjust = 0),
           plot.subtitle     = element_text(size = rel(small_rel),
                                            hjust = 0, vjust = 1,
                                            margin = margin(b = half_line * small_rel)),
           plot.caption      = element_text(size = rel(small_rel),
                                            hjust = 1, vjust = 1,
                                            margin = margin(t = half_line * small_rel)),
           plot.margin       = margin(half_line, font_size, half_line, half_line),
           
           complete = TRUE
         )
     }
   })
  # Perform PCA
  dist_df <- dist(data, method = method)
  pca <- prcomp(data, scale. = TRUE)
  pca_df <- as.data.frame(pca$x)
  
  # Get variance explained by each PC
  per_var_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100

  # Add group variable if provided
  if (!is.null(group_var)) {
    pca_df$group <- group_var
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
  return(pca_plot)
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

args <- commandArgs(trailingOnly = TRUE)
data <- args[1]
group_var <- args[2]
method <- args[3]
palette <- args[4]
alpha <- args[5]
lg.position <- args[6]
plot.centroids <- args[7]
plot_title <- args[8]