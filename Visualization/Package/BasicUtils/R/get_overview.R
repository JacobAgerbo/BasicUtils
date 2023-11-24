#' Generate descriptive plots of a provided dataset
#' Author: Jacob Agerbo Rasmussen
#' 
#' The `get_overview` function is used for generating descriptive plots of a provided dataset. 
#' If `group_var` is specified, the function generates four plots: a ridgeline plot, an empirical cumulative distribution function (ECDF) plot, a Q-Q plot, and a distribution plot of individual features. 
#' If `group_var` is not specified, the function generates a distribution plot of individual features.
#'
#' @param data The dataset to be analyzed
#' @param sample_data The sample information dataset
#' @param group_var The variable used for grouping the data
#' @return A list of plots
#' @import tidyverse
#' @import reshape2
#' @import cowplot
#' @import ggpubr
#' @import wesanderson
#' @import hilldiv
#' @export
get_overview <- function(data = data, sample_data = sample_data, group_var=NULL){
  
  # Suppress warnings and load required packages
  suppressWarnings({ 
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(reshape2))
    suppressPackageStartupMessages(library(cowplot))
    suppressPackageStartupMessages(library(ggpubr))
    suppressPackageStartupMessages(library(wesanderson))
    suppressPackageStartupMessages(library(hilldiv))
    red <- "\033[31m"
    reset <- "\033[0m"
    green <- "\033[32m"
    
    # Check if group variable is specified
    if (!is.null(group_var)) {
      
      # make colouring
      no_colors <- length(unique(sample_data[,group_var]))
      if (no_colors > 10) {
        palette <- distinct_palette(no_colors, pal = "kelly", add = "grey90")
        } else {
        palette <- c("#0B775E","#35274A","#F2300F","#D69C4E","#046C9A","#ECCBAE","#E1BD6D","#EABE94","#ABDDDE","#000000")
      }
      
      # Ridgeline plot
      ridges <- data %>% 
        t() %>%
        as_tibble() %>%
        dplyr::mutate(group = sample_data[,group_var], .before = 1)  %>%
        reshape2::melt() %>%
        ggplot(aes(x = value, y = group, fill = group)) +
        ggridges::geom_density_ridges() +
        labs(title = "Ridgeline Plot") +
        scale_fill_manual(values = palette) +
        ggridges::theme_ridges() +
        theme(legend.position = "bottom") +
        ggtitle(paste("Distribution of numeric values","between groups", sep = "\n"))
      
      # Function to describe ridgeline plots
      describe_ridgeline_plots <- function() {
        cat(paste(green,"A ridgeline plot",reset,"\n","Is a visualization that provides a compact representation of the distribution of a numeric variable across multiple groups. \n It consists of stacked density curves or histograms, allowing for easy comparison between groups. \n", sep=""))
      }
      
      # ECDF plot
      ecdf_plot <- data %>% 
        t() %>%
        as_tibble() %>%
        dplyr::mutate(group = sample_data[,group_var], .before = 1)  %>%
        reshape2::melt() %>%
        ggplot(aes(x = value, color = group)) +
        stat_ecdf(size=1.25) +
        labs(title = "Empirical Cumulative Distribution Function (ECDF) Plot") +
        scale_color_manual(values=palette) +
        theme(legend.position = "bottom") +
        ggtitle(paste("Distribution of numeric values","between groups", sep = "\n"))
      
      # Function to describe ECDF plots
      describe_ecdf_plots <- function() {
        cat(paste(green,"An ECDF plot",reset,"\n","Is a non-parametric way to describe the distribution of a numeric variable. \n The x-axis represents the range of values in the data and the y-axis represents the proportion of values that are less than or equal to each value on the x-axis. \n", sep=""))
      }
      
      # Q-Q plot
      qq_plot <- data %>% 
        t() %>%
        as_tibble() %>%
        dplyr::mutate(group = sample_data[,group_var], .before = 1)  %>%
        reshape2::melt() %>%
        ggplot(aes(sample=value, color=group)) +
        stat_qq(size=1.25) +
        labs(title = "Q-Q Plot") +
        scale_color_manual(values=palette) +
        theme(legend.position = "bottom") +
        ggtitle(paste("Distribution of numeric values","between groups", sep = "\n"))
      
      # Function to describe Q-Q plots
      describe_qq_plots <- function() {
        cat(paste(green,"A Q-Q plot",reset,"\n","Is a graphical tool used to assess whether a dataset follows a particular theoretical distribution, such as the normal distribution. \n If the data follow the theoretical distribution perfectly, the points in the Q-Q plot should fall along a straight line. \n Departures from this line indicate deviations from the theoretical distribution. \n", sep=""))
      }
      
      # Distribution plot of individual features
      dist_plot <- data %>% 
        t() %>%
        as_tibble() %>%
        dplyr::mutate(feature = rownames(.)) %>%
        reshape2::melt(id.vars=c("feature")) %>%
        ggplot(aes(x=value)) +
        geom_density(size=1, alpha=0.7) +
        labs(title = "Distribution Plot of Individual Features") +
        facet_grid(feature ~ .) +
        theme(legend.position = "none") +
        ggtitle(paste("Distribution of individual features","in the dataset", sep = "\n"))
      
      # Function to describe distribution plots of individual features
      describe_dist_plots <- function() {
        cat(paste(green,"A distribution plot of individual features",reset,"\n","Shows the distribution of each feature in the dataset. \n It is useful for identifying outliers and assessing whether the data follow any particular distribution. \n", sep=""))
      }
      
      # Return list of plots and descriptions
      return(list(ridges, describe_ridgeline_plots(), ecdf_plot, describe_ecdf_plots(), qq_plot, describe_qq_plots(), dist_plot, describe_dist_plots()))
    } else {
      
      # Distribution plot of individual features
      dist_plot <- data %>% 
        t() %>%
        as_tibble() %>%
        dplyr::mutate(feature = rownames(.)) %>%
        reshape2::melt(id.vars=c("feature")) %>%
        ggplot(aes(x=value)) +
        geom_density(size=1, alpha=0.7) +
        labs(title = "Distribution Plot of Individual Features") +
        facet_grid(feature ~ .) +
        theme(legend.position = "none") +
        ggtitle(paste("Distribution of individual features","in the dataset", sep = "\n"))
      
      # Function to describe distribution plots of individual features
      describe_dist_plots <- function() {
        cat(paste(green,"A distribution plot of individual features",reset,"\n","Shows the distribution of each feature in the dataset. \n It is useful for identifying outliers and assessing whether the data follow any particular distribution. \n", sep=""))
      }
      
      # Return list of plots and descriptions
      return(list(dist_plot, describe_dist_plots()))
    }
  })
}