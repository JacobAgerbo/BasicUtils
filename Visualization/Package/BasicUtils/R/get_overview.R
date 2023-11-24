#!/usr/bin/env Rscript

source("https://raw.githubusercontent.com/JacobAgerbo/Basic_Utils/main/Visualization/make_test_data.R")

data_list <- make_test_data(100)
data <- data_list$data
sample_data <- data_list$sample_data


get_overview <- function(data = data, sample_data = sample_data, group_var=NULL){
  suppressWarnings({ 
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(doParallel))
    suppressPackageStartupMessages(library(future))
    suppressPackageStartupMessages(library(reshape2))
    suppressPackageStartupMessages(library(cowplot))
    suppressPackageStartupMessages(library(ggpubr))
    suppressPackageStartupMessages(library(wesanderson))
    suppressPackageStartupMessages(library(hilldiv))
    red <- "\033[31m"
    reset <- "\033[0m"
    green <- "\033[32m"
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
        cat(paste(green,"A ridgeline plot",reset,"\n","Is a visualization that provides a compact representation of the distribution of a numeric variable across multiple groups. \n It consists of stacked density curves or histograms, allowing for easy comparison between groups. \n Ridgeline plots are useful for visualizing changes in distributions over space or time and identifying differences in \n central tendency, spread, and multimodality between groups.", sep = ""))
      }
      

      
      # Empirical Cumulative Distribution Function
      ECDF <- data %>% 
        t() %>%
        as_tibble() %>%
        dplyr::mutate(group = sample_data[,group_var], .before = 1)  %>%
        reshape2::melt() %>%
        ggplot(aes(x = value)) +
        stat_ecdf() +
        facet_wrap(~group, scales = "free", nrow = 1) +
        labs(title = "Empirical Cumulative Distribution Function") +
        ggridges::theme_ridges() 
      
      # Function to describe ECDF
      describe_ecdf <- function() {
        cat(paste("\033[32m", "The Empirical Cumulative Distribution Function (ECDF)", "\033[0m", "\n", "is a non-parametric way to describe the distribution of a numeric variable. \nIt shows the proportion or percentage of observations that are less than or equal to a given value. \nThe ECDF is useful for visualizing the shape of the distribution, including its central tendency, spread, and skewness. \nIt can also be used to compare two or more distributions or to test hypotheses about the underlying distribution.", sep = ""))
      }
      
      # Q-Q plot
      QQplot <- data %>% 
        t() %>%
        as_tibble() %>%
        dplyr::mutate(group = sample_data[,group_var], .before = 1)  %>%
        reshape2::melt() %>%
        ggplot(aes(sample = value)) +
        stat_qq() +
        facet_wrap(~group, scales = "free", nrow = 1) +
        labs(title = "Q-Q Plot of Variables") +
        ggridges::theme_ridges()
      
      # Function to describe Q-Q plot
      describe_qq_plot <- function() {
        cat(paste("\033[32m", "A Q-Q plot", "\033[0m", "\n", "is a graphical tool used to assess whether a dataset follows a particular theoretical distribution, such as the normal distribution. \nIt compares the quantiles of the observed data against the quantiles expected from the theoretical distribution. \nIf the data points fall approximately along a straight line, it suggests that the data can be reasonably modeled by the theoretical distribution. \nQ-Q plots are useful for identifying deviations from the assumed distribution, detecting outliers, and assessing the goodness-of-fit of statistical models.", sep = ""))
      }
      
      
      # Shapiro test
      lshap <- apply(data,1, shapiro.test)
      shapiro_test <- data.frame(do.call(rbind, lshap))
      shap_plot <- shapiro_test %>%
        rownames_to_column(var = "feature") %>%
        select(feature, p.value) %>%
        dplyr::mutate(p.value = as.numeric(p.value)) %>%
        dplyr::mutate(p.value = round(p.value,3)) %>%
        dplyr::filter(p.value < 0.15) %>%
        dplyr::mutate(Distribution = ifelse(p.value > 0.05, "Normal", "Not Normal")) %>%
        ggplot(aes(x = reorder(feature,p.value), y = p.value, fill = Distribution)) + geom_bar(stat = "identity") + 
        geom_hline(yintercept = 0.05) + 
        ggridges::theme_ridges() + 
        xlab("Features") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              legend.position = "bottom") +
        scale_fill_manual(values = palette) + 
        ggtitle("Distribution of individual Features")
      
      # Function to describe Shapiro-Wilk test
      describe_shapiro_wilk_test <- function() {
        cat(paste("\033[32m", "The Shapiro-Wilk test", "\033[0m", "\n", "is a statistical test used to determine whether a dataset follows a normal distribution. \nIt is based on the comparison between the observed data and the expected values from a normal distribution. \nThe test provides a p-value that indicates the level of evidence against the null hypothesis of normality. \nIf the p-value is below a chosen significance level (e.g., 0.05), it suggests that the data significantly deviate from a normal distribution. \nThe Shapiro-Wilk test is commonly used in exploratory data analysis and as a preliminary step in many statistical analyses.", sep = ""))
      }
      
      
      plot <- cowplot::plot_grid(ridges,ECDF,QQplot,shap_plot, ncol = 2)
      
      # Call the function to describe ridgeline plots
      describe_ridgeline_plots()
      cat(paste("\n"))
      cat(paste("\n"))
      # Call the function to describe ECDF
      describe_ecdf()
      cat(paste("\n"))
      cat(paste("\n"))
      # Call the function to describe Q-Q plot
      describe_qq_plot()
      cat(paste("\n"))
      cat(paste("\n"))
      # Call the function to describe Shapiro-Wilk test
      describe_shapiro_wilk_test()
      cat(paste("\n"))
      cat(paste("\n"))
      
    } else {
      
      
      {cat(paste(red,"\n","Whoops!",reset, sep = ""))
      cat(paste(reset,"\n","Maybe add a `group_var` from your sample information dataset", sep = ""))
      cat(paste(reset,"\n","then you will get more descriptive plots","\n", sep = ""))}
      cat(paste("\n"))
      cat(paste("\n"))
      # Shapiro test
      lshap <- apply(data,1, shapiro.test)
      shapiro_test <- data.frame(do.call(rbind, lshap))
      shap_plot <- shapiro_test %>%
        rownames_to_column(var = "feature") %>%
        select(feature, p.value) %>%
        dplyr::mutate(p.value = as.numeric(p.value)) %>%
        dplyr::mutate(p.value = round(p.value,3)) %>%
        dplyr::filter(p.value < 0.15) %>%
        dplyr::mutate(Distribution = ifelse(p.value > 0.05, "Normal", "Not Normal")) %>%
        ggplot(aes(x = reorder(feature,p.value), y = p.value, fill = Distribution)) + geom_bar(stat = "identity") + 
        geom_hline(yintercept = 0.05) + 
        ggridges::theme_ridges() + 
        xlab("Features") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              legend.position = "bottom") +
        scale_fill_manual(values = palette) +
        ggtitle("Distribution of individual Features")
      
      # Function to describe Shapiro-Wilk test
      describe_shapiro_wilk_test <- function(){
        cat(paste("\033[32m", "The Shapiro-Wilk test", "\033[0m", "\n", "is a statistical test used to determine whether a dataset follows a normal distribution. \nIt is based on the comparison between the observed data and the expected values from a normal distribution. \nThe test provides a p-value that indicates the level of evidence against the null hypothesis of normality. \nIf the p-value is below a chosen significance level (e.g., 0.05), it suggests that the data significantly deviate from a normal distribution. \nThe Shapiro-Wilk test is commonly used in exploratory data analysis and as a preliminary step in many statistical analyses.", sep = ""))
      }
      
      # Call the function to describe Shapiro-Wilk test 
      describe_shapiro_wilk_test()
      
      plot <- shap_plot
      
    }
    
    return(plot)
    
  } 
  )
}
  

#get_overview(data = data, sample_data = sample_data, group_var = "O_Group")
