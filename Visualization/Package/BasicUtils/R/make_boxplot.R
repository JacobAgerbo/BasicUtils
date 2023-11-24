#' BasicUtils - Generate a boxplot
#' 
#' The `generate_boxplot` function is used for generating a boxplot of a provided dataset. 
#' If `sample_data` and `group_var` are specified, the function generates a grouped boxplot showing the distribution of values in each group. 
#' If `violin` is set to TRUE, the function generates a violin plot instead of a boxplot. 
#' If `jitter` is set to TRUE, the function adds jitter points to the plot. 
#' If `test` is set to TRUE, the function performs a statistical test (Levene's test) and adds p-values and significance labels to the plot.
#'
#' @param data The dataset to be analyzed
#' @param sample_data The sample information dataset
#' @param group_var The grouping variable to be used for generating a grouped boxplot
#' @param plot_title The title of the plot
#' @param palette The color palette to be used for the plot
#' @param violin A logical value indicating whether to generate a violin plot instead of a boxplot
#' @param jitter A logical value indicating whether to add jitter points to the plot
#' @param test A logical value indicating whether to perform a statistical test and add p-values and significance labels to the plot
#' @return A boxplot or violin plot
#' @import ggpubr
#' @import rstatix
#' @import tidyverse
#' @import reshape2
#' @import RColorBrewer
#' @export
#' @examples
#' # Generate a simple boxplot
#' generate_boxplot(data, group_var = group_var, jitter = TRUE,
#'                violin = FALSE, test = TRUE)

# Define a function to generate a boxplot
generate_boxplot <- function(data,
                            sample_data = NULL,
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
    #data$group <- sample_data[,group_var]
    data <- data %>%
    t() %>%
    as_tibble() %>%
    mutate(group = sample_data[,group_var], .before = 1)
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
      geom_jitter(aes(fill = group), alpha = 0.15)
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

