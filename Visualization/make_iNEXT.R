#!/usr/bin/env Rscript
# R Script for Generating iNEXT's sample-size-based richness  Plot
# Author: Jacob Agerbo Rasmussen

# Define function to create iNEXT figure
do_iNEXT <- function(data, sample_data = NULL, group_var = NULL) {
  #Load dependencies
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(iNEXT))
  # Calculate iNEXT diversity
  # Add group variable if provided
  suppressWarnings({if (!is.null(group_var)) {
    n <- sample_data %>%
      select(all_of(group_var)) %>%
      table %>%
      as_tibble()
    
    # Aggregate rows of iNEXT.table by summing values across each row
    colnames(data) <- sample_data$Group
    data.agg <- as.matrix(aggregate(t(data), list(rownames(t(data))), FUN = "sum", na.rm = TRUE))
    data.agg <- t(data.agg)
    colnames(data.agg) <- data.agg[1,]
    data.agg <- data.agg[-1,]
    # Set the class of iNEXT.table as "numeric"
    class(data.agg) <- "numeric"
    # make incidence-based
    data.agg <- ifelse(data.agg > 0,1,0)
    # Add the count (n) as the first row in iNEXT.table
    iNEXT.table <- rbind(n$n, data.agg)
    # Perform iNEXT analysis on the relevant columns of iNEXT.table
    out.inc <- iNEXT(iNEXT.table, q = 0, datatype = "incidence_freq")
    
  } else {
    data.agg <- as.data.frame(colSums(data))
    data.agg <- ifelse(data.agg > 0,1,0)
    n <- length(data)
    iNEXT.table <- rbind(n, data.agg)
    out.inc <- iNEXT(iNEXT.table, q = 0, datatype = "incidence_freq")
  }
  
  # Create iNEXT plot type 1 with reordered facets
  plot <- ggiNEXT(out.inc, type = 1, facet.var = "Assemblage") + 
    # Set the theme to minimal
    theme_minimal() +
    # Set the label for the y-axis
    ylab("Richness") +
    # Add a dashed line representing the intercept and slope
    # Set the title and caption for the plot
    labs(caption = "Sample-size-based R/E curves") +
    # Customize the theme settings
    theme(
      plot.title = element_markdown(
        family = "Helvetica",
        face = "bold", size = 16, margin = margin(12, 0, 12, 0)
      ),
      plot.caption = element_markdown(
        size = 7, color = "grey50", margin = margin(12, 0, 6, 0)
      ),
      legend.position = "none"
    )
  return(plot)
  })
}




#data <- data.frame(matrix(sample(0:75, 800, replace = TRUE), ncol = 20))
#data[data < 70] <- 0

#sample_data <- data.frame("Sample" = paste("Sample",c(20),sep = "_"), 
                          "Group" = rep(c("Group A", "Group B"), each = 10))

#colnames(data) <- sample_data$Sample

#do_iNEXT(data = data, sample_data = sample_data, group_var = "Group")

##