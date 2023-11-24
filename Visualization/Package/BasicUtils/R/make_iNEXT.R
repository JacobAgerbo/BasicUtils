#' BasicUtils - Generate iNEXT Sample-Size-Based Richness Plot
#' 
#' The `do_iNEXT` function generates a sample-size-based richness plot using the iNEXT package. 
#' The function takes the following parameters:
#' 
#' - `data`: The dataset containing the observations.
#' - `sample_data`: The dataset containing the sample data (default: NULL).
#' - `group_var`: The variable used to group the samples (default: NULL).
#' - `table`: A logical value indicating whether to return the extrapolated data as a table (default: FALSE).
#' - `threads`: The number of parallel threads to use (default: 2).
#' 
#' The function loads the required packages (`tidyverse`, `iNEXT`, `ggtext`, `ggpubr`, `doParallel`, and `future`) and calculates iNEXT diversity for the provided data. 
#' If a group variable is provided, the function aggregates the data by summing values across each row and performs iNEXT analysis on the relevant columns of the aggregated data. 
#' If no group variable is provided, the function sums the values across the columns and performs iNEXT analysis on the resulting data.
#'
#' The function then generates an iNEXT plot type 1 with reordered facets, sets the theme to minimal, and adds a dashed line representing the intercept and slope. 
#' The plot title, y-axis label, and caption are customized using the `ggtext` package. 
#'
#' If the `table` parameter is set to `TRUE`, the function returns both the plot and an extrapolated data table as a list. Otherwise, it returns only the plot.
#' 
#' Finally, the function stops the parallel cluster.
#'
#' Example usage:
#' ```R
#' # Load data
#' data(iris)
#'
#' # Generate iNEXT plot for iris dataset
#' do_iNEXT(iris[,1:4])
#' ```
#'
#' @param data The dataset containing the observations
#' @param sample_data The dataset containing the sample data (default: NULL)
#' @param group_var The variable used to group the samples (default: NULL)
#' @param table A logical value indicating whether to return the extrapolated data as a table (default: FALSE)
#' @param threads The number of parallel threads to use (default: 2)
#' @import tidyverse
#' @import iNEXT
#' @import ggtext
#' @import ggpubr
#' @import doParallel
#' @import future



# Define function to create iNEXT figure
do_iNEXT <- function(data, sample_data = NULL, group_var = NULL, table=FALSE, threads=2) {
  #Load dependencies
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(iNEXT))
  suppressPackageStartupMessages(library(ggtext))
  suppressPackageStartupMessages(library(ggpubr))
  suppressPackageStartupMessages(library(doParallel))
  suppressPackageStartupMessages(library(future))

  # get threads  
    cl <- makePSOCKcluster(threads)
    registerDoParallel(cl)

  # Calculate iNEXT diversity
  # Add group variable if provided
  suppressWarnings({if (!is.null(group_var)) {
    n <- sample_data %>%
      select(all_of(group_var)) %>%
      table %>%
      as_tibble()
    
    # Aggregate rows of iNEXT.table by summing values across each row
    colnames(data) <- sample_data[,group_var]
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
    
    Extrapolated_data <- data.frame(
      Assemblage = out.inc[["iNextEst"]][["size_based"]][["Assemblage"]],
      qD = out.inc[["iNextEst"]][["size_based"]][["qD"]],
      Method = out.inc[["iNextEst"]][["size_based"]][["Method"]]) %>%
      filter(Method != "Rarefaction") %>%
      group_by(Assemblage, Method) %>%
      summarise(mean_qD = mean(qD),
                SD = sd(qD))
    
  } else {
    data.agg <- as.data.frame(colSums(data))
    data.agg <- ifelse(data.agg > 0,1,0)
    n <- length(data)
    iNEXT.table <- rbind(n, data.agg)
    out.inc <- iNEXT(iNEXT.table, q = 0, datatype = "incidence_freq")
    
    Extrapolated_data <- data.frame(
      Assemblage = out.inc[["iNextEst"]][["size_based"]][["Assemblage"]],
      qD = out.inc[["iNextEst"]][["size_based"]][["qD"]],
      Method = out.inc[["iNextEst"]][["size_based"]][["Method"]]) %>%
      filter(Method != "Rarefaction") %>%
      group_by(Assemblage, Method) %>%
      summarise(mean_qD = mean(qD),
                SD = sd(qD))
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

  })

  if (isTRUE(table)) {
    # Code to be executed if condition is true
    return(list(plot = plot, table = Extrapolated_data))
  } else {
    # Code to be executed if condition is false
    return(list(plot = plot))
  }
 ## stop parallels
  stopCluster(cl) 
}