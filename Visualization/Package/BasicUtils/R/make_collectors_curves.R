#' Generate Collector's Curves Plot
#' Author: Jacob Agerbo Rasmussen
#' The `make_collectors_curves` function generates a collector's curves plot based on the provided data. 
#' It calculates rarefaction curves by collecting unique observations for a specific group. 
#' The function takes the following parameters:
#' 
#' - `data`: The dataset containing the observations.
#' - `iterations`: The number of iterations to perform for each group (default: 5).
#' - `threads`: The number of parallel threads to use (default: 2).
#' 
#' The function utilizes the `tidyverse`, `reshape2`, `cowplot`, `doParallel`, and `future` packages.
#' It also defines a custom theme called `theme_ridges` for the plot.
#' 
#' Example usage:
#' ```R
#' # Create a sample data frame and group variable
#' set.seed(1234)
#' data <- matrix(rpois(1000, 5), nrow = 100)
#' colnames(data) <- paste0("Group", 1:ncol(data))
#' group_var <- rep(colnames(data), each = nrow(data))
#'
#' # Generate plots for the sample data frame and group variable
#' generate_plots(data, group_var, iterations = 50)
#' ```
#'
#' @param data The dataset containing the observations
#' @param iterations The number of iterations to perform for each group (default: 5)
#' @param threads The number of parallel threads to use (default: 2)
#' @return A collector's curves plot
#' @import tidyverse
#' @import reshape2
#' @import cowplot
#' @import doParallel
#' @import future
#' @export

# Define the function
make_collectors_curves <- function(data, interations = 5, threads = 2){
  suppressWarnings({ 
    suppressPackageStartupMessages(library(tidyverse, quietly = TRUE,
            warn.conflicts = FALSE))
    suppressPackageStartupMessages(library(reshape2))
    suppressPackageStartupMessages(library(cowplot))
    suppressPackageStartupMessages(library(doParallel))
    suppressPackageStartupMessages(library(future))

      # get threads  
    cl <- makePSOCKcluster(threads)
    registerDoParallel(cl)

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
  # Reshape the data to long format
  data_long <- melt(t(data))
  
  # Define a function to collect data for a specific group
  collect <- function(data, group) {
    data %>%
      filter(Var1 == group) %>%
      uncount(value) %>%
      sample_n(n()) %>%
      mutate(observation = row_number()) %>%
      arrange(Var2, observation) %>%
      group_by(Var2) %>%
      mutate(distinct = row_number() == 1) %>%
      ungroup() %>%
      arrange(observation) %>%
      mutate(s = cumsum(distinct)) %>%
      select(observation, s) 
  }
  
  # Remove rows with NA values from data_long
  data_long <- na.omit(data_long)
  
  # Create an empty list for storing plots
  plot_list <- list()
  
  # Select 4 random samples
  random_sample_list <- slice_sample(data_long, n = 4, replace = TRUE) %>%
    select(Var1) 
  random_sample_list <- random_sample_list$Var1 %>%
    as.character()
  
  # Iterate over the selected random samples
  for (i in random_sample_list){
    sample <- i
    
    # Generate curves for collecting ARGs
    collect_curves <- map_dfr(1:interations, ~collect(data_long, sample), .id = "iteration")
    
    # Calculate rarefaction curve
    rarefraction_curve <- collect_curves %>%
      group_by(observation) %>%
      summarize(
        r = mean(s)
      )
    
    # Create and store the plot for the current sample
    i.plot <- collect_curves %>%
      ggplot(aes(x = observation, y = s, group = iteration)) + 
      geom_line(color = "grey80") +
      theme_ridges() +
      geom_line(data = rarefraction_curve, 
                aes(x = observation, y = r), linewidth = 1, inherit.aes = FALSE) +
      xlab("Observations") +
      ylab("Unique Observations") +
      ggtitle(i) +
      labs(color = group_var)
    
    plot_list[[i]] <- i.plot
  }
  
  # Return the list of plots
  plot <- cowplot::plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]], ncol = 2)
  
  return(plot)
  

 ## stop parallels
  suppressWarnings({
    suppressMessages({
      stopCluster(cl)})
    })
    
}

# Example usage
# Create a sample data frame and group variable
#set.seed(1234)
#data <- matrix(rpois(1000, 5), nrow = 100)
#colnames(data) <- paste0("Group", 1:ncol(data))
#group_var <- rep(colnames(data), each = nrow(data))

# Generate plots for the sample data frame and group variable
#generate_plots(data, group_var, interations = 50)