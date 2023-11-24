#' BasicUtils - Generate PCA Plot
#'
#' The `generate_pca_plot` function is used to generate a PCA (Principal Component Analysis) plot. It takes several parameters to customize the plot, such as the data, sample information, group variable, color palette, transparency level, legend position, and more.
#'
#' The function performs the following steps:
#' 1. Loads the necessary packages, including `tidyverse`, which provides a collection of packages for data manipulation and visualization.
#' 2. Performs PCA on the input data using the specified method and scaling.
#' 3. Extracts the principal component scores and stores them in a data frame.
#' 4. Calculates the variance explained by each principal component.
#' 5. If a group variable is provided, adds it to the data frame for coloring the points on the PCA plot.
#' 6. Creates a PCA plot using `ggplot`, a popular plotting package in R.
#' 7. Adds visual elements to the plot, such as points, labels, titles, color palettes, and legends.
#' 8. If `plot.centroids` is set to `TRUE`, calculates centroids for each group and adds them to the PCA plot.
#' 9. Returns the PCA plot as the output.
#'
#' Example usage:
#' ```R
#' # Load data
#' data <- read.csv("data.csv")
#' sample_data <- read.csv("sample_info.csv")
#'
#' # Generate PCA plot
#' generate_pca_plot(data, sample_data = sample_data, group_var = "group", method = "euclidean", palette = "Dark2", alpha = 0.75, lg.position = "bottom", plot.centroids = FALSE, plot_title = "PCA Plot", scale = TRUE)
#' ```
#'
#' @param data A matrix or data frame containing the data for the PCA analysis.
#' @param sample_data A data frame containing sample information such as group variable.
#' @param group_var The name of the group variable in `sample_data` to be used for coloring the points on the PCA plot.
#' @param method The method used for calculating the distance matrix. The default is "euclidean".
#' @param palette The color palette to be used for coloring the points on the PCA plot. The default is "Dark2".
#' @param alpha The transparency level of the points on the PCA plot. The default is 0.75.
#' @param lg.position The position of the legend on the PCA plot. The default is "bottom".
#' @param plot.centroids A logical value indicating whether to plot centroids for each group on the PCA plot. The default is FALSE.
#' @param plot_title The title of the PCA plot. The default is "PCA Plot".
#' @param scale A logical value indicating whether to scale the data before performing PCA. The default is TRUE.
#'
#' @return The PCA plot as an output.
#'
#' @import tidyverse
#' @import Rcolorbrewer
#' @import ggpubr
#' 
#' @examples
#' # Load data
#' data <- read.csv("data.csv")
#' sample_data <- read.csv("sample_info.csv")
#'
#' # Generate PCA plot
#' generate_pca_plot(data, sample_data = sample_data, group_var = "group", method = "euclidean", palette = "Dark2", alpha = 0.75, lg.position = "bottom", plot.centroids = FALSE, plot_title = "PCA Plot", scale = TRUE)
#'



# Function to generate PCA plot
 generate_pca_plot <- function(data, sample_data = NULL, group_var = NULL, method = "euclidean", 
                               palette = "Dark2", alpha = 0.75, lg.position = "bottom", plot.centroids = FALSE,
                               plot_title = "PCA Plot", scale = TRUE) {
                              
   suppressWarnings({ 
     library(tidyverse, quietly = TRUE,
             warn.conflicts = FALSE)
     theme_ridges <- source("https://raw.githubusercontent.com/wilkelab/ggridges/master/R/theme.R")
   })
  # Perform PCA
  dist_df <- dist(t(data), method = method)
  pca <- prcomp(dist_df, scale. = scale)
  pca_df <- as.data.frame(pca$x)
  
  # Get variance explained by each PC
  per_var_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100

  # Add group variable if provided
  if (!is.null(group_var)) {
    pca_df <- pca_df %>%
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
                                      size = 12.5) +
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