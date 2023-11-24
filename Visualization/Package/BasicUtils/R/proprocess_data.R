#' Preprocess Data
#'
#' The `preprocess_data` function is used to preprocess the input data before further analysis. It performs several preprocessing steps, including data cleaning, normalization, and batch effect correction using ComBat.
#'
#' The function performs the following steps:
#' 1. Loads the necessary packages, including `sva`, `tidyverse`, and `cowplot`.
#' 2. Stores a copy of the original data in `prior_data`.
#' 3. If `clean_data` is set to `TRUE`, removes missing or invalid values from the data and removes duplicate samples or features.
#' 4. If `sum_scaling` is set to `TRUE`, performs total sum scaling normalization on the data.
#' 5. If `magic_norm` is set to `TRUE`, applies inverse normal transformation (INT) to normalize continuous phenotype data.
#' 6. If `correct_batch` is set to `TRUE`, identifies batch information from `sample_data` and performs batch effect correction using ComBat.
#' 7. Generates outlier detection plots for raw and processed data if applicable.
#' 8. Generates PCA plots for raw and batch-corrected data if applicable.
#' 9. Combines the outlier detection plots and PCA plots into a single plot.
#' 10. Returns the combined plot, processed data, and raw data as output.
#'
#' Example usage:
#' ```R
#' # Preprocess data
#' preprocessed <- preprocess_data(data, sample_data = sample_data, batch = "Group", clean_data = TRUE, sum_scaling = TRUE, magic_norm = FALSE, correct_batch = TRUE)
#'
#' # Access the combined plot
#' preprocessed$plot
#'
#' # Access the processed data
#' preprocessed$processed_data
#'
#' # Access the raw data
#' preprocessed$raw_data
#' ```
#'
#' @import sva
#' @import tidyverse
#' @import cowplot
#'
#' @param data A matrix or data frame containing the input data to be preprocessed.
#' @param sample_data A data frame containing sample information, such as batch information.
#' @param batch The name of the column in `sample_data` that represents the batch information.
#' @param clean_data A logical value indicating whether to perform data cleaning (removing missing or invalid values and duplicate samples or features). Default is `TRUE`.
#' @param sum_scaling A logical value indicating whether to perform total sum scaling normalization. Default is `TRUE`.
#' @param magic_norm A logical value indicating whether to perform inverse normal transformation (INT) for normalization of continuous phenotype data. Default is `FALSE`.
#' @param correct_batch A logical value indicating whether to perform batch effect correction using ComBat. Default is `TRUE`.
#'
#' @return A list containing the combined plot, processed data, and raw data.
#' 
preprocess_data <- function(data, 
                            sample_data = NULL, 
                            batch = "Group", 
                            clean_data = TRUE, 
                            sum_scaling = TRUE, 
                            magic_norm = FALSE, 
                            correct_batch = TRUE) {
  
  suppressPackageStartupMessages(library(sva))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(cowplot))
  theme_ridges <- source("https://raw.githubusercontent.com/wilkelab/ggridges/master/R/theme.R")
  prior_data <- data
  # Data cleaning
  if (clean_data) {
    # Remove missing or invalid values
    data <- na.omit(data)
    # Remove duplicate samples or features
    data <- unique(data)
    
    Outliers_raw_data <- prior_data %>%
      reshape2::melt() %>%
      ggplot(aes(x = "", y = value)) +
      geom_boxplot(fill = "orange", color = "black") +
      labs(title = "Outlier detection of raw data", x = "", y = "Data") +
      theme_ridges()
    Outliers_processed_data <- data %>%
      reshape2::melt() %>%
      ggplot(aes(x = "", y = value)) +
      geom_boxplot(fill = "lightblue", color = "black") +
      labs(title = "Outlier detection of processed data", x = "", y = "Data") +
      theme_ridges()
    
    boxes <- cowplot::plot_grid(Outliers_raw_data,Outliers_processed_data)
  }
  
  # Normalization
  if (sum_scaling) {
    # Perform normalization steps here (e.g., total sum scaling, log transformation, etc.)
    # Replace the following line with your specific normalization method
    tss <- function (abund){
      if (is.null(dim(abund)) == TRUE) {
        abund.norm <- abund/sum(abund)
      }
      if (is.null(dim(abund)) == FALSE) {
        abund.norm <- sweep(abund, 2, colSums(abund), FUN = "/")
      }
      return(abund.norm)}
    data <- tss(data)
    
    Outliers_raw_data <- prior_data %>%
      reshape2::melt() %>%
      ggplot(aes(x = "", y = value)) +
      geom_boxplot(fill = "orange", color = "black") +
      labs(title = "Outlier detection of raw data", x = "", y = "Data") +
      theme_ridges()
    Outliers_processed_data <- data %>%
      reshape2::melt() %>%
      ggplot(aes(x = "", y = value)) +
      geom_boxplot(fill = "lightblue", color = "black") +
      labs(title = "Outlier detection of processed data", x = "", y = "Data") +
      theme_ridges()
    
    boxes <- cowplot::plot_grid(Outliers_raw_data,Outliers_processed_data)
  }
  
  if (magic_norm) {
    # inverse normal transformation (INT) suggested by Anders Albrechtsen for
    # normalisation of continuous phenotype data for association tests, like GWAS
    
    # The inverse normal transformation is a statistical method used to transform a non-normal distribution to a normal distribution. 
    # It involves converting the original data to a set of quantiles, 
    # which are then transformed to their corresponding z-scores (i.e., the number of standard deviations from the mean). 
    # This transformation results in a distribution that is approximately normal, 
    # which can be useful for certain types of statistical analyses that assume normality. 
    # The inverse normal transformation is often used in genetics research to normalize continuous phenotype data 
    # for genome-wide association studies (GWAS).
    
    qtrans <- function(x) {
      k <- !is.na(x)  # Identify non-missing values
      ran <- rank(x[k])  # Rank the non-missing values
      y <- qnorm((1:sum(k) - 0.5) / sum(k))  # Apply the inverse normal transformation
      x[k] <- y[ran]  # Replace the original values with the transformed values
      x  # Return the modified vector
    }
    
    # df: samples as rows, features as columns, values are abundances
    data <- as.data.frame(apply(t(data), MARGIN = 2, FUN = qtrans))
    data <- t(data) %>%
      as_tibble()
    
    Outliers_raw_data <- prior_data %>%
      reshape2::melt() %>%
      ggplot(aes(x = "", y = value)) +
      geom_boxplot(fill = "orange", color = "black") +
      labs(title = "Outlier detection of raw data", x = "", y = "Data") +
      theme_ridges()
    Outliers_processed_data <- data %>%
      reshape2::melt() %>%
      ggplot(aes(x = "", y = value)) +
      geom_boxplot(fill = "orange", color = "black") +
      labs(title = "Outlier detection of processed data", x = "", y = "Data") +
      theme_ridges()
    
    boxes <- cowplot::plot_grid(Outliers_raw_data,Outliers_processed_data)
  }
  
  # Batch effect correction using ComBat
  if (correct_batch) {
    # Identify batch information in the data (if available)
    batch_info <- sample_data[,batch] 
    
    # Perform batch effect correction using ComBat
    data_corrected <- sva::ComBat(dat = data, batch = batch_info)
    
    # Replace the original expression data with the corrected data
    data <- data_corrected
    
    
    source("https://raw.githubusercontent.com/JacobAgerbo/Basic_Utils/main/Visualization/make_PCA.R")
    batch_raw_data <- generate_pca_plot(data = prior_data, 
                                        sample_data = sample_data, 
                                        group_var = batch, 
                                        method = "euclidean", 
                                        palette = "YlOrRd", alpha = 0.75, 
                                        lg.position = "bottom", 
                                        plot.centroids = FALSE,
                                        plot_title = "Raw", scale = TRUE)  
    batch_processed_data <- generate_pca_plot(data = data, 
                                              sample_data = sample_data, 
                                              group_var = batch, 
                                              method = "euclidean", 
                                              palette = "YlOrRd", alpha = 0.75, 
                                              lg.position = "bottom", 
                                              plot.centroids = FALSE,
                                              plot_title = "Batch-corrected", scale = TRUE)  
    batches <- cowplot::plot_grid(batch_raw_data,batch_processed_data)
  }
  
  if (exists("batches") && exists("boxes")) {
    # Your code here
    plot <- cowplot::plot_grid(boxes,batches, nrow = 2)
  } else if (exists("batches")) {
    plot <- batches
  } else if (exists("boxes")) {
    plot <- boxes
  } else {
    cat("Something went wrong, either plots for processing or batch correction was found.")
  }

     return(list(plot = plot, processed_data = data, raw_data = prior_data))  
}


#processed_data <- preprocess_data(data = abundance_data,
#                sample_data = sample_data,
#                batch = "Batch",
#                clean_data = TRUE,
#                sum_scaling = TRUE,
#                correct_batch = TRUE,
#                magic_norm = TRUE)
