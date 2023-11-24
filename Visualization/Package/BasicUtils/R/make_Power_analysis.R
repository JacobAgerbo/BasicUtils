#' Perform Power Analysis for Sample Size and Effect Size
#' Author: Jacob Agerbo Rasmussen
#' The `do_Power` function is used to perform a power analysis for sample size and effect size. It calculates the required sample size for a given effect size, or the required effect size for a given sample size, based on the desired significance level and power.
#'
#' The function performs the following steps:
#' 1. Loads the necessary packages, including `pwr`, `tidyverse`, `scales`, and `cowplot`.
#' 2. Checks if the maximum sample size is less than 20, and displays a warning message if it is.
#' 3. Defines the effect size categories (low, medium, and large).
#' 4. Calculates the required sample sizes for different effect sizes using the `pwr.t.test` function.
#' 5. Creates a scatter plot of sample size versus effect size, with a logarithmic x-axis.
#' 6. Adds horizontal lines and text annotations to indicate the effect size categories.
#' 7. Calculates the required effect sizes for different sample sizes using the `pwr.t.test` function.
#' 8. Creates a second scatter plot of effect size versus sample size, with a logarithmic x-axis.
#' 9. Adds horizontal lines and text annotations to indicate the effect size categories.
#' 10. Combines both scatter plots into a single plot using `cowplot::plot_grid`.
#' 11. Returns the combined plot as the output.
#'
#' Example usage:
#' ```R
#' # Perform power analysis
#' do_Power(max_sample_size = 100, min_effect_size = 0.01, sig.level = 0.05, power = 0.8)
#' ```
#'
#' @param max_sample_size The maximum sample size to consider in the power analysis.
#' @param min_effect_size The minimum effect size to consider in the power analysis.
#' @param sig.level The desired significance level (Type I error rate) for the power analysis.
#' @param power The desired power (1 - Type II error rate) for the power analysis.
#'
#' @import pwr
#' @import tidyverse
#' @import scales
#' @import cowplot
#'
#' @return A combined scatter plot of sample size versus effect size and effect size versus sample size.
do_Power <- function(max_sample_size = 100, min_effect_size = 0.01, sig.level = 0.05, power = 0.8){
  
  # Load the necessary packages
  suppressPackageStartupMessages(library(pwr))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(scales))
  suppressPackageStartupMessages(library(cowplot))
  if (max_sample_size < 20){
    cat("Sample size should be bigger than 20, since i don't believe in effects size above 1.")
  } else {
  
    suppressWarnings({
      suppressMessages({
        # Define the effect size categories
        low_eff_size <- 0.2
        med_eff_size <- 0.5
        large_eff_size <- 0.8
        
        sample_size = max_sample_size
        
        effect_sizes <- c(seq(min_effect_size,1, by = 0.01))
        results_effect_sizes <- list()
        
        for (effect_size in effect_sizes){
          result <- pwr.t.test(d = effect_size, n = NULL, sig.level = sig.level, power = power)    
          results_effect_sizes[[paste("effect_size_", effect_size, sep = "")]] <- result
        }
        
        sublist_names <- names(results_effect_sizes)
        pwr_results <- data.frame("sample_size" = unlist(lapply(sublist_names, function(name) results_effect_sizes[[name]][["n"]])),
                                  "effect_size" = unlist(lapply(sublist_names, function(name) results_effect_sizes[[name]][["d"]])))
        
        annotate_x <- mean(pwr_results$sample_size)
        
        sample_size_plot <- pwr_results %>%
          ggplot(aes(x = sample_size, y = effect_size)) +
          geom_point() +
          scale_x_log10(labels = scales::label_number(scale = 10^0)) +
          theme_minimal() +
          xlab("Sample Size") +
          ylab("Effect Size") +
          geom_hline(yintercept = low_eff_size, linetype = "dashed", color = "grey50") +
          annotate("text", x = annotate_x, y = low_eff_size, vjust = -0.15, label = "Low Effect Size", color = "black") +
          geom_hline(yintercept = med_eff_size, linetype = "dashed", color = "grey50") +
          annotate("text", x = annotate_x, y = med_eff_size, vjust = -0.15, label = "Medium Effect Size", color = "black") +
          geom_hline(yintercept = large_eff_size, linetype = "dashed", color = "grey50") +
          annotate("text", x = annotate_x, y = large_eff_size, vjust = -0.15, label = "Large Effect Size", color = "black") +
          ggtitle(paste("Range of samples size to perform study with mimimum effect size (",min_effect_size, ")" ))
        
        sample_sizes <- c(seq(20,max_sample_size, by = 1))
        results_samples_sizes <- list()
        
        for (sample_size in sample_sizes){
          result <- pwr.t.test(d = NULL, n = sample_size, sig.level = sig.level, power = power)    
          results_samples_sizes[[paste("sample_size_", sample_size, sep = "")]] <- result
        }
        
        sublist_names <- names(results_samples_sizes)
        pwr_results <- data.frame("sample_size" = unlist(lapply(sublist_names, function(name) results_samples_sizes[[name]][["n"]])),
                                  "effect_size" = unlist(lapply(sublist_names, function(name) results_samples_sizes[[name]][["d"]])))
        
        
        
        
        annotate_x2 <- median(pwr_results$sample_size)
        
        # Create the ggplot with logarithmic x-axis and effect size categories
        effect_results <- pwr_results %>%
          ggplot(aes(x = sample_size, y = effect_size)) +
          geom_point() +
          scale_x_log10(labels = scales::label_number(scale = 10^0)) +
          theme_minimal() +
          xlab("Sample Size") +
          ylab("Effect Size") +
          geom_hline(yintercept = low_eff_size, linetype = "dashed", color = "grey50") +
          annotate("text", x = annotate_x2, y = low_eff_size, vjust = -0.15, label = "Requires Low Effect Size", color = "black") +
          geom_hline(yintercept = med_eff_size, linetype = "dashed", color = "grey50") +
          annotate("text", x = annotate_x2, y = med_eff_size, vjust = -0.15, label = "Requires Medium Effect Size", color = "black") +
          geom_hline(yintercept = large_eff_size, linetype = "dashed", color = "grey50") +
          annotate("text", x = annotate_x2, y = large_eff_size, vjust = -0.15, label = "Requires Large Effect Size", color = "black") +
          ggtitle(paste("Range of required effect size to perform study with maximum: ",max_sample_size, " samples" )) + scale_x_reverse()
        
        
        plot <- cowplot::plot_grid(sample_size_plot, effect_results)    
        
        return(plot)
        
      })
    })
  }
} 