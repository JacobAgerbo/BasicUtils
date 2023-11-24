#!/usr/bin/env Rscript
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