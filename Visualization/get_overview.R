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
    
    if (!is.null(group_var)) {
      data$group <- group_var
    }
    
    suppressWarnings({
      data.melt <- melt(data)  
    })
    
    
}



# Ridgeline plot
data %>% 
  t() %>%
  as_tibble() %>%
  dplyr::mutate(group = sample_data[,group], .before = 1)  %>%
  reshape2::melt() %>%
  ggplot(aes(x = value, y = group, fill = group)) +
  ggridges::geom_density_ridges() +
  labs(title = "Ridgeline Plot") +
  scale_fill_brewer(palette = palette) +
  ggridges::theme_ridges()

data %>% 
  t() %>%
  as_tibble() %>%
  dplyr::mutate(group = sample_data[,group], .before = 1)  %>%
  reshape2::melt() %>%
  ggplot(aes(x = value)) +
  stat_ecdf() +
  facet_wrap(~group, scales = "free", nrow = 1) +
  labs(title = "Empirical Cumulative Distribution Function") +
  ggridges::theme_ridges()

# Q-Q plot
data %>% 
  t() %>%
  as_tibble() %>%
  dplyr::mutate(group = sample_data[,group], .before = 1)  %>%
  reshape2::melt() %>%
  ggplot(aes(sample = value)) +
  stat_qq() +
  facet_wrap(~variable, scales = "free") +
  labs(title = "Q-Q Plot of Variables") +
  ggridges::theme_ridges()



# Perform Shapiro-Wilk test
shapiro_test <- data.frame()
shapiro_test <- data %>% 
  t() %>%
  as_tibble() %>%
  dplyr::select(F1) %>%
  mutate(feature = as.numeric(F1)) %>%
  .$feature %>%
  shapiro.test()


lshap <- lapply(as.data.frame(t(data)), shapiro.test)
lres <- sapply(lshap, `[`, c("statistic","p.value"))

unlist(lapply(lshap, function(x) x$p.value))

# Check if the data is normally distributed
if (shapiro_test$p.value >= 0.05) {
  message("The '", variable, "' variable is normally distributed.")
  normal_data <- data
  non_normal_data <- NULL
} else {
  message("The '", variable, "' variable is not normally distributed.")
  normal_data <- data %>% 
    t() %>%
    as_tibble() %>%
    filter(!.data[[variable]] %in% shapiro_test$estimate)
  non_normal_data <- data %>% 
    t() %>%
    as_tibble() %>%
    filter(.data[[variable]] %in% shapiro_test$estimate)
}

# Plot the distributions
ggplot(data, aes(x = !!sym(variable))) +
  geom_density(fill = "lightblue", alpha = 0.7) +
  geom_density(data = normal_data, fill = "darkblue", alpha = 0.7) +
  geom_density(data = non_normal_data, fill = "red", alpha = 0.7) +
  labs(title = "Distribution of Variable: ", subtitle = variable)
