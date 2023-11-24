#' BasicUtils - Perform Linear Mixed-Effects Modeling and Correlation Analysis
#' 
#' The `make_lme_cor` function is used to perform linear mixed-effects modeling and correlation analysis. 
#' It takes the following parameters:
#'
#' - `data`: A matrix or data frame containing the data for the analysis.
#' - `sample_data`: A data frame containing sample information such as treatment and group.
#' - `feature`: The name of the feature (dependent variable) to be analyzed.
#' - `exp_var`: The name of the explanatory variable.
#' - `random_var`: The name of the random effect variable.
#' - `datatype`: The type of data to be analyzed. It can be "logcpm", "relabu", or "counts".
#' - `tolerance`: The tolerance value used for computing adjusted correlation.
#'
#' The function performs the following steps:
#'
#' 1. Loads required packages (`tidyverse`, `performance`, `ggpubr`, `animalcules`, `lmerTest`, `broom.mixed`, and `ggplot2`).
#' 2. Prepares the data by transforming it based on the specified data type.
#' 3. Removes rows with zero counts, if any.
#' 4. Selects the relevant columns from the data frame.
#' 5. Creates a scatter plot with a linear regression line for each group.
#' 6. Estimates a mixed-effects model using the `lmerTest` package.
#' 7. Retrieves the sign of the coefficient and the adjusted correlation.
#' 8. Computes the p-value using `anova` with Kenward-Roger degrees of freedom adjustment.
#' 9. Returns the plot and statistics as a list.
#'
#' Example usage:
#' ```R
#' # Load data
#' data <- read.csv("data.csv")
#' sample_data <- read.csv("sample_info.csv")
#'
#' # Perform linear mixed-effects modeling and correlation analysis
#' make_lme_cor(data, sample_data, "feature", "exp_var", "random_var", datatype = "logcpm", tolerance = 0.1)
#' ```
#'
#' @param data A matrix or data frame containing the data for the analysis
#' @param sample_data A data frame containing sample information such as treatment and group
#' @param feature The name of the feature (dependent variable) to be analyzed
#' @param exp_var The name of the explanatory variable
#' @param random_var The name of the random effect variable
#' @param datatype The type of data to be analyzed (default: "logcpm")
#' @param tolerance The tolerance value used for computing adjusted correlation
#'
#' @return A list containing the scatter plot and statistics
#' 
#' @import tidyverse
#' @import performance
#' @import ggpubr
#' @import animalcules
#' @import lmerTest
#' @import broom.mixed
#' @import ggplot2
#'

make_lme_cor <- function(data, 
                    sample_data,
                    feature, 
                    exp_var, 
                    random_var, 
                    datatype = c("logcpm", "relabu", "counts"), 
                    tolerance){
  
  
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(performance))
  suppressPackageStartupMessages(library(ggpubr))
  suppressPackageStartupMessages(library(animalcules))
  suppressPackageStartupMessages(library(lmerTest))
  suppressPackageStartupMessages(library(broom.mixed))
  suppressPackageStartupMessages(library(ggplot2))
  
  sam_table = sample_data
  exp_var = exp_var
  random_var = random_var
  datatype = match.arg(datatype)
  feature = feature
  tolerance = tolerance
  
  
  suppressWarnings({
  if (datatype == "relabu"){
    lab = "Relative Abundance"
  } else if (datatype != "relabu"){
    if (datatype == "logcpm"){
      lab = "Log Counts Per Million"
    } else {
      lab = "Counts"
    }
  }
  })
  ## ggplot theme
  theme_shannon <- function(){
    # theme minimal creates white background and removes ticks on axes ----
    theme_minimal() +
      theme(
        # removes grid lines from plot ----
        panel.grid = element_blank(),
        # moves legend to top instead of side ----
        legend.position = "top",
        # removes title from legend, often superfluous ----
        legend.title = element_blank(),
        # creates the light gray box around the plot ----
        panel.background = element_rect(color = "#F2F2F2"),
        # creates the gray background boxes for faceting labels ----
        strip.background = element_rect(
          color = "#F2F2F2",
          fill = "#F2F2F2"
        ),
        # if using facet grid, this rotates the y text to be more readable ----
        strip.text.y = element_text(angle = 0),
        # this produces a fully left justified title ----
        plot.title.position = "plot"
      )
  }
  
  microbe <- data
  sam_table <- sam_table
  counts_table <- data
  
  suppressWarnings({
  df <- counts_table %>%   
    {
      if (datatype == "relabu") {
        counts_to_relabu(.)
      }
      else if (datatype == "logcpm") {
        counts_to_logcpm(.)
      }
      else {
        .
      }
    } %>% {
      if (sum(base::rowSums(as.matrix(.)) == 0) > 0) {
        . <- .[-which(base::rowSums(as.matrix(.)) == 0), 
        ]
      }
      else {
        .
      }
    } %>%  t() %>%  
    as_tibble()
  })
  
  suppressWarnings({
    suppressMessages({
      rownames(df) <- colnames(counts_table)
    })
  })
  #
  
  df <- df %>% 
    dplyr::select(feature) %>% 
    # duplicate species variable for coloring & grouping ---
    dplyr::mutate(exp_var = as.numeric(sam_table[,exp_var])) %>% 
    dplyr::mutate(categories = sam_table[,random_var]) %>%
    dplyr::mutate(random_var = sam_table[,random_var]) %>%
    drop_na()
  
  colnames(df) <- c("feature", colnames(df[2:4]))
  level_df <- c(levels(as.factor(df$categories)))
  
  suppressWarnings({
    suppressMessages({
      plot <- df %>% 
    # add stack for all species to be analyzed together ----
  bind_rows(df %>% mutate(categories = "All")) %>% 
    # now examine by 3 species plus all ----
  group_by(categories) %>% 
    nest() %>%
    # within each group, compute base n and correlation ----
  mutate(
    base_n = map_int(data, nrow),
    corr = map(data, ~ cor.test(x = .x$exp_var, y = .x$feature) %>% broom::tidy())
  ) %>% 
    ungroup() %>% 
    # bring results back to raw data ----
  unnest(c(data, corr)) %>% 
    mutate(
      # create ordered facet label for plotting ----
      categories = fct_relevel(categories, level_df, "All"),
      corr_label =  glue::glue("{categories}\nn = {base_n}\n r = {scales::number(estimate, 0.01)}"),
      corr_label = fct_reorder(as.character(corr_label), as.numeric(categories))
    ) %>% 
    # create scatter plots ----
  ggplot(aes(x = exp_var, y = feature)) +
    geom_point(aes(color = random_var), alpha = 0.85, show.legend = FALSE) +
    geom_smooth(method = "lm", color = "darkgray", se = FALSE) +
    facet_wrap(. ~ corr_label, nrow = 1) +
    scale_color_brewer(palette = "Set2") + 
    theme_shannon() + 
    xlab(paste(exp_var)) + 
    ylab(paste(lab," of ",feature))
    })
  })
  
  # estimate mixed model ----
  mixed_model <- lmerTest::lmer(feature ~ exp_var + (1 | random_var), df)
  
  # retrieve sign of coefficient ----
  coef_sign <- mixed_model %>% 
    broom.mixed::tidy() %>% 
    filter(term == "exp_var") %>% 
    pull(estimate) %>% 
    sign()
  
  # retrieve r2 measure ----
  r2_by_group <- performance::r2_nakagawa(mixed_model, by_group = TRUE, tolerance = tolerance)$R2[1]
  # compute adjusted correlation ----
  adj_corr <- coef_sign * sqrt(r2_by_group)
  aov <- anova(mixed_model, ddf="Kenward-Roger")
  aov$`Pr(>F)`
  stats <- c(adj_corr,aov$`Pr(>F)`)
  names(stats) <- c("Random Effect Adjusted Correlation", "p-value")
  suppressWarnings({
    suppressMessages({
      return(list(plot = plot, stats = stats))  
    })
  })
  
}

# Make test data
#set.seed(1234)
#data1 <- matrix(rpois(10000, 50), nrow = 100)
#data2 <- matrix(rpois(5000, 4), nrow = 50)
#data31 <- matrix(rpois(50, 1000), nrow = 1)
#data32 <- matrix(rpois(50, 1000), nrow = 1)
#data33 <- matrix(rpois(50, 200), nrow = 1)
#data41 <- matrix(rpois(50, 100), nrow = 1)
#data42 <- matrix(rpois(50, 50), nrow = 1)
#data43 <- matrix(rpois(50, 1000), nrow = 1)
#data3 <- cbind(data31,data41)
#data4 <- cbind(data32,data42)
#data5 <- cbind(data33,data43)
#data <- rbind(data3,data4,data5, data1, data2)
#rm(data1, data2, data3, data4,data5)
#rm(data31,data32,data33,data41,data42,data43)
#colnames(data) <- c(paste0("Sample", 1:ncol(data)))
#rownames(data) <- c(paste0("F", 1:nrow(data)))
# Set the length of the vector
#length <- 100
# Generate random values between 0 and 1
#random_values <- runif(length)
# Create the increasing vector
#increasing_vector <- cumsum(random_values)

#sample_data <- data.frame("Sample"= colnames(data),
#                          "Treatment" = increasing_vector,
#                          "Group" = rep(c("GroupA", "GroupB", "GroupC", "GroupD", "GroupE"), each = 20))


### Test Function
#lme_cor(data=data, 
#                    sam_table=sample_data,
#                    feature="F1", 
#                    exp_var = "Treatment", 
#                    random_var = "Group", 
#                    datatype = "counts", 
#                    tolerance = 0.01)


