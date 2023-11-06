#!/usr/bin/env Rscript
# R Script for Generating biomarkers, using random forrest or GLMnet (or both)
# Author: Jacob Agerbo Rasmussen
# Load necessary packages

# Define the function
find_biomarkers <- function(data,
                    sample_data,
                    feature,
                    exp_var, 
                    random_var, 
                    datatype = c("logcpm", "relabu", "counts"), 
                    method = c("GLM", "RF", "both"),
                    top_biomarker=0.1,
                    prevalence_tolerance=0.01){
  
  exp_var = exp_var
  random_var = random_var
  datatype = match.arg(datatype)
  method = match.arg(method)
  feature = feature
  top_biomarker=top_biomarker

  # Set dependencies
  suppressWarnings({ 
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(reshape2))
    suppressPackageStartupMessages(library(cowplot))
    suppressPackageStartupMessages(library(ggpubr))
    suppressPackageStartupMessages(library(animalcules))
    suppressPackageStartupMessages(library(plotROC))
    suppressPackageStartupMessages(library(MASS))
    suppressPackageStartupMessages(library(caret))
    suppressPackageStartupMessages(library(wesanderson))
    suppressPackageStartupMessages(library(hilldiv))
  })
  
  if (datatype == "relabu"){
    lab = "Relative Abundance"
  } else if (datatype != "relabu"){
    if (datatype == "logcpm"){
      lab = "Log Counts Per Million"
    } else {
      lab = "Counts"
    }
  }
  
  ## ggplot theme
  theme_ridges <- source("https://raw.githubusercontent.com/wilkelab/ggridges/master/R/theme.R")
  
  #
  df <- data %>%
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
  
  rownames(df) <- colnames(data)  
  
  df_prevalence <- df  
  
  # add target variable
  df[, "y"] <- 
    sample_data %>% dplyr::pull(exp_var)
  # set up classification model prameters
  fitControl <- caret::trainControl(method = "repeatedcv", 
                                    number = 3, 
                                    repeats = 99, 
                                    classProbs = TRUE, 
                                    summaryFunction = twoClassSummary, 
                                    sampling = "smote", 
                                    savePredictions = TRUE)
  #  
  
  if (method == "GLM"){
    method.cap = "Generalized linear models (GLMs)"
    # Fit model
    model_fit_glm <- 
      caret::train(y ~ ., data = logcpm_table, method = "glmnet", 
                   tuneLength = 5, trControl = fitControl, metric = "ROC")
    model.plot <- ggplot(model_fit_glm) + theme_bw()
    
    # Estimate importance
    percent_top_biomarker = top_biomarker
    # GLMs
    importance_df_glm <- caret::varImp(model_fit_glm)$importance %>% 
      base::as.data.frame() %>% 
      rownames_to_column() %>% 
      dplyr::rename(importance = Overall) %>% 
      dplyr::rename(biomarker = rowname) %>% 
      dplyr::arrange(importance) %>% 
      dplyr::filter(importance > 
                      quantile(importance, 1 - percent_top_biomarker)) %>% 
      dplyr::mutate(biomarker = forcats::fct_inorder(biomarker))
    importance_df_glm$biomarker <- paste("GLM", importance_df_glm$biomarker)
    importance_df_glm$Type <- rep("Generalised Linear Model", length(rownames(importance_df_glm)))
    
    
    importance_df <- importance_df_glm
    #write.csv(importance_df, file = "biomarkers_importance_PRO.csv")
    #
    importance_plot <- ggplot2::ggplot() + 
      geom_col(data = importance_df, aes(x = reorder(biomarker, importance), y = importance)) + 
      coord_flip() + theme_minimal() + xlab(paste("Top ",percent_top_biomarker,"% Biomarkers", sep = "")) + 
      ylab("Biomarker Importance")
    #
    prob_pred <- as.numeric(model_fit_glm$pred$obs)
    prob_pred[prob_pred == 1] <- 0
    prob_pred[prob_pred == 2] <- 1
    #
    df_roc_glm <- data.frame(m = model_fit_glm$pred[, 
                                                    which(colnames(model_fit_glm$pred) == levels(model_fit_glm$pred$obs)[2])], 
                             d = prob_pred, stringsAsFactors = FALSE)
    #
    df_roc_glm$Type <- rep("GLMnet", length(df_roc_glm))
    df_roc <- df_roc_glm
    #
    #
    g <- ggplot(df_roc, aes(m = m, d = d, colour = Type)) + 
      geom_roc(n.cuts = 0) + coord_equal() + 
      style_roc() +
      ggtitle(method.cap)
    #
    roc_plot <- g + annotate("text", 
                             x = 0.75, 
                             y = 0.35, 
                             label = paste("AUC - GLM =", 
                                           round((calc_auc(g))$AUC[1],4)))
    
    
    # Prevalence plot
    df_prevalence <- t(data) %>%
      as_tibble()
    Groups <- unique(sample_data[,exp_var])
    Group_A <- df_prevalence %>%
      hilldiv::tss() %>%
      mutate(Sample = sample_data$Sample, .before = 1)  %>%
      filter(Sample %in% sample_data$Sample[sample_data$Group == Groups[1]]) %>%
      reshape2::melt() %>%
      mutate(presence= ifelse(value > prevalence_tolerance, 1,0)) %>%
      group_by(variable) %>%
      summarise(
        prevalence = sum(presence)
      )
    
    Group_B <- df %>%
      hilldiv::tss() %>%
      mutate(Sample = sample_data$Sample, .before = 1)  %>%
      filter(Sample %in% sample_data$Sample[sample_data$Group == Groups[2]]) %>%
      reshape2::melt() %>%
      mutate(presence= ifelse(value > prevalence_tolerance, 1,0)) %>%
      group_by(variable) %>%
      summarise(
        prevalence = sum(presence)
      )
    
    prevalence_df <- full_join(Group_A,Group_B, by="variable")
    colnames(prevalence_df) <- c("Feature", "x", "y")
    
    prevalence_df <- prevalence_df %>%
      dplyr::mutate(
        deviation = x-y,
        deviation = sqrt(deviation*deviation) # had to make some kind of wierd trans to make all positive integers.
      )
    
    prevalence_plot  <- ggplot(data=prevalence_df, aes(x=x, y=y, label=Feature, color = deviation)) +
      xlab(paste("Prevalence (", Groups[1],")", sep = "")) + ylab(paste("Prevalence (", Groups[2],")", sep = "")) +
      geom_point(alpha = 0.65, size = 3) +
      ggrepel::geom_text_repel() +
      theme_minimal() +
      geom_abline(slope = 1, intercept = 0) +
      xlim(0,length(sample_data$Sample[sample_data$Group == Groups[1]])) +
      ylim(0,length(sample_data$Sample[sample_data$Group == Groups[2]])) +
      theme_ridges() +
      scale_color_gradientn(colors = c(wesanderson::wes_palette(21, name = "Zissou1", type = "continuous")))
    
    # Combine final plots
    plot1 <- cowplot::plot_grid(importance_plot, roc_plot)
    out_plot <- cowplot::plot_grid(plot1,
                                   prevalence_plot, 
                                   ncol = 1,
                                   labels = "AUTO")
    
    return(list(model = model.plot, table = importance_df, plot = out_plot))
    
  } else if (method == "RF"){
    method.cap = "Random Forrest (RF)"
    
    model_fit_rf <- 
      caret::train(y ~ ., data = logcpm_table, method = "ranger", 
                   trControl = fitControl, tuneLength = 5, 
                   metric = "ROC", importance = "impurity")
    
    model.plot <- ggplot(model_fit_rf) + theme_bw()
    
    ## Pick a dired % to investigate
    percent_top_biomarker = top_biomarker
    
    # RF
    importance_df_rf <- caret::varImp(model_fit_rf)$importance %>% 
      base::as.data.frame() %>% 
      rownames_to_column() %>% 
      dplyr::rename(importance = Overall) %>% 
      dplyr::rename(biomarker = rowname) %>% 
      dplyr::arrange(importance) %>% 
      dplyr::filter(importance > 
                      quantile(importance, 1 - percent_top_biomarker)) %>% 
      dplyr::mutate(biomarker = forcats::fct_inorder(biomarker))
    importance_df_rf$biomarker <- paste("RF", importance_df_rf$biomarker)
    importance_df_rf$Type <- rep("Random Forrest", length(rownames(importance_df_rf)))
    #
    importance_df <- rbind(importance_df_glm,importance_df_rf)
    #write.csv(importance_df, file = "biomarkers_importance_PRO.csv")
    #
    importance_plot <- ggplot2::ggplot() + 
      geom_col(data = importance_df, aes(x = reorder(biomarker, importance), y = importance)) + 
      coord_flip() + theme_minimal() + xlab(paste("Top ",percent_top_biomarker,"% Biomarkers", sep = "")) + 
      ylab("Biomarker Importance")
    
    prob_pred <- as.numeric(model_fit_rf$pred$obs)
    prob_pred[prob_pred == 1] <- 0
    prob_pred[prob_pred == 2] <- 1
    #
    df_roc_rf <- data.frame(m = model_fit_rf$pred[, 
                                                  which(colnames(model_fit_rf$pred) == levels(model_fit_rf$pred$obs)[2])], 
                            d = prob_pred, stringsAsFactors = FALSE)
    #
    df_roc_rf$Type <- rep("Random Forrest", length(df_roc_rf))
    #
    #
    df_roc <- df_roc_rf
    #
    g <- ggplot(df_roc, aes(m = m, d = d, colour = Type)) + 
      geom_roc(n.cuts = 0) + coord_equal() + 
      style_roc() +
      ggtitle(method.cap)
    #
    roc_plot <- g + annotate("text", x = 0.75, y = 0.25,
                             label = paste("AUC - RF =", round((calc_auc(g))$AUC[2], 
                                                               4)))
    
    # Prevalence plot
    df_prevalence <- t(data) %>%
      as_tibble()
    Groups <- unique(sample_data[,exp_var])
    Group_A <- df_prevalence %>%
      hilldiv::tss() %>%
      mutate(Sample = sample_data$Sample, .before = 1)  %>%
      filter(Sample %in% sample_data$Sample[sample_data$Group == Groups[1]]) %>%
      reshape2::melt() %>%
      mutate(presence= ifelse(value > prevalence_tolerance, 1,0)) %>%
      group_by(variable) %>%
      summarise(
        prevalence = sum(presence)
      )
    
    Group_B <- df %>%
      hilldiv::tss() %>%
      mutate(Sample = sample_data$Sample, .before = 1)  %>%
      filter(Sample %in% sample_data$Sample[sample_data$Group == Groups[2]]) %>%
      reshape2::melt() %>%
      mutate(presence= ifelse(value > prevalence_tolerance, 1,0)) %>%
      group_by(variable) %>%
      summarise(
        prevalence = sum(presence)
      )
    
    prevalence_df <- full_join(Group_A,Group_B, by="variable")
    colnames(prevalence_df) <- c("Feature", "x", "y")
    
    
    prevalence_df <- prevalence_df %>%
      dplyr::mutate(
        deviation = x-y,
        deviation = sqrt(deviation*deviation) # had to make some kind of wierd trans to make all positive integers.
      )
    
    prevalence_plot <- ggplot(data=prevalence_df, aes(x=x, y=y, label=Feature, color = deviation)) +
      xlab(paste("Prevalence (", Groups[1],")", sep = "")) + ylab(paste("Prevalence (", Groups[2],")", sep = "")) +
      geom_point(alpha = 0.65, size = 3) +
      ggrepel::geom_text_repel() +
      theme_minimal() +
      geom_abline(slope = 1, intercept = 0) +
      xlim(0,length(sample_data$Sample[sample_data$Group == Groups[1]])) +
      ylim(0,length(sample_data$Sample[sample_data$Group == Groups[2]])) +
      theme_ridges() +
      scale_color_gradientn(colors = c(wesanderson::wes_palette(21, name = "Zissou1", type = "continuous")))
    
    # Combine final plots
    plot1 <- cowplot::plot_grid(importance_plot, roc_plot)
    out_plot <- cowplot::plot_grid(plot1,
                                   prevalence_plot, 
                                   ncol = 1,
                                   labels = "AUTO")
    
    return(list(model = model.plot, table = importance_df, plot = out_plot))
    
    
  } else {
    method.cap = "Generalized linear models (GLMs) & Random Forrest (RF)"
    
    # Both
    model_fit_glm <- 
      caret::train(y ~ ., data = logcpm_table, method = "glmnet", 
                   tuneLength = 5, trControl = fitControl, metric = "ROC")
    
    model_fit_rf <- 
      caret::train(y ~ ., data = logcpm_table, method = "ranger", 
                   trControl = fitControl, tuneLength = 5, 
                   metric = "ROC", importance = "impurity")
    
    #
    GLM.plot <- ggplot(model_fit_glm) + theme_bw()
    RF.plot <- ggplot(model_fit_rf) + theme_bw()
    
    # Estimate importance
    ## Pick a dired % to investigate
    percent_top_biomarker = top_biomarker
    
    # GLMs
    importance_df_glm <- caret::varImp(model_fit_glm)$importance %>% 
      base::as.data.frame() %>% 
      rownames_to_column() %>% 
      dplyr::rename(importance = Overall) %>% 
      dplyr::rename(biomarker = rowname) %>% 
      dplyr::arrange(importance) %>% 
      dplyr::filter(importance > 
                      quantile(importance, 1 - percent_top_biomarker)) %>% 
      dplyr::mutate(biomarker = forcats::fct_inorder(biomarker))
    importance_df_glm$biomarker <- paste("GLM", importance_df_glm$biomarker)
    importance_df_glm$Type <- rep("Generalised Linear Model", length(rownames(importance_df_glm)))
    
    # RF
    importance_df_rf <- caret::varImp(model_fit_rf)$importance %>% 
      base::as.data.frame() %>% 
      rownames_to_column() %>% 
      dplyr::rename(importance = Overall) %>% 
      dplyr::rename(biomarker = rowname) %>% 
      dplyr::arrange(importance) %>% 
      dplyr::filter(importance > 
                      quantile(importance, 1 - percent_top_biomarker)) %>% 
      dplyr::mutate(biomarker = forcats::fct_inorder(biomarker))
    importance_df_rf$biomarker <- paste("RF", importance_df_rf$biomarker)
    importance_df_rf$Type <- rep("Random Forrest", length(rownames(importance_df_rf)))
    #
    importance_df <- rbind(importance_df_rf, importance_df_glm)
    #
    #
    importance_plot <- ggplot2::ggplot() + 
      geom_col(data = importance_df, aes(x = reorder(biomarker, importance), y = importance, colour = Type, fill = Type)) + 
      coord_flip() + theme_minimal() + xlab("Top 5% Biomarkers") + ylab("Biomarker Importance") + 
      scale_fill_manual(values = wes_palette("Darjeeling1", n =2)) + 
      scale_color_manual(values = wes_palette("Darjeeling1", n =2))
    
    prob_pred <- as.numeric(model_fit_rf$pred$obs)
    prob_pred[prob_pred == 1] <- 0
    prob_pred[prob_pred == 2] <- 1
    #
    df_roc_rf <- data.frame(m = model_fit_rf$pred[, 
                                                  which(colnames(model_fit_rf$pred) == levels(model_fit_rf$pred$obs)[2])], 
                            d = prob_pred, stringsAsFactors = FALSE)
    #
    df_roc_rf$Type <- rep("Random Forrest", length(df_roc_rf))
    #
    prob_pred <- as.numeric(model_fit_glm$pred$obs)
    prob_pred[prob_pred == 1] <- 0
    prob_pred[prob_pred == 2] <- 1
    #
    df_roc_glm <- data.frame(m = model_fit_glm$pred[, 
                                                    which(colnames(model_fit_glm$pred) == levels(model_fit_glm$pred$obs)[2])], 
                             d = prob_pred, stringsAsFactors = FALSE)
    #
    df_roc_glm$Type <- rep("GLMnet", length(df_roc_glm))
    #
    df_roc <- rbind(df_roc_glm,df_roc_rf)
    #
    g <- ggplot(df_roc, aes(m = m, d = d, colour = Type)) + 
      geom_roc(n.cuts = 0) + coord_equal() + 
      style_roc() + scale_fill_manual(values = wes_palette("Darjeeling1", n =2)) + scale_color_manual(values = wes_palette("Darjeeling1", n =2))
    #
    roc_plot <- g + annotate("text", x = 0.75, y = 0.35, 
                             label = paste("AUC - GLM =", round((calc_auc(g))$AUC[1], 
                                                                4))) + annotate("text", x = 0.75, y = 0.25, 
                                                                                label = paste("AUC - RF =", round((calc_auc(g))$AUC[2], 
                                                                                                                  4)))
    
    # Prevalence plot
    df_prevalence <- t(data) %>%
      as_tibble()
    Groups <- unique(sample_data[,exp_var])
    Group_A <- df_prevalence %>%
      hilldiv::tss() %>%
      mutate(Sample = sample_data$Sample, .before = 1)  %>%
      filter(Sample %in% sample_data$Sample[sample_data$Group == Groups[1]]) %>%
      reshape2::melt() %>%
      mutate(presence= ifelse(value > prevalence_tolerance, 1,0)) %>%
      group_by(variable) %>%
      summarise(
        prevalence = sum(presence)
      )
    
    Group_B <- df %>%
      hilldiv::tss() %>%
      mutate(Sample = sample_data$Sample, .before = 1)  %>%
      filter(Sample %in% sample_data$Sample[sample_data$Group == Groups[2]]) %>%
      reshape2::melt() %>%
      mutate(presence= ifelse(value > prevalence_tolerance, 1,0)) %>%
      group_by(variable) %>%
      summarise(
        prevalence = sum(presence)
      )
    
    prevalence_df <- full_join(Group_A,Group_B, by="variable")
    colnames(prevalence_df) <- c("Feature", "x", "y")
    
    
    prevalence_df <- prevalence_df %>%
      dplyr::mutate(
        deviation = x-y,
        deviation = sqrt(deviation*deviation) # had to make some kind of wierd trans to make all positive integers.
      )
    
    prevalence_plot <- ggplot(data=prevalence_df, aes(x=x, y=y, label=Feature, color = deviation)) +
      xlab(paste("Prevalence (", Groups[1],")", sep = "")) + ylab(paste("Prevalence (", Groups[2],")", sep = "")) +
      geom_point(alpha = 0.65, size = 3) +
      ggrepel::geom_text_repel() +
      theme_minimal() +
      geom_abline(slope = 1, intercept = 0) +
      xlim(0,length(sample_data$Sample[sample_data$Group == Groups[1]])) +
      ylim(0,length(sample_data$Sample[sample_data$Group == Groups[2]])) +
      theme_ridges() +
      scale_color_gradientn(colors = c(wesanderson::wes_palette(21, name = "Zissou1", type = "continuous")))
    
    # Combine final plots
    plot1 <- cowplot::plot_grid(importance_plot, roc_plot)
    out_plot <- cowplot::plot_grid(plot1,
                                   prevalence_plot, 
                                   ncol = 1,
                                   labels = "AUTO")
    
  }
}


# Make test data
set.seed(1234)
data1 <- matrix(rpois(10000, 50), nrow = 100)
data2 <- matrix(rpois(5000, 4), nrow = 50)
data3 <- matrix(rpois(50, 1000), nrow = 1)
data4 <- matrix(rpois(50, 10), nrow = 1)
data3 <- cbind(data3,data4)
data <- rbind(data3, data1, data2)
rm(data1, data2, data3, data4)

colnames(data) <- c(paste0("Sample", 1:ncol(data)))
sample_data <- data.frame("Sample"= colnames(data),
                          "Group" = rep(c("Group A", "Group B"), each = 50))

df <- t(data) %>%
  as_tibble()

  