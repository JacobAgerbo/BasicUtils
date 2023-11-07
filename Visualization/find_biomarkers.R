#!/usr/bin/env Rscript
# R Script for Generating biomarkers, using random forrest or GLMnet (or both)
# Author: Jacob Agerbo Rasmussen
# Load necessary packages

# Define the function
find_biomarkers <- function(data,
                            sample_data,
                            exp_var,
                            datatype = c("logcpm", "relabu", "counts"), 
                            method = c("GLM", "RF", "both"),
                            top_biomarker=0.1,
                            prevalence_tolerance=0.01){
  
  exp_var = exp_var
  datatype = match.arg(datatype)
  method = match.arg(method)
  top_biomarker=top_biomarker
  
  suppressWarnings({
  print("Starting finding biomarkers")
  start_time <- Sys.time()  # Record the start time
    
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
    ## ggplot theme
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
  
  if (datatype == "relabu"){
    lab = "Relative Abundance"
  } else if (datatype != "relabu"){
    if (datatype == "logcpm"){
      lab = "Log Counts Per Million"
    } else {
      lab = "Counts"
    }
  }
  
  
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
      caret::train(y ~ ., data = df, method = "glmnet", 
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
    
    
    # Print the estimated remaining time
    print("Still working...")
    
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
    
    Group_B <- df_prevalence %>%
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
    prevalence_df$Feature <- ifelse(log(prevalence_df$deviation) > 2.5, prevalence_df$Feature, "")
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
      caret::train(y ~ ., data = df, method = "ranger", 
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
    importance_df <- importance_df_rf
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
    
    
    # Print the estimated remaining time
    print("Still working...")

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
    
    Group_B <- df_prevalence %>%
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
    prevalence_df$Feature <- ifelse(log(prevalence_df$deviation) > 2.5, prevalence_df$Feature, "")
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
      caret::train(y ~ ., data = df, method = "glmnet", 
                   tuneLength = 5, trControl = fitControl, metric = "ROC")
    
    model_fit_rf <- 
      caret::train(y ~ ., data = df, method = "ranger", 
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
    
    
    # Print the estimated remaining time
    print("Still working...")

    
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
    
    Group_B <- df_prevalence %>%
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
    prevalence_df$Feature <- ifelse(log(prevalence_df$deviation) > 2.5, prevalence_df$Feature, "")
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
  
})
  end_time <- Sys.time()  # Record the end time
  elapsed_time <- end_time - start_time  # Calculate the elapsed time
  # Print the estimated remaining time
  print("Job is done")
  print("This job took:")
  print(elapsed_time)  
}


# Make test data
set.seed(1234)
data1 <- matrix(rpois(10000, 50), nrow = 100)
data2 <- matrix(rpois(5000, 4), nrow = 50)
data31 <- matrix(rpois(50, 1000), nrow = 1)
data32 <- matrix(rpois(50, 1000), nrow = 1)
data33 <- matrix(rpois(50, 200), nrow = 1)
data41 <- matrix(rpois(50, 100), nrow = 1)
data42 <- matrix(rpois(50, 50), nrow = 1)
data43 <- matrix(rpois(50, 1000), nrow = 1)
data3 <- cbind(data31,data41)
data4 <- cbind(data32,data42)
data5 <- cbind(data33,data43)
data <- rbind(data3,data4,data5, data1, data2)
rm(data1, data2, data3, data4,data5)
rm(data31,data32,data33,data41,data42,data43)
colnames(data) <- c(paste0("Sample", 1:ncol(data)))
sample_data <- data.frame("Sample"= colnames(data),
                          "Group" = rep(c("GroupA", "GroupB"), each = 50))

df <- t(data) %>%
  as_tibble()

find_biomarkers(data = data,
                            sample_data = sample_data,
                            exp_var = "Group", 
                            datatype = "counts", 
                            method = "RF",
                            top_biomarker=0.1,
                            prevalence_tolerance=0.01)
   
  