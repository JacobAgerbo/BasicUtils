#!/usr/bin/env Rscript
# R Script for Generating barplots, with grouping and taxonomic information
# Author: Jacob Agerbo Rasmussen
# Load necessary packages


make_barplot <- function(data = data,
                         sample_data = NULL,
                         tax_data = NULL,
                         taxa = NULL,
                         group = NULL,
                         grouping = FALSE){
                           suppressWarnings({
                             suppressMessages({
                               
                               # add grouping informations
                               data = data
                               sample_data = sample_data
                               tax_data = tax_data
                               taxa = taxa
                               # set dependencies
                               # Define colors using ANSI escape codes
                               red <- "\033[31m"
                               green <- "\033[32m"
                               reset <- "\033[0m"
                               suppressPackageStartupMessages(library(tidyverse))
                               suppressPackageStartupMessages(library(ggpubr))
                               suppressPackageStartupMessages(library(ggplot2))
                               suppressPackageStartupMessages(library(hilldiv))
                               # set colors
                               source("https://raw.githubusercontent.com/david-barnett/microViz/11cdd447c4e7d59a75f7194b249f4ddb1b942867/R/distinct_palette.R")
                               palette <- distinct_palette(10, pal = "kelly", add = "grey90")
                               if (!is.null(tax_data)) {
                                 if (is.null(taxa)) {
                                   cat(paste(red, "Whoops!", reset, "Please add a taxonomic group, using the taxa parameter!\n"))
                                 } else {
                                   plot <- data %>% 
                                     hilldiv::tss() %>%
                                     as_tibble() %>%
                                     dplyr::mutate("taxa" = tax_data[,taxa]) %>%
                                     reshape2::melt() %>%
                                     group_by(taxa,variable) %>%
                                     reframe(
                                       value = sum(value)
                                     ) %>% 
                                     ggplot(aes(x = variable, y = value, fill = taxa)) + 
                                     geom_bar(stat = "identity") +
                                     ylab("Relative abundance") +
                                     xlab("Samples") +
                                     theme_minimal() +
                                     labs(fill = paste(taxa))
                                 }
                                 
                                 if (isTRUE(grouping)){
                                   
                                   plot <- data %>% 
                                     hilldiv::tss() %>%
                                     as_tibble() %>%
                                     dplyr::mutate("taxa" = tax_data[,taxa]) %>%
                                     reshape2::melt() %>%
                                     group_by(taxa,variable) %>%
                                     reframe(
                                       value = sum(value)) %>%
                                    pivot_wider(names_from = taxa, id_cols = variable) %>%
                                     dplyr::mutate(group = sample_data[,group]) %>% 
                                   rename(Sample = variable) %>%
                                   reshape2::melt() %>%
                                     group_by(group,variable) %>%
                                     reframe(
                                       value = mean(value)) %>%
                                     ggplot(aes(x = group, y = value, fill = variable)) + 
                                     geom_bar(stat = "identity") +
                                     ylab("Mean relative abundance") +
                                     xlab(paste(group)) +
                                     theme_minimal() +
                                     labs(fill = paste(taxa))
                                 }
                                 
                               } else {
                                 
                                 data %>% 
                                   hilldiv::tss() %>%
                                   reshape2::melt() %>%
                                   ggplot(aes(x = Var2, y = value, fill = Var1)) + geom_bar(stat = "identity")
                                 
                               }

          plot <- plot + scale_color_manual(values = palette)

          return(plot)
                                       
        })
    })
}


# Generate test data
#set.seed(1234)
#data1 <- matrix(rpois(4700, 50), nrow = 47)
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
#rm(data1, data2, data3, data4,data5)
#rm(data31,data32,data33,data41,data42,data43)
#colnames(data) <- c(paste0("Sample", 1:ncol(data)))
#rownames(data) <- c(paste0("F", 1:nrow(data)))

#sample_data <- data.frame("Sample"= colnames(data),
#                          "Group" = rep(c("GroupA", "GroupB", "GroupC", "GroupD", "GroupE"), each = 20))


#tax_data <- data.frame("Phylum"= rep(c("PhyA", "PhyB", "PhyC", "PhyD", "PhyE"), each = 10),
#                       "Order"= rep(c("OrdA", "OrdB", "OrdC", "OrdD", "OrdE"), each = 10),
#                       "Class"= rep(c("ClassA", "ClassB", "ClassC", "ClassD", "ClassE"), each = 10),
#                       "Family"= rep(c("FamA", "FamB", "FamC", "FamD", "FamE"), each = 10),
#                       "Genus"= rep(c("GenA", "GenB", "GenC", "GenD", "GenE", "GenF", "GenG", "GenH", "GenI", "GenJ"), each = 10))


#make_barplot(data = data,
#                         sample_data = sample_data,
#                         tax_data = tax_data,
#                         taxa = "Genus",
#                         grouping = TRUE, group = "Group")
