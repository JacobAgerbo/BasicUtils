#' BasicUtils - Generate a barplot of a provided dataset
#' 
#' The `make_barplot` function is used for generating a barplot of a provided dataset. 
#' If `tax_data` and `taxa` are specified, the function generates a stacked barplot showing the relative abundance of each taxonomic group in each sample. 
#' If `grouping` is set to TRUE, the function generates a grouped barplot showing the mean relative abundance of each taxonomic group in each sample group. 
#' If `tax_data` is not specified, the function generates a simple barplot showing the values of each variable in each sample.
#'
#' @param data The dataset to be analyzed
#' @param sample_data The sample information dataset
#' @param tax_data The taxonomic information dataset
#' @param taxa The taxonomic group to be plotted
#' @param group The grouping variable to be used for generating a grouped barplot
#' @param grouping A logical value indicating whether to generate a grouped barplot
#' @return A ggplot barplot
#' @import tidyverse 
#' @import ggpubr
#' @import ggplot2
#' @import hilldiv
#' 
#' @export

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
                                 
                                 plot <- data %>% 
                                   hilldiv::tss() %>%
                                   reshape2::melt() %>%
                                   ggplot(aes(x = Var2, y = value, fill = Var1)) + geom_bar(stat = "identity")
                                 
                               }

          plot <- plot + scale_fill_manual(values = palette) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


          return(plot)
                                       
        })
    })
}