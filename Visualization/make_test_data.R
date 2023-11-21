#!/usr/bin/env Rscript
# Make test data

make_test_data <- function(no_samples=100){
  set.seed(1234)
  data1 <- matrix(rpois(no_samples*100, 50), nrow = 100)
  data31 <- matrix(rpois(no_samples/2, 1000), nrow = 1)
  data32 <- matrix(rpois(no_samples/2, 1000), nrow = 1)
  data33 <- matrix(rpois(no_samples/2, 200), nrow = 1)
  data41 <- matrix(rpois(no_samples/2, 100), nrow = 1)
  data42 <- matrix(rpois(no_samples/2, 50), nrow = 1)
  data43 <- matrix(rpois(no_samples/2, 1000), nrow = 1)
  data3 <- cbind(data31,data41)
  data4 <- cbind(data32,data42)
  data5 <- cbind(data33,data43)
  data <- rbind(data3,data4,data5, data1)
  rm(data1, data3, data4,data5)
  rm(data31,data32,data33,data41,data42,data43)
  colnames(data) <- c(paste0("Sample", 1:ncol(data)))
  rownames(data) <- c(paste0("F", 1:nrow(data)))
  sample_data <- data.frame("Sample"= colnames(data),
                            "Group" = rep(c("GroupA", "GroupB","GroupC","GroupD","GroupE"), each = no_samples/5),
                            "Batch" = rep(c("Batch1", "Batch2","Batch3","Batch4","Batch5","Batch6","Batch7","Batch8","Batch9","Batch10"), each = no_samples/10),
                            "O_Group" = rep(c("GroupA", "GroupB"), each = no_samples/2),
                            "Treatment" = c(seq(0, 1, length.out = no_samples) + rnorm(no_samples, mean = 0, sd = 0.05)))
  
  list <- list("data" = data,
               "sample_data" = sample_data)
  
  # Define colors using ANSI escape codes
  red <- "\033[31m"
  green <- "\033[32m"
  reset <- "\033[0m"
  
  cat(paste(green,"\n","Here is some test data for you","\n",reset, sep = ""))
  
  return(list)
  
}


#test <-make_test_data(100)
