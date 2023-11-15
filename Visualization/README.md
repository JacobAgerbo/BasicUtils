# Basic R functions for visualization of multi-variate analysis

*Currently this haven't been made to an R-package, but we are on it*

This repository contains a basic R scripts (functions) that can be used for analyzing multi-variate analysis from bacteria and host-microbe interactions using omics data. 

## Visualization Features
The R functions is thoroughly documented, providing explanations and examples for each function to facilitate its usage and customization. 
Main functions includes preprocessing multi-variate data from omics, such as filtering, normalization, and transformation. Furthermore, it offers various statistical analysis methods for identifying important feature in the dataset. Key functions are listed below.

- **Data overview**: A function getting an overview of features in dataset. 
  
- **Preprocessing**: A function for preprocessing omics data, such as quality control, filtering, normalization, and transformation.

- **Collector's curves**: A function for investigating sufficient data collection.
  
- **Sample-size based extrapolation of richness in data**: The script is thoroughly documented, providing explanations and examples for each function to facilitate its usage and customization.

- **Boxplots and statistical comparison of groups**: A function statistical and visual comparison of features between groups

- **Composition of features, using barplots**: A function for generating the barplot, which always are used in microbiome studies. Here you can group samples per group, if you like.

- **Linaer mixed effect models**: If you know, you know.

- **Principal component analysis (PCA)**: If you know, you know.


## Tutorial
For reproducability purposes I generated a function for generating sandbox data. The output is a list of two dataset. A abundance matrix and a sample information sheet. 

```{R}
# Get function
source("https://raw.githubusercontent.com/JacobAgerbo/Basic_Utils/main/Visualization/make_test_data.R")

# Make data
data_list <- make_test_data(100)
# Export from list to datasets
data <- data_list$data
sample_data <- data_list$sample_data
```

### Data overview

```

```

### Preprocessing

```

```

### Collector's curves

```
# get function
source("https://raw.githubusercontent.com/JacobAgerbo/Basic_Utils/main/Visualization/make_collectors_curves.R")

# generate collector's curves
make_collectors_curves(data = data, interations = 20)
```
![alt text](http://url/to/img.png)
### Sample-size based extrapolation of richness in data
This is based on the iNEXT *R package*

```

```