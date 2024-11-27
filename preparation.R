################################################################################
# This script is used to clean the envrionment and import all necessary packages
################################################################################


# Specify the list of packages that need to be imported ########################
list.of.packages <- c("heatmaply",
                      "dplyr",
                      "DT", 
                      "ggplot2", 
                      "tidyverse", 
                      "tximport", 
                      "tximportData",
                      "comprehenr",
                      "stringr",
                      "AnnotationDbi",
                      "vroom",
                      "plotly",
                      "heatmaply",
                      "DT",
                      "janitor",
                      "knitr",
                      "formattable",
                      "triwise",
                      "hash"
)

bioc.packages <- c("PCAtools", 
                   "tximport",
                   "DESeq2")

# Check if packages are installed ##############################################

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
new.bioc.packages <- bioc.packages[!(bioc.packages %in% installed.packages()[,"Package"])]

# Install new packages #########################################################
if(length(new.packages)) install.packages(new.packages)
if(length(new.bioc.packages)) BiocManager::install(new.bioc.packages)

# Import all packages ##########################################################
package_loading_Status <- lapply(list.of.packages, 
                                 require, 
                                 character.only = TRUE)

package_loading_Status_bioc <- lapply(bioc.packages, 
                                      require, 
                                      character.only = TRUE)


# Load the input sample sheet ##################################################


# EOF ##########################################################################


