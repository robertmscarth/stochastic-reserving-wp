# R script to read in the bootstrap data and create the graphs for the paper

# Set the wd to the top level of the git repo checkout
setwd("C:/Users/rs09/Working Party/stochastic-reserving-wp/")

# load libraries
library("tidyverse")

# set file paths - these are all relative to the git repo
resultsFilePath = "Bootstrapping/Excel/Results/"
outputsFilePath = "/Paper/Paper2/Graphs/"

# Set the list of filenames for the results .csv files
#resultsFileNames = c( "Example1_Mack" , "Example2_Incurred_Mack" )
resultsFileName = "Example1_Mack"

# Read in the results from the .csv file
filePath = paste( resultsFilePath , resultsFileName , ".csv" , sep="" )
results = read_csv(filePath)

# Calculate the same summay columns as in the summarise spreadsheet

