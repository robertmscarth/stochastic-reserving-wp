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
results = results %>% mutate(
  Ult_1 = Ult_1_2007 + Ult_1_2008 + Ult_1_2009 + Ult_1_2010 + Ult_1_2011 + Ult_1_2012 + Ult_1_2013 + Ult_1_2014 + Ult_1_2015,
  Ult_2 = Ult_2_2008 + Ult_2_2009 + Ult_2_2010 + Ult_2_2011 + Ult_2_2012 + Ult_2_2013 + Ult_2_2014 + Ult_2_2015,
  Ult_3 = Ult_3_2009 + Ult_3_2010 + Ult_3_2011 + Ult_3_2012 + Ult_3_2013 + Ult_3_2014 + Ult_3_2015,
  Ult_4 = Ult_4_2010 + Ult_4_2011 + Ult_4_2012 + Ult_4_2013 + Ult_4_2014 + Ult_4_2015,
  Ult_5 = Ult_5_2011 + Ult_5_2012 + Ult_5_2013 + Ult_5_2014 + Ult_5_2015,
  Ult_6 = Ult_6_2012 + Ult_6_2013 + Ult_6_2014 + Ult_6_2015,
  Ult_7 = Ult_7_2013 + Ult_7_2014 + Ult_7_2015,
  Ult_8 = Ult_8_2014 + Ult_8_2015,
  Ult_9 = Ult_9_2015,
  Res_1 = Res_1_2007 + Res_1_2008 + Res_1_2009 + Res_1_2010 + Res_1_2011 + Res_1_2012 + Res_1_2013 + Res_1_2014 + Res_1_2015,
  Res_2 = Res_2_2008 + Res_2_2009 + Res_2_2010 + Res_2_2011 + Res_2_2012 + Res_2_2013 + Res_2_2014 + Res_2_2015,
  Res_3 = Res_3_2009 + Res_3_2010 + Res_3_2011 + Res_3_2012 + Res_3_2013 + Res_3_2014 + Res_3_2015,
  Res_4 = Res_4_2010 + Res_4_2011 + Res_4_2012 + Res_4_2013 + Res_4_2014 + Res_4_2015,
  Res_5 = Res_5_2011 + Res_5_2012 + Res_5_2013 + Res_5_2014 + Res_5_2015,
  Res_6 = Res_6_2012 + Res_6_2013 + Res_6_2014 + Res_6_2015,
  Res_7 = Res_7_2013 + Res_7_2014 + Res_7_2015,
  Res_8 = Res_8_2014 + Res_8_2015,
  Res_9 = Res_9_2015,
  CDR_1 = CDR_1_2007 + CDR_1_2008 + CDR_1_2009 + CDR_1_2010 + CDR_1_2011 + CDR_1_2012 + CDR_1_2013 + CDR_1_2014 + CDR_1_2015,
  CDR_2 = CDR_2_2008 + CDR_2_2009 + CDR_2_2010 + CDR_2_2011 + CDR_2_2012 + CDR_2_2013 + CDR_2_2014 + CDR_2_2015,
  CDR_3 = CDR_3_2009 + CDR_3_2010 + CDR_3_2011 + CDR_3_2012 + CDR_3_2013 + CDR_3_2014 + CDR_3_2015,
  CDR_4 = CDR_4_2010 + CDR_4_2011 + CDR_4_2012 + CDR_4_2013 + CDR_4_2014 + CDR_4_2015,
  CDR_5 = CDR_5_2011 + CDR_5_2012 + CDR_5_2013 + CDR_5_2014 + CDR_5_2015,
  CDR_6 = CDR_6_2012 + CDR_6_2013 + CDR_6_2014 + CDR_6_2015,
  CDR_7 = CDR_7_2013 + CDR_7_2014 + CDR_7_2015,
  CDR_8 = CDR_8_2014 + CDR_8_2015,
  CDR_9 = CDR_9_2015
  )

# Calculate statistics for all the variables
# First defined a function which calculates the statistics of interest
SummaryStatisticsFunc = function(x) c(
  "Mean"= mean(x,na.rm=TRUE),
  "Stand dev" = sd(x),
  "CoeffofVariation" = sd(x)/mean(x,na.rm=TRUE)
)
# Apply the function to get the statistics
SummaryStatistics = results %>% sapply( SummaryStatisticsFunc ) %>% as_tibble

##
t = SummaryStatistics %>% select( CDR_1_2007 , CDR_1_2008 , CDR_1_2009 , CDR_1_2010 , CDR_1_2011 , CDR_1_2012 , CDR_1_2013 , CDR_1_2014 , CDR_1_2015 )
