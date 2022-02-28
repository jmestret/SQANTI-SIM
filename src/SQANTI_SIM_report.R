#######################################
#                                     #
#    SQANTI-SIM report generation     #
#                                     #
#######################################

# Author: Jorge Meste
# Last modified: 28/02/2022 by Jorge Mestre


#######################################
#                                     #
#      PACKAGES AND LIBRARIES         #
#                                     #
#######################################

suppressMessages(library(gridExtra))
suppressMessages(library(knitr))
suppressMessages(library(rmarkdown))
suppressMessages(library(ggplot2))


#######################################
#                                     #
#             FUNCTIONS               #
#                                     #
#######################################




#######################################
#                                     #
#                MAIN                 #
#                                     #
#######################################

# -------------------- Argument parser
args <- commandArgs(trailingOnly = TRUE)
class.file <- args[1] # classification file SQANTI3
junc.file <- args[2] # junctions file SQANTI3
del.file <- args[3] # deleted transcripts file preparatory step
expr.file <- args[4] # sim counts file sim step








