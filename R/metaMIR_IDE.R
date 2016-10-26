#instructions for use:
# ensure downloaded files are in the same directory as this script
# metaMIR_analytical.R
# metaMIR_subs.R
# betafit_parameters_mv_geo.RData
# refposneg.RData
#
# list of genes to be analyzed should be at least 5 genes long. Up to
# 15 genes will be simultaneously analyzed. Longer lists will be partitioned
# according to the nearest neighbours (according to miRNA co-regulation)
# for a provided input file <input_name>.txt the following output files are
# generated:
# <input_name>_cluster_<date_time>.txt
#     - tab-delimited file with the clustering results for genelists longer than
#       the maximum size, showing the frequency of gene-pair occurrence in the
#       resulting lists
# <input_name>_rep_<date_time>.txt
#     - tab-delimited file describing the percentage representation of each
#       provided gene in the resulting lists, when more genes than the max
#       cluster size are provded.
# <input_name>_params_<date_time>.txt
#     - the parameters applied during execution.
# <input_name>_output_<date_time>.txt
#     - resulting output of the script, with miRNA, number of genes predicted
#       to be simultaneously targeted, the standardized score for the
#       corresponding combination, an aggregate score, incorporating standard
#       score and group size, the list from which the combination was derived,
#       and the genes predicted to be targeted.
# output files are written to a "results" folder under the current
# folder, which will be created if it does not exist.

# set flag to indicate script not being run in command line mode
comline <- FALSE

# Run the entire script. The result is the creation of a function called
# batch_metaMIR, which takes the following arguments:
# inputfile - name of file to be analyzed
# mincomb   - minimum gene combination size to include in results, default 5
# maxclust  - maximum gene combination size. Defaults to 15, however can be
#             lowered depending on system resources, if run times per
#             analysis become long, or the system crashes
# Zthresh   - minimum threshold for individual combination scores
#             The default of 2.3 excludes 90% of 5 gene combinations that
#             were detected in random searches
#
# For the analysis, the list of genes is provided via a plain text file, one gene
# symbol per line below, the function "metaMIR" can be run on an individual file
# according to:
#
#      batch_metaMIR("<filename>")
#
# where <filename> is the name of the text file with the gene list, enclosed in double
# quotation marks as above. Alternatively, a set of gene lists can be provided by
# specifying them manually and looping through using the provided for loop:
#
# input_files=c("<file1.txt>","<file2.txt>","<file3.txt>")
# for (zz in 1:length(input_files)){
#   metaMIR(input_files[zz])
# }
#
# or if filenames consist of a pattern,
#
# input_files = list.files(pattern = "<matching pattern common to your filenames>")
#
#
# check for installation status of required packages
pkgs = c("data.table","stringr")
if(length(new.pkgs <- setdiff(pkgs, rownames(installed.packages())))>0) install.packages(new.pkgs)
rm(pkgs,new.pkgs)
suppressPackageStartupMessages(library(data.table,quietly=T))
suppressPackageStartupMessages(library(stringr,quietly=T))
# load reference database of predictions and beta dist fitting parameters
refspace.env <- new.env()
load("refposneg.RData",envir=refspace.env)
load("betafit_mv.RData",envir=refspace.env)

# load parameters for score standardization. Contains the mean and variance
# calculated from the shape1, shape2 fitting parameters determined for each miRNA
# source("metaMIR_analytical_GH.R")
source("metaMIR_analytical.R")
######################################################################
# change parameters here if desired
mincomb  = 5
maxclust = 15
Zthresh  = 1.1

# batch_metaMIR(inputfile,mincomb,maxclust,Zthresh)
inventory <- list.files(pattern="example",full.names = T)

for (zz in 1:length(inventory)){
batch_metaMIR(inventory[[zz]],mincomb = mincomb,maxclust = maxclust,Zthresh = Zthresh)
}
