
# instructions for use:

# ensure downloaded files are in the same directory as this script
# metaMIR_analytical.R
# metaMIR_subs.R
# betafit_mv.RData
# refposneg.RData

# Run the entire script. The result is the creation of a function called
# batch_metaMIR, which takes the following arguments:

# inputfile - name of file to be analyzed, which is a
#             list of genes (one per line) provided via a plain text file
# mincomb   - minimum gene combination size to include in results (default=5)
# maxclust  - maximum gene combination size (default=15)
# Zthresh   - minimum threshold for individual combination scores (default=1)
# outPath   - (optional) the path for the output files

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

scriptPath <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE);
        needle <- "--file=";
        match <- grep(needle, cmdArgs);
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(dirname(sub(needle, "", cmdArgs[match]))));
        } else {
		stop("run with argument --file=metaMIR.R");
	}
}


#################
# read commandline args
args <- commandArgs(trailingOnly=TRUE)

# ensure enough arguments
if (length(args)<4 ) {
	stop("no arguments given!\n\n\t-> use <inputfile> <mincomb> <maxclust> <zthresh> [optOutPath]");
}
outPath <- "."
if (length(args)>=5) {
	outPath <- args[5];
}


# get the absolute path of this script
metaMirPath <- scriptPath();


# check for installation status of required packages
pkgs <- c("data.table","stringr","reshape2")
if(length(new.pkgs <- setdiff(pkgs, rownames(installed.packages(lib.loc=c(.libPaths(),metaMirPath)))))>0) install.packages(new.pkgs, lib=metaMirPath, repos="http://cran.us.r-project.org", destdir=metaMirPath)
rm(pkgs,new.pkgs)

suppressPackageStartupMessages(library(data.table,quietly=T,lib.loc=c(.libPaths(),metaMirPath)))
suppressPackageStartupMessages(library(stringr,quietly=T,lib.loc=c(.libPaths(),metaMirPath)))
suppressPackageStartupMessages(library(reshape2,quietly=T,lib.loc=c(.libPaths(),metaMirPath)))


# load reference database of predictions and beta dist fitting parameters
refspace.env <- new.env()
load( paste(metaMirPath,"refposneg.RData",sep="/"), envir=refspace.env)
load( paste(metaMirPath,"betafit_mv.RData",sep="/"), envir=refspace.env)

# # load parameters for score standardization. Contains the mean and variance
# # calculated from the shape1, shape2 fitting parameters determined for each miRNA
source( paste(metaMirPath,"metaMIR_analytical.R",sep="/") )




#################
# call metaMIR using the command line arguments
batch_metaMIR( inputfile=args[1], mincomb=as.numeric(args[2]), maxclust=as.numeric(args[3]), Zthresh=as.numeric(args[4]), outPath=outPath);


