# About
This folder contains the source code for metaMIR.
The code can be run at the command line, or in an IDE such as RStudio.
Execution of the main file creates a function called batch_metaMIR taking the following arguments:
* inputfile - list of genes (HGNC format) to be analyzed
* mincomb - the minimum combination size to be considered when perfoming the combinatorial analysis (default 5)
* maxclust - the maximum number of genes to analyze simultaneously (default 15)
* Zthresh - the minimum score per miRNA:gene-cluster pair to consider (default 1.1)
Further description is contained in the main batch file.

# Usage
To run metaMIR, edit or execute the main file (one of):
* metaMIR_cl.R - for execution at the command line with a list of arguments
* metaMIR_IDE.R - for execution in an IDE

Dependency files that need to appear in the same directory are:
* metaMIR_analytical.R
* metaMIR_subs.R
* refnegpos.RData (separate download needed, see below)
* betafit_mv.RData

Gene input lists must be provided in text files, one gene per line in 
[HGNC](http://www.genenames.org/) format. Default operation is to create a 
'Results' subfolder under the current folder into which results will be written.

# Contents
This GitHub repository contains the script files and the binary fitting parameters R-object.
* betafit_mv.RData
* metaMIR_cl.R
* metaMIR_IDE.R
* metaMIR_subs.R
* metaMIR_analytical.R

# Additional files
For the script to run, the `refnegpos.RData` binary R object file is also 
required. This file is too large to host in this repository, but is available 
for free download [here](http://www.bioinf.uni-freiburg.de/Software/metaMIR/miRNA_predictions.tar.gz)
