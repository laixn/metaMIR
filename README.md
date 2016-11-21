# **metaMIR** 
**metaMIR** is a framework to predict in human interactions between 
microRNAs (miRNA) and clusters of genes. The user provides a set of genes to be targeted, 
and optionally genes not to be targeted. The analysis is performed to identify 
miRNAs that may simultaneously interact with a number of genes.

The first part of the analysis makes use of individual prediction scores of 
interaction between one miRNA and one miRNA. These scores are derived from 
combination of previously established prediction algorithms. Using 
experimentally-validated interactions for training and testing, a machine 
learning approach is used to integrate the results of the established algorithm 
into a new data set. The resulting scores are based on class prediction 
probabilities; that is, the likelihood that a miRNA targets a target gene.

The individual scores are used in a combinatorial analysis. All combinations of
genes (from a minimum cluster size up to inclusion of all genes) are generated 
and analyzed for miRNAs which will interact with them. The **metaMIR** algorithm 
will return miRNA candidates predicted to co-regulate groups of genes.

A webserver version for ad hoc usage is part of the 
[Freiburg RNA tools webserver](http://rna.informatik.uni-freiburg.de/metaMIR/), hosted on 
servers of the 
[Bioinformatics Group, Technical Faculty, University of Freiburg, Germany](http://www.bioinf.uni-freiburg.de).

## Usage and Installation
In order to run **metaMIR** locally you have to

* download the source files (`R` folder content of this github repo or of [release file](https://github.com/rnagear/metaMIR/releases))
* download the precompiled reference data [`miRNA_predictions.tar.gz`](http://www.bioinf.uni-freiburg.de/Software/metaMIR/miRNA_predictions.tar.gz)
* store and decompress both downloads in one directory
* call one the main scripts (for details see README in `R` folder)
  * command line : `R --file=metaMIR_cl.R --args <ARGUMENTS_SEE_R_SCRIPT>`
  * Rstudio : open and run `metaMIR_IDE.R`

