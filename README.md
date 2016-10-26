# metaMIR 

metaMIR is a microRNA (miRNA) framework to predict interactions in human between miRNAs and clusters of genes. The user provides a set of genes to be targeted, and optionally genes not to be targeted. The analysis is performed to identify miRNAs which may simultaneously interact with a number of genes.

The first part of the analysis makes use of individual prediction scores of interaction between one miRNA and one miRNA. These scores are derived from combination of previously established prediction algorithms. Using experimentally-validated interactions for training and testing, a machine learning approach is used to integrate the results of the established algorithm into a new data set. The resulting scores are based on class prediction probabilities; that is, the likelihood that a miRNA targets a target gene.

The individual scores are used in a combinatorial analysis. All combinations of genes (from a minimum cluster size up to inclusion of all genes) are generated and analyzed for miRNAs which will interact with them. The metaMIR algorithm will return miRNA candidates predicted to co-regulate groups of genes.
A webserver version of the script in operation can be found [here](http://rna.informatik.uni-freiburg.de/metaMIR/), hosted on servers of the Department of Bioinformatics, University of Freiburg, Germany.
