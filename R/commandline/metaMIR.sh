#!/bin/bash

# load correct R version
source /usr/local/tools/R/3.2.1/iuc/package_r_3_2_1/e46a7803f17b/env.sh;

# setup tmp dir for R package installation
export TMPDIR=/scratch/rna/bisge001/Software/metaMIR/1.0.0/bin/tmp;
mkdir -p $TMPDIR;

# call metaMIR
R --slave --vanilla --file=/scratch/rna/bisge001/Software/metaMIR/1.0.0/bin/metaMIR.R --args "$@"

