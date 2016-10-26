# Useful routines
# the "contingency matrix" is a matrix of miRNA:gene interaction scores
# using built-in data.table functionality, a table is generated where
# all miRNAs are collected which have at least one interaction with a gene
# in the provided gene list. Final step is to cast this long-format object
# into wide format
build_contingency_matrix <- function(genelist,refmat){
  gene_matrix <- data.table(GeneName = genelist,key="GeneName")
  setkey(refmat,GeneName)
  gene_matrix <- gene_matrix[refmat,nomatch=0]
  dcast.data.table(gene_matrix,GeneName ~ miRID,value.var="Score")
}

# auxilliary function to generate the combinations for the mask matrix build
# function takes a character list of genenames and the minimum combination
# size to consider. Using sapply, all possible combinations are generated
# from the minimum combination size, up to the 'combination' including all genes
# check is performed in calling script to ensure minimum combination size
# is less than or equal to the length of the gene list
generate_combos <- function(genelist, mincomb){
  test <- sapply(c(mincomb:length(genelist)),function(x){
    as.list(data.frame(combn(genelist,x),stringsAsFactors = F))})
  Reduce(append,test)
}

# convert raw scores from the analysis using the parameters of the beta dist fit
std_score <- function(score1,miR,NGen,bmat){
  if(score1==0 | is.na(score1)) return(0)
  params <- get(bmat,pos=refspace.env)[miR,NGen,]
  fit.mnp <- params[1]
  fit.vrp <- params[2]
  (score1 - fit.mnp)/sqrt(fit.vrp)
}

