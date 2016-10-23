# Useful routines
build_contingency_matrix <- function(genelist,refmat){
  gene_matrix <- data.table(GeneName = genelist,key="GeneName")
  setkey(refmat,GeneName)
  gene_matrix <- gene_matrix[refmat,nomatch=0]
  dcast.data.table(gene_matrix,GeneName ~ miRID,value.var="Score")
}

# auxilliary function to generate the combinations for the mask matrix build
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

