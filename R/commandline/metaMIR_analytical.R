
# get the absolute path of this script 
metaMirPath <- normalizePath(dirname(sys.frame(1)$ofile));



batch_metaMIR <- function(inputfile, mincomb=5, maxclust=15, Zthresh=1, outPath=NA){

print(paste("input=",inputfile,mincomb,maxclust,Zthresh,outPath));

  if (missing(inputfile)){
    stop("Missing input file of genes to analyze.")
  }
  if(Zthresh < 0.75){
    warning("Zthreshold < 1 - results may be noisy.")
  }
  if(maxclust > 25){
    print(paste("Attempting to set maximum cluster size too high",maxclust,"Resetting to 15."))
    maxclust <- 15
  } else if(maxclust < 5){
    print(paste("Attempting to set maximum cluster size too low",maxclust,"Resetting to 10."))
    maxclust <- 10
  }
  # calculate geometric mean of a numeric vector
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }

  # generate the geometric mean from row matrix input for all combinations
  xtrgm <- function(x,mincomb){
    x <- x[1,!is.na(x),drop=F]
    if (length(x) < mincomb) return(NULL)
    if (length(x) == mincomb){
      combos <- paste(colnames(x),collapse=",")
      scores <- gm_mean(x)
      groups <- length(x)
      data.table(miRNA=rownames(x),Score=scores,Group=groups,Combo=combos)
    } else {
      combos <- generate_combos(colnames(x),mincomb)
      scores <- vapply(generate_combos(x,mincomb),FUN = gm_mean,FUN.VALUE = numeric(1))
      groups <- vapply(combos,length,FUN.VALUE = numeric(1))
      combos <- vapply(combos,paste, collapse=",",FUN.VALUE=character(1))
      tmp <- data.table(miRNA=rownames(x),Score=scores,Group=groups,Combo=combos)
      setorder(tmp,-Group,-Score)
      tmp[,head(.SD,1),by=Group]
    }
  }

  # establish subroutines necessary for score calulation and matrix creation
  source( paste(metaMirPath,"metaMIR_subs.R",sep="/") );

  # if list is longer than the defined max cluster size, partition the list
  # finding the (maxclust - 1) nearest neighbours
  partgenes <- function(inputgenes, anatype,mat){
    print(paste0("Partitioning ",anatype," genelist"))
    cm <- build_contingency_matrix(inputgenes,get(mat))
    cm2 <- as.matrix(cm[,-1,with=F])
    rownames(cm2) <- cm[,GeneName]
    dcm <- dist(cm2)
    adjm <- apply(cm2,1,function(x){apply(cm2,1,function(y){length((x+y)[!is.na(x+y)])})})
    prod <- as.matrix(dcm)/sqrt(adjm)
    dcm.l <- melt(data.frame(var1=rownames(prod),prod,check.names=F),id.vars="var1")
    colnames(dcm.l) <- c("var1","var2","value")
    dcm.l$var1 <- as.character(dcm.l$var1)
    dcm.l$var2 <- as.character(dcm.l$var2)
    dcm.l <- dcm.l[dcm.l$var1!=dcm.l$var2,]
    dcm.dt <- data.table(dcm.l,key=c("var1","var2"))
    setorder(dcm.dt,var1,value)
    if(length(inputgenes)>maxclust){
      dcm.dt <- dcm.dt[,head(.SD,(maxclust-1)),by=.(var1)]
      glists <- lapply(unique(dcm.dt[,var1]),function(x){unique(unlist(dcm.dt[var1==x,c("var1","var2"),with=F]))})
      glists <- lapply(glists,function(x){x[order(x)]})
      tmp <- lapply(glists,function(x){
        data.table(expand.grid(Gene1=unlist(x),Gene2=unlist(x),stringsAsFactors = F))})
      tmp2 <- data.table(matrix(unlist(tmp),ncol=2,byrow=T,dimnames = list(NULL,c("Gene1","Gene2"))))
      tmp2 <- tmp2[Gene1!=Gene2][,"Norm":=.N,by=.(Gene1,Gene2)][,"Inter":="pp"]
      tmp2 <- unique(tmp2)
      tmp2 <- tmp2[,.(Gene1,Inter,Gene2,Norm=Norm/min(Norm))]
	  if (is.na(outPath)) {
		  fname1 <- paste0("Results/",gsub(".txt","",ifn),"_cluster_",substring(anatype,1,3),gsub(" ","-",gsub(":","_",Sys.time())),".txt")
	  } else {
		  fname1 <- paste(outPath,paste(ifn,"cluster",substring(anatype,1,3),"txt",sep="."),sep="/")
	  }
      write.table(tmp2,file=fname1,col.names=T,row.names=F,quote=F,sep="\t")
      glists <- unique(glists)
      names(glists) <- paste0("L",sprintf("%02d",seq(1,length(glists),1)))
      rep_data <- melt(summary(as.factor(unlist(glists))))
      rep_data <- data.frame(GeneName=rownames(rep_data),Rep=round(rep_data$value*100/length(glists),1))
      rep_data <- rep_data[order(rep_data$Rep,decreasing=T),]
      fname2 <- gsub("_cluster_","_rep_",fname1)
	  if (is.na(outPath)) {
      	fname2 <- gsub("_cluster_","_rep_",fname1)
	  } else {
      	fname2 <- gsub("cluster","rep",fname1)
	  }
	  write.table(rep_data,file=fname2,sep="\t",quote=F,col.names=T,row.names=F)
      return(glists)
    } else {
      # take the lowest 25th percentile of distances
      # these will be duplicated, so take every other row
      dcmout <- dcm.dt[,.(Gene1=var1,Inter="pp",Gene2=var2,Norm=value*100/max(value))]
#       fname1 <- paste0("Results/",gsub(".txt","",ifn),"_cluster_",substring(anatype,1,3),gsub(" ","-",gsub(":","_",Sys.time())),".txt")
#       write.table(dcmout,file=fname1,col.names=T,row.names=F,quote=F,sep="\t")
      list(L1=inputgenes)
    }

  }

  metaMIRg <- function(inputgenes,state,mincomb=5){
    if(state=="pos"){
      pmat <- "posmat"
      bmat <- "beta_mvp"
      mircounts <- mircountsp
    } else {
      pmat <- "negmat"
      bmat <- "beta_mvn"
      mircounts <- mircountsn
    }
    mincomb <- as.numeric(mincomb);
    cmp <- build_contingency_matrix(inputgenes,get(pmat))
    # the generated contingency matrix is a data.table. Convert to matrix
    # and transpose, assigning appropriate rownames
    cmp2 <- t(as.matrix(cmp[,-1,with=F]))
    selv <- apply(cmp2,1,function(x){length(x[!is.na(x)]) >= mincomb})
    cmp2 <- cmp2[selv,,drop=F]
    if(is.null(cmp2) | nrow(cmp2)==0) return(NULL)
    colnames(cmp2) <- sort(inputgenes)
    outlist <- list()
    for (i in 1:nrow(cmp2)){
      outlist[[i]] <- xtrgm(cmp2[i,,drop=F],mincomb)
    }
    outmat <- Reduce(rbind,outlist)
    if (is.null(outmat)) return(NULL)

    outmat[,Score:=max(Score),by=.(miRNA,Group)]
    setkey(outmat,miRNA,Score,Group)
    outmat <- unique(outmat)
    #     outmat[,Score := mapply(std_score,Score,miRNA,Group)][,Score:=max(Score),by=miRNA]
    outmat[,Score := mapply(std_score,Score,miRNA,Group,bmat)]
    outmat[,adjScore:=Group*Score]
    outmat[,Score:=sum(adjScore)/sum(Group),by=miRNA][,adjScore:=NULL]
    outmat[,Score:=max(Score),by=.(miRNA,Group)]
    setkey(outmat,miRNA,Score,Group)
    outmat <- unique(outmat)
    outmat <- merge(outmat,mircounts,by="miRNA")
    outmat[,Prob:=(choose(miRN,Group)*choose((totNG-miRN),(maxclust-Group)))/(choose(totNG,maxclust))]
    outmat <- outmat[Score>0,newPS:=log10(Score/Prob)]
    outmat[,c("Score","newPS","Prob","miRN"):=.(newPS,NULL,NULL,NULL)]
    outmat <- outmat[!is.na(Score)]
    setorder(outmat,miRNA,-Score,-Group)
    outmat <- outmat[,head(.SD,1),by=miRNA]
    setcolorder(outmat,c("Combo","miRNA","Score","Group"))
    outmat[Score>=Zthresh]
  }

  # create master list of all available genes
  refDB_genes <- unique(get("posDB",pos=refspace.env)[,GeneName])
  posDB <- get("posDB",pos=refspace.env)
  negDB <- get("negDB",pos=refspace.env)
  # get counts of targeted genes per miRNA
  mircountsp <- posDB[,.N,by=miRID]
  mircountsn <- negDB[,.N,by=miRID]
  setnames(mircountsp,c("miRNA","miRN"))
  setnames(mircountsn,c("miRNA","miRN"))
  # get total number of genes in (positive) database
  totNG <- length(refDB_genes)

  # function begins here. For construction, manually set the arguments
  # method = 0  ("general")
  # method = 1  ("proportion")
#   method=0
  # check presence of Results subfolder
  if(!("Results" %in% list.dirs(full.names=F))) dir.create("Results")

  # read in the genelist and perform checks regarding whether there are genes
  # to be tracked, or for negative analysis
  genelist.raw <- readLines(inputfile,warn = F)
  genelist.mat <- data.table(src=genelist.raw,kern=NA,pos=NA,neg=NA)
  genelist.mat[,c("kern","pos","neg"):=
                 .(str_detect(src,"\\*"),!str_detect(src,"^\\-"),str_detect(src,"^\\-"))]
  genelist.mat[,src:=str_replace(str_replace(str_trim(src,side = "both"),"\\*",""),"^\\-","")]
  genelist.mat <- unique(genelist.mat)
  genelist.src <- genelist.mat[,src]
  if (length(genelist.src)<maxclust) maxclust <- length(genelist.src)
  if (length(setdiff(genelist.src,refDB_genes))!=0) {
    print("Supplied genes not all in table. Missing ID(s):")
    print(paste(setdiff(genelist.src,refDB_genes),collapse=","))
    print("Seek alternate name/symbol, or no prediction information available")
    return(NULL)
  } else {
    print("All entered genes ok.")
  }
  if (any(genelist.mat[,kern])){
    genelist.kern <- genelist.mat[kern==T,src]
    method=1
  } else {
    genelist.kern <- NULL
  }

  # set posflag/negflag to track whether the length of the input list
  # is longer than the absolute minimum of 4 and the minimum combo size
  if (any(genelist.mat[,pos])){
    genelist.pos <- genelist.mat[pos==T,src]
    posflag <- (length(genelist.pos) > max(c(4,mincomb)))
  } else {
    posflag <- F
    genelist.pos <- NULL
  }
  if (any(genelist.mat[,neg])){
    genelist.neg <- genelist.mat[neg==T,src]
    negflag <- (length(genelist.neg) > max(c(4,mincomb)))
  } else {
    negflag <- F
    genelist.neg <- NULL
  }
  # to speed subsequent calculations, subset the refDB to the provided gene set
  submat <- rbind(get("posDB",pos=refspace.env)[GeneName %in% genelist.src],
                  get("negDB",pos=refspace.env)[GeneName %in% genelist.src])

  # control structure to ensure assay size large enough
  if((is.null(genelist.pos) & !negflag) | (is.null(genelist.neg) & !posflag)){
    cat("Insufficient gene input. Need minimum 5 genes or minimum combination selected.")
    return(NULL)
  }
  if(!posflag & !negflag){
    cat("Insufficient gene input. Need minimum 5 genes.")
    return(NULL)
  }

  # add further checks on the length of gene lists entered
  inlenp <- length(genelist.pos)
  inlenn <- length(genelist.neg)
  # strip directory information from input filename
  ifn <- unlist(strsplit(inputfile,split="/"))
  ifn <- ifn[[length(ifn)]]
  glistpos <- glistneg <- list()
  # define gene lists to be analyzed, partitioning if necessary
  if(posflag) {
    posman <- ifelse(inlenp < 7,T,F)
    glistpos <- partgenes(genelist.pos,anatype = "positive",mat = "submat")
  } else if(inlenp > 0 & inlenp < 7){
    posman <- T
  } else posman <- F

#   if(posflag & inlenp > maxclust) {
#     posman <- F
#     glistpos <- partgenes(genelist.pos,anatype = "positive",mat = "submat")
#   } else if(posflag & inlenp <= maxclust){
#     posman <- F
#     glistpos <- partgenes(genelist.pos,anatype = "positive",mat = "submat")
#   } else if(any(genelist.mat[,pos])){
#     posman <- T
#   } else posman <- F



  if(negflag){
    negman <- ifelse(inlenn < 7,T,F)
    glistneg <- partgenes(genelist.neg,anatype = "negative",mat = "submat")
  } else if(inlenn > 0 & inlenn < 7){
    negman <- T
  } else negman <- F

#   if(negflag & inlenn > maxclust){
#     negman <- F
#     glistneg <- partgenes(genelist.neg,anatype = "negative",mat = "submat")
#   } else if(negflag & inlenp <= maxclust){
#     negman <- F
#     glistneg <- partgenes(genelist.neg,anatype = "negative", mat = "submat")
#   } else if(any(genelist.mat[,neg])){
#     negman <- T
#   } else negman <- F
  # add control statement to check status of manual flags.
  # both cannot be true, or one true, other non-existent
  if (posman & negman){
    negman <- F
  } else if (posman & (is.null(genelist.neg) | !negflag)) {
    posman <- F
  } else if (negman & (is.null(genelist.pos) | !posflag)){
    negman <- F
  }


  posmat <- submat[Score > 0.5]
  negmat <- copy(submat[Score < -0.4])[,Score:=-Score]
  # analyze positive gene list if necessary
  if(!is.null(genelist.pos) & !posman){
    # run metaMIR on pos
    elaps.timep <- system.time({
      final_tablep <- data.table(matrix(nrow=0,ncol=5))
      setnames(final_tablep,c("Combo","miRNA","Score","Group","List"))
      print(paste0("Analyzing ",length(glistpos)," ",ifelse(length(glistpos)>1,"lists","list")," from ",length(genelist.pos)," pos genes entered."))
      pb = txtProgressBar(min = 0, max = length(glistpos), initial = 0,style=3)
      for (i in 1:length(glistpos)){
        outtab <- metaMIRg(glistpos[[i]],state = "pos",mincomb)
        if(!is.null(outtab)){
          outtab[,List := names(glistpos)[i]]
          final_tablep <- do.call(rbind,list(final_tablep,outtab))
        }
        setTxtProgressBar(pb,i)
      }
    })
    if (nrow(final_tablep)==0) {print("No positive miRNAs above threshold")
      return(NULL)}

    if(negman){
      # manually scan through pos list from final_tablep
      elaps.timep <- elaps.timep + system.time({
        mirs <- final_tablep[,miRNA]
        genes <- lapply(mirs,function(x){
          intersect(negmat[miRID==x,GeneName],genelist.neg)})
        Scores <- mapply(function(x,y){gm_mean(negmat[GeneName %in% x & miRID==y,Score])},
                         genes,mirs)
        Grps <- sapply(genes,length)
        genes <- sapply(genes,paste,collapse=",")
        final_tablep[,c("NegGroup","NegScore","NegCombo"):=
                       .(Grps,Scores,genes)]
        final_tablep <- final_tablep[NegGroup > 1]
        final_tablep[,NegScore:=mapply(std_score,NegScore,miRNA,NegGroup,bmat = "beta_mvn")]
        final_tablep <- merge(final_tablep,mircountsn,by="miRNA")
        final_tablep[,Prob:=(choose(miRN,NegGroup)*choose((totNG-miRN),(inlenn-NegGroup))/choose(totNG,inlenn))]
        final_tablep <- final_tablep[NegScore>0]
        final_tablep[,newPS := -log10(NegScore/Prob)]
        final_tablep[,NegScore := round(newPS,3)][,c("newPS","Prob","miRN"):=NULL]
        final_tablep <- final_tablep[NegScore > 0]
        final_tablep[,finScore:=Score+NegScore]
      })
    } else {
      final_tablep[,finScore:=Score]
    }
    final_tablep[,MeanAg:=mean(finScore),by=miRNA][,MaxAg:=max(finScore),by=miRNA]
    final_tablep <- final_tablep[finScore==MaxAg][,"MaxAg":=NULL]
    final_tablep[,c("Score","MeanAg","finScore"):=.(round(Score,3),round(MeanAg,2),round(finScore,3))]
    # setcolorder(final_tablep,c("miRNA","Group","Score","MeanAg","List","Combo"))
    setkeyv(final_tablep,setdiff(names(final_tablep),"List"))
    final_tablep <- unique(final_tablep)
  } else elaps.timep <- system.time({NULL})

  # analyze negative gene list if necessary
  if(!is.null(genelist.neg) & !negman){
    # run metaMIR on neg
    elaps.timen <- system.time({
      final_tablen <- data.table(matrix(nrow=0,ncol=5))
      setnames(final_tablen,c("Combo","miRNA","Score","Group","List"))
      print(paste0("Analyzing ",length(glistneg)," ",ifelse(length(glistneg)>1,"lists","list")," from ",length(genelist.neg)," neg genes entered."))
      pb = txtProgressBar(min = 0, max = length(glistneg), initial = 0,style=3)
      for (i in 1:length(glistneg)){
        outtab <- metaMIRg(glistneg[[i]],state = "neg",mincomb)
        if(!is.null(outtab)){
          outtab[,List := names(glistneg)[i]]
          final_tablen <- do.call(rbind,list(final_tablen,outtab))
        }
        setTxtProgressBar(pb,i)
      }
    })
    if (nrow(final_tablen)==0) {print("No negative miRNAs above threshold")
      return(NULL)}
    if(posman){
      # manually scan through neg list from final_tablep
      mirs <- final_tablen[,miRNA]
      genes <- lapply(mirs,function(x){
        intersect(posmat[miRID==x,GeneName],genelist.pos)})
      Scores <- mapply(function(x,y){gm_mean(posmat[GeneName %in% x & miRID==y,Score])},
                       genes,mirs)
      Grps <- sapply(genes,length)
      genes <- sapply(genes,paste,collapse=",")
      final_tablen[,c("PosGroup","PosScore","PosCombo"):=
                     .(Grps,Scores,genes)]
      final_tablen <- final_tablen[PosGroup > 1]
      final_tablen[,PosScore:=mapply(std_score,PosScore,miRNA,PosGroup,bmat = "beta_mvp")]
      final_tablen <- merge(final_tablen,mircountsp,by="miRNA")
      final_tablen[,Prob:=(choose(miRN,PosGroup)*choose((totNG-miRN),(inlenp-PosGroup))/choose(totNG,inlenp))]
      final_tablen <- final_tablen[PosScore>0]
      final_tablen[,newPS:=-log10(PosScore/Prob)]
      final_tablen[,PosScore:=round(newPS,3)][,c("newPS","Prob","miRN"):=NULL]
      final_tablen <- final_tablen[PosScore > 0]
      final_tablen[,finScore:=Score+PosScore]
    } else {
      final_tablen[,finScore:=Score]
    }
    final_tablen[,MeanAg:=mean(finScore),by=miRNA][,MaxAg:=max(finScore),by=miRNA]
    final_tablen <- final_tablen[finScore==MaxAg][,"MaxAg":=NULL]
    final_tablen[,c("Score","MeanAg","finScore"):=.(round(Score,3),round(MeanAg,2),round(finScore,3))]
    # setcolorder(final_tablep,c("miRNA","Group","Score","MeanAg","List","Combo"))
    setkeyv(final_tablen,setdiff(names(final_tablen),"List"))
    final_tablen <- unique(final_tablen)
  } else elaps.timen <- system.time({NULL})
  if (exists("final_tablep") & exists("final_tablen")){
    final_table <- merge(final_tablep,final_tablen,by="miRNA")
    final_table[,c("finScore.x","finScore.y"):=NULL][,FinScore:=Score.x+Score.y]
    setnames(final_table,c("miRNA","PosCombo","PosScore","PosGroup","PosList","PosAg",
                           "NegCombo","NegScore","NegGroup","NegList","NegAg","FinScore"))
    setcolorder(final_table,c(1,12,2:11))
    setorder(final_table,-PosScore)
  } else if(exists("final_tablep")){
    if(negman){
      setcolorder(final_tablep,c(1,9,2,3,4,10,5,8,7,6))
      setnames(final_tablep,c("miRNA","FinScore","PosCombo","PosScore","PosGroup",
                              "PosAg","PosList","NegCombo","NegScore","NegGroup"))
      setorder(final_tablep,-PosScore)
      final_table <- final_tablep
    } else {
      final_table <- final_tablep[,finScore:=NULL]
      setcolorder(final_table,c(2,6,1,3,4,5))
      setnames(final_table,c("miRNA","FinScore","PosCombo","PosScore","PosGroup","PosList"))
    }
  } else {
    if(posman){
      setcolorder(final_tablen,c(1,9,2,3,4,10,5,8,7,6))
      setnames(final_tablen,c("miRNA","FinScore","NegCombo","NegScore","NegGroup",
                              "NegAg","NegList","PosCombo","PosScore","PosGroup"))
      setorder(final_tablen,-NegScore)
      final_table <- final_tablen

    } else {
      final_table <- final_tablen[,finScore:=NULL]
      setcolorder(final_table,c(2,6,1,3,4,5))
      setnames(final_table,c("miRNA","FinScore","NegCombo","NegScore","NegGroup","NegList"))

    }
  }

  if (nrow(final_table)==0){
    print("No scores detected above threshold with current parameters.")
    print(paste0("File analyzed: ",inputfile))
    print(paste0("Genes found: ",paste(genelist.src,collapse=", ")))
    print(paste0("Number of Genes: ",length(genelist.src)))
    print(paste0("Minimum combination to consider: ",mincomb))
    print(paste0("Minimum score threshold: ",Zthresh))
#     print(paste0("Proportion method detected: ",(method[2]=="1")))
    return(NULL)
  }

  if (!is.null(genelist.kern)){
    print("Assigning proportions of defined core")
    props <- rep(0,nrow(final_table))
    for (i in 1:nrow(final_table)){
      fincomb <- character()
      if("PosCombo" %in% names(final_table)) {
        tmp <- unlist(strsplit(final_table[i,PosCombo],split=","))
        fincomb <- c(fincomb,tmp)
      }
      if("NegCombo" %in% names(final_table)) {
        tmp <- unlist(strsplit(final_table[i,NegCombo],split=","))
        fincomb <- c(fincomb,tmp)
      }
      fincomb <- unique(fincomb[fincomb!=""])
      props[i] <- length(intersect(fincomb,genelist.kern))/length(genelist.kern)
    }
    final_table[,Prop:=round(props,3)]
    final_table <- final_table[props>=0.2]
  }
  if ("PosScore" %in% names(final_table)){
    setorder(final_table,-PosScore)
  } else {
    setorder(final_table,-NegScore)
  }
  #   setorder(final_table,-FinScore)

  print("Writing results")
  # trim directory information from input filename
  ifn <- unlist(strsplit(inputfile,split="/"))
  ifn <- ifn[[length(ifn)]]

  # write output file, tab delimited with header of run information
  if(!is.null(genelist.kern)){
    trackcore <- paste(genelist.kern,collapse=", ")
  } else {
    trackcore <- "none"
  }
  if(is.null(genelist.pos)){
    genelist.pos <- "none"
  }
  if(is.null(genelist.neg)){
    genelist.neg <- "none"
  }
  # if(!is.null(genelist.part)){
  #   trackpart <- paste(genelist.part,collapse=",")
  # } else {
  #   trackpart <- "none"
  # }
  elaps.time <- elaps.timen + elaps.timep
  heading <- rep("##",14)
  heading[1] <- "#### metaMIR output ####"
  heading[2] <- "# Run Parameters:"
  heading[3] <- paste0("# Max allowed cluster size: ",maxclust)
  heading[4] <- paste0("# Minimum analyzed cluster: ",mincomb)
  heading[c(5,10,12)] <- "# "
  heading[6] <- paste0("# Input list: ",paste(genelist.src,collapse=", "))
  heading[7] <- paste0("# Core of genes tracked: ",trackcore)
  heading[8] <- paste0("# Genes for positive interaction: ",paste(genelist.pos,collapse=", "))
  heading[9] <- paste0("# Genes for negative interaction: ",paste(genelist.neg,collapse=", "))
  if(elaps.time[3]>3600){
    elaps.time <- round(elaps.time[3]/3600,1)
    un <- "hr"
  } else if(elaps.time[3]>60){
    elaps.time <- round(elaps.time[3]/60,2)
    un <- "min"
  } else {
    elaps.time <- round(elaps.time[3],2)
    un <- "sec"
  }
  heading[11] <- paste0("# Elapsed time: ",elaps.time," ",un)
  heading[13] <- "# Analyzed gene lists:"

  if(genelist.pos[[1]]!="none" & !posman){
    glout1 <- rep("#",length(glistpos)+2)
    glout1[[1]] <- "# Positive list(s)"
    for (i in 1:length(glistpos)){
      glout1[i+1] <- paste0("# ",names(glistpos)[i],": ",paste(glistpos[[i]],collapse=", "))
    }
  } else if(genelist.pos[[1]]!="none" & posman){
    glout1 <- rep("#",2)
    glout1[[1]] <- "# Positive list(s)"
    glout1[[2]] <- paste0("# L01: ",paste(genelist.pos,collapse=", "))
  }
  if(genelist.neg[[1]]!="none" & !negman){
    glout2 <- rep("#",length(glistneg)+2)
    glout2[[1]] <- "# Negative list(s)"
    for (i in 1:length(glistneg)){
      glout2[i+1] <- paste0("# ",names(glistneg)[i],": ",paste(glistneg[[i]],collapse=", "))

    }
  } else if(genelist.neg[[1]] != "none" & negman){
    glout2 <- rep("#",2)
    glout2[[1]] <- "# Negative list(s)"
    glout2[[2]] <- paste0("# L01: ",paste(genelist.neg,collapse=", "))

  }
  if(exists("glout1")&exists("glout2")){
    glout <- c(glout1,glout2)
  } else if (exists("glout1") & !exists("glout2")){
    glout <- glout1
  } else glout <- glout2

  # setup output path
  if (is.na(outPath)) {
	  fname1 <- paste0("Results/",gsub(".txt","",ifn),"_params_",gsub(" ","-",gsub(":","_",Sys.time())),".txt")
	  fname2 <- gsub("_params_","_output_",fname1)
  } else {
	  fname1 <- paste(outPath,paste(ifn,"params.txt",sep="."),sep="/")
	  fname2 <- gsub("params","output",fname1)
  }
  write.table(c(heading,glout),file=fname1,quote=F,col.names=F,row.names=F)
  write.table(final_table,file=fname2,col.names=T,row.names=F,quote=F,sep="\t")
}
