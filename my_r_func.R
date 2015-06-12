# Version info: R 3.2.0, bioconductor 3.1
# Install packages
source("http://bioconductor.org/biocLite.R")
biocLite()
# GEOquery, limma 
biocLite("GEOquery")
biocLite("limma")
library(GEOquery)
library(limma)

ParseGSE <- function(f){
  # extract expression values, probes and patients information from GEO SOFT.GZ files.
  #
  # Arvgs:
  #  f: GES number (e.x."GSE31023")
  #
  # Returns:
  #  list which contains:
  #    "oexp": original expression values matrix (prob.ID vs. patients)
  #    "pat": patients information matrix (info.col vs. patients.ID)
  #    "prob": probes information matrix (prob.ID vs. info.col)
  
  gse <- getGEO(f,GSEMatrix=FALSE)
  gpls <- GPLList(gse)
  gsms <- GSMList(gse)

  # Extract values into matrix 
  exp <- do.call("cbind", lapply(gsms, function(x) as.numeric(Table(x)$VALUE)))  
  rownames(exp) <- Table(gsms[[1]])$ID
  
  # Extract patients information
  pat.info <- do.call("cbind", lapply(gsms, function(x) unlist(Meta(x))))
  # Usage: unlist() --> list to matrix (value vs. names)\
  
  # Extract probes information
  prob.info <- Table(gpls[[1]])
  
  return(list(oexp=exp, prob=prob.info,pat=pat.info))
}

GroupPat <- function(pat){
  # Divide patients into different groups (control vs treatment)
  #
  # Args:
  #  pat: matrix which contains patients information (ex. Meta(gsm))
  #
  # Returns:
  #  transposed original matrix with one more colum of "group" contains "G0" or "G1"
  #  matrix (patient.ID vs. information.categories)
  
  # Print all possible useful information to decide group.
  ii <- 0
  index <- vector()
  for(i in 1:nrow(pat)){
    f <- levels(factor(pat[i,]))  # summarize all items in one catalog of information.
    if(length(f)>1 & length(f) < ncol(pat)){
      ii <- ii+1
      index <- c(index,i)
      print(paste(ii,"-",f))
    }
  }
  
  # Choose and decide groups.
  iput <- as.numeric(readline(prompt="Which group contains infomation to devide patients to control or treatment:\t"))
  iput <- index[iput]
  cat("typing the name for each group\n")
  cat("'G0' for control, 'G1' for disease, leaving blank will be not used.:\n")
  decision <- as.character(levels(factor(pat[iput,])))  # all items in selected catalog of information 
  iput.name <- character()
  for(i in decision){
    iput.name <- c(iput.name,readline(prompt=paste(i,":\t")))
  }
  
  # Add colum contain grouping informaiton.
  pat <- rbind(pat,group = NA)
  for(i in 1:length(decision)){
    pat[nrow(pat),which(pat[iput,]==decision[i])] =iput.name[i]
  }
  
  return(t(pat))
  # Usage: t(matrix) : transpose matrix. 
}

PreProc <- function(exp){
  # normalizate microarray data with quantiles normalize
  # TODO: deal with NA values
  # 
  # Args:
  #  exp: expression value matrix (probs.ID vs. patients)
  #
  # Returns:
  #  not decided
  
  # Normalization
  print("normalization")
  exp <- normalizeQuantiles(exp)
  cat("Done!\n")

  # TODO: deal with NA values

  return(exp)
}

MatrixTest <- function(x,y,
                       manner=1){
  # perform t.test and wilcox.text to every row of matrix
  #
  # Args:
  #  if manner == 1 (default): x is matrix 1, and y is matrix 2 for test(x,y)
  #  if manner == 2 : x is total matrix, y is binary factor vector for test(x~y)
  #
  # Returns:
  #  matrix probes vs. categories of test result
  #    test result: t.test and wilcox.test's statstic, p-value and adjust p-value
  
  d1 <- matrix()
  d2 <- matrix()
  
  # determine the data matrix for test
  if (manner==1){
    d1 <- x
    d2 <- y
  }else if(manner ==2){
    if (length(y)==ncol(x)){
      sp <- split(1:ncol(x),y)  # sp=list, name is factor, and value is colum index number
      d1 <- x[,sp[[1]]]
      d2 <- x[,sp[[2]]]
    }else print("length of y is not match for colum of x")
  }else {print("wrong manner")}

  # test whether varance equal
  var.equ <- FALSE
  var.p <- var.test(d1,d2)$p.value
  if(var.p <= 0.05) var.equal <- TRUE
  
  # Perform test
  output <- matrix(ncol=4)
  
  cat("* calculating number of value\n\n")
  avl.d1 <- rowSums(!is.na(d1))
  avl.d2 <- rowSums(!is.na(d2))

  # Available value should be more than 3
  cat("* calculating t.test and wilcox.test\n")
  pd <- txtProgressBar(min=0,max=nrow(d1),style=3)
  for (i in 1:nrow(d1)){
    setTxtProgressBar(pd,i)
    if(avl.d1[i]<3 | avl.d2[i]<3){
      output <- rbind(output,c(NA,NA,NA,NA))
    }else{
      t <- t.test(d1[i,],d2[i,],var.equal=var.equ)
      w <- wilcox.test(d1[i,],d2[i,],var.equal=var.equ)
      output <- rbind(output,c(t$p.value,t$statistic,w$p.value,w$statistic))
    }
  }
  close(pd)
  cat("\n")

  # Add adjust p value
  cat("* Adjusting p value with FDR\n\n")
  output <- cbind(output[,1:2],fdr.t.p=p.adjust(output[,1],method="BH"),
                  output[,3:4],dfr.t.p=p.adjust(output[,3],method="BH"))

  # Add rowname and colname
  colnames(output) <- c("t.p-value","t.statistic","t.adjust-p","w.p-value","w.statistic","w.adjust-p")
  output <- output[-1,]
  rownames(output) <- rownames(d1)
  return(output)
}


SelectThres <- function(m,p.col,
                       thres=0.05){
  # Extract elements which p-values is under threshord
  #
  # Args:
  #  m: data matrix (elements vs. value categories)
  #  p.col: colum name which stand for the p-value.
  #  thres: threshord setted, usually 0.05 or 0.01 or 0.001 
  #
  # Returns:
  #  matrix which is subsection of input matrix, contains elements under thresord (elements vs. value categories)

  output<- m[which(m[,colname]<=thres),]
  return(output)
}

ProbsToGene <- function(lst){
  # select common chars among each list'element
  #
  # Args:
  #  lst: list with each element contains a list with indivdual exptression matrx and related probes infomation matrix
  #  lst[[i]]: one set of data, lst[[i]][[1]]: expression matrix in set, lst[[i]][[2]]: probes matrix
  #
  # Returns:
  #  list with each elements contains a expression matrix but with last column(colnames: "geneID") stand for gene identifier
  #  if input list contains names, output list will contains same names.

  # Choose category to stand for gene
  gene.ide <- character()
  probe.ide <- character()
  comm <- lapply(lst,function(x) colnames(x[[2]]))  # extract probe info categories
  comm <- Reduce(intersect,comm) # choose common categories
  print(comm)
  if(length(comm)>1){
    iput <- readline("choose common elements will be used follow as gene identifier. (input index number)\t")
    iput.2 <- readline("choose common elements stand for probes ID. (input index number)\t")
    gene.ide <- comm[as.numeric(iput)]
    probe.ide <- comm[as.numeric(iput.2)]
  }else{
    print("only one common coloums, something wrong")
  }

  # Add gene identifier to exp matrix
  output <- lapply(lst,function(x){
    data.frame(x[[1]],"geneID"=as.character(x[[2]][match(rownames(x[[1]]),x[[2]][,probe.ide]),gene.ide]))  # cbind cannot bind numeric and char together.
  })
  return(output)
}

o.list <- list("colon"=colon$exp$exp,"liver"=liver$exp$exp)
p.lst <- list(colon=list(colon$exp$exp,colon$exp$prob),liver=list(liver$exp$exp,liver$exp$prob))
o.list <- ProbsToGene(p.lst) 
a<-SelectSpec(o.list, "geneID")
print(class(a))
print(a[1:10,])
print(colnames(a))
print(nrow(a[which(a[,"num.cancer"]==1),]))
print(nrow(a[which(a[,"num.cancer"]==2),]))
print(nrow(a[which(a[,"num.cancer"]==3),]))
print(nrow(a[which(a[,"num.cancer"]==4),]))
print(nrow(a[which(a[,"num.cancer"]==5),]))

SelectSpec <- function(dats,identi){
  # Classify specific genes which exist in only one caner, or two cancers, three cancers....
  #
  # Args:
  #  dats: list of expression data.frame
  #    data.frames: probsID vs. patientsID, with last colums(colnames: "geneID") is geneID (e.x. GB_ACC)
  #  identi: coloum name which point out the coloums stand for the genes ID number.
  #         It's should be same with last coloums of data.frame in dats.
  #
  # Returns:
  #  list with "mat" contains matrix geneID vs. the gene in each matrix and total num matrix containing gene
  #            "sum" contains statistic information summary matrix 

  # Combine all elements into one vector
  identi.lst <- lapply(dats,function(x) x[,identi])
  identi.all <- factor()
  for (i in identi.lst) identi.all <- factor(c(as.character(identi.all),as.character(i)))
  # usage: link two factor: factor(c(as.chararcter(),as.character()))
  identi.all <- unique(identi.all)  # remove duplicated elem
  identi.all <- identi.all[!identi.all %in% ""]  # remove blank elem like ""

  # Match total identifier list to each matrix data
  mat <- do.call("cbind",lapply(dats,function(x) match(identi.all, x[,identi])))
  mat <- data.frame("geneID"=identi.all,"num.cancer"=rowSums(!is.na(mat)),mat)
   
  # summarize specific genes numbers and other infomation.
  summ <- data.frame()
  cancer.1 <- mat[which(mat[,"num.cancer"]==1),]
  for(i in 1:length(dats)){
    dat.name <- names(dats[i])
    gene.total <- nrow(dats[i])
    gene.1 <- nrow(cancer.1[which(!is.na(cancer.1[,dat.name])),])
    summ<- rbind(summ,c(dat.name,gene.total,gene.1))
  }
  print(summ)

  #return(output)
}

Main_Select <- function(obj.list){
  mat.list <- list() 
  # select genes with p<0.05
  for(i in 1:length(obj.list)){
    n <- names(obj.list[i])
    exp <- obj.list[[i]]$test
    mat.list[[n]] <- SelectGene(SelectGene(exp,"t.adjust-p",0.05),"w.adjust-p",0.05)
  }
 # add gene info 
  for(i in 1:length(obj.list)){
    mat.p <- obj.list[[i]]$exp$prob
    mat.e <- mat.list[[i]]
    mat.list[[i]] <-cbind(mat.e,geneID=as.character(mat.p[,gene.col]))
  }
  output <- CommAndSpec(mat.list,"geneID")
  print(nrow(output))
  print(nrow(output[which(output[,"cancer_num"]==1),]))
  print(nrow(output[which(output[,"cancer_num"]==2),]))
  print(nrow(output[which(output[,"cancer_num"]==3),]))
  print(nrow(output[which(output[,"cancer_num"]==4),]))
  print(nrow(output[which(output[,"cancer_num"]==5),]))

  t<- output[which(output[,"cancer_num"]==1),]
  print(nrow(t[which(!is.na(t[,"colon"])),]))
  print(nrow(t[which(!is.na(t[,"liver"])),]))
  print(nrow(t[which(!is.na(t[,"panc"])),]))
  print(nrow(t[which(!is.na(t[,"renal"])),]))
  print(nrow(t[which(!is.na(t[,"lung"])),]))
  print(nrow(t[which(!is.na(t[,"breast"])),]))

  return(output)
}
Step1 <- function(f){
  # integrated function to anaysis GSE data
  #
  # Args:
  #  f: GES number
  #
  # Returns:
  #  list which contains: 
  #     "oexp" for normalized expression data
  #     "pat" for patients information with divided for groups(control vs treat) 
  #     "prob" for probes information
  #     "test" for t.test and wilcox.test result

  obj <- ParseGSE(f)
  obj$pat <- GroupPat(obj$pat)  # Update obj.pat with adding group information
  obj$oexp <- PreProc(obj$oexp)  # Update ojb.exp with normalization.
  obj$test <- MatrixTest(obj$oexp,obj$pat[,"group"],manner=2)
  obj$sig <- SelectThres(SelectThres(obj$oexp,"t.adjust-p",thres=0.001),"w.adjust-p",thres=0.001)  # Select significant genes with t.test and wilcox test
  #o.list <- list("colon"=colon$exp$exp,"liver"=liver,"panc"=panc,"renal"=renal,"lung"=lung,"breast"=breast)
  o.list <- list("colon"=colon$exp$exp,"liver"=liver$exp$exp,"panc"=panc,"renal"=renal,"lung"=lung,"breast"=breast)
  p.lst <- list(colon=list(colon$exp$exp,colon$exp$prob),liver=list(liver$exp$exp,liver$exp$prob))
  vprob <- list(colnames(colon$prob),colnames(liver$prob),colnames(panc$prob),colnames(renal$prob),colnames(lung$prob),colnames(breast$prob))

  gene.ident <- CommChar(intersect,vprob)


  return(obj)
}

elected <- Main_Select(o.list)
