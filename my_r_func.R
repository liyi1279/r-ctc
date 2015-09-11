# Version info: R 3.2.0, bioconductor 3.1
# Install packages
source("http://bioconductor.org/biocLite.R")
biocLite()
# GEOquery, limma 
biocLite("GEOquery")
biocLite("limma")
library(GEOquery)
library(limma)
install.packages("rSelRF")
install.packages("varSelRF")

ParseGSE <- function(f){
  # extract expression values, probes and patients information from GEO SOFT.GZ files.
  #
  # Arvgs:
  #  f: GES number (e.x."GSE31023")
  #
  # Returns:
  #  list which contains:
  #    "exp": original expression values matrix (prob.ID vs. patients)
  #    "pat": patients information matrix (info.col vs. patients.ID)
  #    "prob": probes information matrix (prob.ID vs. info.col)
  
  gse <- getGEO(f,GSEMatrix=FALSE)
  #gse <- getGEO(file=f)
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
  
  return(list(exp=exp, prob=prob.info,pat=pat.info))
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

PreProc <- function(exp,prob){
  # normalizate microarray data with quantiles normalize
  # TODO: deal with NA values
  # 
  # Args:
  #  exp: expression value matrix (probs.ID vs. patients)
  #
  # Returns:
  #  same matrix but with value is normalized. 
  
  # Mapping Probe with EntrezGeneID
  idents <- colnames(prob)
  for (i in 1:length(idents)){
      idents[i] <- paste(i,"-",idents[i])
  }
  cat("\n")
  print(idents)
  cat("\n")
  iput <- readline("choose which coloum stands for gene ID. (input index number)\n")
  identifier <- colnames(prob)[as.numeric(iput)]
  print(identifier)
  geneList <- prob[,identifier]
  print(head(geneList))
  # Integrate rows with same ID.
  exp <- aggregate(exp,by=list(EntrezID=geneList),FUN=mean)
  rownames(exp) <- exp[,1]
  exp <- exp[-1,-1]                    # Discard blank row, and ID col. 

  # TODO: deal with NA values
  cat(paste("sum(is.na(exp)= ",sum(is.na(exp))))
  cat("\n")
  cat(paste("nrow(exp)= ",nrow(exp)))
  cat("\n")
  exp <- na.omit(exp)
  cat(paste("nrow(na.omit(exp))= ",nrow(exp)))


  # Normalization
  print("normalization")
  exp <- normalizeQuantiles(exp)
  cat("Done!\n")

  return(exp)
}

MatrixTest <- function(x,y, manner=1){
  # perform t.test and wilcox.text to every row of matrix
  #
  # Args:
  #  if manner == 1 (default): x is matrix 1, and y is matrix 2 for test(x,y)
  #  if manner == 2 : x is total matrix, y is binary factor vector for test(x~y)
  #
  # Returns:
  #  matrix probes vs. categories of test result
  #    test result: t.test and wilcox.test's statstic, p-value and adjust p-value
  
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


  # shapiro.test & var.test
  print("Confirming distribution parameter and varence equance")
  pTotal <- cbind(d1,d2)
  d1d2 <- ncol(d1)
  dT <- ncol(d1)+ncol(d2)
  ShaVar <<- do.call("cbind",apply(pTotal,1,function(x){
                         list(shapiro.test(x[1:d1d2])$p.value,
                         shapiro.test(x[(d1d2+1):dT])$p.value,
                         var.test(x[1:d1d2],x[(d1d2+1):dT])$p.value)
                    }))
  ShaVar <<- t(ShaVar)
  colnames(ShaVar) <- c("d1.sha.p","d2.sha.p","var.p")

  # Perform test
  output <- matrix(ncol=3)
  
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
      testMethod=c()
      var.equ=FALSE
      if (ShaVar[i,3]<0.05) var.equ=TRUE
      if (ShaVar[i,1]>0.05 & ShaVar[i,2]>0.05){
        t <- t.test(as.numeric(d1[i,]),as.numeric(d2[i,]),var.equal=var.equ)
        testMethod=1
      }else{
        t <- wilcox.test(as.numeric(d1[i,]),as.numeric(d2[i,]),var.equal=var.equ)
        testMethod=0
      }
      output <- rbind(output,c(t$p.value,t$statistic,testMethod))
    }
  }
  close(pd)
  cat("\n")

  # Add adjust p value
  cat("* Adjusting p value with FDR\n\n")
  output <- cbind(output,fdr.t.p=p.adjust(output[,1],method="BH"))

  # Add rowname and colname
  colnames(output) <- c("P.value","Statistic","1-ttest/0-wilcox","FDR.adjust.p")
  output <- output[-1,]
  rownames(output) <- rownames(d1)

  # Calculate FoldChange
  mean.d1 <- apply(d1,1, function(x) mean(x))
  mean.d2 <- apply(d2,1, function(x) mean(x))
  fc <- foldchange(mean.d1,mean.d2)
  output <- cbind(output,"FoldChange"=fc)
  return(output)
}

SelectThres <- function(m,p.col, thres=0.05){
  # Extract elements which p-values is under threshord
  #
  # Args:
  #  m: data matrix (elements vs. value categories)
  #  p.col: colum name which stand for the p-value.
  #  thres: threshord setted, usually 0.05 or 0.01 or 0.001 
  #
  # Returns:
  #  matrix which is subsection of input matrix, contains elements under thresord (elements vs. value categories)

  output<- m[which(m[,p.col]<=thres),]
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
  comm <- lapply(lst,function(x) colpames(x[[2]]))  # extract probe info categories
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
   
  # summarize specific genes numbers for each cancers.
  spec.gene <- mat[which(mat[,"num.cancer"]==1),]  # all genes which expresed in only one cancer
  for(i in 1:length(dats)){
    dat.name <- names(dats[i])  # cancer name
    gene.total <- nrow(dats[[i]])  # total gene num
    gene.1 <- nrow(spec.gene[which(!is.na(spec.gene[,dat.name])),])  # specific gene num
    rat <- round(gene.1 / gene.total,4)*100
    s <- c(dat.name,gene.total,gene.1,rat)
    ifelse(!exists("summ"),summ<-s, summ <- rbind(summ,s))
  }
  rownames(summ) <- summ[,1]
  summ <- summ[,-1]
  colnames(summ) <- c("total","specific gene","percent(%)")

  # Add common gene numbers for 1 cancers, 2 cancers, 3 cancers,...all cancers. 
  gene.num <- vector()
  for (i in 1:length(dats)) gene.num <- c(gene.num,nrow(mat[which(mat[,"num.cancer"]==i),]))
  t <- cbind(c(1:length(dats)),gene.num,round(gene.num/length(identi.all),4)*100)
  rownames(t) <- c(1:length(dats))
  colnames(t) <- c("# cancer", "com.gene.num","%.in.all")

  return(list(mat=mat,summ=summ,cancsum=t))
}

S2MultiCanc <- function(dat.lst){
  # Extract specific significant genes for each cancer.
  #
  # Args:
  #  dat.lst: list of objects returned from S1Parse
  #  list format: list(cancer.name=cancer.obj, cancer.name=cancer.obj)
  #
  # Returns:
  #  list with "mat" contains matrix geneID vs. the gene in each matrix and total num matrix containing gene
  #            "sum" contains statistic information summary matrix 
 
  cat("\n1. add gene ID to expression matrix\n")
  Sys.sleep(1.5)
  temp.lst <- lapply(dat.lst,function(x) list("expression"=x$sig,"probes"=x$prob))
  temp.geneid <- ProbsToGene(temp.lst)
  cat("\n2. Select specific genes in signifiant matrix\n")
  Sys.sleep(1.5)
  output <- SelectSpec(temp.geneid,"geneID")
  return(output)
 }

MainPrimary <- function(){
  # input all constants and use functions porperaly, MainPrimary() will do all the thing.
  #
  # Args:
  # 
  #
  # Returns:
  #  no returns, but this funciton will create global varabile to contains MA data and last result.
  #  global vars: "breast", "lung", "colon", "panc", "renal", "liver"
  #  results vars: inputed "v"

  cat("\n1. Parse breast cancer\n")
  Sys.sleep(1.5)
  breast <<- S1Parse("GSE15852")
  cat("\n1. Parse lung cancer\n")
  Sys.sleep(1.5)
  lung <<- S1Parse("GSE10072")
  cat("\n1. Parse colon cancer\n")
  Sys.sleep(1.5)
  colon <<- S1Parse("GSE6988")
  cat("\n1. Parse panc cancer\n")
  Sys.sleep(1.5)
  panc <<- S1Parse("GSE15471")
  cat("\n1. Parse renal cancer\n")
  Sys.sleep(1.5)
  renal <<- S1Parse("GSE15641")
  cat("\n1. Parse liver cancer\n")
  Sys.sleep(1.5)
  liver <<- S1Parse("GSE25097")

  #canc.lst <- list("renal"=renal, "panc"=panc)
  canc.lst <- list("breast"=breast, "lung"=lung, "colon"=colon, "renal"=renal, "panc"=panc, "liver"=liver)
  
  PriMarkers <<- S2MultiCanc(canc.lst)
}

MainPrimary <- function(){
  # input all constants and use functions porperaly, MainPrimary() will do all the thing.
  #
  # Args:
  # 
  #
  # Returns:
  #  no returns, but this funciton will create global varabile to contains MA data and last result.
  #  global vars: "breast", "lung", "colon", "panc", "renal", "liver"
  #  results vars: inputed "v"

  cat("\n1. Parse breast cancer\n")
  Sys.sleep(1.5)
  breast <<- S1Parse("GSE15852")
  cat("\n1. Parse lung cancer\n")
  Sys.sleep(1.5)
  lung <<- S1Parse("GSE10072")
  cat("\n1. Parse colon cancer\n")
  Sys.sleep(1.5)
  colon <<- S1Parse("GSE6988")
  cat("\n1. Parse panc cancer\n")
  Sys.sleep(1.5)
  panc <<- S1Parse("GSE15471")
  cat("\n1. Parse renal cancer\n")
  Sys.sleep(1.5)
  renal <<- S1Parse("GSE15641")
  cat("\n1. Parse liver cancer\n")
  Sys.sleep(1.5)
  liver <<- S1Parse("GSE25097")

  #canc.lst <- list("renal"=renal, "panc"=panc)
  canc.lst <- list("breast"=breast, "lung"=lung, "colon"=colon, "renal"=renal, "panc"=panc, "liver"=liver)
  
  PriMarkers <<- S2MultiCanc(canc.lst)
}

########## Code Sending to R ##########
normal <- function(){
    load("~/GSEDataSets.RData")
    a <- ParseGSE("~/r-ctc/GSE21986_family.soft.gz")
    aa <- GroupPat(a$pat)
    
    a$prob[100:110,1:5]
    a$prob[a$prob[,'GENE_SYMBOL']=='RGR',1:5]
    a$prob[a$prob[,'GENE_SYMBOL']=='Rgr',1:5]
    a$pat <- aa
    aaa <- PreProc(a$exp,a$prob)           # GB _ACC as gene ID 
    head(aaa)
    head(a$pat)
    colnames(aaa) <- a$pat[,ncol(a$pat)]
    dim(aaa)
    a$expQnMap <- aaa
    head(a$expQnMap)
    GSE21986 <- list(expQnMap=a$expQnMap)
    
    head(aaa)
    head(a$expQnMap)
    sum(is.na(aaa))
    
    b <- ParseGSE("~/r-ctc/GSE28303_family.soft.gz")
    bb <- GroupPat(b$pat)
    b$pat <- bb
    bbb <- PreProc(b$exp,b$prob)           # GB _ACC as gene ID 
    colnames(bbb) <- b$pat[,ncol(b$pat)]
    b$expQnMap <- bbb
    GSE28303 <- b
    
    c <- ParseGSE("~/r-ctc/GSE28303_family.soft.gz")
    cc <- GroupPat(c$pat)
    c$pat <- cc
    ccc <- PreProc(c$exp,c$prob)           # GB _ACC as gene ID 
    colnames(ccc) <- c$pat[,ncol(c$pat)]
    c$expQnMap <- ccc
    GSE28303 <- c
}
##### change table style #####
Change <- function(){
    g1 <- GSE14671
    head(g1$expQnMap)
    s <- rbind(response=colnames(g1$expQnMap),label=c(rep(0,12),rep(1,24),rep(0,6),rep(1,17)),g1$expQnMap)
    s <- s[-3,]
    rownames(g1$data2)
    m <- match(rownames(g1$data2),g1$prob[,'ENTREZ_GENE_ID'])
    gene <- g1$prob[m,'Gene Symbol']
    ss <- cbind(geneSymbol=gene,s)
    g1$data2 <- ss
    head(g1$data2)
    GSE14671 <- g1
    head(GSE14671$data2)
    
    g2 <- GSE2535
    head(g2$expQnMap)
    s <- rbind(response=colnames(g2$expQnMap),label=c(rep(0,5),rep(1,8),rep(0,7),rep(1,8)),g2$expQnMap)  
    head(s)
    s <- s[-3,]
    head(s)
    m <- match(rownames(s),g2$prob[,'ENTREZ_GENE_ID'])
    gene <- g2$prob[m,'Gene Symbol']
    ss <- cbind(geneSymbol=gene,s)
    g2$data2 <- ss
    head(g2$data2)
    GSE2535 <- g2
    head(GSE2535$data2)
     
    g3 <- GSE33224
    head(g3$expQnMap)
    s <- rbind(response=colnames(g3$expQnMap),label=c(rep(0,8),rep(1,6),rep(0,8),rep(1,6)),g3$expQnMap)  
    head(s)
    s <- s[-3,]
    head(s)
    m <- match(rownames(s),g3$prob[,'GENE'])
    g3$prob[100:110,1:10]
    gene <- g3$prob[m,'GENE_SYMBOL']
    ss <- cbind(geneSymbol=gene,s)
    g3$data2 <- ss
    head(g3$data2)
    GSE33224 <- g3
    head(GSE33224$data2)
 
    f1 <- e                            # PareseGSE obj 
    f2 <- ee                           # pat 
    f3 <- eee                          # geneid to rownames 
    head(f3)
    print(f2)
    f4 <- rbind(response=f2[,'group'],label=c(rep(1,5),rep(0,5)),f3)  
    head(f4)

    colnames(f1$prob)
    f1$prob[100:110,1:10]
    f5 <- match(rownames(f4),f1$prob[,'ENTREZ_GENE_ID']) # add ENTREZ_GENE_ID 
    print(f5[1:50])                    # should be started with NA, NA, ..., ... 
    f6 <- f1$prob[f5,'Gene Symbol']
    print(f6[1:50])                    # should be started with NA, NA, ..., ... 

    f7 <- cbind(geneSymbol=f6,f4)
    head(f7)

    GSE7114 <- list(expQnMap=f3, data2=f7) 
    GSE21986 <- list(expQnMap=f3,data2=f4) 

    head(GSE21986$expQnMap)
    head(GSE7114$expQnMap)
    head(GSE21986$data2)
    head(GSE7114$data2)

    save(GSE33224,GSE2535,GSE14671,file="dataForQC.RData")
# new data
    save(GSE7114,GSE21986,file="dataForQC_2.RData")
    
}

##### GSE62121 #####
NewData <- function(){
    d <- read.csv(file="GSE62121.csv")
    dd <- d[,c(2,6:ncol(d))]
    colnames(dd) <- c('ID','R','NR-IMA','NR-NILO',rep('GNF',5))
    ddn <- as.factor(dd[,1])
    ddd <- aggregate(dd[,2:ncol(dd)],by=list(GENEID=ddn),FUN=mean)
    rownames(ddd) <- ddd[,1]
    ddd <- ddd[,-1]                    # Discard blank row, and ID col. 
    
    # TODO: deal with NA values
    exp <- na.omit(exp)
    
    # Normalization
    print("normalization")
    exp <- normalizeQuantiles(exp)
    cat("Done!\n")
    
    
    length(ddn)
    length(unique(ddn))
    head(d)
    head(dd)
    head(ddd)
    
    #   exp <- aggregate(exp,by=list(EntrezID=geneList),FUN=mean)

    e <- ParseGSE("GSE7114")
    ee <- GroupPat(e$pat)
    eee <- PreProc(e$exp,e$prob)
    colnames(eee) <- ee[,'group']
    head(eee)
    
    be <- GSE7114
    bg <- GSE21986
    head(GSE21986$data2)
    s <- rownames(GSE21986$data2)
    head(s)
    rownames(GSE21986$data2)
    s <- s[-3]
    head(s)
    length(s)
    nrow(GSE21986$data2)
    head(GSE21986$data2)

    GSE21986$data2 <- cbind(s,GSE21986$data2)


# Read from raw.tar
    untar("GSE28303_RAW.tar",exdir="28303")
    cels <- list.files("28303/",pattern="CEL")
    sapply(paste("7117",cels,sep="/"),gunzip)
    cels <- list.files("7117/",pattern="CEL")
    biocLite("affy")
    biocLite("simpleaffy")
    library(simpleaffy)
    raw.data=ReadAffy(verbose=TRUE,filenames='7117/GSM162954.CEL')

} 

##### modifie two-channel #####
two <- function(){
h <- GSE33224
h1 <- h$pat
h2 <- h$expQnMap
head(h2)
colnames(h2) <- c(11,12,21,22,31,32,41,42,51,52,61,62,71,72,13,14,23,24,33,34,43,44,53,54,63,64,73,74)
h2[1:30,c('11','12','13','14','51','52','53','54')]
h2[1:30,c('51','52','53','54')]
h3 <- h2
colnames(h3) <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,1,1,2,2,3,3,4,4,5,5,6,6,7,7)
head(h3)
h33 <- t(h3)
head(t(h33))
h333 <- aggregate(h33,by=list(rownames(h33)),FUN=mean)
head(t(h333[,-1]))
h3333 <- t(h333[,-1])
head(h3333)
colnames(h3333) <- c(rep('NR',4),rep('R',3))

h4 <- h2
colnames(h4) <- c(11,11,21,21,31,31,41,41,51,51,61,61,71,71,12,12,22,22,32,32,42,42,52,52,62,62,72,72)
head(h4)
h44 <- t(h4)
head(t(h44))
h444 <- aggregate(h44,by=list(rownames(h44)),FUN=mean)
head(t(h444))
head(t(h444[,-1]))
h4444 <- t(h444[,-1])
head(h4444)
colnames(h4444) <- c(rep('NR',8),rep('R',6))
GSE33224$data7 <- h3333
GSE33224$data14 <- h4444
h5 <- rbind(response=colnames(GSE33224$data7),GSE33224$data7)
as.factor('nr','r')
h5 <- rbind(factor(nr,nr,nr,nr,r,r,r),GSE33224$data7)
head(GSE33224$data7)
head(GSE33224$data14)
head(GSE33224$data2)
head(h5)
}
code <- function(){
save(GSE33224,file="GSE33224_7Sample_14Sample.RData")


write.csv(GSE14671$expQnMap,file="GSE33224_expQnMap.csv")

write.csv(GSE2535$data2,file="GSE2535_data2.csv")
write.csv(GSE2535$expQnMap,file="GSE2535_expQnMap.csv")

write.csv(GSE14671$data2,file="GSE14671_data2.csv")
write.csv(GSE14671$expQnMap,file="GSE14671_expQnMap.csv")

write.csv(GSE7114$data2,file="GSE7114_data2.csv")
write.csv(GSE7114$expQnMap,file="GSE7114_expQnMap.csv")

write.csv(GSE21986$data2,file="GSE21986_data2.csv")
write.csv(GSE21986$expQnMap,file="GSE21986_expQnMap.csv")

write.csv(GSE33224$data2,file="GSE33224_data2_28sample.csv")
write.csv(GSE33224$data7,file="GSE33224_data2_7sample.csv")
write.csv(GSE33224$data14,file="GSE33224_data2_14sample.csv")
write.csv(GSE21986$expQnMap,file="GSE21986_expQnMap.csv")

# 100130541
GSE33224$prob[100:110,1:5]
nrow(GSE33224$prob)
head(GSE33224$data7 )
GSE33224$prob[GSE33224$prob[,'GENE_SYMBOL']=='MBNL3',1:10]
GSE33224$data2[rownames(GSE33224$data2)==100130541,1:5]

##### Meta-analysis object #####
load("MetaAnalysis-example.RData")
install.packages("MetaDE_1.0.5.tar.gz",repos=NULL,type="source")
biocLite('impute')
install.packages('MetaDE')
install.packages('MetaQC')
library(MetaDE)
library(MetaQC)

MetaDE.Res1$raw.data)
names(MetaDE.Res1)
class(MetaDE.Res1$ind.stat)
class(MetaDE.Res1$ind.p)
class(MetaDE.Res1$raw.data)
head(MetaDE.Res1$ind.p)
class(MetaDE.Res1$meta.analysis)
names(MetaDE.Res1$meta.analysis)
head(MetaDE.Res1$meta.analysis$pval)
head(MetaDE.Res1$meta.analysis$FDR)

head(MetaDE.Res1$meta.analysis$AW)

count.DEnumber(MetaDE.Res1, p.cut=c(0.001,0.005),q.cut=c(0.01,0.05))
count.DEnumber
 <- MetaDE.Res1
ss <- s$meta.analysis
sss <- ss[ss[,'Fisher']>0.001,]

head(s$ind.p)
head(ss$pval)
class(s)
table.p <- matrix(NA, 2, 1005)
table.p[,1:10]

pm <- cbind(result$ind.p, result$meta.analysis$pval)
qm <- cbind(apply(pm, 2, function(x) p.adjust(x, method = "BH")))
}

commonParse <- function(f){
  # integrated function to anaysis GSE data
  #
  # Args:
  #  f: GES number
  #
  # Returns:
  #  list which contains: 
  #     "exp" for normalized expression data
  #     "pat" for patients information with divided for groups(control vs treat) 
  #     "prob" for probes information
  #     "test" for t.test and wilcox.test result
  cat("\n1. extract information from GSE\n")
  Sys.sleep(1.5)
  f <- idCTCBreast
  ctcBreast <- ParseGSE(f)
  obj <- ctcBreast
  obj$pat <- GroupPat(obj$pat)  # Update obj.pat with adding group information
  obj$exp_pre <- PreProc(obj$exp,obj$pro)  # Use ENTREZ_ID to map, (or use Gene_symbol)
  #obj$exp_pre <- PreProc(obj$exp,tmp.prob)  # Use ENTREZ_ID to map, (or use Gene_symbol)
  #obj$test <- MatrixTest(obj$exp_pre,obj$pat[,"group"],manner=2)  # "group" contain only 0 and 1
  d1Index <- which(obj$pat[,"group"]=='0')
  d2Index <- which(obj$pat[,"group"]=='1')
  d1 <- obj$exp_pre[,d1Index] 
  d2 <- obj$exp_pre[,d2Index] 
  obj$test <- MatrixTest(d1,d2)
  obj$sig <- SelectThres(obj$test,"FDR.adjust.p") # default a=0.05
  #obj$sig <- SelectThres(obj$test,"P.value") # default a=0.05
  ctcBreast <- obj
}

sMatrixSelect <- function(m,cname=c(),v=c(),notNA=FALSE,matchRname=FALSE){
# Return Rows which contains elements in vectors(v) in column(cname) of matrix(m)
    seq.2 <- v
    if(matchRname==TRUE){seq.1 <- rownames(m)}else{seq.1 <- m[,cname]}
    if(notNA==FALSE){rule1 <- seq.1 %in% seq.2}else{rule1 <- TRUE}
    rule2 <- !is.na(seq.1)
    oput <- m[rule1 & rule2,]
    return(oput)
}



##### End of defination of functions #####


##### 0911 known marker list #####
kKnownMarkers <- list()
for(i in colnames(tmp.a)){
    mark <- c()
    aaa <- sMatrixSelect(tmp.a,i,notNA=TRUE)
    mark <- as.character(aaa[,1])
    kKnownMarkers[[i]] <- mark
}

seq1 <- rownames(ctcColon$sigSpec.colon)
seq2 <- kKnownMarkers$colon
match(seq1,seq2)
#  [1] NA NA NA NA NA NA NA NA NA NA NA NA NA NA
sMatrixSelect(ctcColon$test,seq2,matchRname=TRUE)
#      P.value Statistic 1-ttest/0-wilcox FDR.adjust.p
sMatrixSelect(colon$test,seq2,matchRname=TRUE)
#      P.value Statistic 1-ttest/0-wilcox FDR.adjust.p

seq3 <- rownames(ctcBreast$sigSpec.breast)
seq4 <- kKnownMarkers$breast
match(seq3,seq4)
# [1] NA NA NA NA NA
sMatrixSelect(ctcBreast$test,seq4,matchRname=TRUE)



##### 0905 statistics info ######
obj <- panc
ctcObj <- ctcPanc
# tissue DEGs
dim(obj$sig)
# [1] 14770     4
# [1] 2582    4
# [1] 3855    4

# Speci in One tissue
dim(obj$sigSpec) 
# [1] 2730    8
# [1] 525   8
# [1] 1444    8

# One UP
sum(obj$sigSpec[,2]>0) 
# [1] 2282
# [1] 282
# 687

# One DOWN
sum(obj$sigSpec[,2]<0)
# [1] 448
# [1] 243
# [1] 757

# one and ctc
sum(!is.na(match(rownames(ctcObj$sig),rownames(obj$sigSpec))))
# [1] 379
# [1] 13
# [1] 28

# one and ctc and same direction
dim(ctcObj$sigSpec.panc)
# [1] 203   6
# [1] 5 6
# [1] 14  7

# same UP
sum(ctcObj$sigSpec.panc[,2]>0)
# [1] 176
# [1] 2
# [1] 7

# same DOWN
sum(ctcObj$sigSpec.panc[,2]<0)
# [1] 27
# [1] 3
# [1] 7

tmp.a <- na.omit(match(rownames(ctcObj$sig),rownames(obj$sigSpec)))
match(rownames(ctcObj$sigSpec.panc),rownames(obj$sigSpec[tmp.a,]))
tmp.b<-is.na(match(rownames(obj$sigSpec[tmp.a,]),rownames(ctcObj$sigSpec.panc)))
# diff - UP
sum(obj$sigSpec[tmp.a,][tmp.b,2]>0)
# [1] 154
# [1] 2
# [1] 5

# diff - DOWN
sum(obj$sigSpec[tmp.a,][tmp.b,2]<0)
# [1] 22
# [1] 6
# [1] 9

# ctc DEGs
dim(ctcObj$sig) 
# [1] 1769    4
# [1] 736   5
# [1] 387   4

# ctc UP
sum(ctcObj$sig[,2]>0)  
# [1] 946
# [1] 487
# [1] 172

# ctc DOWN
sum(ctcObj$sig[,2]<0)
# [1] 727
# [1] 169
# [1] 163

sum(ctcObj$sig[,2]==0)                  # ctc DOWN 

##### 0824 cancer primary marker in CTC #####
idCTCBreast <- "GSE55470"
tmp.m <- match(rownames(ctcBreast$exp),ctcBreast$prob[,'ID'])
tmp.prob <- ctcBreast$prob[m,]
tmp.d<- apply(d1,2,function(x) as.numeric(x))
tmp.a <- c()
for(i in 1:nrow(d1)){
    print(i)
   tmp.a <- c(tmp.a,shapiro.test(as.numeric(d2[i,]))$p.value)
}
print(tmp.a)
dim(ctcBreast$exp)      # Original:         29377, 67
dim(ctcBreast$exp_pre)  # No Duplicate/NA:  20818, 67
dim(ctcBreast$sig)      # Significant:      736, 5
sum(ctcBreast$pat[,'group']=='0')   # Normal Sample: 6
sum(ctcBreast$pat[,'group']=='1')   # Cancer Sample: 7

match(rownames(breast$sigSpec),rownames(ctcBreast$test))
tmp.b <- ctcBreast$test[match(rownames(breast$sigSpec),rownames(ctcBreast$test)),c('FDR.adjust.p','Statistic')]
tmp.c <- cbind(breast$sigSpec[,c(1,2,3,4)],tmp.b)
tmp.d <- na.omit(tmp.c)
tmp.e <- tmp.d[-which(tmp.d[,5]<=0.05&tmp.d[,2]*tmp.d[,6]>0),] # Del same DEGs in array 2 
tmp.f <- tmp.d[which(tmp.d[,5]<=0.05&tmp.d[,2]*tmp.d[,6]>0),] # Del same DEGs in array 2 
# dim(tmp.e): 517 6     517 is no DEG in ctcBreast
# dim(tmp.d): 522 6     522 tested in ctcBreast
# dim(tmp.f): 5 6        203 DEGs same in panc specific DEGs
ctcBreast$sigSpec.breast <- tmp.f


idCTCPanc <- "GSE18670"
tmp.m <- match(rownames(ctcPanc$exp),ctcPanc$prob[,'ID'])
tmp.prob <- ctcPanc$prob[m,]
tmp.d<- apply(d1,2,function(x) as.numeric(x))
tmp.a <- c()
for(i in 1:nrow(d1)){
    print(i)
   tmp.a <- c(tmp.a,shapiro.test(d2[i,])$p.value)
}
print(tmp.a)
dim(ctcPanc$exp)      # Original:         8152, 9
dim(ctcPanc$exp_pre)  # No Duplicate/NA:  5625, 9
dim(ctcPanc$sig)      # Significant:      1769, 4
sum(ctcPanc$pat[,'group']=='0')   # Normal Sample: 6
sum(ctcPanc$pat[,'group']=='1')   # Cancer Sample: 6

match(rownames(panc$sigSpec),rownames(ctcPanc$test))
tmp.b <- ctcPanc$test[match(rownames(panc$sigSpec),rownames(ctcPanc$test)),c('FDR.adjust.p','Statistic')]
tmp.c <- cbind(panc$sigSpec[,c(1,2,3,4)],tmp.b)
tmp.d <- na.omit(tmp.c)
tmp.e <- tmp.d[-which(tmp.d[,5]<=0.05&tmp.d[,2]*tmp.d[,6]>0),] # Del same DEGs in array 2 
tmp.f <- tmp.d[which(tmp.d[,5]<=0.05&tmp.d[,2]*tmp.d[,6]>0),] # Del same DEGs in array 2 
# dim(tmp.e): 947 6     947 is no DEG in ctcPanc
# dim(tmp.d): 1150 6     1150 tested in ctcColon
# dim(tmp.f): 203 6        203 DEGs same in panc specific DEGs
ctcPanc$sigSpec.panc <- tmp.f

idCTCColon <- "GSE31023"
ctcColon$prob <- ctcColon$prob[-34184,]
dim(ctcColon$exp)      # Original:         34183, 9
dim(ctcColon$exp_pre)  # No Duplicate/NA:  3868, 9
dim(ctcColon$sig)      # Significant:      387, 4
sum(ctcColon$pat[,'group']=='0')   # Normal Sample: 3
sum(ctcColon$pat[,'group']=='1')   # Cancer Sample: 6

match(rownames(colon$sigSpec),rownames(ctcColon$test))
tmp.b <- ctcColon$test[match(rownames(colon$sigSpec),rownames(ctcColon$test)),c('P.value','Statistic')]
tmp.c <- cbind(colon$sigSpec[,c(1,2,3,4)],tmp.b)
tmp.d <- na.omit(tmp.c)
tmp.e <- tmp.d[-which(tmp.d[,5]<=0.05&tmp.d[,2]*tmp.d[,6]>0),] # Del same DEGs in array 2 
tmp.f <- tmp.d[which(tmp.d[,5]<=0.05&tmp.d[,2]*tmp.d[,6]>0),] # Del same DEGs in array 2 
# dim(tmp.e): 231 6     231 is no DEG in ctcColon
# dim(tmp.d): 245 6     245 tested in ctcColon
# dim(tmp.f): 14 6        14 DEGs same in colon specific DEGs
ctcColon$sigSpec.colon <- tmp.f

# Writing
write.csv(ctcColon$sigSpec.colon,file="ctcColon_oput.csv")
write.csv(ctcPanc$sigSpec.panc,file="ctcPanc_oput.csv")

##### 0806 cancer primary marker #####
idColon <- "GSE6988"
dim(colon$exp)      # Original:         117664, 123
dim(colon$exp_pre)  # No Duplicate/NA:  8475, 123
dim(colon$sig)      # Significant:      3855, 4
sum(colon$pat[,'group']=='0')   # Normal Sample: 28
sum(colon$pat[,'group']=='1')   # Cancer Sample: 53
# Random Forest
G0 <- rownames(colon$pat[colon$pat[,'group']=='0',])
G1 <- rownames(colon$pat[colon$pat[,'group']=='1',])
x <- colon$exp_pre[match(rownames(colon$sig),rownames(colon$exp_pre)),] 
x <- x[,c(match(G1,colnames(x)),match(G0,colnames(x)))]
x <- t(x)
cl <- factor(colon$pat[match(rownames(x),rownames(colon$pat)),'group'])
# SelVarRF
SelRf <- colon$selrf <- varSelRF(x,cl,ntree=10001)
colon$selrf.res <- SelRf$selec.history # 14 genes with 0 OOB, 11 genes with 1.23% OOB
colon$selrf.sel <- SelRf$selected.vars
colon$selrf.mod <- SelRf$rf.model
colon$selrf.imp <- SelRf$initialImportances
colon$selrf.imp <- colon$selrf.imp[order(colon$selrf.imp[,1],decreasing=T),]
colon$selrf.imp.ord  <- cbind(colon$selrf.imp,1:length(colon$selrf.imp))
tmp <- colon$selrf.imp.ord[colon$selrf.sel,]
colnames(tmp) <- c("Gini.Imp","Order")
colon$selrf.sel <- tmp[order(tmp[,2]),]
colon$selrf.sel
colon$selrf.res[,-2]

pdf(file="RFplot.pdf")
varImpPlot(colon$selrf,type=2,n.var=30,scale=F,main="Varible Importance(Gini) for top 30 predictors")
dev.off()



idBreast <- "GSE15852"
dim(breast$exp)      # Original:         22283, 86
dim(breast$exp_pre)  # No Duplicate/NA:  13211, 86
dim(breast$sig)      # Significant:      2582, 4
sum(breast$pat[,'group']=='0')  # Normal Sample: 43

sum(breast$pat[,'group']=='1')  # Cancer Sample: 43

idPanc <- "GSE15471"
dim(panc$exp)      # Original:         54675, 86
dim(panc$exp_pre)  # No Duplicate/NA:  21049, 86
dim(panc$sig)      # Significant:      14770, 4
sum(panc$pat[,'group']=='0')  # Normal Sample: 39
sum(panc$pat[,'group']=='1')  # Cancer Sample: 39



# Cross-validation
## colon.sel
### in breast
G0 <- rownames(breast$pat[breast$pat[,'group']=='0',])
G1 <- rownames(breast$pat[breast$pat[,'group']=='1',])
x <- breast$exp_pre[match(rownames(colon$selrf.sel),rownames(breast$exp_pre)),] 
x <- x[,c(match(G1,colnames(x)),match(G0,colnames(x)))]
x <- t(x)
cl <- factor(breast$pat[match(rownames(x),rownames(breast$pat)),'group'])
rf.b <- breast$rf.colonSel<- randomForest(x,cl,ntree=10001,importance=T)
breast$rf.colonSel$confusion
#    0  1 class.error
# 0 31 12   0.2790698
# 1 17 26   0.3953488
selrf.b <- varSelRF(x,cl)
selrf.b$selec.history[,-2][1,]
#   Number.Variables       OOB    sd.OOB
# 1               14 0.3372093 0.0509787

### in colon
G0 <- rownames(colon$pat[colon$pat[,'group']=='0',])
G1 <- rownames(colon$pat[colon$pat[,'group']=='1',])
x <- colon$exp_pre[match(rownames(colon$selrf.sel),rownames(colon$exp_pre)),] 
x <- x[,c(match(G1,colnames(x)),match(G0,colnames(x)))]
x <- t(x)
cl <- factor(colon$pat[match(rownames(x),rownames(colon$pat)),'group'])
rf.c <- colon$rf.colonSel<- randomForest(x,cl,ntree=10001,importance=T)
colon$rf.colonSel$confusion
#    0  1 class.error
# 0 27  1  0.03571429
# 1  0 53  0.00000000
selrf.c <- varSelRF(x,cl)
selrf.c$selec.history[,-2][1,]
#   Number.Variables OOB sd.OOB
# 1               14   0      0


### in panc
G0 <- rownames(panc$pat[panc$pat[,'group']=='0',])
G1 <- rownames(panc$pat[panc$pat[,'group']=='1',])
x <- panc$exp_pre[match(rownames(colon$selrf.sel),rownames(panc$exp_pre)),] 
x <- x[,c(match(G1,colnames(x)),match(G0,colnames(x)))]
x <- t(x)
cl <- factor(panc$pat[match(rownames(x),rownames(panc$pat)),'group'])
rf.p <- panc$rf.colonSel<- randomForest(x,cl,ntree=10001,importance=T)
panc$rf.colonSel$confusion
#    0  1 class.error
# 0 34  5  0.12820513
# 1  3 36  0.07692308
selrf.c <- varSelRF(x,cl)
selrf.c$selec.history[,-2][1,]
#   Number.Variables       OOB     sd.OOB
# 1               14 0.1025641 0.03435201

SelVarRf <- function(dat,pat,sig){
# dat is matrix of expression value, "ParseObj$test": row=gene, colon=sample
# pat is sample information, "ParseObj$pat".
# sig is sig genes, "rownames(ParseObj$sig)".
    G0 <- rownames(pat[pat[,'group']=='0',]) # Normal samples ID
    G1 <- rownames(pat[pat[,'group']=='1',]) # Disease samples ID
    x <- dat[match(sig,rownames(dat)),]  # Keep signifcant gene's info 
    x <- x[,c(match(G0,colnames(x)),match(G1,colnames(x)))] # Keep related Samples ID
    x <- t(x)                          # For random forest
    cl <- factor(pat[match(rownames(x),rownames(pat)),'group'])
    rf <- varSelRF(x,cl,ntree=10001)
    return(rf)
}

IterCalOOB <- function(dat,oput,pat,n){
# dat: "ParseObj$exp_pre"
# gened: "oput.3"
# n: number of the end.
    tmp.oob <- c()
    G0 <- rownames(pat[pat[,'group']=='0',]) # Normal samples ID
    G1 <- rownames(pat[pat[,'group']=='1',]) # Disease samples ID
    for(i in 1:n){
        x <- dat[match(rownames(oput)[1:i],rownames(dat)),] 
        x <- x[,c(match(G1,colnames(x)),match(G0,colnames(x)))]
        x <- t(x)
        cl <- factor(pat[match(rownames(x),rownames(pat)),'group'])
        tmp.rf <- varSelRF(x,cl,ntree=10001,vars.drop.frac=0.9)
        tmp.oob <- c(tmp.oob,tmp.rf$selec.history[1,"OOB"])
        print(i)
    }
    return(tmp.oob)
}

DEGinOneSet <- function(dat.1,dat.2,dat.3){
# dat: "ParseObj"
    p.2 <- dat.2$test[match(rownames(dat.1$sig),rownames(dat.2$test)),'FDR.adjust.p']
    t.2 <- dat.2$test[match(rownames(dat.1$sig),rownames(dat.2$test)),'Statistic']
    p.3 <- dat.3$test[match(rownames(dat.1$sig),rownames(dat.3$test)),'FDR.adjust.p']
    t.3 <- dat.3$test[match(rownames(dat.1$sig),rownames(dat.3$test)),'Statistic']
    oput <- cbind(dat.1$sig,fdr.p.2=p.2,stat.2=t.2,fdr.p.3=p.3,stat.3=t.3)
    oput <- na.omit(oput)              # Del gene not detect in other array 
    oput <- oput[-which(oput[,'fdr.p.2']<=0.05&oput[,'stat.2']*oput[,'Statistic']>0),] # Del same DEGs in array 2 
    oput <- oput[-which(oput[,'fdr.p.3']<=0.05&oput[,'stat.3']*oput[,'Statistic']>0),] # Del same DEGs in array 3 
    return(oput)
}




# find Specific DEGs
colon$sigSpec <- tmp.a <- DEGinOneSet(colon,breast,panc)
breast$sigSpec <- tmp.a <- DEGinOneSet(breast,colon,panc)
panc$sigSpec <- tmp.a <- DEGinOneSet(panc,breast,colon)

## common part
dis$rf <- tmp.rf <- SelVarRf(dis$exp_pre,dis$pat,rownames(dis$sigSpec))
oput <- dis$sigSpec
oput.1 <- cbind(oput[order(oput[,4]),],p.order=c(1:nrow(oput))) # Order by FDR.p 
oput.2 <- cbind(oput.1,imp=tmp.rf$initialImportances[rownames(oput.1),1]) # Add Importance
oput.3 <- cbind(oput.2[order(oput.2[,"imp"],decreasing=T),],imp.order=c(1:nrow(oput.2)))
tmp.rf$selec.history[,-2]
tmp.oob <- IterCalOOB(dis$exp_pre,oput.3,dis$pat,25)
oput.4 <- cbind(oput.3,imp.OOB=c(tmp.oob,rep(NA,dim(oput.3)[1]-length(tmp.oob))))
head(oput.4,30)
tmp.rf$selec.history[,-2]
oput.5 <- cbind(oput.4[,-c(5,6,7,8)],colon$test[match(rownames(oput.4),rownames(colon$test)),])
oput.6 <- cbind(oput.5,breast$test[match(rownames(oput.4),rownames(breast$test)),])
head(oput.6,30)
dis$oput <- oput.6

## colon
dis <- colon
# go to common part
colon <- dis

## breast
dis <- breast
# common part
breast <- dis

## panc
dis <- panc
# common part
panc <- dis

# Writing
write.csv(colon$oput,file="colon_oput.csv")
write.csv(breast$oput,file="breas_oput.csv")
write.csv(panc$oput,file="panc_oput.csv")


genes <- c('SCL16A3','GMEB2','FLT3LG','RPS6KA1')
genes <- rownames(ctcColon$sigSpec.colon)
tmp.n <- names(kEntrezID[match(genes,kEntrezID)])
write.csv(cbind(genes,tmp.n),file="GeneEntrezID_oput.csv")
###### Information Index ######
#####k. constant #####
# kKnownMarkers: known marker list
length(names(kKnownMarkers))
# [1] 25
length(kKnownMarkers$colon)
# [1] 24
length(kKnownMarkers$breast)
# [1] 33
length(kKnownMarkers$pancreas)
# [1] 8

names(colon)
#  [1] "exp"           "prob"          "pat"           "exp_pre"       "test"         
#  [6] "sig"           "rf"            "rf2"           "rf3"           "rf4"          
# [11] "selrf.res"     "selrf.mod"     "selrf.imp"     "selrf.sel"     "selrf.imp.ord"
# [16] "selrf"         "rf.colonSel"   "oput"          "oob.2"         "SpecDEG"      
# [21] "sigSpec"      
nrow(colon$sig)
# [1] 3855

names(ctcColon)
# [1] "exp"           "prob"          "pat"           "exp_pre"       "test"         
# [6] "sig"           "sigSpec.colon"
nrow(ctcColon$sig)
# [1] 387

names(breast)
#  [1] "exp"           "prob"          "pat"           "exp_pre"       "test"         
#  [6] "sig"           "selrf.mod"     "selrf.imp"     "selrf.imp.ord" "rf.colonSel"  
# [11] "rf"            "oput"          "SpecDEG"       "sigSpec"      

names(ctcBreast)
# [1] "exp"            "prob"           "pat"            "exp_pre"        "test"          
# [6] "sig"            "sigSpec.breast"

names(panc)
#  [1] "exp"           "prob"          "pat"           "exp_pre"       "test"         
#  [6] "sig"           "selrf"         "selrf.res"     "selrf.sel"     "selrf.mod"    
# [11] "selrf.imp"     "selrf.imp.ord" "rf.colonSel"   "SpecDEG"       "sigSpec"      
# [16] "rf"            "oput"         

names(ctcPanc)
# [1] "exp"          "prob"         "pat"          "exp_pre"      "test"        
# [6] "sig"          "sigSpec.panc"
