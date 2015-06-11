## to be add to git
## Version info: R 3.2.0, bioconductor 3.1
## Usage: 

## install packages
## biocLite()
source("http://bioconductor.org/biocLite.R")
biocLite()

## GEOquery, limma 
biocLite("GEOquery")
biocLite("limma")

## read into R
library(GEOquery)
library(limma)

## prase GSE files
# add G0 as disease, G1 as normal
ParseGSE <- function(f){
  # recive GSEnumber , and return a list with "exp"=expression data, "prob"=probe infomation, "pat"=patients information
  # "norml"=is.normalized, "deal.NA"=is dealed with NA, "group"=is pat grouped as control or treated for test.
  gse <- getGEO(f,GSEMatrix=FALSE)
  gpls <- GPLList(gse) #list of array: gpls[[1]]==> first array
  gsms <- GSMList(gse) #list of patients: gsms[[1]] ==> first patient
  
  # validation-1,2: if gse's platform_id == gpl's geo_accession / gsm's platform_id
  # validation-3: if probes in each patients is same.
  g1<-Meta(gse)$platform_id # common platform id
  pro1<-Table(gsms[[1]])$ID #icommon probes
  
  for(i in 1:length(gpls)){
    g2 <- Meta(gpls[[i]])$geo_accession
    pro2 <- Table(gpls[[i]])$ID
    ifelse(g1==g2,NA,print(paste("Waring: gpl(array):",i,"is wrong and should be checked!!!!!!!!!"))
    )
  }
  cat("\n") # print can not use "\n"
 
  for(i in 1:length(gsms)){
    g2<-Meta(gsms[[i]])$platform_id
    pro2 <- Table(gsms[[i]])$ID
    ifelse(g1==g2,NA,print(paste("Waring: gsm(array):",i,"is wrong and should be checked!!!!!!!!!"))
    )
    ifelse(pro1==pro2,NA,print(paste("Waring: probes of patient:",i,"is wrong and should be checked!!!!!!!!!"))
    )
  }
  
  # Extract values into matrix 
  ## should be changed. rowname is wroong..................................
  exp <- do.call("cbind", lapply(gsms, function(x) as.numeric(Table(x)$VALUE)))  
  rownames(exp) <- pro1
  
  # Extract patients infomation
  #pat_info <- do.call("rbind", lapply(gsms, function(x) Meta(x)$characteristics_ch1))
  #pat_info <- do.call("rbind", lapply(gsms, function(x) Meta(x)[1:5]))
  #pat_info <- do.call("rbind", lapply(gsms, function(x) Meta(x)[1:length(Meta(x))]))
  pat_info <- do.call("cbind", lapply(gsms, function(x) unlist(Meta(x))))
  #pat_info <- matrix(unlist(), ncol = 10, byrow = TRUE)
    # Extract probes info
  p_info <- Table(gpls[[1]])
  return(list(exp=exp, prob=p_info,pat=pat_info))
}

GroupPat <- function(m.pat){
  ind <- vector()
  ii <- 0
  for(i in 1:nrow(m.pat)){
    f <- levels(factor(m.pat[i,]))
    if(length(f)>1 & length(f) < ncol(m.pat)){
      ii <- ii+1
      ind <- c(ind,i)
      print(paste(ii,"-",f))
    }
  }
  inp <- as.numeric(readline(prompt="choose which col contains infomation to identify groups:\t"))
  inp <- ind[inp]
  cat("input the name for each groups, 'G0' for control, 'G1' for disease, no_input for not use:\n")
  #fac.1 <- as.character(levels(factor(m.pat[inp,]))[1])
  #fac.2 <- as.character(levels(factor(m.pat[inp,]))[2])
  fac <- as.character(levels(factor(m.pat[inp,])))
  g.name <- character()
  for(i in fac){
    g.name <- c(g.name,readline(prompt=paste(i,":\t")))
  }
  #g.name.1 <- readline(prompt=paste(fac.1,":\t"))
  #g.name.2 <- readline(prompt=paste(fac.2,":\t"))
  m.pat <- rbind(m.pat,group = NA)
  for(i in 1:length(fac)){
    m.pat[nrow(m.pat),which(m.pat[inp,]==fac[i])] =g.name[i]
  }
  #m.pat[nrow(m.pat),which(m.pat[inp,]==fac.1)] =g.name.1
  #m.pat[nrow(m.pat),which(m.pat[inp,]==fac.2)] =g.name.2
  return(t(m.pat))
}

PreProc <- function(lm){
# recive list obj from parseGSE/groupPat fun.
  exp <- lm$exp
  m0 <- matrix()
  m1 <- matrix()
  print("normalization")
  lm$exp <- normalizeQuantiles(exp)
  cat("Done!\n")
  cname.0 <- rownames(lm$pat[which(lm$pat[,"group"]=="G0"),])
  cname.1 <- rownames(lm$pat[which(lm$pat[,"group"]=="G1"),])
  #lm$g0 <- exp[,which(colnames(exp)==cname.0)]
  #lm$g1 <- exp[,which(colnames(exp)==cname.1)]
  lm$g0 <- exp[,match(cname.0, colnames(exp))]
  lm$g1 <- exp[,match(cname.1, colnames(exp))]
  return(lm)
}
  
MatrixTest <- function(d1,d2){
  var.p <- var.test(d1,d2)$p.value
  var.equ <- FALSE
  num.d1 <- vector()
  num.d2 <- vector()
  if(var.p <= 0.05) var.equal <- TRUE

  output <- matrix(ncol=4)
  total <- nrow(d1)
  cat("calculate number of value\n")
  pd <- txtProgressBar(min=0,max=total,style=3)
  # count value number
  for (i in 1:nrow(d1)){
    setTxtProgressBar(pd,i)
    num.d1 <- c(num.d1,sum(!is.na(d1[i,])))
    num.d2 <- c(num.d2,sum(!is.na(d2[i,])))
  }
  close(pd)
  # for t.test and wil.test
  cat("calculate t.test and wilcox.test\n")
  for (i in 1:nrow(d1)){
    setTxtProgressBar(pd,i)
    if(num.d1[i]<3 | num.d2[i]<3){
      output <- rbind(output,c(NA,NA,NA,NA))
    }else{
      t <- t.test(d1[i,],d2[i,],var.equal=var.equ)
      w <- wilcox.test(d1[i,],d2[i,],var.equal=var.equ)
      output <- rbind(output,c(t$p.value,t$statistic,w$p.value,w$statistic))
    }
  }
  close(pd)
  cat("\n")

  # add adjust p value
  output <- cbind(output[,1:2],fdr.t.p=p.adjust(output[,1],method="BH"),
                  output[,3:4],dfr.t.p=p.adjust(output[,3],method="BH"))

  # add rowname and colname
  colnames(output) <- c("t.p-value","t.statistic","t.adjust-p","w.p-value","w.statistic","w.adjust-p")
  output <- output[-1,]
  rownames(output) <- rownames(d1)
  return(output)
}

Main <- function(gse){
  dat <- ParseGSE(gse)
  dat$pat <- GroupPat(dat$pat)
  dat.proc <- PreProc(dat)
  dat.test <- MatrixTest(dat.proc$g0,dat.proc$g1)
  return(list(exp=dat.proc, test=dat.test))
}

SelectGene <- function(m,colname,thres){
  output<- m[which(m[,colname]<=thres),]
  return(output)
}

CommAndSpec <- function(m.list,colname){
  probs <- vector()
  for (i in m.list){
    probs <- c(probs,i[,colname])
  }
  # remove duplicated
  probs <- unique(probs)
  # remove ""
  probs <- probs[!probs %in% ""]

  output <- cbind(1:length(probs),probs)
  colna <- vector()
  for(i in 1:length(m.list)){
    colna <- c(colna,names(m.list[i]))
    mat <- match(probs,m.list[[i]][,colname])
    output <- cbind(output,mat)
  }
  num <- vector()
  num <- rowSums(!is.na(output[,-1:-2]))
  output <- cbind(output,num)
  colnames(output)<- c("Index","Probes", colna,"cancer_num")
  return(output[,-1])
}
# small sample
a <- renal$test
a_p <- renal$exp$prob
b <- lung$test
c <- panc$test
b_p <- lung$exp$prob
c_p <- panc$exp$prob

a <- cbind(a,"geneID"=as.character(a_p[match(rownames(a),a_p[,"ID"]),"ENTREZ_GENE_ID"]))
b <- cbind(b,"geneID"=as.character(b_p[match(rownames(b),b_p[,"ID"]),"ENTREZ_GENE_ID"]))
c <- cbind(c,"geneID"=as.character(c_p[match(rownames(c),c_p[,"ID"]),"ENTREZ_GENE_ID"]))

a.sel <- SelectGene(a,"t.adjust-p",0.05)
b.sel <- SelectGene(b,"t.adjust-p",0.05)
c.sel <- SelectGene(c,"t.adjust-p",0.05)

a.sel <- SelectGene(a.sel,"w.adjust-p",0.05)
b.sel <- SelectGene(b.sel,"w.adjust-p",0.05)
c.sel <- SelectGene(c.sel,"w.adjust-p",0.05)
list.sample <- list("renal"=a.sel,"lung"=b.sel,"panc"=c.sel)

res <- CommAndSpec(list.sample,"geneID")

o.list <- list("colon"=colon,"liver"=liver,"panc"=panc,"renal"=renal,"lung"=lung,"breast"=breast)
o.list <- list("colon"=colon,"liver"=liver)
Main_Select <- function(obj.list){
  mat.list <- list() 
  # select genes with p<0.05
  for(i in 1:length(obj.list)){
    n <- names(obj.list[i])
    exp <- obj.list[[i]]$test
    mat.list[[n]] <- SelectGene(SelectGene(exp,"t.adjust-p",0.05),"w.adjust-p",0.05)
  }
  # identify common gene indientfier
  p_info <- list()
  gene.col <- character()
  for(i in 1:length(obj.list)){
    p_info[[i]]<- colnames(obj.list[[i]]$exp$prob)
  }
  comm.pinfo <- Reduce(intersect,p_info)
  print(comm.pinfo)
  inp <- readline("choose which is used as gene indentifer(GB_ACC:Gene Bank ID)\t")
  gene.col <- comm.pinfo[as.numeric(inp)]
  print(gene.col)
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
selected <- Main_Select(o.list)
