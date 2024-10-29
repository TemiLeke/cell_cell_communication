
.myPackageInstaller=function(){
  
  #Bash command
  ##sudo apt-get update
  ##sudo apt-get install libgmp3-dev
  install.packages('BiocManager')
  packagesRequired=c("ggbeeswarm","CrossClustering","parallel","colorRamps","ggplot2","PRROC","pROC","reshape","reshape2","Hmisc","patchwork","dplyr","ClusterR","hues","effsize","ggwordcloud")
  
  if(("XML" %in% rownames(installed.packages())) == FALSE){
    if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
    require(devtools)
    install_version("XML", version = "3.98-1.17", repos = "http://cran.us.r-project.org")
  }
  
  
  for(i in 1:length(packagesRequired)){
    if((packagesRequired[i] %in% rownames(installed.packages())) == FALSE) {install.packages(packagesRequired[i])}
  }
  
  packagesRequired=c("Biobase","qvalue","lfa","multtest","limma","biomaRt","GenomicFeatures","EnsDb.Hsapiens.v75","EnsDb.Hsapiens.v86","monocle","org.Hs.eg.db","DropletUtils","scran","AUCell", "RcisTarget","GENIE3","SCENIC","zoo", "mixtools", "rbokeh","DT", "NMF", "pheatmap", "R2HTML", "Rtsne","doMC", "doRNG","MAST","EnsDb.Mmusculus.v79","gskb","SC3",'PCAtools')
  
  for(i in 1:length(packagesRequired)){
    if((packagesRequired[i] %in% rownames(installed.packages())) == FALSE) {BiocManager::install(packagesRequired[i], update = FALSE,ask =FALSE)}
  }
  
  packagesRequired=c("Seurat","jackstraw","caret")
  for(i in 1:length(packagesRequired)){
    if((packagesRequired[i] %in% rownames(installed.packages())) == FALSE) {install.packages(packagesRequired[i])}
  }
  
  if(("SCENIC" %in% rownames(installed.packages())) == FALSE){
    if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
    devtools::install_github("aertslab/SCENIC") 
  }
  
  if(("harmony" %in% rownames(installed.packages())) == FALSE){
    if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
      devtools::install_github("immunogenomics/harmony")
  }
  if(("ggraph" %in% rownames(installed.packages())) == FALSE){
    install.packages("ggraph")
  }
  
  if(("ape" %in% rownames(installed.packages())) == FALSE){
    install.packages("ape")
  }
  if(("seriation" %in% rownames(installed.packages())) == FALSE){
    install.packages("seriation")
  }
  
  
}


#library(rliger)
library(plyr)
library(dplyr)
library(Seurat)
library(patchwork)
library(future)
library(ggplot2)
#library(gskb)
library(scran)


#path="~/myBucket/Alec_BCell_Unfiltered/Brain/SCP_ChoroidPlexus_83k";mtx="sn_rna_counts.mtx";genes="sn_rna_genes.tsv";barcodes=NULL;addBarcodes=T
.myRead10X=function(path,mtx=NULL,genes=NULL,barcodes=NULL,addBarcodes=T,allow_nonInteger=F){
  require(scran)
  require(Matrix)
  require(DropletUtils)
  #	Read 10X Genomics Matrix Exchange Format files created by CellRanger
  if(is.null(mtx) || is.null(genes) || is.null(barcodes)) {
    files <- dir(path)
    if(is.null(mtx)) {
      if(sum(grepl("matrix.mtx", files))) mtx <- files[grepl("matrix.mtx",files)]
      
      if(is.null(mtx)) {
        stop("Can't find matrix.mtx file, please specify filename explicitly")
        } else {
          if(length(mtx)>1){
            print("Multiple matrix files was found! only the first matrix file will be considered")
          }
        mtx=mtx[1]
      }
    }
    if(is.null(genes)) {
      if(sum(grepl("genes.tsv", files))>0) genes <- files[grepl("genes.tsv", files)]
      if(sum(grepl("features.tsv", files))>0) genes <- files[grepl("features.tsv", files)]
      
      if(is.null(genes)){
        stop("Can't find genes.tsv or features.tsv, please specify filename explicitly")
      } else {
        genes=genes[1]
      }
      
    }
    if(is.null(barcodes)&addBarcodes) {
      if(sum(grepl("barcodes.tsv", files))>0) barcodes <- files[grepl("barcodes.tsv", files)]
      if(!is.null(barcodes)){
        barcodes=barcodes[1]
      }
    }
  }
  
  #	Add path
  mtx <- file.path(path,mtx)
  genes <- file.path(path,genes)
  if(!is.null(barcodes)) barcodes <- file.path(path,barcodes)
  
  #	Fetch header info for checking
  skip_lines=3
  N <- scan(mtx,skip=skip_lines-1,sep=" ",nmax=3,quiet=TRUE)
  ngenes <- N[1]
  ncells <- N[2]
  nmtx <- N[3]
  
  if(ngenes<1000){
    skip_lines=2
    N <- scan(mtx,skip=skip_lines-1,sep=" ",nmax=3,quiet=TRUE)
    ngenes <- N[1]
    ncells <- N[2]
    nmtx <- N[3]
  }
  
  
  #	Read gene Ids
  Genes <- read.table(genes,header=FALSE,comment.char="",sep="\t",colClasses="character")
  if(nrow(Genes) != ngenes) stop("Number of feature IDs doesn't agree with header information in mtx file")
  colnames(Genes)[1] <- "Gene_id"
  if(ncol(Genes)>1){
    names(Genes)[2] <- "Symbol"
    if(ncol(Genes) > 2L) names(Genes)[3] <- "Type"
  }
  
  
  #	Read mtx file of counts
  if(allow_nonInteger){
    m <- read.table(mtx,skip=skip_lines,header=FALSE,comment.char="",sep=" ",nrows=nmtx)
  } else {
    m <- read.table(mtx,skip=skip_lines,header=FALSE,comment.char="",sep=" ",colClasses="integer",nrows=nmtx)
  }
  
  
  y=sparseMatrix(m[,1],m[,2], x=m[,3], dims=c(ngenes,ncells))
  
  #	Optionally read barcodes
  if(is.null(barcodes)) {
    Samples <- NULL
  } else {
    Barcodes <- scan(barcodes,what="",quiet=TRUE)
    if(length(Barcodes) != ncells) stop("Number of barcodes doesn't agree with header information in mtx file")
    Samples <- data.frame(Barcode=Barcodes)
    Samples$addCol="alaki"
  }
  
  
  cellData=unique(m$V2)
  Samples=Samples[cellData,]
  y=y[,cellData]
  ind=which(colSums(y)>50)
  Samples=Samples[ind,]
  y=y[,ind]
  
  ind=which(rowSums(y>0)>0)
  y=y[ind,,drop=F]
  Genes=Genes[ind,,drop=F]
  Samples=Samples[,-which(colnames(Samples)=="addCol")]
  if(is.null(nrow(Samples))){
    Samples=data.frame(Barcode=Samples,stringsAsFactors = F)
  }
  if(class(Genes)==class("")){
    Genes=data.frame(Gene_id=Genes,stringsAsFactors = F)
  }
  row.names(y)=Genes$Gene_id
  colnames(y)=Samples$Barcode
  res=SingleCellExperiment(assays = list(counts = y),colData = Samples,rowData=Genes)
  
  return(res)
}

.myRead10X_h5=function (filename,organism,nGene_thr=NULL) {
  library(Matrix)
  library(sparseMatrixStats)
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  
  infile = hdf5r::H5File$new(filename = filename, mode = "r")
  genomes = names(x = infile)
  output = list()
  if (hdf5r::existsGroup(infile, "matrix")) {
    feature_name_slot <- "features/name"
    feature_id_slot <- "features/id"
    
  }else {
    feature_name_slot = "gene_names"
    feature_id_slot = "genes"
  }
  for (genome in genomes) {
    counts <- infile[[paste0(genome, "/data")]]
    indices <- infile[[paste0(genome, "/indices")]]
    indptr <- infile[[paste0(genome, "/indptr")]]
    shp <- infile[[paste0(genome, "/shape")]]
    features_id <- infile[[paste0(genome, "/", feature_id_slot)]][]
    features_name <- infile[[paste0(genome, "/", feature_name_slot)]][]
    barcodes <- infile[[paste0(genome, "/barcodes")]]
    sparse.mat <- sparseMatrix(i = indices[] + 1, p = indptr[], 
                               x = as.numeric(x = counts[]), dims = shp[], giveCsparse = FALSE)
    
    rownames(x = sparse.mat) <- features_id
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = "dgCMatrix")
    if(!is.null(nGene_thr)){
      nGene=sparse.mat
      nGene@x=rep(1,length(nGene@x))
      nGene=colSums(nGene)
      sparse.mat=sparse.mat[,which(nGene>=nGene_thr)]
    }
    
    fd=data.frame(id=features_id,gene_name=features_name,stringsAsFactors = F)
    row.names(fd)=fd$id
    res=.myExpSetCreatorFn(inputExpData=sparse.mat,organism=organism,minExpCells=0,inputPdata=NULL,inputFdata=fd,addExtraAnno=T,server=T,redownload_files=F,ncores=5)
    output[[genome]] <- res
  }
  infile$close_all()
  if (length(x = output) == 1) {
    return(output[[genome]])
  }
  else {
    return(output)
  }
}


.extracBindDetailFn=function(x1,x2,batchNames){
  require(Matrix)
  
  tmp1=setdiff(row.names(x1),row.names(x2))
  if(class(x1)=="SingleCellExperiment"){
    fd1=rowData(x1)
    expr1=counts(x1)
    pd1=as.data.frame(colData(x1))
  } else if(class(x1)=="Seurat") {
    fd1=as.data.frame(x1@assays$RNA@meta.features)
    expr1=x1@assays$RNA@counts
    pd1=as.data.frame(x1@meta.data)
  } else {
    stop("Unrecognized dataset class!")
  }
  
  if(class(x2)=="SingleCellExperiment"){
    fd2=rowData(x2)
    expr2=counts(x2)
    pd2=as.data.frame(colData(x2))
  } else if(class(x2)=="Seurat") {
    fd2=as.data.frame(x2@assays$RNA@meta.features)
    expr2=x2@assays$RNA@counts
    pd2=as.data.frame(x2@meta.data)
  } else {
    stop("Unrecognized dataset class!")
  }
  
  tmpAnno1=fd1[which(row.names(x1) %in% tmp1),]
  
  tmpMat1=sparseMatrix(i=NULL,j=NULL,dims = c(length(tmp1),ncol(x2)))
  row.names(tmpMat1)=tmp1
  
  x2_c=rbind(expr2,tmpMat1)
  
  x2_c_pd=as.data.frame(pd2)
  x2_c_fd=as.data.frame(plyr::rbind.fill(as.data.frame(fd2),as.data.frame(tmpAnno1)))
  x2_c=SingleCellExperiment(assays = list(counts = x2_c),colData = x2_c_pd,rowData=x2_c_fd)
  
  tmp2=setdiff(row.names(x2),row.names(x1))
  tmpAnno2=fd2[which(row.names(x2) %in% tmp2),]
  
  tmpMat2=sparseMatrix(i=NULL,j=NULL,dims = c(length(tmp2),ncol(x1)))
  row.names(tmpMat2)=tmp2
  
  x1_c=rbind(expr1,tmpMat2)
  
  x1_c_pd=as.data.frame(pd1)
  x1_c_fd=plyr::rbind.fill(as.data.frame(fd1),as.data.frame(tmpAnno2))
  x1_c=SingleCellExperiment(assays = list(counts = x1_c),colData = x1_c_pd,rowData=x1_c_fd)
  
  x2_c=x2_c[match(row.names(x1_c),row.names(x2_c)),]
  
  if(all(row.names(x2_c)==row.names(x1_c))){
    
    if(batchNames[1]!=""){
      row.names(colData(x1_c))=paste0(batchNames[1],"_",colnames(x1_c))
      colnames(x1_c)=paste0(batchNames[1],"_",colnames(x1_c))
      x1_c$anno_batch=batchNames[1]
    }
    
    if(batchNames[2]!=""){
      colnames(x2_c)=paste0(batchNames[2],"_",colnames(x2_c))
      row.names(colData(x2_c))=paste0(batchNames[2],"_",colnames(x2_c))
      x2_c$anno_batch=batchNames[2]
    }
    
    #tst1=counts(x1_c)
    #tst2=counts(x2_c)
    #tst1=summary(tst1)
    #tst2=summary(tst2)
    #tst2$j=tst2$j+max(tst1$j)
    #tst=rbind(tst1,tst2)
    #tst=sparseMatrix(i = tst,j = tstt,x = tsttt)
    
    x_m_exp=cbind(counts(x1_c),counts(x2_c))
    pd_m_exp=plyr::rbind.fill(as.data.frame(colData(x1_c)),as.data.frame(colData(x2_c)))
    fd=as.data.frame(rowData(x1_c))
    x_m=SingleCellExperiment(assays=list(counts=x_m_exp),colData=pd_m_exp,rowData=fd)
  } else {
    stop("Error in the merging!")
  }
  return(x_m)
}

#inputList=.tmp;batchNames=NULL;verbose=F
.mycBindFn=function(inputList,batchNames=NULL,verbose_level=1){
  #cbinds multiple singleCellExpression datasets with differring number of rows.
  #inputList: the list of datasets to be merged
  #batchNames: the batch name to be assigned to each dataset. length(batchNames)==length(inputList)
  res_m=""
  if(!is.null(batchNames)){
    if(length(inputList)==1){
      res_m=inputList[[1]]
    } else if(length(inputList)==2){
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = batchNames[1:2])
    } else if(length(inputList)>2){
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = batchNames[1:2])
      for(i in 3:length(inputList)){
        res_m=.extracBindDetailFn(x1=res_m,x2=inputList[[i]],batchNames=c("",batchNames[i]))
      }
    }
  } else {
    if(length(inputList)==1){
      res_m=inputList[[1]]
    } else if(length(inputList)==2){
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = c("",""))
    } else if(length(inputList)>2){
      #x1=inputList[[1]];x2=inputList[[2]];batchNames = c("","")
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = c("",""))
      for(i in 3:length(inputList)){
        
        #x1=res_m;x2=inputList[[i]];batchNames=c("","")
        res_m=.extracBindDetailFn(x1=res_m,x2=inputList[[i]],batchNames=c("",""))
        if(verbose_level>1){
          print(paste("dataset:",i,"; nrow:",nrow(res_m),"; ncol",ncol(res_m)))
        }
      }
    }
  }
  
  if(verbose_level>=1){
    print("batch information is in the anno_batch variable")
  }
  
  return(res_m)
  
}

.fData=function(inputData){
  fd=as.data.frame(rowData(inputData))
  return(fd)
}

.pData=function(inputData){
  pd=""
  if(class(inputData)=="Seurat"){
    pd=res=inputData[[]]
  } else {
    if(class(inputData)=="SingleCellExperiment"){
      pd=as.data.frame(colData(inputData))
    } else {
      pd=pData(inputData)
    }
    
  }
  
  return(pd)
}

.extraMitoGenes=function(organism,redownload_files=T){
  #organism: Human, Mouse
  
  if(tolower(organism)=="human"){
    gns=.extraHumanGeneAnnoAdderFn()
  } else if (tolower(organism)=="mouse"){
    gns=.extraMouseGeneAnnoAdderFn()
  } else if(tolower(organism)=="macaque"){
    gns=.extraMacaqueGeneAnnoAdderFn()
  } else {
    stop("Wrong organism name!")
  }
  
  gns=gns[gns$seqnames=="MT",]
  gns2=unlist(gns$entrezid)
  gns$entrezgene_id=gns2
  gns$ensembl_gene_id=gns$gene_id
  return(gns)
}

.extrarRNAGenes=function(organism,redownload_files=T){
  #organism: Human, Mouse
  
  if(tolower(organism)=="human"){
    gns=.extraHumanGeneAnnoAdderFn()
  } else if (tolower(organism)=="mouse"){
    gns=.extraMouseGeneAnnoAdderFn()
  } else if(tolower(organism)=="macaque"){
    gns=.extraMacaqueGeneAnnoAdderFn()
  } else {
    stop("Wrong organism name!")
  }
  
  gns=gns[which(gns$gene_biotype=="rRNA"),]
  gns2=unlist(gns$entrezid)
  gns$entrezgene_id=gns2
  gns$ensembl_gene_id=gns$gene_id
  return(gns)
}

.extraRibosomalProteinGenes=function(organism,redownload_files=T){
  #organism: Human, Mouse
  
  if(tolower(organism)=="human"){
    gns=.extraHumanGeneAnnoAdderFn()
  } else if (tolower(organism)=="mouse"){
    gns=.extraMouseGeneAnnoAdderFn()
  } else if(tolower(organism)=="macaque"){
    gns=.extraMacaqueGeneAnnoAdderFn()
  } else {
    stop("Wrong organism name!")
  }
  
  gns=gns[which(grepl("^RP[LS]",gns$gene_short_name,ignore.case = T)),]
  gns2=unlist(lapply(gns$entrezid,function(x) x[1]))
  gns$entrezgene_id=gns2
  gns$ensembl_gene_id=gns$gene_id
  return(gns)
}


.extraOXPHOSGenes=function(organism,inputGeneset=NULL,server=F){
  #organism: Human, Mouse
  
  if(tolower(organism)=="human"){
    gns=.extraHumanGeneAnnoAdderFn()
  } else if (tolower(organism)=="mouse"){
    gns=.extraMouseGeneAnnoAdderFn()
  } else if(tolower(organism)=="macaque"){
    gns=.extraMacaqueGeneAnnoAdderFn()
  } else {
    stop("Wrong organism name!")
  }
  
  if(is.null(inputGeneset)){
      if(!file.exists("../data/gazestani_serverFiles/OXPHOS_gene_symbols.txt")){
        stop("Error! file ../data/gazestani_serverFiles/OXPHOS_gene_symbols.txt is missing")
      }
      
      oxphosGenes=read.table("../data/gazestani_serverFiles/OXPHOS_gene_symbols.txt",sep="\t",stringsAsFactors = F)
    
  }
  
  if(sum(colnames(gns)=="gene_short_name")>0){
    gns=gns[which(toupper(gns$gene_short_name) %in% toupper(oxphosGenes[,1])),]
  } else {
    gns=gns[which(toupper(gns$symbol) %in% toupper(oxphosGenes[,1])),]
  }
  
  
  return(gns)
}

.extraIEGGenes=function(organism,server=F){
  #organism: Human, Mouse
  
  gns=.extraMouseGeneAnnoAdderFn()
  
  {
    #library(googleCloudStorageR)
    stop("Error! file gazestani_serverFiles/IEG_gene_symbols is missing")
    
    IEGenes=read.table("../data/gazestani_serverFiles/IEG_gene_symbols",sep="\t",stringsAsFactors = F)
  }
  
  gns=gns[tolower(gns$gene_name) %in% tolower(IEGenes$V1),]
  
  if(tolower(organism)=="human"){
    
    {
      #library(googleCloudStorageR)
      if(!file.exists("../data/gazestani_serverFiles/ortholog_mapping.rda")){
        stop("Error! file ../data/gazestani_serverFiles/ortholog_mapping.rda is missing")
        #gcs_get_object("vgazesta/data/gazestani_serverFiles/orthologsFeb3/ortholog_mapping.rda", saveToDisk = "../data/gazestani_serverFiles/ortholog_mapping.rda",overwrite=T)
      }
      
      load("../data/gazestani_serverFiles/ortholog_mapping.rda")
    }
    
    
    
    mapping=mapping[mapping$mouse %in% gns$ensembl_gene_id,]
    gnsHuman=.extraHumanGeneAnnoAdderFn()
    gns=gnsHuman[gnsHuman$ensembl_gene_id %in% mapping$human,]
  } else if (tolower(organism)=="mouse"){
    
  } else if(tolower(organism)=="macaque"){
    gns=.extraMacaqueGeneAnnoAdderFn()
    gns=gns[which(tolower(gns$symbol) %in% tolower(IEGenes$V1)|tolower(gns$gene_name) %in% tolower(IEGenes$V1)),]
  } else {
    stop("Wrong organism name!")
  }
  
  gns$ensembl_gene_id=gns$gene_id
  return(gns)
}

.extraCellCycleGenes=function(organism){
  if(tolower(organism)=="human"){
    gns=.extraHumanGeneAnnoAdderFn()
  } else if (tolower(organism)=="mouse"){
    gns=.extraMouseGeneAnnoAdderFn()
  } else if(tolower(organism)=="macaque"){
    gns=.extraMacaqueGeneAnnoAdderFn()
  } else {
    stop("Wrong organism name!")
  }
  
  cc.s=Seurat::cc.genes.updated.2019$s.genes
  cc.g2m=Seurat::cc.genes.updated.2019$g2m.genes
  
  cc.g2m[cc.g2m=="PIMREG"]="FAM64A"
  cc.g2m[cc.g2m=="JPT1"]="HN1"
  
  cc.s=gns[toupper(gns$symbol) %in% toupper(cc.s),]
  cc.g2m=gns[toupper(gns$symbol) %in% toupper(cc.g2m),]
  
  cc.s$ensemble_gene_id=row.names(cc.s)
  cc.g2m$ensemble_gene_id=row.names(cc.g2m)
  
  return(list(cc.s=cc.s,cc.g2m=cc.g2m))
}

.matrixExtraction=function (object) {
  #This function extracts numeric matrix from various data formats
  y <- list()
  if (is(object, "list")) {
    if (is(object, "EList")) {
      y$exprs <- as.matrix(object$E)
      y$Amean <- base::rowMeans(y$exprs, na.rm = TRUE)
    }
    else {
      if (is(object, "EListRaw")) 
        stop("EListRaw object: please normalize first")
      if (is(object, "RGList")) 
        stop("RGList object: please normalize first")
      y$printer <- object$printer
      if (is.null(object$M)) 
        stop("data object isn't of a recognized data class")
      y$exprs <- as.matrix(object$M)
      if (!is.null(object$A)) 
        y$Amean <- base::rowMeans(as.matrix(object$A), na.rm = TRUE)
    }
    y$weights <- object$weights
    y$probes <- object$genes
    y$design <- object$design
    
  }
  else {
    if (is(object, "ExpressionSet")) {
      if (!requireNamespace("Biobase", quietly = TRUE)) 
        stop("Biobase package required but is not available")
      y$exprs <- Biobase::exprs(object)
    }
    else {
      if (is(object, "PLMset")) {
        y$exprs <- object@chip.coefs
        if (length(y$exprs) == 0) 
          stop("chip.coefs has length zero")
      }
      else {
        if (is(object, "marrayNorm")) {
          y$exprs <- object@maM
        }
        else {
          if (is(object, "eSet")) {
            if (!requireNamespace("Biobase", quietly = TRUE)) 
              stop("Biobase package required but is not available")
            y$exprs <- Biobase::assayDataElement(object, 
                                                 "exprs")
            if (length(object@featureData@data)) 
              y$probes <- object@featureData@data
            y$Amean <- base::rowMeans(y$exprs, na.rm = TRUE)
            if ("weights" %in% Biobase::assayDataElementNames(object)) 
              y$weights <- Biobase::assayDataElement(object, 
                                                     "weights")
          }
          else {
            if (is.vector(object)) 
              y$exprs <- matrix(object, nrow = 1)
            else {
              if(is(object,"Seurat")){
                y$exprs=as.matrix(GetAssayData(object = object, assay = DefaultAssay(object = object), slot = "counts"))
              } else {
                if(is(object,"SingleCellExperiment")){
                  y$exprs <- data.matrix(counts(object))
                } else {
                  if(is(object,"data.frame")){
                    y$exprs <- as.matrix(object)
                  } else {
                    y$exprs <- data.matrix(object)
                  }
                }
                
              }
            }
          }
        }
      }
    }
  }
  return(y$exprs)
}

.extraHVGvstFn=function (object, loess.span = 0.3, mean.function = Seurat:::FastExpMean, dispersion.function = Seurat:::FastLogVMR, ...) {
  
  if (!inherits(x = object, "Matrix")) {
    object <- as(object = as.matrix(x = object), Class = "Matrix")
  }
  if (!inherits(x = object, what = "dgCMatrix")) {
    object <- as(object = object, Class = "dgCMatrix")
  }
  {
    clip.max <- sqrt(x = ncol(x = object))
    
    hvf.info <- data.frame(mean = rowMeans(x = object))
    hvf.info$variance <- Seurat:::SparseRowVar2(mat = object, mu = hvf.info$mean, 
                                                display_progress = F)
    hvf.info$variance.expected <- 0
    hvf.info$variance.standardized <- 0
    not.const <- hvf.info$variance > 0
    fit <- loess(formula = log10(x = variance) ~ log10(x = mean), 
                 data = hvf.info[not.const, ], span = loess.span)
    hvf.info$variance.expected[not.const] <- 10^fit$fitted
    hvf.info$variance.standardized <- Seurat:::SparseRowVarStd(mat = object, 
                                                               mu = hvf.info$mean, sd = sqrt(hvf.info$variance.expected), 
                                                               vmax = clip.max, display_progress = F)
    colnames(x = hvf.info) <- paste0("vst.", colnames(x = hvf.info))
  }
  
  return(hvf.info)
}

.extraMitoPctFn=function(inputData,organism,inputMTgenes=NULL,inputGeneName=NULL,redownload_files=T,recalculate_nUMI=T){
  
  rwNames=tolower(row.names(inputData))

  x=.extraMitoGenes(organism = organism,redownload_files=redownload_files)

  if(sum(grepl("\\.",rwNames))>0){
    rwNames=strsplit(rwNames,"\\.")
    rwNames=unlist(lapply(rwNames,function(x)x[1]))
  }

  if(is.null(inputGeneName)){
    tmpCols=c("ensembl_gene_id","entrezgene_id",'gene_name')
    
    slCounts=0
    slCol=""
    
    for(i in tmpCols){
      tmp=sum(rwNames %in% tolower(x[,i]))
      if(tmp>slCounts){
        slCounts=tmp
        slCol=i
      }
    }
    if(slCol!=""){
      inputMTgenes=tolower(x[,slCol])
    } else {
      inputMTgenes=""
    }
    
  } else {
    if(inputGeneName=="ensembl_gene_id"){
      inputMTgenes=tolower(x$ensembl_gene_id)
    } else {
      if(inputGeneName=="entrezgene_id"){
        x=.extraMitochondrialGenes()
        inputMTgenes=tolower(x$entrezgene_id)
      }
    }
  }
  
  print(paste("Number of MT genes in the dataset:",length(which(rwNames %in% inputMTgenes)),"/",sum(x$gene_biotype=="protein_coding")))
  
  inputMTgenes=inputMTgenes[inputMTgenes!=""]
  if(length(inputMTgenes)>0){
    tmpColSums=c()
    if(sum(colnames(colData(inputData))=='QC_Gene_total_count')==0|recalculate_nUMI){
      if(ncol(inputData)>10000){
        for(i in seq(1,ncol(inputData),10000)){
          tmpColSums=c(tmpColSums,apply(counts(inputData)[,i:min(i+10000-1,ncol(inputData))],2,sum))
        }
      } else {
        tmpColSums=apply(counts(inputData),2,sum)
      }
    } else {
      tmpColSums=inputData$QC_Gene_total_count
    }
    
    
    
    res=colSums(x = counts(inputData)[which(rwNames %in% inputMTgenes), , drop = FALSE])/tmpColSums
    res=res*100
  } else {
    res=rep(NA,ncol(inputData))
  }

  return(res)
}

.extrarRNAPctFn=function(inputData,organism,inputrRNAgenes=NULL,inputGeneName=NULL,redownload_files=T,recalculate_nUMI=T){
  
  rwNames=tolower(row.names(inputData))
  x=.extrarRNAGenes(organism = organism,redownload_files=redownload_files)
  if(sum(grepl("\\.",rwNames))>0){
    rwNames=strsplit(rwNames,"\\.")
    rwNames=unlist(lapply(rwNames,function(x)x[1]))
  }
  if(is.null(inputGeneName)){
    tmpCols=c("ensembl_gene_id","entrezgene_id",'gene_name')
    
    slCounts=0
    slCol=""
    
    for(i in tmpCols){
      tmp=sum(rwNames %in% tolower(x[,i]))
      if(tmp>slCounts){
        slCounts=tmp
        slCol=i
      }
    }
    if(slCol!=""){
      inputrRNAgenes=tolower(x[,slCol])
    } else {
      inputrRNAgenes=""
    }
    
  } else {
    if(inputGeneName=="ensembl_gene_id"){
      inputrRNAgenes=tolower(x$ensembl_gene_id)
    } else {
      if(inputGeneName=="entrezgene_id"){
        x=.extrarRNAGenes()
        inputrRNAgenes=tolower(x$entrezgene_id)
      }
    }
  }
  
  print(paste("Number of rRNA genes in the dataset:",length(which(rwNames %in% inputrRNAgenes)),"/",sum(x$gene_biotype=="rRNA")))
  
  inputrRNAgenes=inputrRNAgenes[inputrRNAgenes!=""]
  if(length(inputrRNAgenes)>0){
    tmpColSums=c()
    if(sum(colnames(colData(inputData))=='QC_Gene_total_count')==0|recalculate_nUMI){
      if(ncol(inputData)>10000){
        for(i in seq(1,ncol(inputData),10000)){
          tmpColSums=c(tmpColSums,apply(counts(inputData)[,i:min(i+10000-1,ncol(inputData))],2,sum))
        }
      } else {
        tmpColSums=apply(counts(inputData),2,sum)
      }
    } else {
      tmpColSums=inputData$QC_Gene_total_count
    }
    
    
    
    res=colSums(x = counts(inputData)[which(toupper(rwNames) %in% toupper(inputrRNAgenes)), , drop = FALSE])/tmpColSums
    res=res*100
  } else {
    res=rep(NA,ncol(inputData))
  }
  
  return(res)
}

.extraRPctFn=function(inputData,organism,inputRPgenes=NULL,inputGeneName=NULL,redownload_files=T,recalculate_nUMI=T){
  
  rwNames=tolower(row.names(inputData))
  x=.extraRibosomalProteinGenes(organism = organism,redownload_files=redownload_files)
  if(sum(grepl("\\.",rwNames))>0){
    rwNames=strsplit(rwNames,"\\.")
    rwNames=unlist(lapply(rwNames,function(x)x[1]))
  }
  if(is.null(inputGeneName)){
    tmpCols=c("ensembl_gene_id","entrezgene_id",'gene_name')
    
    slCounts=0
    slCol=""
    
    for(i in tmpCols){
      tmp=sum(rwNames %in% tolower(x[,i]))
      if(tmp>slCounts){
        slCounts=tmp
        slCol=i
      }
    }
    if(slCol!=""){
      inputRPgenes=tolower(x[,slCol])
    } else {
      inputRPgenes=""
    }
    
  } else {
    if(inputGeneName=="ensembl_gene_id"){
      inputRPgenes=tolower(x$ensembl_gene_id)
    } else {
      if(inputGeneName=="entrezgene_id"){
        x=.extraRibosomalProteinGenes()
        inputRPgenes=tolower(x$entrezgene_id)
      }
    }
  }
  
  print(paste("Number of ribosomal protein genes in the dataset:",length(which(rwNames %in% inputRPgenes)),"/",nrow(x)))
  
  inputRPgenes=inputRPgenes[inputRPgenes!=""]
  if(length(inputRPgenes)>0){
    tmpColSums=c()
    if(sum(colnames(colData(inputData))=='QC_Gene_total_count')==0|recalculate_nUMI){
      if(ncol(inputData)>10000){
        for(i in seq(1,ncol(inputData),10000)){
          tmpColSums=c(tmpColSums,apply(counts(inputData)[,i:min(i+10000-1,ncol(inputData))],2,sum))
        }
      } else {
        tmpColSums=apply(counts(inputData),2,sum)
      }
    } else {
      tmpColSums=inputData$QC_Gene_total_count
    }
    
    
    
    res=colSums(x = counts(inputData)[which(toupper(rwNames) %in% toupper(inputRPgenes)), , drop = FALSE])/tmpColSums
    res=res*100
  } else {
    res=rep(NA,ncol(inputData))
  }
  
  return(res)
}

.extraOXPHOSPctFn=function(inputData,organism,inputGeneset=NULL,inputGeneName_col=NULL,server=T,recalculate_nUMI=T,standardize_gene_names=T){
  
  rwNames=tolower(row.names(inputData))
  x=.extraOXPHOSGenes(organism = organism,server = server,inputGeneset=inputGeneset)
  if(sum(grepl("\\.",rwNames))>0){
    rwNames=strsplit(rwNames,"\\.")
    rwNames=unlist(lapply(rwNames,function(x)x[1]))
  }
  if(is.null(inputGeneName_col)){
    tmpCols=c("ensembl_gene_id",'gene_short_name',"entrezid","gene_name")
    
    slCounts=0
    slCol=""
    
    for(i in tmpCols){
      tmp=sum(rwNames %in% tolower(x[,i]))
      if(tmp>slCounts){
        slCounts=tmp
        slCol=i
      }
    }
    if(slCol!=""){
      inputOXPHOSgenes=tolower(x[,slCol])
    } else {
      inputOXPHOSgenes=""
    }
    
  } else {
    if(inputGeneName=="ensembl_gene_id"){
      inputOXPHOSgenes=tolower(x$ensembl_gene_id)
    } else {
      if(inputGeneName=="entrezgene_id"){
        x=.extraMitochondrialGenes()
        inputOXPHOSgenes=tolower(x$entrezgene_id)
      }
    }
  }
  
  
  if(!standardize_gene_names){
    inputOXPHOSgenes=unique(tolower(inputGeneset))
  }
  inputOXPHOSgenes=inputOXPHOSgenes[inputOXPHOSgenes!=""]
  if(is.null(inputGeneset)){
    print(paste("Number of OXPHOS genes in the dataset:",length(which(rwNames %in% inputOXPHOSgenes)),"/",length(inputOXPHOSgenes)))
  } else {
    print(paste("Number of genes in the dataset:",length(which(rwNames %in% inputOXPHOSgenes)),"/",length(inputOXPHOSgenes)))
  }
  
  if(length(inputOXPHOSgenes)>0){
    if(sum(colnames(colData(inputData))=="QC_Gene_total_count")==0|recalculate_nUMI){
      tmpColSums=c()
      if(ncol(inputData)>10000){
        for(i in seq(1,ncol(inputData),10000)){
          tmpColSums=c(tmpColSums,apply(counts(inputData)[,i:min(i+10000-1,ncol(inputData))],2,sum))
        }
      } else {
        tmpColSums=apply(counts(inputData),2,sum)
      }
    } else {
      tmpColSums=inputData$QC_Gene_total_count
    }
    
    
    
    res=colSums(x = counts(inputData)[which(rwNames %in% inputOXPHOSgenes), , drop = FALSE])/tmpColSums
    res=res*100
  } else {
    res=rep(NA,ncol(inputData))
  }
  
  return(res)
}

.extraSignaturePCAscoreFn=function(inputData,organism,inputGeneset=NULL,inputGeneName_col=NULL,server=T,align=T,standardize_gene_names=T,renormalize_data=T){
  
  if(sum(class(inputData) %in% c("matrix","Matrix","dgCMatrix"))==0){
    stop("inputData should be in (sparse) matrix format!")
  }
  
  rwNames=tolower(row.names(inputData))
  inputOXPHOSgenes=""
  if(standardize_gene_names){
    x=.extraOXPHOSGenes(organism = organism,server = server,inputGeneset=inputGeneset)
    if(sum(grepl("\\.",rwNames))>0){
      rwNames=strsplit(rwNames,"\\.")
      rwNames=unlist(lapply(rwNames,function(x)x[1]))
    }
    if(is.null(inputGeneName_col)){
      tmpCols=c("ensembl_gene_id",'gene_short_name',"entrezid","gene_name")
      
      slCounts=0
      slCol=""
      
      for(i in tmpCols){
        tmp=sum(rwNames %in% tolower(x[,i]))
        if(tmp>slCounts){
          slCounts=tmp
          slCol=i
        }
      }
      if(slCol!=""){
        inputOXPHOSgenes=tolower(x[,slCol])
      } else {
        inputOXPHOSgenes=""
      }
      
    } else {
      if(inputGeneName=="ensembl_gene_id"){
        inputOXPHOSgenes=tolower(x$ensembl_gene_id)
      } else {
        if(inputGeneName=="entrezgene_id"){
          x=.extraMitochondrialGenes()
          inputOXPHOSgenes=tolower(x$entrezgene_id)
        }
      }
    }
  }
  
  
  
  if(!standardize_gene_names){
    inputOXPHOSgenes=unique(tolower(inputGeneset))
  }
  inputOXPHOSgenes=inputOXPHOSgenes[inputOXPHOSgenes!=""]
  if(is.null(inputGeneset)){
    print(paste("Number of OXPHOS genes in the dataset:",length(which(rwNames %in% inputOXPHOSgenes)),"/",length(inputOXPHOSgenes)))
  } else {
    print(paste("Number of genes in the dataset:",length(which(rwNames %in% inputOXPHOSgenes)),"/",length(inputOXPHOSgenes)))
  }
  
  
  if(length(inputOXPHOSgenes)==1){
    return(inputData[inputOXPHOSgenes,])
  }
  
  if(length(inputOXPHOSgenes)>0){
    if(renormalize_data){
      inputData=Seurat::NormalizeData(inputData,verbose=F)
    }
    
    inputData=Seurat::ScaleData(inputData[rwNames %in% inputOXPHOSgenes,])
    res=Seurat::RunPCA(inputData,npcs=2)
    res=res@cell.embeddings[,1]
    if(align){
      mean_exp=colSums(inputData)
      if(cor(mean_exp,res)<0){
        res=res*(-1)
      }
    }
  } else {
    res=rep(NA,ncol(inputData))
  }
  
  return(res)
}

.extraIEGPctFn=function(inputData,organism,inputIEGgenes=NULL,inputGeneName=NULL,server=F,recalculate_nUMI=T){
  
  rwNames=tolower(row.names(inputData))
  x=.extraIEGGenes(organism = organism,server = server)
  if(sum(colnames(x)=="entrezid")>0){
    colnames(x)[colnames(x)=="entrezid"]="entrezgene_id"
  }
  if(sum(grepl("\\.",rwNames))>0){
    rwNames=strsplit(rwNames,"\\.")
    rwNames=unlist(lapply(rwNames,function(x)x[1]))
  }
  if(is.null(inputGeneName)){
    tmpCols=c("ensembl_gene_id","entrezgene_id",'gene_name')
    
    slCounts=0
    slCol=""
    
    for(i in tmpCols){
      tmp=sum(rwNames %in% tolower(x[,i]))
      if(tmp>slCounts){
        slCounts=tmp
        slCol=i
      }
    }
    if(slCol!=""){
      inputIEGgenes=tolower(x[,slCol])
    } else {
      inputIEGgenes=""
    }
    
  } else {
    if(inputGeneName=="ensembl_gene_id"){
      inputIEGgenes=tolower(x$ensembl_gene_id)
    } else {
      if(inputGeneName=="entrezgene_id"){
        x=.extraMitochondrialGenes()
        inputIEGgenes=tolower(x$entrezgene_id)
      }
    }
  }
  
  print(paste("Number of IEG genes in the dataset:",length(which(rwNames %in% inputIEGgenes)),"/",sum(x$gene_biotype=="protein_coding")))
  
  inputIEGgenes=inputIEGgenes[inputIEGgenes!=""]
  if(length(inputIEGgenes)>0){
    if(sum(colnames(colData(inputData))=="QC_Gene_total_count")==0|recalculate_nUMI){
      tmpColSums=c()
      if(ncol(inputData)>10000){
        for(i in seq(1,ncol(inputData),10000)){
          tmpColSums=c(tmpColSums,apply(counts(inputData)[,i:min(i+10000-1,ncol(inputData))],2,sum))
        }
      } else {
        tmpColSums=apply(counts(inputData),2,sum)
      }
    } else {
      tmpColSums=inputData$QC_Gene_total_count
    }
    
    
    
    res=colSums(x = counts(inputData)[which(rwNames %in% inputIEGgenes), , drop = FALSE])/tmpColSums
    res=res*100
  } else {
    rep(NA,ncol(inputData))
  }
  
  return(res)
}

.extraExprsSeurat=function(inputData,dataType="counts"){
  #dataType: counts, data (normalized values)
  inputData=GetAssayData(object = inputData, assay = DefaultAssay(object = inputData), slot = dataType)
  return(inputData)
}

.extraExport2CellDataSetFn=function(inputData,assay="counts"){
  
  fd <- new("AnnotatedDataFrame", data = as.data.frame(rowData(inputData)))
  pd <- new("AnnotatedDataFrame", data = as.data.frame(colData(inputData)))
  
  res <- newCellDataSet(assays(inputData)[[assay]], phenoData = pd, featureData = fd)
  return(res)
}

.extraExport2ExpressionSetFn=function(counts,pd=NULL,fd=NULL){
  
  if(is.null(fd)){
    fd =data.frame(geneName=row.names(counts),stringsAsFactors = F)
    row.names(fd)=row.names(counts)
  }
  if(is.null(pd)){
    pd =data.frame(sample=colnames(counts),stringsAsFactors = F)
    row.names(pd)=colnames(counts)
  }
  
  library(Matrix)
  
  counts <- as(counts, "sparseMatrix")  
  
  res=SingleCellExperiment(assays = list(counts = counts),colData = pd,rowData=fd)
  return(res)
}

.extraScaleFn=function(inputData){
  #scale function from Seurat
  #it clips the positive values at 10, but doesn't do the clipping on the negative ones!
  if(ncol(inputData)>70000){
    res=NULL
    if(nrow(inputData)>10000){
      for(i in seq(1,nrow(inputData),10000)){
        tmp=t(scale(t(inputData[i:min((i+10000-1),nrow(inputData)),])))
        tmp[which(tmp>10)]=10
        tmp[is.na(tmp)]=0
        res=rbind(res,tmp)
      }
    } else {
      tmp=t(scale(t(inputData)))
      tmp[which(tmp>10)]=10
      tmp[is.na(tmp)]=0
      res=tmp
    }
    
  } else {
    res=Seurat::ScaleData(inputData,do.center=T,do.scale=T,verbose=F)
  }
  
  return(res)
}

.extraPCApermFn=function (mat, max.rank = 100, ..., niters = 50, transposed = FALSE, BSPARAM = BiocSingular:::ExactParam(), BPPARAM = BiocParallel::SerialParam()){
  if (!transposed) {
    mat <- t(mat)
  }
  
  original <- pca(mat, rank = max.rank, ..., transposed = TRUE, 
                  BSPARAM = BSPARAM)
  original.s2 <- original$variance
  pcg.states <- PCAtools:::.setup_pcg_state(niters)
  permuted <- bpmapply(FUN = PCAtools:::.parallel_PA, seed = pcg.states$seeds[[1]], 
                       stream = pcg.states$streams[[1]], MoreArgs = list(mat = mat, 
                                                                         ..., max.rank = max.rank, BSPARAM = BSPARAM), BPPARAM = BPPARAM, 
                       SIMPLIFY = FALSE, USE.NAMES = FALSE)
  permutations <- do.call(cbind, permuted)
  prop <- rowMeans(permutations >= original.s2)
  prop=data.frame(PC=1:length(prop),permPval=prop)
  
  return(list(original = original, n = prop))
}

.extraExport2SeuratFn=function(inputData,project="scRNA"){
  require(scran)
  nbt=Seurat::CreateSeuratObject(counts=counts(inputData),
                                 project = "SeuratProject",
                                 assay = "RNA",
                                 min.cells = 0,
                                 min.features = 0,
                                 names.field = 1,
                                 names.delim = "-",
                                 meta.data = .pData(inputData))

  if(ncol(.fData(inputData = inputData))>0){
    nbt@assays$RNA@meta.data=.fData(inputData = inputData)
  }
  
  # Take all genes in > 3 cells, all cells with > 1k genes, use an expression threshold of 1
  # Cell type is encoded in the second _ field, will be stored in nbt@ident and also placed in the "orig.ident" field of object@data.info
  
  return(nbt)
}

.extraPCAfn=function (object, npcs = 50,seed.use = 12345,reduction.key = "PC_",weight.by.var=T, approx = TRUE,findElbowPoint=T, ...) {
  
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  object=.matrixExtraction(object)
  features.var <- apply(X = object, MARGIN = 1, FUN = var)
  object =object[features.var > 0,]
  total.variance=sum(features.var)
  
  if(all(object>=0)){
    print("Warning! variables are not scaled!")
  }
  
  npcs <- min(npcs, nrow(x = object) - 1)
  if (approx) {
    pca.results <- irlba::prcomp_irlba(t(x = object), n = npcs,...)
  } else {
    pca.results <- prcomp(x = t(object), rank. = npcs,...)
  }
  feature.loadings <- pca.results$rotation
  sdev <- pca.results$sdev
  if (weight.by.var) {
    cell.embeddings <- pca.results$x %*% diag(pca.results$sdev[1:npcs]^2)
  }
  else {
    cell.embeddings <- pca.results$x
  }
  
  rownames(x = feature.loadings) <- rownames(x = object)
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:npcs)
  rownames(x = cell.embeddings) <- colnames(x = object)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  
  elbowPoint=NULL
  if(findElbowPoint){
    elbowPoint=PCAtools::findElbowPoint(sdev^2/total.variance*100)
  }
  
  reduction.data <- list(embeddings = cell.embeddings, loadings = feature.loadings,
                         stdev = sdev,total.variance = total.variance,var.explained=sdev^2/total.variance*100,elbowPoint=elbowPoint)
  
  return(reduction.data)
}

.extraJackStrawFn=function(inputPCA,inputData, 
                           num.replicate = 100, prop.freq = 0.01, maxit = 1000) {
  require(future.apply)
  my.lapply <- future.apply::future_lapply
  my.sapply <- future.apply::future_sapply
  
  dims <- ncol(inputPCA$embeddings)
  
  
  loadings <- inputPCA$loadings
  
  reduc.features <- rownames(x = loadings)
  if (length(x = reduc.features) < 3) {
    stop("Too few features")
  }
  #Num. of genes to permute at each iteration
  if (length(x = reduc.features) * prop.freq < 3) {
    warning("Number of variable genes given ", prop.freq, 
            " as the prop.freq is low. Consider including more variable genes and/or increasing prop.freq. ", 
            "Continuing with 3 genes in every random sampling.")
  }
  
  
  #randomly shuffles prop.freq% of rows and performs PCA. Returns the gene loadings
  fake.vals.raw <- my.lapply(X = 1:num.replicate, FUN = Seurat:::JackRandom, 
                             scaled.data = .matrixExtraction(inputData), prop.use = prop.freq, r1.use = 1, 
                             r2.use = dims, rev.pca = F, weight.by.var = F, 
                             maxit = maxit)
  
  fake.vals <- sapply(X = 1:dims, FUN = function(x) {
    return(as.numeric(x = unlist(x = lapply(X = 1:num.replicate, 
                                            FUN = function(y) {
                                              return(fake.vals.raw[[y]][, x])
                                            }))))
  })
  fake.vals <- as.matrix(x = fake.vals)
  jackStraw.empP <- as.matrix(my.sapply(X = 1:dims, FUN = function(x) {
    return(unlist(x = lapply(X = abs(loadings[, x]), FUN = Seurat:::EmpiricalP, 
                             nullval = abs(fake.vals[, x]))))
  }))
  colnames(x = jackStraw.empP) <- paste0("PC", 1:ncol(x = jackStraw.empP))
  
  score.thresh = 1e-05
  
  pAll <- jackStraw.empP
  pAll <- as.data.frame(pAll)
  pAll$Contig <- rownames(x = pAll)
  score.df <- NULL
  for (i in 1:ncol(jackStraw.empP)) {
    pc.score <- suppressWarnings(prop.test(x = c(length(x = which(x = pAll[, 
                                                                           i] <= score.thresh)), floor(x = nrow(x = pAll) * 
                                                                                                         score.thresh)), n = c(nrow(pAll), nrow(pAll)))$p.val)
    if (length(x = which(x = pAll[, i] <= score.thresh)) == 
        0) {
      pc.score <- 1
    }
    if (is.null(x = score.df)) {
      score.df <- data.frame(PC = paste0("PC", i), Score = pc.score)
    }
    else {
      score.df <- rbind(score.df, data.frame(PC = paste0("PC", 
                                                         i), Score = pc.score))
    }
  }
  score.df$PC <- 1:ncol(pAll)
  score.df <- as.matrix(score.df)
  
  
  return(list(empP=jackStraw.empP,PCscores=score.df,fake.vals=fake.vals))
}

.netMagicFn=function(inputPCAembeddings,inputExpData,n.adaptiveKernel=5,nn.eps=0,nPropIter=5){
  #nPropIter: number of propagation iterations
  #input: PCA res
  n.cells <- nrow(inputPCAembeddings)
  nn.ranked <- RANN::nn2(data = inputPCAembeddings, k = n.adaptiveKernel*6, eps = nn.eps)
  
  dists=nn.ranked$nn.dists
  affinities=t(apply(dists,1,function(x) exp((-1)*(x/x[n.adaptiveKernel])^2)))
  affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
  affinities2=affinities-affinitiesThr
  affinities[affinities2<=0]=0
  
  affCounts=apply(affinities,2,function(x) sum(x==0))
  affinities=affinities[,-which(affCounts==nrow(affinities))]
  nn.ranked$nn.idx=nn.ranked$nn.idx[,-which(affCounts==nrow(affinities))]
  j <- as.numeric(t(nn.ranked$nn.idx))
  i <- ((1:length(j)) - 1)%/%ncol(affinities) + 1
  k=as.numeric(t(affinities))
  
  
  graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(inputPCAembeddings), nrow(inputPCAembeddings)))
  graph=(graph+t(graph))/2
  graph <- Matrix::Diagonal(x = 1 / (rowSums(graph))) %*% graph
  
  if(class(graph)!="dgCMatrix"){
    graph=as(graph,"dgCMatrix")
  }
  
  rownames(graph) <- rownames(inputPCAembeddings)
  colnames(graph) <- rownames(inputPCAembeddings)
  
  diffRes=data.frame(iteration=0,rsq=0)
  tmp=1
  Dmat=t(inputExpData)
  imputedData=Dmat
  slNew=T
  for(i in 1:nPropIter){
    print(i)
    if(all(tmp==1)){
      tmp2=graph
    } else {
      tmp2=tmp%*%graph
    }
    imputedData2=tmp2 %*% Dmat
    score1=sum((imputedData2-imputedData)^2)
    score2=sum((imputedData- mean(imputedData))^2)
    
    score=1-score1/score2
    diffRes=rbind(diffRes,data.frame(iteration=i,rsq=score))
    
    if(median(score)>0.95&slNew&i>min(2,nPropIter-1)){
      slM=tmp2
      slImputed=imputedData2
      slRound=i
      slNew=F
    }
    tmp=tmp2
    imputedData=imputedData2
  }
  diffRes=diffRes[-1,]
  
  return(list(affinityGraph=graph,propagatedGraph=slM,numIter=slRound,imputedData=imputedData,R_sq=diffRes))
}

.myNetCombination=function(inputNet,order=T){
  #This function returns a name for each interaction in the input network
  
  inputNet$Fnode=as.character(inputNet$Fnode)
  inputNet$Snode=as.character(inputNet$Snode)
  if(order){
    tst1=apply(inputNet[,c("Fnode","Snode")],1,min)
    tst2=apply(inputNet[,c("Fnode","Snode")],1,max)
    res=paste(tst1,tst2,sep=".")
  } else {
    res=paste(inputNet$Fnode,inputNet$Snode,sep=".")
  }
  return(res)
}

.netGenie3=function (exprMat, inputTFs=NULL,weightThr=0.001, nParts = 10, nCores=5,rand.seed=12345) {
  #input: log Normalized data
  genesSplit <- suppressWarnings(split(sort(rownames(exprMat)), 
                                       1:nParts))
  weightMatrices <- list()
  for (i in 1:length(genesSplit)) {
    set.seed(seed = rand.seed)
    weightMatrix <- GENIE3::GENIE3(exprMat, regulators = inputTFs, 
                                   nCores = nCores, targets = genesSplit[[i]])
    weightMatrices[[i]] <- weightMatrix
  }
  linkList_list <- list()
  for (i in 1:length(genesSplit)) {
    weightMatrix <- weightMatrices[[i]]
    linkList_list[[i]] <- GENIE3::getLinkList(weightMatrix, 
                                              threshold = weightThr)
  }
  rm(weightMatrices)
  linkList <- do.call(rbind, linkList_list)
  colnames(linkList) <- c("TF", "Target", "weight")
  linkList <- linkList[order(linkList[, "weight"], decreasing = TRUE), ]
  
  mapper=data.frame(name=row.names(exprMat),indx=1:nrow(exprMat))
  mapper=mapper[match(linkList$TF,mapper$name),]
  linkList$TF=mapper$indx
  
  mapper=data.frame(name=row.names(exprMat),indx=1:nrow(exprMat))
  mapper=mapper[match(linkList$Target,mapper$name),]
  linkList$Target=mapper$indx
  
  graph <- sparseMatrix(i = linkList$TF, j = linkList$Target, x = linkList$weight,dims = c(nrow(exprMat), nrow(exprMat)))
  if(is.null(inputTFs)){
    graph=(graph+t(graph))/2
  }
  
  row.names(graph)=row.names(exprMat)
  colnames(graph)=row.names(exprMat)
  
  return(graph)
}

.reductionPCAFn=function(inputData,npcs=50,sigPCmethod=NULL,weight.by.var=T,approx=T){
  #method for identification of the signaificant PCs: horn, jackstraw, and elbow. elbow is done automatically for each run!
  
  pcaRes=.extraPCAfn(object = inputData,npcs = npcs,weight.by.var=weight.by.var,approx = approx)
  
  
  horn=NULL
  jackstraw=NULL
  
  if(sum(sigPCmethod=="horn")>0){
    horn <- .extraPCApermFn(inputData)
    horn=horn$n
  }
  
  
  if(sum(sigPCmethod=="jackstraw")>0){
    jackstraw=.extraJackStrawFn(inputPCA = pcaRes,inputData)
  }
  
  return(list(pca=pcaRes,horn=horn,jackstraw=jackstraw,elbowPoint=pcaRes$pca$elbowPoint))
}

.netEuclidean=function (inputPCAscores, input.is.distance.matrix = FALSE, k.param = 20, compute.SNN = TRUE, 
                        jaccard.indx.thr = 1/15, nn.eps = 0, verbose = TRUE) {
  #inputPCAscores: PCA results
  
  n.cells <- nrow(x = inputPCAscores)
  
  if (!input.is.distance.matrix) {
    if (verbose) {
      message("Computing nearest neighbor graph")
    }
    nn.ranked <- RANN::nn2(data = inputPCAscores, k = k.param, eps = nn.eps)
    nn.ranked <- nn.ranked$nn.idx
  }
  else {
    if (verbose) {
      message("Building SNN based on a provided distance matrix")
    }
    knn.mat <- matrix(data = 0, ncol = k.param, nrow = n.cells)
    for (i in 1:n.cells) {
      knn.mat[i, ] <- order(inputPCAscores[i, ])[1:k.param]
    }
    nn.ranked <- knn.mat[, 1:k.param]
  }
  j <- as.numeric(t(nn.ranked))
  i <- ((1:length(j)) - 1)%/%k.param + 1
  nn.matrix <- sparseMatrix(i = i, j = j, x = 1,dims = c(nrow(inputPCAscores), nrow(inputPCAscores)))
  rownames(nn.matrix) <- rownames(inputPCAscores)
  colnames(nn.matrix) <- rownames(inputPCAscores)
  neighbor.graphs <- list(nn = nn.matrix)
  if (compute.SNN) {
    if (verbose) {
      message("Computing SNN")
    }
    snn.matrix <- Seurat:::ComputeSNN(nn_ranked = nn.ranked, prune = jaccard.indx.thr)
    rownames(x = snn.matrix) <- rownames(inputPCAscores)
    colnames(x = snn.matrix) <- rownames(inputPCAscores)
    snn.matrix <- as.Graph(x = snn.matrix)
    
    neighbor.graphs[["snn"]] <- snn.matrix
  }
  return(neighbor.graphs)
}

.netFindClusters=function (inputGraph, algorithm = 1,resolution = 0.8,group.singletons = TRUE,modularity.fxn = 1,verbose=F) {
  #algorithm:
  #1: louvain
  #2: Louvain algorithm with multilevel refinement
  #3: SLM
  
  if (all(as.matrix(inputGraph[1:min(1000,nrow(inputGraph)),])==round(as.matrix(inputGraph[1:min(1000,nrow(inputGraph)),])))){
    warning("Input is not an SNN graph")
  }
  
  if(class(inputGraph)=="matrix"){
    inputGraph=as(inputGraph, "dgTMatrix")
  }
  
  if(class(inputGraph)=="dgTMatrix"){
    inputGraph=Seurat::as.Graph(inputGraph)
    DefaultAssay(object = inputGraph) <- "RNA"
  }
  
  slAlgorithm="Louvain"
  if(algorithm==2){
    slAlgorithm="Louvain algorithm with multilevel refinement"
  }
  if(algorithm==3){
    slAlgorithm="SLM"
  }
  print(paste("clustering algorithm:",slAlgorithm))
  if (nbrOfWorkers() > 1) {
    clustering.results <- future.apply::future_lapply(X = resolution, FUN = function(r) {
      if (algorithm %in% c(1:3)) {
        ids <- Seurat:::RunModularityClustering(SNN = inputGraph, 
                                                modularity = modularity.fxn, resolution = r, 
                                                algorithm = algorithm, n.start = 10, n.iter = 10, 
                                                random.seed = 0, print.output = F, 
                                                temp.file.location = NULL, edge.file.name = NULL)
      } else {
        stop("algorithm not recognised, please specify as an integer or string")
      }
      names(x = ids) <- colnames(inputGraph)
      ids <- Seurat:::GroupSingletons(ids = ids, SNN = inputGraph, verbose = F)
      results <- list(factor(x = ids))
      names(x = results) <- paste0("res.", r)
      return(results)
    })
    clustering.results <- as.data.frame(x = clustering.results)
  } else {
    clustering.results <- data.frame(row.names = colnames(inputGraph))
    for (r in resolution) {
      if (algorithm %in% c(1:3)) {
        ids <- Seurat:::RunModularityClustering(SNN = inputGraph, 
                                                modularity = modularity.fxn, resolution = r, 
                                                algorithm = algorithm, n.start = 10, n.iter = 10, 
                                                random.seed = 0, print.output = FALSE, 
                                                temp.file.location = NULL, edge.file.name = NULL)
      } else {
        stop("algorithm not recognised, please specify as an integer or string")
      }
      names(x = ids) <- colnames(inputGraph)
      ids <- Seurat:::GroupSingletons(ids = ids, SNN = inputGraph, group.singletons = group.singletons, 
                                      verbose = F)
      clustering.results[, paste0("res.", r)] <- factor(x = ids)
    }
  }
  return(clustering.results)
}

.findMarkersfunction=function (inputLogNormData, inputCountData, cellGroup1 = NULL, cellGroup2 = NULL, 
                               logfc.threshold = 0.25, test.use = "wilcox", latent.vars.df = NULL,
                               min.pct = 0.1, min.diff.pct = -Inf, verbose = TRUE, only.upregulated = FALSE, 
                               max.cells.per.ident = Inf, 
                               min.cells.feature = 3, min.cells.group = 3, ...) {
  "%||%" <- function(a, b) {
    if (!is.null(a)) a else b
  }
  
  data.slot ="data"
  
  if(!is.null(cellClusters1)){
    ident.1 <- cellGroup1
  } else {
    stop("provide the list of cells in the first group")
  }
  
  
  if (is.null(cellGroup2)) {
    ident.2 <- setdiff(colnames(inputLogNormData), ident.1)
  } else {
    ident.2=cellGroup2
  }
  inputLogNormData=inputLogNormData[,colnames(inputLogNormData) %in% c(ident.1,ident.2)]
  if (!is.null(latent.vars.df)) {
    if(class(latent.vars.df)!="data.frame"){
      stop("latent.vars is supposed to be a data.frame with row.names as the cell names")
    } else {
      latent.vars=latent.vars.df[row.names(latent.vars.df) %in% colnames(inputLogNormData),]
      if(nrow(latent.vars)!=ncol(inputLogNormData)){
        stop("some cells are not represented in the latent variable data.frame")
      }
    }
    
  }
  de.results <- Seurat::FindMarkers(object = inputLogNormData, counts = inputCountData, cells.1 = ident.1, cells.2 = ident.2, 
                                    features = features, reduction = reduction, logfc.threshold = logfc.threshold, 
                                    test.use = test.use, min.pct = min.pct, min.diff.pct = min.diff.pct, 
                                    verbose = TRUE, only.pos = only.up, max.cells.per.ident = max.cells.per.ident, 
                                    random.seed = 12345, latent.vars = latent.vars, 
                                    min.cells.feature = min.cells.feature, min.cells.group = min.cells.group, 
                                    pseudocount.use = 1)
  return(de.results)
}

.netAUCellFn=function(inputExp,inputGmt,numAUCgenes=round(0.2*nrow(inputExp))){
  #geneset scores are based on the rank statistics. Normalization method shouldn't in theory affect the results
  
  require(AUCell)
  require(GSEABase)
  cells_rankings <- AUCell_buildRankings(inputExp)
  
  #load("~/Documents/data/early_late_hcASD/mapper/ensembl_entrez.rda")
  #geneSets =.readGeneSetData("~/Documents/data/early_late_hcASD/GeneSets/rASDgeneGroups.gmt",geneLabelsInDataset = unique(gsIds$entrezgene_id),nmin=25,nmax=500)
  #gsIds=gsIds[match(colnames(geneSets),gsIds$entrezgene_id),]
  #colnames(geneSets)=gsIds$ensembl_gene_id
  
  cells_AUC <- AUCell_calcAUC(geneSets=inputGmt, cells_rankings,aucMaxRank=numAUCgenes)
  
  scores=cells_AUC@assays@data$AUC
  #cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
  
  return(scores)
  
}

.extraQCAnnoAdderFn=function(inputExpSet,organism,server=F,redownload_files=T,ncores=5){
  
  {
    x=counts(inputExpSet)
    
    if(ncol(x)>10000){
      
      if(ncores==1){
        nGene=c()
        nUMI=c()
        for(i in seq(1,ncol(x),10000)){
          nGene=c(nGene,apply(x[,i:min(i+10000-1,ncol(x))],2,function(x) sum(x>0)))
          nUMI=c(nUMI,apply(x[,i:min(i+10000-1,ncol(x))],2,function(x) sum(x)))
        }
      } else {
        nGene=c()
        nUMI=c()
        nGeneUMI=parallel::mclapply(seq(1,ncol(x),10000),function(i,x){
          tmp_nGene=apply(x[,i:min(i+10000-1,ncol(x))],2,function(x) sum(x>0))
          tmp_nUMI=apply(x[,i:min(i+10000-1,ncol(x))],2,function(x) sum(x))
          return(list(nGene=tmp_nGene,nUMI=tmp_nUMI,i=i))
        },x=x,mc.cores = ncores)
        ilist=c()
        for(i in 1:length(nGeneUMI)){
          nGene=c(nGene,nGeneUMI[[i]]$nGene)
          nUMI=c(nUMI,nGeneUMI[[i]]$nUMI)
          ilist=c(ilist,nGeneUMI[[i]]$i)
        }
        
        if(any(ilist!=seq(1,ncol(x),10000),na.rm = F)){
          stop("Error in calculation of nGene/nUMI!")
        }
      }
      
    } else {
      nGene=apply(x,2,function(x) sum(x>0))
      nUMI=apply(x,2,function(x) sum(x))
    }
    
    inputExpSet$QC_Gene_total_count=nUMI
    inputExpSet$QC_Gene_unique_count=nGene
  }
  
  tmpExpData=counts(inputExpSet)
  
  inputExpSet$QC_MT.pct=.extraMitoPctFn(inputData = inputExpSet,organism = organism,redownload_files=redownload_files,recalculate_nUMI = F)
  inputExpSet$QC_IEG.pct=.extraIEGPctFn(inputData = inputExpSet,organism = organism,server = server,recalculate_nUMI = F)
  inputExpSet$QC_OXPHOS.pct=.extraOXPHOSPctFn(inputData = inputExpSet,organism = organism,server=server,recalculate_nUMI = F)
  inputExpSet$QC_rRNA.pct=.extrarRNAPctFn(inputData = inputExpSet,organism = organism,redownload_files=redownload_files,recalculate_nUMI = F)
  inputExpSet$QC_RibosomalProtein.pct=.extraRPctFn(inputData = inputExpSet,organism = organism,redownload_files=redownload_files,recalculate_nUMI = F)
  
  
  
  tmpMedian=c()
  #tmpMean=c()
  tmpDetectionRate=c()
  #tmpVar=c()
  if(ncores==1|T){
    for(i in seq(1,nrow(tmpExpData),3000)){
      tmpMedian=c(tmpMedian,apply(tmpExpData[i:min(i+3000-1,nrow(tmpExpData)),],1,median))
      #tmpMean=c(tmpMean,apply(tmpExpData[i:min(i+3000-1,nrow(tmpExpData)),],1,mean))
      #tmpVar=c(tmpVar,apply(tmpExpData[i:min(i+3000-1,nrow(tmpExpData)),],1,var))
      tmpDetectionRate=c(tmpDetectionRate,apply(tmpExpData[i:min(i+3000-1,nrow(tmpExpData)),],1,function(x) sum(x>0)/length(x)))
    }
    
  } else {
    MedianDetectionRate=parallel::mclapply(seq(1,nrow(tmpExpData),3000),function(i,tmpExpData){
      tmpMedian=apply(tmpExpData[i:min(i+3000-1,nrow(tmpExpData)),],1,median)
      #tmpMean=c(tmpMean,apply(tmpExpData[i:min(i+3000-1,nrow(tmpExpData)),],1,mean))
      #tmpVar=c(tmpVar,apply(tmpExpData[i:min(i+3000-1,nrow(tmpExpData)),],1,var))
      tmpDetectionRate=apply(tmpExpData[i:min(i+3000-1,nrow(tmpExpData)),],1,function(x) sum(x>0)/length(x))
      return(list(median=tmpMedian,detectionRate=tmpDetectionRate,i=i))
    },tmpExpData=tmpExpData,mc.cores =ncores)
  }
  
  
  names(tmpMedian)=row.names(tmpExpData)
  rowData(inputExpSet)$QC_medianExp=round(tmpMedian,3)
  tmp=data.frame(gene=names(tmpMedian),exp=tmpMedian,stringsAsFactors = F)
  tmp=tmp[order(tmp$exp,decreasing = T),]
  tmp=tmp[which(tmp$exp>=tmp$exp[50]),]
  rowData(inputExpSet)$QC_top50_expressed="No"
  rowData(inputExpSet)$QC_top50_expressed[row.names(inputExpSet) %in% tmp$gene]="Yes"
  tmp=inputExpSet[row.names(inputExpSet) %in% tmp$gene,]
  inputExpSet$QC_top50_pct=colSums(assays(tmp)[["counts"]])/inputExpSet$QC_Gene_total_count*100#  colSums(tmpExpData)*100
  
  
  
  if(sum(colnames(rowData(inputExpSet))=="gene_biotype")==1){
    tmp=inputExpSet[which(rowData(inputExpSet)$gene_biotype=="lincRNA"),]
    
    tmp=colSums(counts(tmp))/inputExpSet$QC_Gene_total_count
    inputExpSet$QC_lncRNA_pct=tmp*100
  }
  
  #rowData(inputExpSet)$QC_meanExp=round(tmpMean,3)
  
  rowData(inputExpSet)$QC_detectionRate=round(tmpDetectionRate,3)
  
  #rowData(inputExpSet)$QC_varExp=round(tmpVar,3)
  
  #tst=cor(tmpExpData,method = "spearman")
  
  #colData(inputExpSet)$QC_raw_data_mean_cor=apply(tst,1,mean)
  #colData(inputExpSet)$QC_raw_data_median_cor=apply(tst,1,median)
  
  
  return(inputExpSet)
}



.extraHumanGeneAnnoAdderFn=function(inputGeneNames=NULL,server=T,redownload_files=T){
  #require(EnsDb.Hsapiens.v75)
  library(ensembldb)
  library(EnsDb.Hsapiens.v86)

  require(EnsDb.Hsapiens.v86)
  
  if(!dir.exists("../data/gazestani_serverFiles")){
    dir.create("../data/gazestani_serverFiles",recursive = T)
  }
  
  gns <- as.data.frame(ensembldb::genes(EnsDb.Hsapiens.v86))
  gns$gene_short_name=gns$gene_name
  gns$symbol=toupper(gns$symbol)
  gns$ensembl_gene_id=row.names(gns)
  tst=unique(as.character(gns$gene_biotype))
  tst=setdiff(tst,c("protein_coding","lincRNA"))
  tst=factor(as.character(gns$gene_biotype),levels=c("protein_coding","lincRNA",tst))
  gns=gns[order(tst,decreasing = F),]
  if(!is.null(inputGeneNames)){
    rwNames=toupper(inputGeneNames)
    psCols=c("gene_short_name","ensembl_gene_id")
    slCounts=0
    slCol=""
    if(sum(grepl("\\.",rwNames)&grepl("^ENS",rwNames))>0){
      rwNames=strsplit(rwNames,"\\.")
      rwNames=unlist(lapply(rwNames,function(x)x[1]))
    }
    {
      
      if(!file.exists("../data/gazestani_serverFiles/human_map_to_ensembl.rda")){
        stop("Error! file ../data/gazestani_serverFiles/human_map_to_ensembl.rda is missing")
      }
      
      load("../data/gazestani_serverFiles/human_map_to_ensembl.rda")
    }
    
    map_to_ensmbl$source=toupper(map_to_ensmbl$source)
    
    if(!file.exists("../data/gazestani_serverFiles/human_mapping_hg19.rda")){
      stop("Error! file ../data/gazestani_serverFiles/human_mapping_hg19.rda is missing")
    }
    
    load("../data/gazestani_serverFiles/human_mapping_hg19.rda")
    human_hg19$source=toupper(human_hg19$source)
    
    if(sum(toupper(rwNames) %in% human_hg19$source) > sum(toupper(rwNames) %in% map_to_ensmbl$source)){
      map_to_ensmbl=rbind(human_hg19,map_to_ensmbl)
      map_to_ensmbl=map_to_ensmbl[!duplicated(map_to_ensmbl$source),]
    } else {
      map_to_ensmbl=rbind(map_to_ensmbl,human_hg19)
      map_to_ensmbl=map_to_ensmbl[!duplicated(map_to_ensmbl$source),]
    }
    
    c_map=data.frame(source=gns$gene_short_name,target=gns$ensembl_gene_id,stringsAsFactors = F)
    c_map=rbind(c_map,data.frame(source=gns$ensembl_gene_id,target=gns$ensembl_gene_id,stringsAsFactors = F))
    c_map=c_map[!is.na(c_map$ensembl_gene_id),]
    map_to_ensmbl=rbind(map_to_ensmbl,c_map)
    map_to_ensmbl=map_to_ensmbl[!duplicated(map_to_ensmbl$source),]
    map_to_ensmbl=merge(map_to_ensmbl,data.frame(source=toupper(rwNames),stringsAsFactors = F),by="source",all.y=T)
    gns=merge(gns,map_to_ensmbl,by.x="ensembl_gene_id",by.y="target",all.y=T)
    
    gns=gns[match(rwNames,gns$source),]
    row.names(gns)=inputGeneNames
    gns$gene_id=inputGeneNames
    gns=gns[,-which(colnames(gns) %in% c("source","target"))]
  }
  
  return(gns)
}
.extraMouseGeneAnnoAdderFn=function(inputGeneNames=NULL,server=T){
  require(EnsDb.Mmusculus.v79)
  
  if(!dir.exists("../data/gazestani_serverFiles")){
    dir.create("../data/gazestani_serverFiles",recursive = T)
  }
  
  gns <- as.data.frame(genes(EnsDb.Mmusculus.v79))
  gns$symbol=toupper(gns$symbol)
  gns$gene_short_name=gns$gene_name
  gns$ensembl_gene_id=row.names(gns)
  tst=unique(as.character(gns$gene_biotype))
  tst=setdiff(tst,c("protein_coding","lincRNA"))
  tst=factor(as.character(gns$gene_biotype),levels=c("protein_coding","lincRNA",tst))
  gns=gns[order(tst,decreasing = F),]
  
  if(!is.null(inputGeneNames)){
    rwNames=toupper(inputGeneNames)
    
    if(sum(grepl("\\.",rwNames)&grepl("^ENS",rwNames))>0){
      rwNames=strsplit(rwNames,"\\.")
      rwNames=unlist(lapply(rwNames,function(x)x[1]))
    }
    {
      #library(googleCloudStorageR)
      if(!file.exists("../data/gazestani_serverFiles/mouse_map_to_ensembl.rda")){
        stop("Error! file ../data/gazestani_serverFiles/mouse_map_to_ensembl.rda is missing")
      }
      
      load("../data/gazestani_serverFiles/mouse_map_to_ensembl.rda")
      
    }
    map_to_ensmbl$source=toupper(map_to_ensmbl$source)
    
    
    c_map=data.frame(source=gns$gene_short_name,target=gns$ensembl_gene_id,stringsAsFactors = F)
    c_map=rbind(c_map,data.frame(source=gns$ensembl_gene_id,target=gns$ensembl_gene_id,stringsAsFactors = F))
    c_map=c_map[!is.na(c_map$ensembl_gene_id),]
    map_to_ensmbl=rbind(map_to_ensmbl,c_map)
    map_to_ensmbl=map_to_ensmbl[!duplicated(map_to_ensmbl$source),]
    map_to_ensmbl=merge(map_to_ensmbl,data.frame(source=toupper(rwNames),stringsAsFactors = F),by="source",all.y=T)
    
    
    gns=merge(gns,map_to_ensmbl,by.x="ensembl_gene_id",by.y="target",all.y=T)
    gns=gns[match(rwNames,gns$source),]
    row.names(gns)=inputGeneNames
    gns=gns[,-which(colnames(gns) %in% c("source","target"))]
  }
  
  
  return(gns)
}

.extraMacaqueGeneAnnoAdderFn=function(inputGeneNames=NULL,server=T){
  library(AnnotationHub)
  hub <- AnnotationHub()
  avail.resources=query(hub, c("fascicularis"))
  avail.resources=data.frame(id=avail.resources$ah_id,desc=avail.resources$description,genome=avail.resources$genome,title=avail.resources$title,stringsAsFactors = F)
  avail.resources=avail.resources[grepl("Ensembl",avail.resources$desc),]
  ensdb <- hub[[avail.resources$id[nrow(avail.resources)]]]
  
  if(!dir.exists("../data/gazestani_serverFiles")){
    dir.create("../data/gazestani_serverFiles",recursive = T)
  }
  
  gns <- as.data.frame(genes(ensdb))
  gns$symbol=toupper(gns$symbol)
  gns$gene_short_name=gns$gene_name
  gns$ensembl_gene_id=row.names(gns)
  
  if(!is.null(inputGeneNames)){
    rwNames=toupper(inputGeneNames)
    
    if(sum(grepl("\\.",rwNames)&grepl("^ENS",rwNames))>0){
      rwNames=strsplit(rwNames,"\\.")
      rwNames=unlist(lapply(rwNames,function(x)x[1]))
    }
    
    sl_cols=c('gene_id','gene_name',"symbol")
    map_to_ensmble=NULL
    for(icol in sl_cols){
      map_to_ensmble=rbind(map_to_ensmble,data.frame(source=gns[,icol],ensembl_gene_id=gns$gene_id,stringsAsFactors = F))
    }
    map_to_ensmble=map_to_ensmble[!duplicated(map_to_ensmble$source),]
    map_to_ensmble$source=toupper(map_to_ensmble$source)
    
    map_to_ensmble=merge(map_to_ensmble,data.frame(source=toupper(rwNames),stringsAsFactors = F),by="source",all.y=T)
    
    
    gns=merge(map_to_ensmble,gns,by.x="ensembl_gene_id",by.y="gene_id",all.y=T)
    gns=gns[match(rwNames,gns$source),]
    row.names(gns)=inputGeneNames
    if(sum(colnames(gns) %in% c("source","target"))>0){
      gns=gns[,-which(colnames(gns) %in% c("source","target"))]
    }
    
  }
  
  gns$gene_short_name=gns$symbol
  return(gns)
}

.plotMeanHGVexplorerFn=function(inputFD,vstScaleRank=2000,vstNotScaleRank=2000,logSizeFactor="auto",logNorm="auto",scTransform="auto",sqrtNorm="auto"){
  
  meanHGVfn=function(inputFD,slCol,slThr,slSign){
    tmpDf=quantile(inputFD$QC_meanExp,seq(0,1,0.02))
    
    inputFD$slVal=cut(inputFD$QC_meanExp,breaks = unique(c(0,quantile(inputFD$QC_meanExp,seq(0,1,0.02)))))
    tmp=strsplit(as.character(inputFD$slVal),",")
    tmp=unlist(lapply(tmp,function(x) x[1]))
    tmp=(gsub("\\(","",tmp))
    tmp=as.numeric(gsub("\\[","",tmp))
    inputFD$slVal=tmp
    tmp=unique(tmp)
    tmp=tmp[order(tmp)]
    
    inputFD$slVal=factor(inputFD$slVal,levels=tmp)
    if(slSign=="gt"){
      tmp=aggregate(inputFD[,slCol]~slVal,data=inputFD,function(x) mean(x[!is.na(x)]>slThr))
    } else {
      tmp=aggregate(inputFD[,slCol]~slVal,data=inputFD,function(x) mean(x[!is.na(x)]<slThr))
    }
    
    colnames(tmp)[2]="Prob"
    
    p=ggplot(tmp,aes(slVal,Prob*100))+geom_bar(stat="identity")+.th+theme(axis.text.x = element_text(angle = 30,vjust=1,hjust=1))+xlab("Mean Exp.")+ylab("1st 2000 genes")
  }
  
  count=0
  if(sum(colnames(inputFD)=="hvg_vst_rank_scaled")>0){
    count=count+1
    inputFD[is.na(inputFD[,"hvg_vst_rank_scaled"]),"hvg_vst_rank_scaled"]=0
    pVstScaled=meanHGVfn(inputFD=inputFD,slCol="hvg_vst_rank_scaled",slThr=vstScaleRank,slSign = "lt")+ggtitle("Seurat VST scaled counts")
    if(count==1){
      p=pVstScaled
    } else {
      p=p+pVstScaled
    }
    
  }
  
  if(sum(colnames(inputFD)=="hvg_vst_rank_notScaled")>0){
    count=count+1
    inputFD[is.na(inputFD[,"hvg_vst_rank_notScaled"]),"hvg_vst_rank_notScaled"]=0
    pVstNotScaled=meanHGVfn(inputFD=inputFD,slCol="hvg_vst_rank_notScaled",slThr=vstNotScaleRank,slSign = "lt")+ggtitle("Seurat VST not scaled counts")
    if(count==1){
      p=pVstNotScaled
    } else {
      p=p+pVstNotScaled
    }
  }
  
  if(sum(colnames(inputFD)=="hvg_mvp.dispersion.std_lognorm")>0){
    inputFD[is.na(inputFD[,"hvg_mvp.dispersion.std_lognorm"]),"hvg_mvp.dispersion.std_lognorm"]=0
    
    if(logNorm=="auto"){
      logNorm=inputFD[,"hvg_mvp.dispersion.std_lognorm"]
      logNorm=logNorm[order(logNorm,decreasing = T)]
      logNorm=logNorm[2000]
    }
    count=count+1
    pLogNorm=meanHGVfn(inputFD=inputFD,slCol="hvg_mvp.dispersion.std_lognorm",slThr=logNorm,slSign = "gt")+ggtitle("Log normalized")
    if(count==1){
      p=pLogNorm
    } else {
      p=p+pLogNorm
    }
    
  }
  
  if(sum(colnames(inputFD)=="hvg_mvp.dispersion.std_sqrtnorm")>0){
    count=count+1
    inputFD[is.na(inputFD[,"hvg_mvp.dispersion.std_sqrtnorm"]),"hvg_mvp.dispersion.std_sqrtnorm"]=0
    
    if(sqrtNorm=="auto"){
      sqrtNorm=inputFD[,"hvg_mvp.dispersion.std_sqrtnorm"]
      sqrtNorm=sqrtNorm[order(sqrtNorm,decreasing = T)]
      sqrtNorm=sqrtNorm[2000]
    }
    pSqrt=meanHGVfn(inputFD=inputFD,slCol="hvg_mvp.dispersion.std_sqrtnorm",slThr=sqrtNorm,slSign = "gt")+ggtitle("Sqrt normalized")
    if(count==1){
      p=pSqrt
    } else {
      p=p+pSqrt
    }
  }
  
  if(sum(colnames(inputFD)=="hvg_mvp.dispersion.std_logSizeFactor")>0){
    count=count+1
    inputFD[is.na(inputFD[,"hvg_mvp.dispersion.std_logSizeFactor"]),"hvg_mvp.dispersion.std_logSizeFactor"]=0
    
    if(logSizeFactor=="auto"){
      logSizeFactor=inputFD[,"hvg_mvp.dispersion.std_logSizeFactor"]
      logSizeFactor=logSizeFactor[order(logSizeFactor,decreasing = T)]
      logSizeFactor=logSizeFactor[2000]
    }
    pSizeFactor=meanHGVfn(inputFD=inputFD,slCol="hvg_mvp.dispersion.std_logSizeFactor",slThr=logSizeFactor,slSign = "gt")+ggtitle("Size factor normalized")
    if(count==1){
      p=pSizeFactor
    } else {
      p=p+pSizeFactor
    }
  }
  
  if(sum(colnames(inputFD)=="scTransform_resid_var")){
    count=count+1
    inputFD[is.na(inputFD[,"scTransform_resid_var"]),"scTransform_resid_var"]=0
    
    if(scTransform=="auto"){
      scTransform=inputFD[,"scTransform_resid_var"]
      scTransform=scTransform[order(scTransform,decreasing = T)]
      scTransform=scTransform[2000]
    }
    pSCtransform=meanHGVfn(inputFD=inputFD,slCol="scTransform_resid_var",slThr=scTransform,slSign = "gt")+ggtitle("ScTransform (NB regression)")
    if(count==1){
      p=pSCtransform
    } else {
      p=p+pSCtransform
    }
  }
  
  myNormPlotsFn=function(inputData){
    
    counter=0
    if(sum(colnames(inputData)=="scTransform_resid_var")>0){
      counter=counter+1
      tmpP1=ggplot(inputData,aes(scTransform_resid_var))+geom_density()+scale_x_log10()+xlab("Gene variance")+ggtitle("ScTransform (NB regression)")
      if(counter==1){
        p=tmpP1
      }else{
        p=p+tmpP1
      }
    }
    
    if(sum(colnames(inputData)=="hvg_vst_score_scaled")>0){
      counter=counter+1
      tmpP2=ggplot(inputData,aes(hvg_vst_score_scaled))+geom_density()+scale_x_log10()+xlab("Gene variance")+ggtitle("Seurat VST scaled counts")
      if(counter==1){
        p=tmpP2
      }else{
        p=p+tmpP2
      }
    }
    
    if(sum(colnames(inputData)=="hvg_vst_score_notScaled")>0){
      counter=counter+1
      tmpP3=ggplot(inputData,aes(hvg_vst_score_notScaled))+geom_density()+scale_x_log10()+xlab("Gene variance")+ggtitle("Seurat VST not scaled counts")
      if(counter==1){
        p=tmpP3
      }else{
        p=p+tmpP3
      }
    }
    
    if(sum(colnames(inputData)=="hvg_mvp.dispersion.std_lognorm")>0){
      counter=counter+1
      tmpP4=ggplot(inputData,aes(hvg_mvp.dispersion.std_lognorm))+geom_density()+xlab("Relative gene variance")+ggtitle("Log normalized")
      if(counter==1){
        p=tmpP4
      }else{
        p=p+tmpP4
      }
    }
    
    if(sum(colnames(inputData)=="hvg_mvp.dispersion.std_logSizeFactor")>0){
      counter=counter+1
      tmpP5=ggplot(inputData,aes(hvg_mvp.dispersion.std_logSizeFactor))+geom_density()+xlab("Relative gene variance")+ggtitle("Size factor normalized")
      if(counter==1){
        p=tmpP5
      }else{
        p=p+tmpP5
      }
    }
    
    if(sum(colnames(inputData)=="hvg_mvp.dispersion.std_sqrtnorm")>0){
      counter=counter+1
      tmpP6=ggplot(inputData,aes(hvg_mvp.dispersion.std_sqrtnorm))+geom_density()+xlab("Relative gene variance")+ggtitle("Sqrt normalized")
      if(counter==1){
        p=tmpP6
      }else{
        p=p+tmpP6
      }
    }
    return(p)
  }
  p2=myNormPlotsFn(inputFD)
  
  return(list(normPlots=p2,HVGplots=p))
}

#inputPCAembeddings=harmony_embeddings
#umap.method = "uwot"; n.neighbors = 30L; 
#n.components = 2L; metric = "cosine"; n.epochs = NULL; learning.rate = 1; 
#min.dist = 0.3; spread = 1; set.op.mix.ratio = 1; local.connectivity = 1L; 
#repulsion.strength = 1; negative.sample.rate = 5; a = NULL; 
#b = NULL; uwot.sgd = FALSE; seed.use = 42; metric.kwds = NULL; 
#angular.rp.forest = FALSE; reduction.key = "UMAP_"; verbose = TRUE

.reductionUMAPFn=function (inputPCAembeddings,testPCAembeddings=NULL, umap.method = "umap-learn", n.neighbors = 30L, 
                           n.components = 2L, metric = "cosine", n.epochs = NULL, learning.rate = 1, 
                           min.dist = 0.3, spread = 1, set.op.mix.ratio = 1, local.connectivity = 1L, 
                           repulsion.strength = 1, negative.sample.rate = 5, a = NULL, 
                           b = NULL, uwot.sgd = FALSE, seed.use = 42, metric.kwds = NULL, 
                           angular.rp.forest = FALSE, reduction.key = "UMAP_", verbose = TRUE, 
                           ...) {
  #adapted from Seurat
  require(Seurat)
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  umap.model=""
  if (umap.method != "umap-learn") {
    warning("The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric", 
            "\nTo use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'", 
            call. = FALSE, immediate. = TRUE)
  }
  umap.output <- switch(EXPR = umap.method, `umap-learn` = {
    require(reticulate)
    if (! reticulate::py_module_available(module = "umap")) {
      stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
    }
    if (!is.null(x = seed.use)) {
      py_set_seed(seed = seed.use)
    }
    if (typeof(x = n.epochs) == "double") {
      n.epochs <- as.integer(x = n.epochs)
    }
    umap_import <- import(module = "umap", delay_load = TRUE)
    umap <- umap_import$UMAP(n_neighbors = as.integer(x = n.neighbors), 
                             n_components = as.integer(x = n.components), metric = metric, 
                             n_epochs = n.epochs, learning_rate = learning.rate, 
                             min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio, 
                             local_connectivity = local.connectivity, repulsion_strength = repulsion.strength, 
                             negative_sample_rate = negative.sample.rate, a = a, 
                             b = b, metric_kwds = metric.kwds, angular_rp_forest = angular.rp.forest, 
                             verbose = verbose)
    #umap$fit_transform(as.matrix(x = inputPCAembeddings))
    
    tst=umap$fit(as.matrix(x = inputPCAembeddings))
    
    umap.output=tst$transform(inputPCAembeddings)
    
    colnames(umap.output) <- paste0(reduction.key, 1:ncol(umap.output))
    rownames(umap.output) <- rownames(inputPCAembeddings)
    
    res=c()
    if(!is.null(testPCAembeddings)){
      umap.test=tst$transform(testPCAembeddings)
      colnames(umap.test) <- paste0(reduction.key, 1:ncol(umap.test))
      rownames(umap.test) <- rownames(testPCAembeddings)
      res=list(embedding=umap.output,model=tst,test=umap.test)
    } else {
      res=list(embedding=umap.output,model=tst)
    }
    
    res
    
  }, uwot = {
    if (metric == "correlation") {
      warning("UWOT does not implement the correlation metric, using cosine instead", 
              call. = FALSE, immediate. = TRUE)
      metric <- "cosine"
    }
    resModel=uwot::umap(X = inputPCAembeddings, n_threads = 1, n_neighbors = as.integer(n.neighbors), 
               n_components = as.integer(n.components), metric = metric, 
               n_epochs = n.epochs, learning_rate = learning.rate, 
               min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio, 
               local_connectivity = local.connectivity, repulsion_strength = repulsion.strength, 
               negative_sample_rate = negative.sample.rate, a = a, 
               b = b, fast_sgd = uwot.sgd, verbose = verbose,ret_model = T)
    
    umap.model=resModel
    umap.output=resModel[["embedding"]]
    
    colnames(umap.output) <- paste0(reduction.key, 1:ncol(umap.output))
    rownames(umap.output) <- rownames(inputPCAembeddings)
    
    res=c()
    if(!is.null(testPCAembeddings)){
      umap.test=uwot::umap_transform(testPCAembeddings,resModel)
      res=list(embedding=umap.output,model=resModel,test=umap.test)
    } else {
      res=list(embedding=umap.output,model=umap.model)
    }
    
    res
  }, stop("Unknown umap method: ", umap.method, call. = FALSE))
  
  
  
  return(umap.output)
}

.plot2D=function (resUMAP,clusters, cells = NULL, cols = NULL, 
                  pt.size = NULL, group.by = NULL, 
                  shapeCol=NULL, order = NULL, label = FALSE, label.size = 4, 
                  repel = FALSE, cells.highlight = NULL, cols.highlight = "#DE2D26", 
                  sizes.highlight = 1, na.value = "grey50", ncol = NULL) 
{
  "%||%" <- function(a, b) {
    if (!is.null(a)) a else b
  }
  
  
  cells <- cells %||% row.names(resUMAP)
  data <- resUMAP[cells,]
  data <- as.data.frame(x = data)
  dims <- colnames(resUMAP)
  
  data$ident <- clusters[row.names(resUMAP) %in% cells , drop = FALSE]
  if (!is.factor(x = data$ident)) {
    data$ident <- factor(data$ident)
  }
  
  if (!is.null(x = shapeCol)) {
    data$shapeCol <- clusterRes[,shapeCol]
  }
  
  plot <- Seurat:::SingleDimPlot(data = data[, c(dims, "ident", 
                                                 shapeCol)], dims = dims, col.by = "ident", cols = cols, 
                                 pt.size = pt.size, shape.by = shapeCol, order = order, 
                                 label = FALSE, cells.highlight = cells.highlight, 
                                 cols.highlight = cols.highlight, sizes.highlight = sizes.highlight, 
                                 na.value = na.value)
  if (label) {
    plot <- LabelClusters(plot = plot, id = x, repel = repel, 
                          size = label.size, split.by = split.by)
  }
  
  
  return(plot)
}

.plotFeatures=function (inputReducedDim,inputFeatureData,inputClusters, dims = c(1, 2), cells = NULL,blend = FALSE, blend.threshold = 0.5,  cols = if (blend) {c("lightgrey", "#ff0000", "#00ff00")} else {c("lightgrey", "blue")}, pt.size = NULL, order = FALSE, min.feature.cutoff = NA, max.feature.cutoff = NA, 
                        reduction = NULL, split.by = NULL, shape.by = NULL, slot = "data", 
                        label = FALSE, label.size = 4, 
                        repel = FALSE, ncol = NULL, coord.fixed = FALSE, by.col = TRUE, 
                        sort.cell = FALSE, combine = TRUE) 
{
  
  require(ggplot2)
  "%||%" <- function(a, b) {
    if (!is.null(a)) a else b
  }
  
  inputReducedDim=as.data.frame(inputReducedDim[,dims])
  
  inputReducedDim$ident=factor(inputClusters)
  
  no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
                    axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold", size = 14, margin = margin(r = 7)))
  
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  if (blend && length(x = features) != 2) {
    stop("Blending feature plots only works with two features")
  }
  if (blend) {
    default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
    cols <- switch(EXPR = as.character(x = length(x = cols)), 
                   `0` = {
                     warning("No colors provided, using default colors", 
                             call. = FALSE, immediate. = TRUE)
                     default.colors
                   }, `1` = {
                     warning("Only one color provided, assuming specified is double-negative and augmenting with default colors", 
                             call. = FALSE, immediate. = TRUE)
                     c(cols, default.colors[2:3])
                   }, `2` = {
                     warning("Only two colors provided, assuming specified are for features and agumenting with '", 
                             default.colors[1], "' for double-negatives", 
                             call. = FALSE, immediate. = TRUE)
                     c(default.colors[1], cols)
                   }, `3` = cols, {
                     warning("More than three colors provided, using only first three", 
                             call. = FALSE, immediate. = TRUE)
                     cols[1:3]
                   })
  }
  if (blend && length(x = cols) != 3) {
    stop("Blending feature plots only works with three colors; first one for negative cells")
  }
  
  cells <- cells %||% row.names(inputReducedDim)
  data <- cbind(inputReducedDim,inputFeatureData)
  
  features <- colnames(inputFeatureData)
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = min(data[, 
                                                           feature]), no = cutoff))
  }, cutoff = min.feature.cutoff, feature = features)
  
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = max(data[, 
                                                           feature]), no = cutoff))
  }, cutoff = max.feature.cutoff, feature = features)
  
  check.lengths <- unique(x = vapply(X = list(features, min.cutoff, 
                                              max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols, 
                                                                            ]$maxcolors, no = length(x = cols))
  
  data[, 4:ncol(x = data)] <- sapply(X = 4:ncol(x = data), 
                                     FUN = function(index) {
                                       data.feature <- as.vector(x = data[, index])
                                       min.use <- Seurat:::SetQuantile(cutoff = min.cutoff[index - 3], data.feature)
                                       max.use <- Seurat:::SetQuantile(cutoff = max.cutoff[index - 3], data.feature)
                                       data.feature[data.feature < min.use] <- min.use
                                       data.feature[data.feature > max.use] <- max.use
                                       if (brewer.gran == 2) {
                                         return(data.feature)
                                       }
                                       data.cut <- if (all(data.feature == 0)) {
                                         0
                                       }
                                       else {
                                         as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature), 
                                                                          breaks = brewer.gran)))
                                       }
                                       return(data.cut)
                                     })
  colnames(x = data)[4:ncol(x = data)] <- features
  rownames(x = data) <- cells
  data$split <- Seurat:::RandomName()
  
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  if (!is.null(x = shape.by)) {
    if(class(shape.by)!="vector"){
      print("shapes are supposed to be in the format of a vector")
    } else {
      data$shape.by <- shape.by
    }
    
  }
  plots <- vector(mode = "list", length = ifelse(test = blend, 
                                                 yes = 4, no = length(x = features) * length(x = levels(x = data$split))))
  xlims <- c(floor(x = min(inputReducedDim[,1])), ceiling(x = max(inputReducedDim[,1])))
  ylims <- c(floor(min(inputReducedDim[,2])), ceiling(x = max(inputReducedDim[,2])))
  if (blend) {
    ncol <- 4
    color.matrix <- Seurat:::BlendMatrix(two.colors = cols[2:3], col.threshold = blend.threshold, 
                                         negative.color = cols[1])
    cols <- cols[2:3]
    colors <- list(color.matrix[, 1], color.matrix[1, ], 
                   as.vector(x = color.matrix))
  }
  for (i in 1:length(levels(data$split))) {
    ident <- levels(data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident, 
                      , drop = FALSE]
    if (blend) {
      features <- features[1:2]
      no.expression <- features[colMeans(x = data.plot[, 
                                                       features]) == 0]
      if (length(x = no.expression) != 0) {
        stop("The following features have no value: ", 
             paste(no.expression, collapse = ", "), call. = FALSE)
      }
      data.plot <- cbind(data.plot[, c(dims, "ident")], 
                         BlendExpression(data = data.plot[, features[1:2]]))
      features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
    }
    for (j in 1:length(x = features)) {
      feature <- features[j]
      if (blend) {
        cols.use <- as.numeric(x = as.character(x = data.plot[, 
                                                              feature])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      } else {
        cols.use <- NULL
      }
      data.single <- data.plot[, c(dims, "ident", feature, 
                                   shape.by)]
      if (sort.cell) {
        data.single <- data.single[order(data.single[, 
                                                     feature]), ]
      }
      plot <- SingleDimPlot(data = data.single, dims = dims, 
                            col.by = feature, order = order, pt.size = pt.size, 
                            cols = cols.use, shape.by = shape.by, label = FALSE) + 
        scale_x_continuous(limits = xlims) + scale_y_continuous(limits = ylims) + 
        theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
      if (label) {
        plot <- LabelClusters(plot = plot, id = "ident", 
                              repel = repel, size = label.size)
      }
      if (length(x = levels(x = data$split)) > 1) {
        plot <- plot + theme(panel.border = element_rect(fill = NA, 
                                                         colour = "black"))
        plot <- plot + if (i == 1) {
          labs(title = feature)
        }
        else {
          labs(title = NULL)
        }
        if (j == length(x = features) && !blend) {
          suppressMessages(expr = plot <- plot + scale_y_continuous(sec.axis = dup_axis(name = ident), 
                                                                    limits = ylims) + no.right)
        }
        if (j != 1) {
          plot <- plot + theme(axis.line.y = element_blank(), 
                               axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
                               axis.title.y.left = element_blank())
        }
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(axis.line.x = element_blank(), 
                               axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
                               axis.title.x = element_blank())
        }
      }
      else {
        plot <- plot + labs(title = feature)
      }
      if (!blend) {
        plot <- plot + guides(color = NULL)
        cols.grad <- cols
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        }
        else if (length(x = cols) > 1) {
          unique.feature.exp <- unique(data.plot[, feature])
          if (length(unique.feature.exp) == 1) {
            warning("All cells have the same value (", 
                    unique.feature.exp, ") of ", feature, ".")
            if (unique.feature.exp == 0) {
              cols.grad <- cols[1]
            }
            else {
              cols.grad <- cols
            }
          }
          plot <- suppressMessages(expr = plot + scale_color_gradientn(colors = cols.grad, 
                                                                       guide = "colorbar"))
        }
      }
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      plot <- plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  if (blend) {
    blend.legend <- BlendMap(color.matrix = color.matrix)
    for (ii in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(x = plots, 
                                              values = list(blend.legend + scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) > 
                                                                                                                                  1, yes = levels(x = data$split)[ii], no = "")), 
                                                                                              expand = c(0, 0)) + labs(x = features[1], y = features[2], 
                                                                                                                       title = if (ii == 1) {
                                                                                                                         paste("Color threshold:", blend.threshold)
                                                                                                                       } else {
                                                                                                                         NULL
                                                                                                                       }) + no.right), after = 4 * ii - 1))
    }
  }
  plots <- Filter(f = Negate(f = is.null), x = plots)
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = features) == 1) {
      ncol <- 1
    }
    if (length(x = features) > 6) {
      ncol <- 3
    }
    if (length(x = features) > 9) {
      ncol <- 4
    }
  }
  ncol <- ifelse(test = is.null(x = split.by) || blend, yes = ncol, 
                 no = length(x = features))
  legend <- if (blend) {
    "none"
  }
  else {
    split.by %iff% "none"
  }
  if (combine) {
    if (by.col && !is.null(x = split.by) && !blend) {
      plots <- lapply(X = plots, FUN = function(x) {
        return(suppressMessages(expr = x + theme_cowplot() + 
                                  ggtitle("") + scale_y_continuous(sec.axis = dup_axis(name = ""), 
                                                                   limits = ylims) + no.right))
      })
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) + 
                 1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(expr = plots[[i]] + 
                                         scale_y_continuous(sec.axis = dup_axis(name = features[[idx]]), 
                                                            limits = ylims) + no.right)
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots)%%length(x = features) == 
                      1)) {
        plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) + 
          theme(plot.title = element_text(hjust = 0.5))
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) + 
            theme(plot.title = element_text(hjust = 0.5))
          idx <- idx + 1
        }
        ncol <- nsplits
        nrow <- 1
      }
      else {
        nrow <- split.by %iff% length(x = levels(x = data$split))
      }
      plots <- plots[c(do.call(what = rbind, args = split(x = 1:length(x = plots), 
                                                          f = ceiling(x = seq_along(along.with = 1:length(x = plots))/length(x = features)))))]
      plots <- wrap_plots(plots, ncol = ncol, nrow = nrow)
      if (!is.null(x = legend) && legend == "none") {
        plots <- plots & NoLegend()
      }
    }
    else {
      plots <- wrap_plots(plots, ncol = ncol, nrow = split.by %iff% 
                            length(x = levels(x = data$split)))
    }
    if (!is.null(x = legend) && legend == "none") {
      plots <- plots & NoLegend()
    }
  }
  return(plots)
}

.genesetMouse=function(genesetName=NULL,geneId="ensemble_gene_id",minGeneSetSize=15,maxGenesetSize=600,server=F){
  #genesetName: mm_GO, mm_pathway, mm_coexpression,, mm_TFperturb, mm_TF, mm_miRNA
  
  if(is.null(genesetName)){
    stop("Select one of the following genesets: enrichR, mm_GO, mm_pathway, mm_coexpression,, mm_TFperturb, mm_TF, mm_miRNA")
  }
  
  if(genesetName=="mm_coexpression"){
    d1 <- scan("http://ge-lab.org/gskb/2-MousePath/MousePath_Co-expression_gmt.gmt", what="", sep="\n", skip=1)
    geneset <- strsplit(d1, "\t")
    names(geneset) <- sapply(geneset, '[[', 1)
    names(geneset)=gsub(" ","_",names(geneset))
  } else if(genesetName=="mm_TFperturb"){
    d1 <- scan("https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=TF_Perturbations_Followed_by_Expression", what="", sep="\n", skip=1)
    geneset <- strsplit(d1, "\t")
    names(geneset) <- sapply(geneset, '[[', 1)
    names(geneset)=gsub(" ","_",names(geneset))
  } else if(genesetName=="enrichR"){
    if(server){
      d1=scan("~/Documents/data/SC/genesets/genesets_enrichR.gmt", what="", sep="\n")
    } else {
      d1=scan("~/Documents/data/SC/genesets/genesets_enrichR.gmt", what="", sep="\n")
      
    }
    
    geneset <- strsplit(d1, "\t")
    names(geneset) <- sapply(geneset, '[[', 1)
    names(geneset)=gsub(" ","_",names(geneset))
  } else {
    data(list=genesetName)
    geneset=get(genesetName)
    rm(list=ls(pattern = genesetName))
  }
  
  
  gns=.extraMouseGeneAnnoAdderFn()
  if(geneId=="symbol"){
    res=geneset
  } else if(geneId=="ensemble_gene_id") {
    gns$symbol=toupper(gns$symbol)
    res=list()
    for(i in 1:length(geneset)){
      res=c(res,list(gns$ensembl_gene_id[which(gns$symbol %in% geneset[[i]])]))
    }
    names(res)=names(geneset)
  } else {
    stop("unrecognized gene id!")
  }
  toKeep=c()
  for(i in 1:length(res)){
    if(length(res[[i]])>minGeneSetSize&length(res[[i]])<maxGenesetSize){
      toKeep=c(toKeep,i)
    }
  }
  res=res[toKeep]
  return(res)
}


.myPvalFn=function(nn_data,nn.thr,anno_data,coord_data){
  
  
  tst_nodes=nn_data$nn.idx
  tst_dist=nn_data$nn.dists
  
  tst_nodes=apply(tst_nodes,2,function(x) anno_data$value[x])
  
  tst_dist=apply(nn_data$nn.dists,1,min)
  
  for(k in unique(as.character(anno_data$value))){
    pvalDf=NULL
    for(i in 0:ncol(nn_data$nn.idx)){
      tmp=sum(dhyper(i:ncol(nn_data$nn.dists),sum(anno_data$value==k,na.rm = T),sum(anno_data$value!=k,na.rm = T),ncol(nn_data$nn.dists)))
      pvalDf=rbind(pvalDf,data.frame(count=i,pval=tmp))
    }
    
    tmp_count=apply(tst_nodes,1,function(x) sum(x==k,na.rm = T))
    tmp_count=pvalDf$pval[match(tmp_count,pvalDf$count)]
    tmp_count[tst_dist>nn.thr]=NA
    coord_data$new=tmp_count
    colnames(coord_data)[colnames(coord_data)=="new"]=k
  }
  
  
  
  
  return(coord_data)
}

.myPvalPermFn=function(X,nn_data,nn.thr,anno_data,coord_data){
  anno_data$value=sample(as.character(anno_data$value))
  
  resEnrichment=.myPvalFn(nn_data = nn_data,nn.thr = nn.thr,anno_data = anno_data,coord_data = coord_data)
  
  return(resEnrichment)
}

.my2dPlot=function(inputSeuratObject,attribute,reduction_name,perm=F,k=0.3){
  require(Seurat)
  require(ggplot2)
  require(RANN)
  
  if(sum(colnames(inputSeuratObject@meta.data)==attribute)==0){
    stop("Specified attribute couldn't be found in this Seurat object")
  }
  
  anno=data.frame(sample=row.names(inputSeuratObject@meta.data),value=as.factor(inputSeuratObject@meta.data[,attribute]))
  if(reduction_name=="umap"){
    if(!is.null(inputSeuratObject@reductions$umap)){
      reduction=inputSeuratObject@reductions$umap
    } else {
      stop("No UMAP data was found in this Seurat object")
    }
    
  } else {
    if(!is.null(inputSeuratObject@reductions$tsne)){
      reduction_name="tsne"
      reduction=inputSeuratObject@reductions$tsne
    } else {
      if(!is.null(inputSeuratObject@reductions$umap)){
        print("Using the UMAP data")
        reduction_name="umap"
        reduction=inputSeuratObject@reductions$umap
      } else {
        stop("No tSNE data was found in this Seurat object")
      }
    }
  }
  
  xval=reduction@cell.embeddings[,1]
  yval=reduction@cell.embeddings[,2]
  
  rescoord=expand.grid(seq(min(xval),max(xval),length.out = 100),seq(min(yval),max(yval),length.out = 100))
  colnames(rescoord)=c("x","y")
  rescoord$status=paste0("point_",1:nrow(rescoord))
  row.names(rescoord)=rescoord$status
  
  
  nn <- RANN::nn2(query = rescoord[,-3],data = data.frame(x=xval,y=yval,stringsAsFactors = F), k = 40, eps = 0)
  
  
  nn_thr=quantile(as.numeric(nn$nn.dists),k)
  
  pval_res=.myPvalFn(nn_data=nn,nn.thr=nn_thr,anno_data=anno,coord_data=rescoord)
  
  if(perm){
    require(parallel)
    set.seed(12345)
    res_perm=mclapply(1:1000,.myPvalPermFn,nn_data=nn,nn.thr=nn_thr,anno_data=anno,coord_data=rescoord,mc.cores = 5)
    res_perm=do.call("rbind",res_perm)
    x=split(res_perm,f=paste0(res_perm$x,"_",res_perm$y))
    
    xname=unlist(lapply(x,function(x) paste0(x[1,1],"_",x[1,2])))
    o=match(paste0(pval_res$x,"_",pval_res$y),xname)
    x=x[o]
    
    emp_pval=mclapply(1:nrow(pval_res),function(y) unlist(lapply(4:ncol(pval_res),function(z) sum(x[[y]][,z]<=pval_res[y,z]))),mc.cores = 5)
    emp_pval=do.call("rbind",emp_pval)
    emp_pval=emp_pval/1000
    
    the_pval=pval_res[,4:ncol(pval_res)]
    emp_pval[the_pval>emp_pval]=the_pval[the_pval>emp_pval]
    pval_res[4:ncol(pval_res)]=emp_pval
    
  }
  
  if(sum(colnames(pval_res)=="")>0){
    colnames(pval_res)[which(colnames(pval_res)=="")]="X"
  }
  pval_res=reshape::melt(pval_res,id.vars=colnames(pval_res)[1:3])
  colnames(pval_res)[colnames(pval_res)=="value"]="pval"
  pval_res$score=-log10(pval_res$pval)
  
  pval_res$score[pval_res$score<(-log10(0.01))]=0
  pval_res$score[pval_res$score<0]=0
  pval_res$score[is.na(pval_res$score)]=(-1)
  p=ggplot(pval_res, aes(x, y,fill = score)) +geom_raster()+
    scale_fill_gradientn(colours = c("white","black","black","yellow","red"),values=c(-1,0,0.5,-log10(0.01),max(pval_res$score))/max(pval_res$score))+
    theme_bw()+theme(panel.grid = element_blank())+facet_wrap(~variable)+xlab(toupper(paste0(reduction_name,"_1")))+ylab(toupper(paste0(reduction_name,"_2")))
  
  return(figure=p)
}

.my2dPlot2=function(inputPCA,batchCol,reductionCols,perm=F,k=0.3){
  require(Seurat)
  require(ggplot2)
  require(RANN)
  
  if(sum(colnames(inputPCA)==batchCol)==0){
    stop("Specified batch column couldn't be found")
  }
  
  if(sum(colnames(inputPCA) %in% reductionCols)!=2){
    stop("Specified reduction columns couldn't be found ")
  }
  
  if(sum(is.na(inputPCA[,batchCol]))>0){
    inputPCA=inputPCA[!is.na(inputPCA[,batchCol]),]
  }
  
  anno=data.frame(sample=row.names(inputPCA),value=as.factor(inputPCA[,batchCol]))
  
  yval=inputPCA[,reductionCols]
  xval=yval[,1]
  yval=yval[,2]
  
  rescoord=expand.grid(seq(min(xval),max(xval),length.out = 100),seq(min(yval),max(yval),length.out = 100))
  colnames(rescoord)=c("x","y")
  rescoord$status=paste0("point_",1:nrow(rescoord))
  row.names(rescoord)=rescoord$status
  
  
  nn <- RANN::nn2(query = rescoord[,-3],data = data.frame(x=xval,y=yval,stringsAsFactors = F), k = 40, eps = 0)
  
  
  nn_thr=quantile(as.numeric(nn$nn.dists),k)
  
  pval_res=.myPvalFn(nn_data=nn,nn.thr=nn_thr,anno_data=anno,coord_data=rescoord)
  
  if(perm){
    require(parallel)
    set.seed(12345)
    res_perm=mclapply(1:1000,.myPvalPermFn,nn_data=nn,nn.thr=nn_thr,anno_data=anno,coord_data=rescoord,mc.cores = 5)
    res_perm=do.call("rbind",res_perm)
    x=split(res_perm,f=paste0(res_perm$x,"_",res_perm$y))
    
    xname=unlist(lapply(x,function(x) paste0(x[1,1],"_",x[1,2])))
    o=match(paste0(pval_res$x,"_",pval_res$y),xname)
    x=x[o]
    
    emp_pval=mclapply(1:nrow(pval_res),function(y) unlist(lapply(4:ncol(pval_res),function(z) sum(x[[y]][,z]<=pval_res[y,z]))),mc.cores = 5)
    emp_pval=do.call("rbind",emp_pval)
    emp_pval=emp_pval/1000
    
    the_pval=pval_res[,4:ncol(pval_res)]
    emp_pval[the_pval>emp_pval]=the_pval[the_pval>emp_pval]
    pval_res[4:ncol(pval_res)]=emp_pval
    
  }
  
  pval_res=reshape::melt(pval_res,id.vars=colnames(pval_res)[c(1,2,3)])
  colnames(pval_res)[colnames(pval_res)=="value"]="pval"
  pval_res$score=-log10(pval_res$pval)
  pval_res$variable=factor(as.character(pval_res$variable),levels=levels(anno$value))
  pval_res$score[pval_res$score<(-log10(0.01))]=0
  pval_res$score[pval_res$score<0]=0
  pval_res$score[is.na(pval_res$score)]=0
  pval_res2=pval_res
  pval_res$score[which(pval_res$score==0)]=NA
  pval_res=pval_res[!is.na(pval_res$score),]
  anno2=inputPCA[,c(reductionCols,batchCol)]
  colnames(anno2)=c("x","y","variable")
  p=ggplot(pval_res, aes(x, y))+geom_point(data=anno2,aes(x,y),size=0.3) +geom_tile(aes(fill = score))+
    scale_fill_gradientn(colours = c("yellow","red"),values=c(-log10(0.01),max(pval_res$score))/max(pval_res$score))+
    theme_bw()+theme(panel.grid = element_blank())+facet_wrap(~variable)+xlab(toupper(reductionCols[1]))+ylab(toupper(reductionCols[2]))
  
  return(figure=p)
}

.my2dPlot_continuous=function(inputPCA,continuous_batch_values,reductionCols){
  require(Seurat)
  require(ggplot2)
  require(RANN)
  
  if(sum(colnames(inputPCA)==continuous_batch_values)==0){
    stop("Specified batch column couldn't be found")
  }
  
  if(sum(colnames(inputPCA) %in% reductionCols)!=2){
    stop("Specified reduction columns couldn't be found ")
  }
  
  if(sum(is.na(inputPCA[,continuous_batch_values]))>0){
    inputPCA=inputPCA[!is.na(inputPCA[,continuous_batch_values]),]
  }
  
  yval=inputPCA[,reductionCols]
  xval=yval[,1]
  yval=yval[,2]
  
  rescoord=expand.grid(seq(min(xval),max(xval),length.out = 100),seq(min(yval),max(yval),length.out = 100))
  xrange=abs(seq(min(xval),max(xval),length.out = 100)[2]-seq(min(xval),max(xval),length.out = 100)[1])/1.5
  yrange=abs((seq(min(yval),max(yval),length.out = 100)[2]-seq(min(yval),max(yval),length.out = 100)[1]))/1.5
  colnames(rescoord)=c("x","y")
  rescoord$status=paste0("point_",1:nrow(rescoord))
  row.names(rescoord)=rescoord$status
  

  
  tmpMeanFn=function(index,inputPCA,reductionCols,xrange,yrange,rescoord){
    value=NA
    tmp=which(abs(inputPCA[,reductionCols[1]]-rescoord[index,"x"])<=xrange&abs(inputPCA[,reductionCols[2]]-rescoord[index,"y"])<=yrange)
    if(!is.null(tmp)){
      if(length(tmp)>0){
        value=mean(inputPCA[tmp,continuous_batch_values],na.rm=T)
      }
    }
    return(value)
  }
  
  tmpRes=lapply(1:nrow(rescoord),tmpMeanFn,inputPCA=inputPCA,reductionCols=reductionCols,xrange=xrange,yrange=yrange,rescoord=rescoord)
  
  anno2=rescoord
  anno2$value=unlist(tmpRes)
  p=ggplot(anno2, aes(x, y))+geom_tile(aes(fill = log2(value)))+
    scale_fill_gradientn(colours = c("black","yellow","orange","red"))
  
  return(figure=p)
}

.my2dPlot_counts=function(inputPCA,batch_values,reductionCols,geneNameList=NULL,geneNameCol="gene_short_name",expData=NULL,ncolumns=6,useEffectiveCount=T,combine_figs=T,ncores=1,...){
  
  if(useEffectiveCount){
    p=.extra2dPlot_counts_effectiveCounts(inputPCA=inputPCA,batch_values=batch_values,reductionCols=reductionCols,geneNameList=geneNameList,geneNameCol=geneNameCol,expData=expData,ncolumns=ncolumns,combine_figs=combine_figs,ncores = ncores,...)
  } else {
    p=.extra2dPlot_counts_original(inputPCA=inputPCA,batch_values=batch_values,reductionCols=reductionCols,geneNameList=geneNameList,geneNameCol=geneNameCol,expData=expData,ncolumns=ncolumns,combine_figs=combine_figs)
  }
  
  return(figure=p)
}

.extra2dPlot_counts_original=function(inputPCA,batch_values,reductionCols,geneNameList=NULL,geneNameCol="gene_short_name",expData=NULL,ncolumns=6,combine_figs=T){
  require(Seurat)
  require(ggplot2)
  require(RANN)
  
  
  expData@assays$RNA@meta.features[,geneNameCol]=toupper(expData@assays$RNA@meta.features[,geneNameCol])
  
  
  geneNameList=toupper(geneNameList)
  
  if(sum(colnames(inputPCA)==batch_values)==0){
    stop("Specified batch column couldn't be found")
  }
  
  if(sum(colnames(inputPCA) %in% reductionCols)!=2){
    stop("Specified reduction columns couldn't be found ")
  }
  
  if(sum(is.na(inputPCA[,batch_values]))>0){
    inputPCA=inputPCA[!is.na(inputPCA[,batch_values]),]
    print("NA values are being discarded")
  }
  
  yval=inputPCA[,reductionCols]
  xval=yval[,1]
  yval=yval[,2]
  
  rescoord=expand.grid(seq(min(xval),max(xval),length.out = 100),seq(min(yval),max(yval),length.out = 100))
  xrange=abs(seq(min(xval),max(xval),length.out = 100)[2]-seq(min(xval),max(xval),length.out = 100)[1])/1.5
  yrange=abs((seq(min(yval),max(yval),length.out = 100)[2]-seq(min(yval),max(yval),length.out = 100)[1]))/1.5
  colnames(rescoord)=c("x","y")
  rescoord$status=paste0("point_",1:nrow(rescoord))
  row.names(rescoord)=rescoord$status
  
  tmpMeanFn=function(index,inputPCA,reductionCols,xrange,yrange,rescoord){
    value=NA
    tmp=which(abs(inputPCA[,reductionCols[1]]-rescoord[index,"x"])<=xrange&abs(inputPCA[,reductionCols[2]]-rescoord[index,"y"])<=yrange)
    if(!is.null(tmp)){
      if(length(tmp)>0){
        value=length(unique(inputPCA[tmp,batch_values],na.rm=T))
      }
    }
    return(value)
  }
  p=""
  
  if(is.null(geneNameList)){
    
    
    tmpRes=lapply(1:nrow(rescoord),tmpMeanFn,inputPCA=inputPCA,reductionCols=reductionCols,xrange=xrange,yrange=yrange,rescoord=rescoord)
    
    anno2=rescoord
    anno2$value=unlist(tmpRes)
    p=ggplot(anno2, aes(x, y))+geom_tile(aes(fill = value))+
      scale_fill_gradientn(colours = c("black",colorRampPalette(c("yellow","orange","red"))(5),"red"))
  } else {
    plots=list()
    
    if(sum(is.na(expData@assays$RNA@meta.features[,geneNameCol]))>0){
      expData@assays$RNA@meta.features[which(is.na(expData@assays$RNA@meta.features[,geneNameCol])),geneNameCol]=row.names(expData@assays$RNA@meta.features)[which(is.na(expData@assays$RNA@meta.features[,geneNameCol]))]
    }
    expData@assays$RNA@meta.features[,geneNameCol]=toupper(expData@assays$RNA@meta.features[,geneNameCol])
    geneNameList=lapply(geneNameList,toupper)
    
    for(i in 1:length(geneNameList)){
      if(length(which(expData@assays$RNA@meta.features[,geneNameCol] %in% unlist(strsplit(geneNameList[[i]],","))))>0){
        tmpExp=expData[which(expData@assays$RNA@meta.features[,geneNameCol] %in% unlist(strsplit(geneNameList[[i]],","))),]
        
        tmpExp=as.matrix(tmpExp@assays$RNA@data)
        tmpExp=apply(tmpExp,2,function(x) sum(x))
        tmpExp=tmpExp>max(quantile(tmpExp,0.9),1)
        tmpPCA=inputPCA[tmpExp,]
        tmpRes=lapply(1:nrow(rescoord),tmpMeanFn,inputPCA=tmpPCA,reductionCols=reductionCols,xrange=xrange,yrange=yrange,rescoord=rescoord)
        anno2=rescoord
        anno2$ds_count=unlist(tmpRes)
        tmpPlot=ggplot(anno2, aes(x, y))+geom_tile(aes(fill = ds_count))+
          scale_fill_gradientn(colours = c("black",colorRampPalette(c("yellow","orange","red"))(5),"red"))+ggtitle(geneNameList[[i]])+theme_bw()+theme(panel.grid = element_blank())
        plots=c(plots,list(tmpPlot))
        names(plots)[length(plots)]=geneNameList[[i]]
      }
      
    }
    
    if(combine_figs){
      p=patchwork::wrap_plots(plots,ncol=ncolumns,nrow=ceiling(length(plots)/ncolumns))
    } else {
      p=plots
    }
  }
  
  
  return(figure=p)
}

#ncolumns=6;exponetial_pwr=8;combine_figs=T
.extra2dPlot_counts_effectiveCounts=function(inputPCA,batch_values,reductionCols,geneNameList=NULL,geneNameCol="gene_short_name",expData=NULL,ncolumns=6,exponetial_pwr=8,combine_figs=T,ncores=1){
  require(Seurat)
  require(ggplot2)
  require(RANN)
  require(cowplot)
  
  
  geneNameList=toupper(geneNameList)
  
  if(sum(colnames(inputPCA)==batch_values)==0){
    stop("Specified batch column couldn't be found")
  }
  
  if(sum(colnames(inputPCA) %in% reductionCols)!=2){
    stop("Specified reduction columns couldn't be found ")
  }
  
  if(sum(is.na(inputPCA[,batch_values]))>0){
    inputPCA=inputPCA[!is.na(inputPCA[,batch_values]),]
    print("NA values are being discarded")
  }
  
  yval=inputPCA[,reductionCols]
  xval=yval[,1]
  yval=yval[,2]
  
  rescoord=expand.grid(seq(min(xval),max(xval),length.out = 100),seq(min(yval),max(yval),length.out = 100))
  xrange=abs(seq(min(xval),max(xval),length.out = 100)[2]-seq(min(xval),max(xval),length.out = 100)[1])/1.5
  yrange=abs((seq(min(yval),max(yval),length.out = 100)[2]-seq(min(yval),max(yval),length.out = 100)[1]))/1.5
  colnames(rescoord)=c("x","y")
  rescoord$status=paste0("point_",1:nrow(rescoord))
  row.names(rescoord)=rescoord$status
  
  tmpMeanFn=function(index,inputPCA,reductionCols,xrange,yrange,rescoord){
    effectiveSize=NULL
    tmp=which(abs(inputPCA[,reductionCols[1]]-rescoord[index,"x"])<=xrange&abs(inputPCA[,reductionCols[2]]-rescoord[index,"y"])<=yrange)
    if(length(tmp)>0){
      if(sum(is.na(tmp))==0){
        df=data.frame(coord_index=index,pd_index=tmp)
        effectiveSize=rbind(effectiveSize,df)
      }
    }
    return(effectiveSize)
  }
  tmpQuantileFn=function(inputExp,geneList=NULL){
    cnames=colnames(inputExp)
    
    if(!is.null(geneList)){
      inputExp@assays$RNA@meta.features[is.na(inputExp@assays$RNA@meta.features[,geneNameCol]),geneNameCol]=row.names(inputExp)[is.na(inputExp@assays$RNA@meta.features[,geneNameCol])]
      if(length(which(toupper(inputExp@assays$RNA@meta.features[,geneNameCol]) %in% toupper(geneList)))>0){
        inputExp=inputExp[which(toupper(inputExp@assays$RNA@meta.features[,geneNameCol]) %in% toupper(geneList)),]
      } else{
        inputExp=inputExp[1,]
        inputExp@assays$RNA@data[1,]=0
      }
    }
    
    updateFn=function(x){
      
      if(length(which(x>0.95*max(x)))>0&max(x)!=0){
        x[which(x>0.95*max(x))]=1
      } else if(length(which(x>0.9*max(x)&x<0.95*max(x)))>0&max(x)!=0){
        x[which(x>=0.9*max(x)&x<0.95*max(x))]=0.965
      }
      return(x)
    }
    
    inputExp=as.matrix(inputExp@assays$RNA@data)
    inputExp2=t(apply(inputExp,1,function(x) ecdf(x)(x)))
    inputExp2[which(inputExp==0)]=0
    inputExp2=t(apply(inputExp2,1,updateFn))
    inputExp2=inputExp2^exponetial_pwr
    inputExp2[which(inputExp==0)]=0
    inputExp2=apply(inputExp2,2,prod)
    names(inputExp2)=cnames
    return(inputExp2)
  }
  
  p=""
  
  indexList=lapply(1:nrow(rescoord),tmpMeanFn,inputPCA=inputPCA,reductionCols=reductionCols,xrange=xrange,yrange=yrange,rescoord=rescoord)
  indexList=do.call("rbind",indexList)
  
  batch_factor=as.data.frame(table(inputPCA[,batch_values]))
  batch_factor$factor=10000/batch_factor$Freq
  
  
  if(is.null(geneNameList)){
    
    tmpRes=indexList
    tmpRes$batch=inputPCA[tmpRes$pd_index,batch_values]
    tmpRes$expWeight=1
    tmpRes$index=tmpRes$coord_index
    batchDf_sl=aggregate(expWeight~batch+index,tmpRes,sum)
    batchDf_sl=merge(batchDf_sl,batch_factor,by.x="batch",by.y="Var1")
    batchDf_sl$batch_probs=batchDf_sl$expWeight*batchDf_sl$factor
    batch_probs=aggregate(batch_probs~index,batchDf_sl,sum)
    batchDf_sl=merge(batchDf_sl,batch_probs,by="index")
    batchDf_sl$batch_probs=batchDf_sl$batch_probs.x/batchDf_sl$batch_probs.y
    batchDf_sl$effectiveSize=batchDf_sl$batch_probs^2
    effectiveSize=aggregate(effectiveSize~index,batchDf_sl,function(x) 1/sum(x))
    anno2=rescoord
    effectiveSize=effectiveSize$effectiveSize[match(1:nrow(anno2),effectiveSize$index)]
    
    anno2$ds_count=effectiveSize
    p=ggplot(anno2, aes(x, y))+geom_tile(aes(fill = ds_count))+
      scale_fill_gradientn(colours = c("black",colorRampPalette(c("yellow","orange","red"))(5),"red"))+ggtitle(geneNameList[[i]])+theme_cowplot()+ theme(plot.title = element_text(hjust = 0.5))
    
  } else {
    plots=list()
    
    
    
    geneNameList=lapply(geneNameList,toupper)
    total_gene_list=c()
    if(class(expData)!=class(list())){
      expData@assays$RNA@meta.features[,geneNameCol]=toupper(expData@assays$RNA@meta.features[,geneNameCol])
      tmpExp_main=expData[which(expData@assays$RNA@meta.features[,geneNameCol] %in% unlist(strsplit(unlist(geneNameList),","))),]
      tmpExp_main=SplitObject(tmpExp_main, split.by = batch_values)
      total_gene_list=unique(expData@assays$RNA@meta.features[,geneNameCol])
    } else {
      tmpExp_main=list()
      for(ik in 1:length(expData)){
        tmp=expData[[ik]]
        if(sum(is.na(tmp@assays$RNA@meta.features[,geneNameCol]))>0){
          tmp@assays$RNA@meta.features[which(is.na(tmp@assays$RNA@meta.features[,geneNameCol])),geneNameCol]=row.names(tmp@assays$RNA@meta.features)[which(is.na(tmp@assays$RNA@meta.features[,geneNameCol]))]
        }
        tmp@assays$RNA@meta.features[,geneNameCol]=toupper(tmp@assays$RNA@meta.features[,geneNameCol])
        tmp=tmp[which(tmp@assays$RNA@meta.features[,geneNameCol] %in% unlist(strsplit(unlist(geneNameList),","))),]
        tmpExp_main=c(tmpExp_main,list(tmp))
        if(length(names(expData))>=ik){
          names(tmpExp_main)[length(tmpExp_main)]=names(expData)[ik]
        }
        total_gene_list=unique(c(total_gene_list,tmp@assays$RNA@meta.features[,geneNameCol]))
      }
    }
    
    
    for(i in 1:length(geneNameList)){
      if(length(which(total_gene_list %in% unlist(strsplit(geneNameList[[i]],","))))>0){
        
        tmpExpList=parallel::mclapply(tmpExp_main,tmpQuantileFn,geneList=unlist(strsplit(geneNameList[[i]],",")),mc.cores = ncores)
        tmpExpList=unlist(unname(tmpExpList))
        tmpExpList=tmpExpList[match(row.names(inputPCA),names(tmpExpList))]
        if(!all(row.names(inputPCA)==names(tmpExpList))){
          stop("Error in the matching names!")
        }
        
        tmpRes=indexList
        tmpRes$batch=inputPCA[tmpRes$pd_index,batch_values]
        tmpRes$expWeight=tmpExpList[tmpRes$pd_index]
        if(sum(tmpRes$expWeight>0)>0){
          tmpRes$index=tmpRes$coord_index
          batchDf_sl=aggregate(expWeight~batch+index,tmpRes,sum)
          batchDf_sl=batchDf_sl[which(batchDf_sl$expWeight>0),]
          batchDf_sl=merge(batchDf_sl,batch_factor,by.x="batch",by.y="Var1")
          batchDf_sl$batch_probs=batchDf_sl$expWeight*batchDf_sl$factor
          batch_probs=aggregate(batch_probs~index,batchDf_sl,sum)
          batchDf_sl=merge(batchDf_sl,batch_probs,by="index")
          batchDf_sl$batch_probs=batchDf_sl$batch_probs.x/batchDf_sl$batch_probs.y
          batchDf_sl$effectiveSize=batchDf_sl$batch_probs^2
          effectiveSize=aggregate(effectiveSize~index,batchDf_sl,function(x) 1/sum(x))
          anno2=rescoord
          effectiveSize=effectiveSize$effectiveSize[match(1:nrow(anno2),effectiveSize$index)]
          
          anno2$ds_count=effectiveSize
          tmpPlot=ggplot(anno2, aes(x, y))+geom_tile(aes(fill = ds_count))+
            scale_fill_gradientn(colours = c(colorRampPalette(c("black","yellow","orange","red"))(8),"red"))+ggtitle(geneNameList[[i]])+theme_cowplot()+ theme(plot.title = element_text(hjust = 0.5))
          plots=c(plots,list(tmpPlot))
          names(plots)[length(plots)]=geneNameList[[i]]
        }
       
      }
      
    }
    
    if(combine_figs){
      p=patchwork::wrap_plots(plots,ncol=ncolumns,nrow=ceiling(length(plots)/ncolumns))
      
    } else {
      p=plots
    }
    
  }
  
  
  return(figure=p)
}

.extra_gmean=function (x, eps = 1,weights=NULL) {
  
  if(is.null(weights)){
    if (inherits(x = x, what = "matrix")) {
      return(exp(rowMeans(log(x + eps))) - eps)
    }
    if (inherits(x = x, what = "dgCMatrix")) {
      
      ret <- sctransform:::row_gmean_dgcmatrix(matrix = x, eps = eps)
      names(ret) <- rownames(x)
      return(ret)
    }
  } else {
    x=as.matrix(x)
    x=log(x+eps)
    x=apply(x,1,function(x) weighted.mean(x,w=weights))
    return(exp(x)-eps)
  }
  
  stop("matrix x needs to be of class matrix or dgCMatrix")
}

.my2dEntropyPlot2=function(inputPCA,batchCol,reductionCols,perm=F,k=0.3){
  require(Seurat)
  require(ggplot2)
  require(RANN)
  
  if(sum(colnames(inputPCA)==batchCol)==0){
    stop("Specified batch column couldn't be found")
  }
  
  if(sum(colnames(inputPCA) %in% reductionCols)!=2){
    stop("Specified reduction columns couldn't be found ")
  }
  
  anno=data.frame(sample=row.names(inputPCA),value=as.factor(inputPCA[,batchCol]))
  
  yval=inputPCA[,reductionCols]
  xval=yval[,1]
  yval=yval[,2]
  
  rescoord=expand.grid(seq(min(xval),max(xval),length.out = 100),seq(min(yval),max(yval),length.out = 100))
  colnames(rescoord)=c("x","y")
  rescoord$status=paste0("point_",1:nrow(rescoord))
  row.names(rescoord)=rescoord$status
  
  
  nn <- RANN::nn2(query = rescoord[,-3],data = data.frame(x=xval,y=yval,stringsAsFactors = F), k = 40, eps = 0)
  
  
  nn_thr=quantile(as.numeric(nn$nn.dists),k)
  
  pval_res=.myPvalFn(nn_data=nn,nn.thr=nn_thr,anno_data=anno,coord_data=rescoord)
  
  if(perm){
    require(parallel)
    set.seed(12345)
    res_perm=mclapply(1:1000,.myPvalPermFn,nn_data=nn,nn.thr=nn_thr,anno_data=anno,coord_data=rescoord,mc.cores = 5)
    res_perm=do.call("rbind",res_perm)
    x=split(res_perm,f=paste0(res_perm$x,"_",res_perm$y))
    
    xname=unlist(lapply(x,function(x) paste0(x[1,1],"_",x[1,2])))
    o=match(paste0(pval_res$x,"_",pval_res$y),xname)
    x=x[o]
    
    emp_pval=mclapply(1:nrow(pval_res),function(y) unlist(lapply(4:ncol(pval_res),function(z) sum(x[[y]][,z]<=pval_res[y,z]))),mc.cores = 5)
    emp_pval=do.call("rbind",emp_pval)
    emp_pval=emp_pval/1000
    
    the_pval=pval_res[,4:ncol(pval_res)]
    emp_pval[the_pval>emp_pval]=the_pval[the_pval>emp_pval]
    pval_res[4:ncol(pval_res)]=emp_pval
    
  }
  
  pval_res=reshape::melt(pval_res,id.vars=colnames(pval_res)[c(1,2,3)])
  colnames(pval_res)[colnames(pval_res)=="value"]="pval"
  pval_res$score=-log10(pval_res$pval)
  
  pval_res$score[pval_res$score<(-log10(0.01))]=0
  pval_res$score[pval_res$score<0]=0
  p=ggplot(pval_res, aes(x, y,fill = score)) +geom_raster()+
    scale_fill_gradientn(colours = c("black","black","yellow","red"),values=c(0,0.5,-log10(0.01),max(pval_res$score))/max(pval_res$score))+
    theme_bw()+theme(panel.grid = element_blank())+facet_wrap(~variable)+xlab(toupper(reductionCols[1]))+ylab(toupper(reductionCols[2]))
  
  return(figure=p)
}

.myExpSetCreatorFn=function(inputExpData,organism,minExpCells=0,inputPdata=NULL,inputFdata=NULL,addExtraAnno=T,server=T,redownload_files=T,ncores=5){
  require(scran)
  require(DropletUtils)
  require(Matrix)
  if(is.null(colnames(inputExpData))){
    colnames(inputExpData)=1:ncol(inputExpData)
  }
  
  if(tolower(class(inputExpData))=="seurat"){
    inputExpData=.sconline.convertSeuratToExpressionSet(object=inputExpData)
  }
  
  if(is(inputExpData,"SingleCellExperiment")){
    tmp=as.data.frame(rowData(inputExpData))
    if(is.null(inputFdata)){
      inputFdata=tmp
    } else {
      inputFdata=cbind(inputFdata,tmp)
    }
    
    tmpPdata=as.data.frame(colData(inputExpData))
    if(is.null(inputPdata)){
      inputPdata=tmpPdata
    } else {
      inputPdata=cbind(inputPdata,tmpPdata)
    }
  }
  
  if(is(inputExpData,"Seurat")){
    tmpPdata=as.data.frame(inputPdata[[]])
    if(is.null(inputPdata)){
      inputPdata=tmpPdata
    } else {
      inputPdata=cbind(inputPdata,tmpPdata)
    }
  }
  
  
  if(!is.null(inputPdata)){
    pd <- inputPdata
  } else {
    pd=data.frame(sample=colnames(inputExpData),stringsAsFactors = F)
    row.names(pd)=colnames(inputExpData)
  }
  
  if(tolower(organism)=="human"){
    tmpName=.extraHumanGeneAnnoAdderFn(row.names(inputExpData),server = server)
  } else if (tolower(organism)=="mouse"){
    tmpName=.extraMouseGeneAnnoAdderFn(row.names(inputExpData),server = server)
  } else if(tolower(organism)=="macaque"){
    tmpName=.extraMacaqueGeneAnnoAdderFn(row.names(inputExpData),server = server)
  } else {
    warning("unknown organism!")
  }
  
  fd=NULL
  if(tolower(organism) %in% c("human","mouse","macaque")){
    if(all(row.names(inputExpData)==row.names(tmpName))){
    fd=tmpName
    
    if(!is.null(inputFdata)){
      fd=cbind(inputFdata,fd)
    }
    if(sum(colnames(fd)=="ID")==0){
      fd$ID=row.names(inputExpData)
    }
    row.names(fd)=row.names(inputExpData)
    
    if(tolower(organism) %in% c("mouse","human","macaque")){
      mycc.genes=.extraCellCycleGenes(organism = organism)
      fd$QC_ssGenes="No"
      fd$QC_g2mGenes="No"
      fd$QC_ssGenes[fd$ensembl_gene_id %in% row.names(mycc.genes$cc.s)]="s_phase"
      fd$QC_g2mGenes[fd$ensembl_gene_id %in% row.names(mycc.genes$cc.g2m)]="g2m_phase"
      
      mymito.genes=.extraMitoGenes(organism=organism)
      fd$QC_mtGenes="No"
      fd$QC_mtGenes[fd$ensembl_gene_id %in% mymito.genes$ensembl_gene_id]="Yes"
      
      myIEG.genes=.extraIEGGenes(organism=organism,server = server)
      fd$QC_IEG_Genes="No"
      fd$QC_IEG_Genes[fd$ensembl_gene_id %in% myIEG.genes$ensembl_gene_id]="Yes"
    }
   
    
  } else {
      print("Error!")
    }
    if(sum(colnames(fd) %in% c("seqnames", "ranges", "strand", "start", "end", "width","element"))>0){
      colnames(fd)[colnames(fd) %in% c("seqnames", "ranges", "strand", "start", "end", "width","element")]=paste0("anno_",colnames(fd)[colnames(fd) %in% c("seqnames", "ranges", "strand", "start", "end", "width","element")])
      
    }
  } else {
    fd=inputFdata
  }
  
  
  
  if(class(inputExpData)=="SingleCellExperiment"){
    res=counts(inputExpData)
  } else {
    if(class(inputExpData)!="dgCMatrix"){
      if(class(inputExpData)=="data.frame"){
        res=Seurat::as.sparse(inputExpData)
      } else {
        res=.matrixExtraction(inputExpData)
      }
      row.names(res)=row.names(inputExpData)
      colnames(res)=colnames(inputExpData)
    } else {
      res=inputExpData
    }
    
  }
  
  colnames(res)=colnames(inputExpData)
  
  res <- as(res, "dgCMatrix")
  
  res=SingleCellExperiment(assays = list(counts = res),colData = pd,rowData=fd)
  rm(inputExpData)
  
  if(minExpCells>0){
    tmpFilter=c()
    
    if(nrow(res)>10000){
      for(i in seq(1,nrow(res),10000)){
        tmpFilter=c(tmpFilter,apply(assays(res)[["counts"]][i:min(i+10000-1,nrow(res)),],1,function(x) sum(x>0)))
      }
    } else {
      tmpFilter=apply(assays(res)[["counts"]],1,function(x) sum(x>0))
    }
    
    
    res=res[which(tmpFilter>=minExpCells),]
  }
  
  
  
  if(addExtraAnno){
    #inputExpSet = res;organism=organism;server=server;redownload_files=redownload_files
    res=.extraQCAnnoAdderFn(inputExpSet = res,organism=organism,server=server,redownload_files=redownload_files,ncores=ncores)
  }
  
  for(i in 1:ncol(colData(res))){
    if(class(colData(res)[,i])=="factor"){
      colData(res)[,i]=droplevels(colData(res)[,i])
    }
  }
  
  return(res)
}

.myRTNgsea=function (listOfRegulons, phenotype,method='gsea1', pAdjustMethod = "BH", nPermutations = 1000, 
                     exponent = 1, orderAbsValue = TRUE, verbose = TRUE) {
  require(RTN)
  
  
  if(method=="gsea1"){
    for (i in names(listOfRegulons)) {
      regs <- listOfRegulons[[i]]
      listOfRegulons[[i]] <- regs[regs %in% names(phenotype)]
    }
    GSEA.results <- RTN:::run.gsea1(listOfRegulons = listOfRegulons, phenotype = phenotype, 
                                    pAdjustMethod = pAdjustMethod, nPermutations = nPermutations, 
                                    exponent = exponent, orderAbsValue = orderAbsValue, 
                                    verbose = verbose)
  } else if(method=='gsea2') {
    for (i in names(listOfRegulons)) {
      regs <- listOfRegulons[[i]]
      listOfRegulons[[i]] <- regs[names(regs) %in% names(phenotype)]
    }
    if(is.null(names(listOfRegulons[[1]]))){
      stop("for two sided gsea, a named numeric vector should be provided!")
    }
    GSEA.results <- RTN:::run.gsea2(listOfRegulonsAndMode = listOfRegulons, phenotype = phenotype, 
                                    pAdjustMethod = pAdjustMethod, nPermutations = nPermutations, 
                                    exponent = exponent, 
                                    verbose = verbose)
    #RTN:::run.gsea2(listOfRegulonsAndMode = listOfRegulonsAndMode, 
    #                phenotype = object@phenotype, pAdjustMethod = object@para$gsea2$pAdjustMethod, 
    #                nPermutations = object@para$gsea2$nPermutations, exponent = object@para$gsea2$exponent, 
    #                verbose = verbose)
    
  }
  
  idx <- GSEA.results$Observed.Score <= 0
  if(sum(idx)>0){
    GSEA.results$Pvalue[idx] <- 1
    GSEA.results$Adjusted.Pvalue[idx] <- 1
  }
  
  
  return(GSEA.results)
}

.myQCfn_archived=function(inputData){
  .th=theme_bw()+theme(panel.grid = element_blank(),axis.text=element_text(color = "black"))
  
  require(ggplot2)
  require(patchwork)
  require(scran)
  
  top50_expFn=function (object, feature_names_to_plot = "gene_short_name") {
    set.seed(12345)
    exprs_mat <- counts(object)
    ave_exprs <- rowSums(exprs_mat)
    oo <- order(ave_exprs, decreasing = TRUE)
    chosen <- which(rowData(object)$QC_top50_expressed=="Yes")
    chosen=intersect(oo,chosen)
    sub_mat <- exprs_mat[chosen, , drop = FALSE]
    sub_ave <- ave_exprs[chosen]
    
    feature_names <- .fData(object)[,feature_names_to_plot]
    
    sub_names <- feature_names[chosen]
    {
      total_exprs <- sum(ave_exprs)
      top_pctage <- 100 * sum(sub_ave)/total_exprs
      sub_mat <- 100 * sweep(sub_mat, 2, colSums(exprs_mat), 
                             "/", check.margin = FALSE)
    }
    ordered_names <- factor(sub_names, rev(sub_names))
    
    sub_mat2=sub_mat[,sample(ncol(sub_mat),min(300,ncol(sub_mat)))]
    df_exprs_by_cell_long <- data.frame(Cell = rep(seq_len(ncol(sub_mat2)), 
                                                   each = nrow(sub_mat2)), Tag = rep(ordered_names, ncol(sub_mat2)), 
                                        value = as.numeric(sub_mat2))
    
    
    plot_most_expressed <- ggplot(df_exprs_by_cell_long, aes(y = Tag, x = value)) + 
      geom_point(alpha = 0.6, shape = 124)
    plot_most_expressed <- plot_most_expressed + xlab("% of total exp. genes") + 
      ylab("") + theme_bw(8) + theme(legend.position = c(1,0), legend.justification = c(1, 0), axis.text = element_text(colour = "black"), 
                                     axis.title = element_text(colour = "black"), title = element_text(colour = "black"))
    
    df_to_plot <- data.frame(Feature = ordered_names)
    df_to_plot$pct_total <- 100 * sub_ave/total_exprs
    
    plot_most_expressed=plot_most_expressed + geom_point(data = df_to_plot, aes(x = pct_total, y = Feature),
                                                         colour = "red",fill="red", shape = 21)
    
    return(plot_most_expressed)
  }
  
  if(sum(colnames(.pData(inputData = inputData))=="QC_MT.pct")>0){
    pMT=ggplot(.pData(inputData),aes(QC_MT.pct))+geom_histogram()+ggtitle("QC - %MT",paste0("MT >5%: ",sum(inputData$QC_MT.pct>5)," cells (",round(sum(inputData$QC_MT.pct>5)/ncol(inputData),2)*100,"%); excluded!"))+ylab("Cell count")+xlab("% MT genes")+geom_vline(xintercept = 5,color="red")+.th
    inputData=inputData[,which(.pData(inputData)$QC_MT.pct<=5)]
  } else {
    print("No Mito gene is annotated in this file!")
    pMT=""
  }
  
  qt995=quantile(inputData$QC_Gene_unique_count,0.995)
  qt005=quantile(inputData$QC_Gene_unique_count,0.005)
  
  p1=ggplot(.fData(inputData),aes(QC_meanExp,QC_varExp))+geom_point(alpha=0.1)+geom_density_2d()+geom_abline(slope = 1,color="red")+scale_x_continuous(trans = 'log2')+scale_y_continuous(trans='log2')+xlab("Mean Expression")+ylab("variance")+.th
  
  x = seq(from = min(log10(.fData(inputData)$QC_meanExp+1e-5)), to = max(log10(.fData(inputData)$QC_meanExp)), length.out = 1000)
  poisson_model <- data.frame(log_mean = x, detection_rate = 1 - dpois(0, lambda = 10^x))
  p2=ggplot(.fData(inputData), aes(log10(QC_meanExp), QC_detectionRate)) + 
    geom_point(alpha=0.3, shape=16) + 
    geom_line(data=poisson_model,aes(log_mean, detection_rate), color='red') +
    theme_bw()+xlab("Log Mean Exp")+ylab("Detection rate")
  
  pMeanVar=p1+p2
  
  br.out <- DropletUtils::barcodeRanks(counts(inputData))
  
  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
         legend=c("knee", "inflection"))
  pKneePoint=recordPlot()
  
  pd=colData(inputData)
  
  tst=pd$QC_Gene_total_count
  tst=tst[order(tst,decreasing = T)]
  tst2=cumsum(tst)
  tst=data.frame(Rank=1:length(tst2),cumulative_sum=tst2/tst2[length(tst2)])
  pCumSum=ggplot(tst,aes(Rank,cumulative_sum))+geom_point()+theme_bw()
  
  thr10th99ile=quantile(inputData$QC_Gene_total_count,0.99)/10
  pTotalCount=ggplot(.pData(inputData),aes("sample",QC_Gene_total_count))+ggtitle("",paste("0.5%ile:",quantile(inputData$QC_Gene_unique_count,0.005),";","99.5%ile:",quantile(inputData$QC_Gene_unique_count,0.995),"; excluded!"))+geom_violin(color="red")+geom_jitter(alpha=0.1,height = 0)+scale_y_continuous(breaks = seq(0,100000,500))+ylab("Unique UMI")+xlab("")+.th
  inputData=inputData[,which(.pData(inputData)$QC_Gene_unique_count<qt995&.pData(inputData)$QC_Gene_unique_count>qt005)]
  pQC1=pMT+pTotalCount
  
  colData(inputData)$QC_pass_kneePoint="No"
  colData(inputData)$QC_pass_kneePoint[which(colData(inputData)$QC_Gene_total_count>=metadata(br.out)$knee)]="Yes"
  colData(inputData)$QC_pass_inflectionPoint="No"
  colData(inputData)$QC_pass_inflectionPoint[which(colData(inputData)$QC_Gene_total_count>=metadata(br.out)$inflection)]="Yes"
  
  top50exp=top50_expFn(inputData)
  
  return(list(data=inputData,QC1plot=pQC1,CumulativeSumPlot=pCumSum,KneePlot=pKneePoint,MeanVarPlot=pMeanVar,top50exp=top50exp,thr10th99ile=thr10th99ile))
  
}

.myDataNorm=function(inputData,Method=c("LogTransform","sqrtTransform","sizeFactor","scTransform"),zTransform=F,verbose=T){
  require(Matrix)
  require(scran)
  if(verbose){
    print("Before normalization, low quality cells and low expressed genes should be removed!")
  }
  
  
  if(max(counts(inputData))<20){
    print("Error! data is supposed to be in the count format!")
  }
  
  
  if(sum(Method %in% c("LogTransform","sqrtTransform"))>0){
    
    dataScaled=Seurat::NormalizeData(counts(inputData),normalization.method="RC",scale.factor=10000,verbose=F)
    
    assay(inputData,"colScaled")=dataScaled
    if(zTransform){
      assay(inputData,"z_colScaled")=.extraScaleFn(dataScaled)
    }
    
    tmp1=log1p(dataScaled)
    assay(inputData,"lognorm")=tmp1
    if(zTransform){
      assay(inputData,"z_lognorm")=.extraScaleFn(tmp1)
    }
    
    if(max(assay(inputData,"lognorm"))>50){
      print("Error in logTransformation!")
    }
    
    tmp2=sqrt(dataScaled)
    assay(inputData,"sqrtnorm")=tmp2
    if(zTransform){
      assay(inputData,"z_sqrtnorm")=.extraScaleFn(tmp2)
    }
    
    
    if(max(assay(inputData,"sqrtnorm"))>60){
      print("Error in sqrt transformation!")
    }
    
  }
  
  if(sum(Method=="sizeFactor")>0){
    #quickCluster() automatically normalize the data by logTransformation (log2)
    #getMethod(quickCluster,signature = "ANY")
    
    clusters <- quickCluster(inputData)
    #Counts+1 are divided by the size factors since sf is directly correlated with the library size
    sf <- computeSumFactors(assay(inputData,"counts"), clusters=clusters)
    
    normData=scater::logNormCounts(inputData,size_factors=sf)
    assay(inputData,"logSizeFactor")=as(assay(normData,"logcounts"), "dgCMatrix")
    if(zTransform){
      assay(inputData,"z_sizeFactor")=.extraScaleFn(assay(normData,"logcounts"))
    }
    
  }
  
  if(sum(Method=="scTransform")>0){
    normData=sctransform::vst(counts(inputData),n_genes = min(nrow(inputData),7000),return_cell_attr = T,return_gene_attr = T,show_progress=F)
    inputData=inputData[row.names(inputData) %in% row.names(normData$y),]
    assay(inputData,"scTransform")=normData$y
    rowData(inputData)$variance=normData$gene_attr$variance
    rowData(inputData)$scTransform_resid_var=normData$gene_attr$residual_variance
    #plot(log10(normData$gene_attr$gmean),normData$gene_attr$residual_variance)
  }
  
  return(inputData)
}

.myHighlyVarGenesFn=function(inputData,minExpThr=0.08,maxExpThr=exp(8)){
  
  hvgVst=.extraHVGvstFn(counts(inputData))
  hvgVst=hvgVst[match(row.names(inputData),row.names(hvgVst)),]
  if(all(row.names(inputData)==row.names(hvgVst))){
    rowData(inputData)$hvg_vst_score_notScaled=round(hvgVst$vst.variance.standardized,3)
    rowData(inputData)$hvg_vst_rank_notScaled=rank(hvgVst$vst.variance.standardized*(-1))
  } else {
    print("Error!")
  }
  
  if(any(names(assays(inputData)) %in% "colScaled")){
    hvgVst=.extraHVGvstFn(assay(inputData,"colScaled"))
    hvgVst=hvgVst[match(row.names(inputData),row.names(hvgVst)),]
    if(all(row.names(inputData)==row.names(hvgVst))){
      rowData(inputData)$hvg_vst_score_scaled=round(hvgVst$vst.variance.standardized,3)
      rowData(inputData)$hvg_vst_rank_scaled=rank(hvgVst$vst.variance.standardized*(-1))
    } else {
      print("Error!")
    }
  }
  
  for(i in c("lognorm","logSizeFactor","sqrtnorm")){
    if(any(names(assays(inputData)) %in% i)){
      binningMethod="equal_width"
      tmpThr=c(minExpThr,maxExpThr)
      if(i=="lognorm"){
        tmpThr=log1p(tmpThr)
      } else {
        if(i=="logSizeFactor"){
          tmpThr=log2(tmpThr+1)
        } else {
          if(i=="sqrtnorm"){
            tmpThr=sqrt(tmpThr)
            binningMethod="equal_frequency"
          }
        }
      }
      
      tmp=assay(inputData,i)
      tmp=Seurat:::FindVariableFeatures.default(tmp,selection.method="mvp",binning.method=binningMethod,verbose=F)
      
      tmp$mvp.dispersion.scaled[(tmp$mvp.mean<tmpThr[1]|tmp$mvp.mean>tmpThr[2])]=1
      tmp=tmp[match(row.names(inputData),row.names(tmp)),]
      rowData(inputData)$new=round(tmp$mvp.dispersion.scaled,3)
      rowData(inputData)$new_rank=rank(tmp$mvp.dispersion.scaled*-1)
      if(all(row.names(inputData)==row.names(tmp))){
        colnames(rowData(inputData))[colnames(rowData(inputData))=="new"]=paste0("hvg_mvp.dispersion.std_",i)
        colnames(rowData(inputData))[colnames(rowData(inputData))=="new_rank"]=paste0("hvg_rank_mvp.dispersion.std_",i)
      } else {
        print(paste("Error in",i))
      }
      rm(tmp)
    }
  }
  
  if(sum(colnames(rowData(inputData))=="scTransform_resid_var")>0){
    rowData(inputData)$scTransform_rank_resid_var=rank(rowData(inputData)$scTransform_resid_var*-1)
  }
  
  return(inputData)
}

.myOneHotFn=function (inputVector) {
  if(length(unique(inputVector))>1){
    if(T){
      #using caret package
      require(caret)
      if(sum(is.na(inputVector))>0){
        inputVector[is.na(inputVector)]="NA"
      }
      formula="~."
      data=data.frame(data=inputVector,stringsAsFactors = F)
      sep = "."
      levelsOnly = FALSE
      fullRank = FALSE
      
      formula <- as.formula(formula)
      if (!is.data.frame(data)) 
        data <- as.data.frame(data, stringsAsFactors = FALSE)
      vars <- all.vars(formula)
      if (any(vars == ".")) {
        vars <- vars[vars != "."]
        vars <- unique(c(vars, colnames(data)))
      }
      isFac <- unlist(lapply(data[, vars, drop = FALSE], is.factor))
      if (sum(isFac) > 0) {
        facVars <- vars[isFac]
        lvls <- lapply(data[, facVars, drop = FALSE], levels)
        if (levelsOnly) {
          tabs <- table(unlist(lvls))
          if (any(tabs > 1)) {
            stop(paste("You requested `levelsOnly = TRUE` but", 
                       "the following levels are not unique", "across predictors:", 
                       paste(names(tabs)[tabs > 1], collapse = ", ")))
          }
        }
      } else {
        facVars <- NULL
        lvls <- NULL
      }
      trms <- attr(model.frame(formula, data), "terms")
      out <- list(call = match.call(), form = formula, vars = vars, 
                  facVars = facVars, lvls = lvls, sep = sep, terms = trms, 
                  levelsOnly = levelsOnly, fullRank = fullRank)
      class(out) <- "dummyVars"
      trsf <- data.frame(predict(out, newdata = data))
      
      tmp=gsub("^data","",colnames(trsf))
      tmp=gsub("^\\.","",tmp)
      colnames(trsf)=tmp
    } else {
      if(length(names(inputVector))==0){
        trsf=reshape2::dcast(data=data.frame(var=1:length(inputVector),outcome=as.character(inputVector)),var ~ outcome, length)
        row.names(trsf)=trsf[,1]
        trsf=trsf[,-1]
      } else {
        inputVector=inputVector[order(inputVector)]
        inputVector=names(inputVector)
        trsf=reshape2::dcast(data=data.frame(var=1:length(inputVector),outcome=as.character(inputVector)),var ~ outcome, length)
        trsf=as.matrix(trsf)
        row.names(trsf)=trsf[,1]
        trsf=trsf[,-1]
      }
      
    }
  } else {
    trsf=matrix(1,ncol=1,nrow=length(inputVector))
    row.names(trsf)=names(inputVector)
    colnames(trsf)=unique(inputVector)
  }
  
  
  
  return(trsf)
}

.myHVGadder=function(inputData,organism,server=server){
  
  if(class(inputData)=="Seurat"){
    x=inputData@assays$RNA@counts
    x1=.myExpSetCreatorFn(inputExpData = x,organism = organism,minExpCells = 3,addExtraAnno=F,server=server)
  } else{
    if(class(inputData)=="SingleCellExperiment"){
      x1=inputData
      inputData=.extraExport2SeuratFn(inputData = inputData)
    } else {
      stop("Wrong input data format!")
      
    }
  }
  
  
  x2=.myDataNorm(inputData = x1)
  x3=.myHighlyVarGenesFn(x2)
  x4=as.data.frame(rowData(x3))
  inputData=inputData[row.names(inputData) %in% row.names(x4),]
  x4=x4[match(row.names(inputData) , row.names(x4)),]
  if(any(row.names(inputData) !=row.names(x4))){
    stop("Error in mapping the names")
  }
  
  x5=inputData@assays$RNA@meta.features
  x4=x4[,setdiff(colnames(x4),colnames(x5))]
  
  inputData@assays$RNA@meta.features=x4
  return(inputData)
}

#inputExpData=tmpExp;inputPCAembeddings=tmp_embeddings; n.adaptiveKernel=5; nPropIter=3;nPCs=30;verbose=T
.myPseudoCellfn=function(inputExpData,inputPCAembeddings=NULL, n.adaptiveKernel=5, nPropIter=3,nPCs=30,verbose=T){
  #nPropIter: number of propagation iterations
  #inputPCAembeddings: PCA res
  #inputExpData: the expression data is expected to be raw counts (singlecell expression set)
  
  if(!("spam" %in% rownames(installed.packages()))){
    install.packages("spam")
  }
  
  if(ncol(inputExpData)<100){
    n.adaptiveKernel=4
  }
  if(ncol(inputExpData)<50){
    n.adaptiveKernel=3
  }
  if(ncol(inputExpData)<25){
    n.adaptiveKernel=2
  }
  
  if(ncol(inputExpData)<18){
    n.adaptiveKernel=1
  }
  
  if(is.null(inputPCAembeddings)){
    tmpData=.extraExport2SeuratFn(inputExpData)
    tmpData = NormalizeData(tmpData)
    tmpData = FindVariableFeatures(tmpData, selection.method = "vst", nfeatures = 2000)
    tmpData <- ScaleData(tmpData)
    tmpData <- RunPCA(tmpData)
    inputPCAembeddings=tmpData@reductions$pca@cell.embeddings[,1:nPCs]
  } else {
    inputPCAembeddings=inputPCAembeddings[,1:nPCs]
  }
  
  n.cells <- nrow(inputPCAembeddings)
  
  minNeighborList=c(6,7,8,9,10,11,12,13,20)
  simThr_coef=c(1,1,1,1.1,1.2,1.5,1.5,1.5,1.5)
  
  
  nn.ranked.org <- RANN::nn2(data = inputPCAembeddings, k = min(ncol(inputExpData),n.adaptiveKernel*12), eps = 0)
  
  nn.ranked=nn.ranked.org
  
  dists=nn.ranked$nn.dists
  affinities=t(apply(dists,1,function(x) exp((-1)*(x/x[n.adaptiveKernel])^2)))
  affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
  affinities2=affinities-affinitiesThr
  affinities[affinities2<=0]=0
  
  affCounts=apply(affinities,2,function(x) sum(x==0))
  affinities=affinities[,-which(affCounts==nrow(affinities))]
  nn.ranked$nn.idx=nn.ranked$nn.idx[,-which(affCounts==nrow(affinities))]
  j <- as.numeric(t(nn.ranked$nn.idx))
  i <- ((1:length(j)) - 1)%/%ncol(affinities) + 1
  k=as.numeric(t(affinities))
  
  
  graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(inputPCAembeddings), nrow(inputPCAembeddings)))
  graph=(graph+t(graph))/2
  graph <- Matrix::Diagonal(x = 1 / (rowSums(graph))) %*% graph
  
  if(class(graph)!="dgCMatrix"){
    graph=as(graph,"dgCMatrix")
  }
  
  rownames(graph) <- rownames(inputPCAembeddings)
  colnames(graph) <- rownames(inputPCAembeddings)
  
  for(i in 1:nPropIter){
    if(verbose){
      #print(i)
    }
    
    if(i==1){
      prop_mat=graph
    } else {
      prop_mat=prop_mat%*%graph
    }
  }
  
  maxScores=spam::diag(prop_mat)
  
  affinityNorm <- Matrix::Diagonal(x = 1 / maxScores) %*% prop_mat
  
  df.org=list()
  for(i in 2:ncol(nn.ranked.org$nn.idx)){
    tmp=data.frame(id1=nn.ranked.org$nn.idx[,1],id2=nn.ranked.org$nn.idx[,i],value=1)
    df.org=c(df.org,list(tmp))
  }
  df.org=do.call("rbind",df.org)
  
  df.org <- sparseMatrix(i = df.org$id1, j = df.org$id2, x = df.org$value,dims = c(nrow(affinityNorm), ncol(affinityNorm)))
  
  df.org=df.org*affinityNorm
  
  df.org=summary(df.org)
  df.org=df.org[order(df.org$x,decreasing = T),]
  
  df=df.org
  
  simThr=quantile(df$x,0.8)
  df=df[df$x>=simThr,]
  
  ids=aggregate(x~i,data=df,sum)
  ids=ids[order(ids$x,decreasing = T),]
  
  df2=ids[sample(nrow(ids),nrow(ids)/10),]
  df2=df[df$i %in% df2$i,]
  df2=df2[!duplicated(df2$j),]
  df2=df2[!(df2$j %in% df2$i),]
  
  df2=as.data.frame(table(df2$i))
  df2=df2[order(df2$Freq,decreasing=T),]
  minNeighborCount=min(6,max(df2$Freq))
  df2=df2[df2$Freq>=minNeighborCount,]
  
  mappingList=NULL
  
  counterMain=0
  while(nrow(df2)>0){
    counterMain=counterMain+1
    if(counterMain %%100==5){
      if(verbose){
        print(paste0(counterMain," : ",nrow(ids)))
      }
      
    }
    
    tmpMapping.main=as.numeric(as.character(df2$Var1))
    tmpMapping.neighbors=df[df$i %in% tmpMapping.main,]
    tmpMapping.neighbors=tmpMapping.neighbors[!duplicated(tmpMapping.neighbors$j),]
    tmpMapping.neighbors=tmpMapping.neighbors[!(tmpMapping.neighbors$j %in% tmpMapping.neighbors$i),]
    tmpMapping.neighbors=tmpMapping.neighbors[order(tmpMapping.neighbors$x,decreasing = T),]
    tmpMapping.neighbors=aggregate(j~i,data=tmpMapping.neighbors,function(x) x[1:minNeighborCount])
    
    for(j in 1:minNeighborCount){
      mappingList=rbind(mappingList,data.frame(Fnode=tmpMapping.neighbors[,1],Snode=tmpMapping.neighbors[,2][,j]))
    }
    
    tmpMapping=c(tmpMapping.neighbors[,1],as.numeric(tmpMapping.neighbors[,2]))
    df=df[!df$i %in% tmpMapping &!df$j %in% tmpMapping,]
    
    ids=aggregate(x~i,data=df,sum)
    ids=ids[order(ids$x,decreasing = T),]
    counter=0
    while(counter<100){
      counter=counter+1
      df2=ids[sample(nrow(ids),nrow(ids)/10),]
      df2=df[df$i %in% df2$i,]
      df2=df2[!duplicated(df2$j),]
      df2=df2[!(df2$j %in% df2$i),]
      
      df2=as.data.frame(table(df2$i))
      df2=df2[order(df2$Freq,decreasing=T),]
      df2=df2[df2$Freq>=minNeighborCount,]
      if(nrow(df2)>0){
        break;
      }
      
    }
    
  }
  
  if(sum(duplicated(c(unique(mappingList$Fnode),mappingList$Snode)))>0){
    stop("error in the mapping file1!")
  }
  
  tst=as.data.frame(table(df$i))
  tst=tst[tst$Freq>=minNeighborCount,]
  if(nrow(tst)>0){
    tst=df[df$i %in% as.numeric(as.character(tst$Var1)),]
    tst=tst[!duplicated(tst$j),]
    tst=tst[!(tst$j %in% tst$i),]
    
    tmpMapping.neighbors=tst[order(tst$x,decreasing = T),]
    tmpMapping.neighbors=aggregate(j~i,data=tmpMapping.neighbors,function(x) x[1:minNeighborCount])
    
    for(j in 1:minNeighborCount){
      mappingList=rbind(mappingList,data.frame(Fnode=tmpMapping.neighbors[,1],Snode=tmpMapping.neighbors[,2][,j]))
    }
    mappingList=mappingList[!is.na(mappingList$Snode),]
    tmpMapping=c(tmpMapping.neighbors[,1],as.numeric(tmpMapping.neighbors[,2]))
    df=df[!df$i %in% tmpMapping &!df$j %in% tmpMapping,]
  }
  df=df.org[!df.org$i %in% c(mappingList$Fnode,mappingList$Snode),]
  
  if(sum(duplicated(c(unique(mappingList$Fnode),mappingList$Snode)))>0){
    stop("error in the mapping file2!")
  }
  
  if(length(minNeighborList)!=length(simThr_coef)){
    stop("minNeighborList doesn't match with the simThr_coef list")
  }
  
  for(iter in 1:length(minNeighborList)){
    if(nrow(df)>0){
      df2=df.org
      df2=df2[df2$x>=simThr*simThr_coef[iter],]
      df2=df2[df2$i %in% c(mappingList$Fnode,mappingList$Snode),]
      df2=df2[order(df2$x,decreasing = T),]
      
      for(i in intersect(unique(df$i),c(df2$j))){
        tmp=df2[df2$j==i,]
        
        tst=mappingList$Fnode[mappingList$Fnode %in% tmp$i|mappingList$Snode %in% tmp$i]
        tst=mappingList[mappingList$Fnode %in% tst,]
        tst=as.data.frame(table(tst$Fnode))
        tst=tst[tst$Freq==min(min(tst$Freq),minNeighborList[iter]),]
        if(nrow(tst)>0){
          tst=as.numeric(as.character(tst$Var1))
          tst2=tst[match(tmp$i,tst)]
          tst2=tst2[!is.na(tst2)]
          if(length(tst2)>0){
            tst=tst2[1]
          } else {
            tst=tst[1]
          }
          mappingList=rbind(mappingList,data.frame(Fnode=tst,Snode=i))
        }
        
        
      }
      
      df=df[!df$i %in% c(mappingList$Fnode,mappingList$Snode),]
    }
    if(sum(duplicated(c(unique(mappingList$Fnode),mappingList$Snode)))>0){
      print("error in the mapping file!")
      print(itr)
    }
  }
  
  tst=as.data.frame(table(df$i[!df$j %in% c(mappingList$Fnode,mappingList$Snode)]))
  
  tst=df[df$i %in% as.numeric(as.character(tst$Var1[tst$Freq>=minNeighborCount])),]
  tst=tst[!(duplicated(tst$j)|tst$j %in% tst$i|tst$j %in% c(mappingList$Fnode,mappingList$Snode)),]
  if(nrow(tst)>(minNeighborCount-1)){
    mappingList=rbind(mappingList,data.frame(Fnode=tst$i,Snode=tst$j))
  }
  df=df[!df$i %in% c(mappingList$Fnode,mappingList$Snode),]
  df=df[!df$j %in% c(mappingList$Fnode,mappingList$Snode),]
  if(sum(duplicated(c(unique(mappingList$Fnode),mappingList$Snode)))>0){
    stop("error in the mapping file3!")
  }
  
  coverage_res=round(length(c(unique(mappingList$Fnode),mappingList$Snode))/nrow(inputPCAembeddings),3)
  if(verbose){
    print(paste("Coverage:",coverage_res))
    
  }
  df=unique(df$i)
  
  
  tstMap1=aggregate(Snode~Fnode,data=mappingList,function(x) x[1:min(length(x),20)])
  tstMap=cbind(tstMap1[,1],tstMap1[,2])
  tstMap=apply(tstMap,1,function(x) list(x[!is.na(x)]))
  
  tstMap=c(tstMap,lapply(df,function(x) x))
  
  inputExpMat=counts(inputExpData)
  tst=lapply(tstMap,function(x) if(length(unlist(x))>1){rowSums(inputExpMat[,unlist(x)])} else {as.numeric(inputExpMat[,unlist(x)])})
  
  tst=do.call("cbind",tst)
  
  if(length(df)>0){
    tstColnames=c(paste0("group_",colnames(inputExpData)[tstMap1[,1]]),paste0("singleton_",colnames(inputExpData)[df]))
    
    colnames(tst)=tstColnames
  } else {
    tstColnames=paste0("group_",colnames(inputExpData)[tstMap1[,1]])
    
    colnames(tst)=tstColnames
  }
  
  
  return(list(mappingRes=mappingList,notMapped=df,collapsedExpData=tst,coverage=coverage_res))
}

.extraHarmony_dataset_integratorFn=function(dsSource,dsTarget,cellType_source,cellType_target,covariates=c("depth_per_gene"),calculate_depth_per_gene=T,source_batch_label=NULL,target_batch_label=NULL,inputVarGenes=NULL,ds_specific_hvg=T,indScaling=F){
  require(Seurat)
  require(ggplot2)
  require(ggalluvial)
  
  print("Combining the datasets...")
  
  if(any(row.names(dsSource)!=row.names(dsTarget))){
    stop("Row names should match between the two datasets")
  }
  
  ortholog_mapped=F
  if(sum(colnames(rowData(dsSource))=="mouse_ensembl_gene_id")>0|sum(colnames(rowData(dsTarget))=="mouse_ensembl_gene_id")>0){
    ortholog_mapped=T
    print("Ortholog mapping was detected")
  }
  
  
  if(sum(colnames(colData(dsSource))=="anno_batch")>0){
    dsSource$org_batch=dsSource$anno_batch
    colData(dsSource)=colData(dsSource)[,-which(colnames(colData(dsSource))=="anno_batch")]
    if(sum(colnames(colData(dsTarget))=="anno_batch")==0){
      dsTarget$org_batch="target"
    }
  }
  
  if(sum(colnames(colData(dsTarget))=="anno_batch")>0){
    dsTarget$org_batch=dsTarget$anno_batch
    colData(dsTarget)=colData(dsTarget)[,-which(colnames(colData(dsTarget))=="anno_batch")]
    if(sum(colnames(colData(dsSource))=="anno_batch")==0){
      dsSource$org_batch="source"
    }
  }
  
  data.list =list(source=dsSource,target=dsTarget)
  res=.mycBindFn(data.list,names(data.list))
  data.list =list(source=.extraExport2SeuratFn(dsSource),target=.extraExport2SeuratFn(dsTarget))

  
  res$batch_label=res$anno_batch
  res$batch_label2=res$anno_batch
  if(sum(colnames(colData(res))=="org_batch")>0){
    res$anno_batch=res$org_batch
  }
  
  if(!is.null(target_batch_label)){
    res$batch_label2[which(res$batch_label=="target")]=colData(res)[which(res$batch_label=="target"),target_batch_label]
  }
  
  if(!is.null(source_batch_label)){
    res$batch_label2[which(res$batch_label=="source")]=colData(res)[which(res$batch_label=="source"),source_batch_label]
  }
  
  
  res=.extraExport2SeuratFn(res)
  

  for (i in 1:length(data.list)) {
    data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
    data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }
  
  res <- NormalizeData(res, normalization.method = "LogNormalize", scale.factor = 10000)
  if(is.null(inputVarGenes)){
    if(ds_specific_hvg){
      var_genes=unique(c(data.list[[1]]@assays$RNA@var.features,data.list[[2]]@assays$RNA@var.features))
      res@assays$RNA@var.features=var_genes
    } else {
      res=FindVariableFeatures(res, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    }
  } else {
    res@assays$RNA@var.features=row.names(res)[which(toupper(row.names(res)) %in% toupper(inputVarGenes))]
  }
  
  if(ortholog_mapped){
    var_genes=as.character(res@assays$RNA@var.features)
    fd=res@assays$RNA@meta.features
    row.names(fd)=row.names(res@assays$RNA@data)
    fd=fd[match(var_genes,row.names(fd)),]
    
    #library(googleCloudStorageR)
    if(!file.exists("../data/gazestani_serverFiles/ortholog_mapping.rda")){
      stop("Error! file ../data/gazestani_serverFiles/ortholog_mapping.rda is missing")
    }
    load("../data/gazestani_serverFiles/ortholog_mapping.rda")
    
    var_genes=row.names(fd)[which(fd$ensembl_gene_id %in% as.character(mapping$human))]
    res@assays$RNA@var.features=var_genes
  }
  
  if(calculate_depth_per_gene){
    res$depth_per_gene=res$QC_Gene_total_count/(res$QC_Gene_unique_count+1)
  }
  
  
  
  if(indScaling){
    dataList=.mySplitObject(res,"batch_label2")
    mySeuratFn2=function(inputData,data){
      inputData@assays$RNA@var.features=data@assays$RNA@var.features
      inputData = ScaleData(inputData, verbose = FALSE,features=data@assays$RNA@var.features)
      return(inputData)
    }
    
    tmpInd=c()
    if(sum(!is.null(covariates))>0){
      for(i in covariates[!is.null(covariates)]){
        tmpInd=c(tmpInd,which(is.na(res@meta.data[,i])))
      }
      if(length(tmpInd)>0){
        print("Removing some cells with the covariate value of NA")
        res=res[,-unique(tmpInd)]
      }
    }
    
    
    dataList=parallel::mclapply(dataList,mySeuratFn2,data=res,mc.cores = length(dataList))
      
      resScaled=dataList[[1]]@assays$RNA@scale.data
      for(i in 2:length(dataList)){
        tmp=dataList[[i]]@assays$RNA@scale.data
        resScaled=cbind(resScaled,dataList[[i]]@assays$RNA@scale.data)
      }
      
      if(!all(colnames(res)==colnames(resScaled))){
        stop("Error in the matching!")
      }
      if(!all(row.names(res)[row.names(res) %in% res@assays$RNA@var.features]==row.names(resScaled))){
        stop("Error in the matching!")
      }
      
      res@assays$RNA@scale.data=resScaled
      
      res = RunPCA(res,features=res@assays$RNA@var.features,verbose = F)
      
  } else {
    res <- ScaleData(res, verbose = FALSE,features=res@assays$RNA@var.features,vars.to.regress =covariates)
    res <- RunPCA(res,verbose = F)
  }
  
  
  return(res)
  
}

.extraAligned_dataset_integratorFn=function(pca_source,pca_target,cellType_source,cellType_target){
  require(Seurat)
  require(ggplot2)
  require(ggalluvial)
  
  print("Combining the datasets...")
  
  
  
  
  return(res)
  
}

.myLabelTransfer_harmony=function(dataset_source,dataset_target,source_label_col,indScaling,return_seurat_obj=T,inputVarGenes=NULL,ds_specific_hvg=T,target_label_col=NULL,source_batch_label_col=NULL,target_batch_label_col=NULL,covariates=NULL,prenatalDataset=F,calculate_depth_per_gene=T,doubletGroups=NULL,nPCs=30,n.adaptiveKernel=5,nPropIter=3){
  
  #dataset_source: dataset whose labels will be transfered
  #dataset_target: dataset for which we are predicting the cell types 
  #source_label_col: the column that specifies the cell type in the source dataset
  #target_label_col: the column that specifies the cell type in the target dataset
  #nPCs: number of PCs for the construction of knn network
  #covariates: covariates to be regressed out from the dataset
  
  if(prenatalDataset){
    covariates=unique(c(covariates,"depth_per_gene"))
  }
  
  if(is.null(target_label_col)){
    dataset_target$target_lbl="unknown"
    target_label_col="target_lbl"
  }
  
  if(class(dataset_source)!="SingleCellExperiment"){
    stop("Source dataset should be in the format of SingleCellExperiment")
  }
  
  if(class(dataset_target)!="SingleCellExperiment"){
    stop("Target dataset should be in the format of SingleCellExperiment")
  }
  
  
  if(!is.null(doubletGroups)){
    dataset_source=.extraDoubletMakerFn(inputData=dataset_source,label_col=source_label_col,sel_labels=doubletGroups)
  }
  
  dataset_source$madeCluster=colData(dataset_source)[,source_label_col]
  dataset_target$madeCluster=colData(dataset_target)[,target_label_col]

  #dsSource=dataset_source;dsTarget=dataset_target;cellType_source=colData(dataset_source)[,source_label_col];cellType_target=colData(dataset_target)[,target_label_col];covariates=covariates;calculate_depth_per_gene=calculate_depth_per_gene;source_batch_label=source_batch_label_col;target_batch_label=target_batch_label_col;indScaling=indScaling
  res_seurat=.extraHarmony_dataset_integratorFn(dsSource=dataset_source,dsTarget=dataset_target,cellType_source=colData(dataset_source)[,source_label_col],indScaling = indScaling,cellType_target=colData(dataset_target)[,target_label_col],covariates=covariates,calculate_depth_per_gene=calculate_depth_per_gene,source_batch_label=source_batch_label_col,target_batch_label=target_batch_label_col,inputVarGenes=inputVarGenes,ds_specific_hvg=ds_specific_hvg)
  
  if(sum(grepl("^inferred_",colnames(res_seurat@meta.data)))>0){
    colnames(res_seurat@meta.data)=gsub("^inferred_","org_inferred_",colnames(res_seurat@meta.data))
  }
  
  meta_data=as.data.frame(res_seurat@meta.data)
  if(sum(is.na(meta_data$madeCluster))>0){
    meta_data$madeCluster[is.na(meta_data$madeCluster)]=""
  }
  

  label_col="madeCluster"
  training_idx=which(res_seurat$batch_label=="source")

  harmony_embeddings <- harmony::HarmonyMatrix(res_seurat@reductions$pca@cell.embeddings[,1:nPCs], meta_data, 'batch_label2', do_pca = FALSE, verbose=FALSE)

  res=.myKnnLabelTransferFn(inputPCAembeddings=harmony_embeddings,meta_data=meta_data,training_idx=training_idx,label_col=label_col,n.adaptiveKernel=n.adaptiveKernel,nPropIter=nPropIter)
  
  resTest=res$test_labels
  
  slCols=colnames(resTest)[which(grepl("^inferred_",colnames(resTest)))]


  plotTestDf=NULL
  for(i in slCols){
    tmp=data.frame(label=resTest[,'madeCluster'],inferred=gsub("inferred_","",i),weight=resTest[,i],stringsAsFactors = F)
    plotTestDf=rbind(plotTestDf,tmp)
  }
  
  plotTestDf2=aggregate(weight~label+inferred,data=plotTestDf,sum)
  
  plotTestDf2=plotTestDf2[order(plotTestDf2$weight,decreasing = T),]
  
  plotTrainingDf=NULL
  for(i in slCols){
    tmp=data.frame(label=res$training_labels[,'madeCluster'],inferred=gsub("inferred_","",i),weight=res$training_labels[,i],stringsAsFactors = F)
    plotTrainingDf=rbind(plotTrainingDf,tmp)
  }
  
  plotTrainingDf2=aggregate(weight~label+inferred,data=plotTrainingDf,sum)
  
  res=c(res,list(dSource=plotTrainingDf2,dTarget=plotTestDf2,source_label_col=source_label_col,harmony_embeddings=harmony_embeddings))
  
  if(return_seurat_obj){
    umap_res=.reductionUMAPFn(inputPCAembeddings = harmony_embeddings,umap.method = "uwot")
    umap_res=umap_res$embedding
    
    umap_res=umap_res[match(row.names(res$combined_labels), row.names(umap_res)),]
    if(sum(colnames(res$combined_labels) %in% c("UMAP_1","UMAP_2"))>0){
      res$combined_labels=res$combined_labels[,!grepl("UMAP",colnames(res$combined_labels))]
    }
    res$combined_labels=cbind(res$combined_labels,umap_res)
    sl_cols=colnames(res$combined_labels)[grepl("^inferred\\_",colnames(res$combined_labels))]
    res$combined_labels$anno=apply(res$combined_labels[,sl_cols],1,function(x) sl_cols[which(x==max(x,na.rm=T))[1]])
    
    arranged_data=Seurat::CreateSeuratObject(counts=res_seurat@assays$RNA@counts,
                                                 project = "SeuratProject",
                                                 assay = "RNA",
                                                 min.cells = 0,
                                                 min.features = 0,
                                                 names.field = 1,
                                                 names.delim = "-",
                                                 meta.data = as.data.frame(res$combined_labels))
    
    if(ncol(res_seurat@assays$RNA@meta.features)>0){
      arranged_data@assays$RNA@meta.features=res_seurat@assays$RNA@meta.features
    }
    
    umapdata=as.matrix(umap_res[,c("UMAP_1","UMAP_2")])
    umap.reduction <- CreateDimReducObject(embeddings = umapdata, 
                                           key = "UMAP_", assay = "RNA", global = TRUE)
    arranged_data[["umap"]]=umap.reduction
    
    if(sum(!grepl("^PC_",colnames(harmony_embeddings)))>0){
      colnames(harmony_embeddings)=gsub("PC_","",colnames(harmony_embeddings))
    }
    
    reduction.data <- CreateDimReducObject(embeddings = harmony_embeddings, assay = "RNA", key = "PC_")
    arranged_data[["pca"]] <- reduction.data
    
    res=c(res,list(seurat_obj=arranged_data))
    
  }
  
  return(res)
}

.myLabelTransfer_aligned=function(pca_source,meta_source,pca_target,meta_target,source_label_col,source_data=NULL,target_data=NULL,target_label_col=NULL,nPCs=NULL,n.adaptiveKernel=5,nPropIter=3,return_seurat_obj=T){
  
  #pca_source: aligned pc space for the source dataset
  #meta_source: meta data for the source dataset
  #pca_target: aligned pc space for the target dataset
  #meta_target: meta data for the target dataset
  #source_label_col: the column that specifies the cell type in the source dataset
  #target_label_col: the column that specifies the cell type in the target dataset
  #nPCs: number of PCs for the construction of knn network; if NULL, set to ncol(pca_source)
  
  if(ncol(pca_source)!=ncol(pca_target)){
    stop("There should be identical number of PCs between pca_source and pca_target")
  }
  
  if(return_seurat_obj){
    if(is.null(source_data)|is.null(target_data)){
      stop("source_data and target_data should be provided!")
    }
  }
  
  if(nrow(pca_source)!=nrow(meta_source)){
    stop("meta data and PC space for the source dataset don't match!")
  }
  
  if(nrow(pca_target)!=nrow(meta_target)){
    stop("meta data and PC space for the target dataset don't match!")
  }
  
  
  if(!all(row.names(pca_source)==row.names(meta_source))){
    stop("row.names don't match between meta data and PC space for the source dataset!")
  }
  
  if(!all(row.names(pca_target)==row.names(meta_target))){
    stop("row.names don't match between meta data and PC space for the target dataset!")
  }
  
  if(is.null(target_label_col)){
    meta_target$target_lbl="unknown"
    target_label_col="target_lbl"
  }
  
  if(sum(colnames(meta_source)==source_label_col)==0){
    stop("source_label_col could not be found!")
  }
  
  meta_source$madeCluster=meta_source[,source_label_col]
  meta_target$madeCluster=meta_target[,target_label_col]
  
  meta_source$batch_label="source"
  meta_target$batch_label="target"
  
  meta_data=plyr::rbind.fill(meta_source,meta_target)
  row.names(meta_data)=c(row.names(meta_source),row.names(meta_target))
  
  if(sum(is.na(meta_data$madeCluster))>0){
    meta_data$madeCluster[is.na(meta_data$madeCluster)]=""
  }
  
  if(sum(grepl("^inferred_",colnames(meta_data)))>0){
    colnames(meta_data)=gsub("^inferred_","org_inferred_",colnames(meta_data))
  }
  
  label_col="madeCluster"
  training_idx=which(meta_data$batch_label=="source")
  
  harmony_embeddings <- rbind(pca_source,pca_target)
  
  row.names(harmony_embeddings)=c(row.names(pca_source),row.names(pca_target))
  
  
  res=.myKnnLabelTransferFn(inputPCAembeddings=harmony_embeddings,meta_data=meta_data,training_idx=training_idx,label_col=label_col,n.adaptiveKernel=n.adaptiveKernel,nPropIter=nPropIter)
  
  resTest=res$test_labels
  
  slCols=colnames(resTest)[which(grepl("^inferred_",colnames(resTest)))]
  
  plotTestDf=NULL
  for(i in slCols){
    tmp=data.frame(label=resTest[,'madeCluster'],inferred=gsub("inferred_","",i),weight=resTest[,i],stringsAsFactors = F)
    plotTestDf=rbind(plotTestDf,tmp)
  }
  
  plotTestDf2=aggregate(weight~label+inferred,data=plotTestDf,sum)
  
  plotTestDf2=plotTestDf2[order(plotTestDf2$weight,decreasing = T),]
  
  plotTrainingDf=NULL
  for(i in slCols){
    tmp=data.frame(label=res$training_labels[,'madeCluster'],inferred=gsub("inferred_","",i),weight=res$training_labels[,i],stringsAsFactors = F)
    plotTrainingDf=rbind(plotTrainingDf,tmp)
  }
  
  plotTrainingDf2=aggregate(weight~label+inferred,data=plotTrainingDf,sum)
  
  
  res=c(res,list(dSource=plotTrainingDf2,dTarget=plotTestDf2,source_label_col=source_label_col,harmony_embeddings=harmony_embeddings))
  
  if(return_seurat_obj){
    umap_res=.reductionUMAPFn(inputPCAembeddings = harmony_embeddings,umap.method = "uwot")
    umap_res=umap_res$embedding
    
    umap_res=umap_res[match(row.names(res$combined_labels), row.names(umap_res)),]
    if(sum(colnames(res$combined_labels) %in% c("UMAP_1","UMAP_2"))>0){
      res$combined_labels=res$combined_labels[,!grepl("UMAP",colnames(res$combined_labels))]
    }
    res$combined_labels=cbind(res$combined_labels,umap_res)
    sl_cols=colnames(res$combined_labels)[grepl("^inferred\\_",colnames(res$combined_labels))]
    res$combined_labels$anno=apply(res$combined_labels[,sl_cols],1,function(x) sl_cols[which(x==max(x,na.rm=T))[1]])
    
    res_seurat=.extraExport2SeuratFn(.mycBindFn(list(source_data,target_data)))
    res_seurat=res_seurat[,match(row.names(res$combined_labels),colnames(res_seurat))]
    arranged_data=Seurat::CreateSeuratObject(counts=res_seurat@assays$RNA@counts,
                                             project = "SeuratProject",
                                             assay = "RNA",
                                             min.cells = 0,
                                             min.features = 0,
                                             names.field = 1,
                                             names.delim = "-",
                                             meta.data = as.data.frame(res$combined_labels))
    
    if(ncol(res_seurat@assays$RNA@meta.features)>0){
      arranged_data@assays$RNA@meta.features=res_seurat@assays$RNA@meta.features
    }
    
    umapdata=as.matrix(umap_res[,c("UMAP_1","UMAP_2")])
    umap.reduction <- CreateDimReducObject(embeddings = umapdata, 
                                           key = "UMAP_", assay = "RNA", global = TRUE)
    arranged_data[["umap"]]=umap.reduction
    
    if(sum(!grepl("^PC_",colnames(harmony_embeddings)))>0){
      colnames(harmony_embeddings)=gsub("PC_","",colnames(harmony_embeddings))
    }
    
    reduction.data <- CreateDimReducObject(embeddings = harmony_embeddings, assay = "RNA", key = "PC_")
    arranged_data[["pca"]] <- reduction.data
    
    res=c(res,list(seurat_obj=arranged_data))
    
  }
  
  return(res)
}

.myKnnLabelTransferFn=function(inputPCAembeddings,meta_data,training_idx,label_col,n.adaptiveKernel=5,nPropIter=3,NNmethod="annoy",n.trees=50){
  
  #training_idx: the row index of training data in the inputPCAembedings
  #label_col: the column in the meta_data that specifies the cell labels (labels of test cells can be set as unknown; the labels of test sets would be excluded from the propagation step)
  
  training_data=inputPCAembeddings[training_idx,]
  training_labels=meta_data[training_idx,]
  
  test_data=inputPCAembeddings[-training_idx,]
  test_labels=meta_data[-training_idx,]
  ref_ind=NULL
  if(NNmethod=="annoy"){
    ref_ind=Seurat:::AnnoyBuildIndex(data = training_data, metric = "euclidean", 
                                 n.trees = n.trees)
    nn.ranked.1=Seurat:::AnnoySearch(index = ref_ind, query = test_data,k=n.adaptiveKernel*12,include.distance = T,search.k = -1)
    
  } else {
    nn.ranked.1 <- RANN::nn2(data = training_data,query = test_data, k = n.adaptiveKernel*12, eps = 0)
  }
  
  
  nn.ranked.1$nn.dists=cbind(rep(0,nrow(nn.ranked.1$nn.dists)),nn.ranked.1$nn.dists)
  nn.ranked.1$nn.idx=cbind(1:nrow(nn.ranked.1$nn.idx)+nrow(training_data),nn.ranked.1$nn.idx)
  nn.ranked.1$nn.dists=nn.ranked.1$nn.dists[,1:(n.adaptiveKernel*12)]
  nn.ranked.1$nn.idx=nn.ranked.1$nn.idx[,1:(n.adaptiveKernel*12)]
  
  if(NNmethod=="annoy"){
    idx=Seurat:::AnnoyBuildIndex(data = test_data, metric = "euclidean", 
                                 n.trees = n.trees)
    nn.ranked.purityIndex=Seurat:::AnnoySearch(index = idx, query = test_data,k=n.adaptiveKernel*3,include.distance = T,search.k = -1)
    
  } else {
    nn.ranked.purityIndex=RANN::nn2(data = test_data, k = n.adaptiveKernel*3, eps = 0)
  }
  
  nn.ranked.purityIndex=nn.ranked.purityIndex$nn.dists[,n.adaptiveKernel*3]
  
  if(NNmethod=="annoy"){
    nn.ranked.2=Seurat:::AnnoySearch(index = ref_ind, query = training_data,k=n.adaptiveKernel*12,include.distance = T,search.k = -1)
    
  } else {
    nn.ranked.2 <- RANN::nn2(data = training_data,query = training_data, k = n.adaptiveKernel*12, eps = 0)
    
  }
  nn.training.purityIndx=nn.ranked.2$nn.dists[,n.adaptiveKernel]
  
  test_labels$status="test_set"
  training_labels$status="training_set"
  labels=rbind(training_labels,test_labels)
  rm(test_labels,training_labels,test_data,training_data)
  
  dists=nn.ranked.1$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/nn.ranked.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
  affinities2=affinities-affinitiesThr
  affinities[affinities2<=0]=0
  #removing the outlier matches
  toRM= which(scale(affinities[,2])<(-3))
  if(length(toRM)>0){
    for(i in 1:ncol(affinities)){
      affinities[toRM,i]=0
    }
  }
  nn.ranked.1$affinities=affinities
  rm(affinities)
  
  dists=nn.ranked.2$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/nn.training.purityIndx[x])^2))
  affinities=do.call("rbind",affinities)
  affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
  affinities2=affinities-affinitiesThr
  affinities[affinities2<=0]=0
  nn.ranked.2$affinities=affinities
  rm(affinities)
  
  nn.ranked.org=nn.ranked.2
  nn.ranked.org$nn.idx=rbind(nn.ranked.2$nn.idx,nn.ranked.1$nn.idx)
  nn.ranked.org$nn.dists=rbind(nn.ranked.2$nn.dists,nn.ranked.1$nn.dists)
  nn.ranked.org$affinities=rbind(nn.ranked.2$affinities,nn.ranked.1$affinities)
  rm(nn.ranked.1,nn.ranked.2)
  
  
  nn.ranked=nn.ranked.org
  affinities=nn.ranked.org$affinities
  
  affCounts=apply(affinities,2,function(x) sum(x==0))
  affinities=affinities[,-which(affCounts==nrow(affinities))]
  nn.ranked$nn.idx=nn.ranked$nn.idx[,-which(affCounts==nrow(affinities))]
  j <- as.numeric(t(nn.ranked$nn.idx))
  i <- ((1:length(j)) - 1)%/%ncol(affinities) + 1
  k=as.numeric(t(affinities))
  
  graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(inputPCAembeddings), nrow(inputPCAembeddings)))
  rm(i,j,k)
  graph <- Matrix::Diagonal(x = 1 / (rowSums(graph)+0.0001)) %*% graph
  
  if(class(graph)!="dgCMatrix"){
    graph=as(graph,"dgCMatrix")
  }
  
  rownames(graph) <- rownames(inputPCAembeddings)
  colnames(graph) <- rownames(inputPCAembeddings)
  
  for(i in 1:nPropIter){
    print(i)
    if(i==1){
      prop_mat=graph
    } else {
      prop_mat=prop_mat%*%graph
    }
  }
  
  
  tmpLables=as.data.frame(labels)[,label_col]
  tmpLables[labels$status=="test_set"]="unknown_toexclude"
  one_hot_labels=.myOneHotFn(tmpLables)
  training_col_labels=unique(as.data.frame(labels)[training_idx,label_col])
  if(length(which(colnames(one_hot_labels)=="unknown_toexclude"))>0){
    one_hot_labels=one_hot_labels[,-which(colnames(one_hot_labels)=="unknown_toexclude")]
  }
  
  
  res=as.matrix(prop_mat %*% as.matrix(one_hot_labels))
  res=round(res,3)
  
  res2=as.data.frame(res)
  colnames(res)=paste0("inferred_",colnames(res))
  
  colnames(res2)=colnames(res)
  res=cbind(labels,res2)
  
  res_testset=res[res$status=="test_set",]
  res_training=res[res$status=="training_set",]
  
  return(list(test_labels=res_testset,training_labels=res_training,combined_labels=res))
}

.myKnnLabelTransferFn_randWalkAffinities=function(inputPCAembeddings,meta_data,training_idx,label_col,n.adaptiveKernel=5,nPropIter=3,NNmethod="annoy",n.trees=50){
  
  #training_idx: the row index of training data in the inputPCAembedings
  #label_col: the column in the meta_data that specifies the cell labels (labels of test cells can be set as unknown; the labels of test sets would be excluded from the propagation step)
  
  training_data=inputPCAembeddings[training_idx,]
  training_labels=meta_data[training_idx,]
  
  test_data=inputPCAembeddings[-training_idx,]
  test_labels=meta_data[-training_idx,]
  ref_ind=NULL
  myL2normFn=function(inputMat){
    prop_mat2=rowSums(inputMat^2)
    prop_mat2=sqrt(prop_mat2)
    res=Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
    return(res)
  }
  
  if(is.null(batchPCAdata)){
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    batchPCAdata=harmony_embeddings
  }
  
  if(sum(argList$singleton.method %in% c("snn","fast"))==0){
    stop("unrecognized singleton.method")
  }
  
  n.neighbors=20
  
  batchPCAdata=as.matrix(batchPCAdata)[,1:argList$nPCs]
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  if(NNmethod=="annoy"){
    if(!file.exists(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))){
      idx=Seurat:::AnnoyBuildIndex(data = batchPCAdata, metric = "euclidean", 
                                   n.trees = n.trees)
      nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = batchPCAdata,k=k.param,include.distance = T,search.k = -1)
      qsave(nn.ranked.1,file=.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
    } else {
      nn.ranked.1=qread(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
    }
  } else {
    nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  }
  
  
  affinities=Matrix::Diagonal(x=1/(nn.ranked.1$nn.dists[,2]+0.000001)) %*% nn.ranked.1$nn.dists
  affinities@x=-1*affinities@x^2
  affinities@x=exp(affinities@x)
  affinities[,1]=affinities[,2]
  j <- as.numeric(t(nn.ranked.1$nn.idx))
  i <- ((1:length(j)) - 1)%/%n.neighbors + 1
  x=as.numeric(t(affinities))
  adj = sparseMatrix(i = i, j = j, x = x, dims = c(nrow(batchPCAdata),nrow(batchPCAdata)))
  rownames(adj) <- row.names(batchPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  adj= Matrix::Diagonal(x=1/rowSums(adj)) %*% adj
  
  
  
  load(.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
  adj_t=t(adj)
  adj_t=Matrix::Diagonal(x=1/rowSums(adj_t)) %*% adj_t
  adj=adj[sl_pseudo$pseudocell,]
  row.names(adj)=sl_pseudo$cluster
  adj=adj %*% adj_t
  
  
  
  adj <- Matrix::Diagonal(x = 1 / rowSums(adj)) %*% adj
  
  prop_mat=centroid.affinities=adj
  
  if(F){
    colMax_vals=c()
    for(i in seq(1,ncol(prop_mat),5000)){
      tmp_max=as.numeric(qlcMatrix::colMax(prop_mat[,i:min(i+4999,ncol(prop_mat))]))
      colMax_vals=c(colMax_vals,as.numeric(tmp_max))
    }
    prop_mat =  prop_mat %*% Matrix::Diagonal(x = 1 / colMax_vals)
    prop_mat <- Matrix::Diagonal(x = 1 / rowSums(prop_mat)) %*% prop_mat
  }
  
  
  prop_mat2=myL2normFn(inputMat=prop_mat)
  c_c_aff=t(prop_mat2)
  c_c_aff=prop_mat2 %*% c_c_aff
  #c_c_aff=c_c_aff %*%c_c_aff %*% c_c_aff
  
  #c_c_aff=Matrix::drop0(c_c_aff,tol=quantile(c_c_aff@x,0.1))
  c_c_aff <- Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff_t=c_c_aff2=c_c_aff
  
  inc_check=T
  inc_count=nrow(prop_mat)
  counter=0
  while((inc_check|counter<10)&counter<15){
    
    counter=counter+1
    
    c_c_aff=c_c_aff %*% c_c_aff_t*0.7+c_c_aff_t *0.3
    
    c_c_aff=Matrix::drop0(c_c_aff,tol=quantile(c_c_aff@x,0.00001))
    c_c_aff <- Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
    
    tst=as.numeric(qlcMatrix::rowMax(c_c_aff))
    tst2=as.numeric(diag(c_c_aff))
    tst=as.numeric(qlcMatrix::rowMax(c_c_aff))
    tst=Matrix::Diagonal(x=1/tst) %*% c_c_aff
    tst=which(tst>0.9,arr.ind = T)
    sl_ps_list=unique(tst[,2])
    #print(paste(counter,":",length(sl_ps_list)))
    if(length(sl_ps_list)/inc_count>0.95&length(sl_ps_list)<2000){
      inc_check=F
    }
    inc_count=length(sl_ps_list)
  }
  
  if(F){
    
    load(.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
    
    length(unique(pd$anno_cellState[pd$sample %in% sl_pseudo$pseudocell]))
    setdiff(pd$anno_cellState,unique(pd$anno_cellState[pd$sample %in% sl_pseudo$pseudocell]))
    
    
    #adj2=diff2
    
    adj2=prop_mat#=res$prop_mat
    #adj2=t(adj_all)
    
    adj2=Matrix::Diagonal(x=1/rowSums(adj2)) %*% adj2
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    pd=pd[colnames(adj2),]
    
    table(is.na(pd$anno_cellState))
    tmp= adj2 %*% as.matrix(.myOneHotFn(pd$anno_cellState))
    
    tmp=apply(tmp,1,function(x){
      y=which(x==max(x))[1]
      tmp2=data.frame(cluster=colnames(tmp)[y],purity=x[y])
      tmp2
    })
    
    tmp=do.call("rbind",tmp)
    if(F){
      tmp2=aggregate(purity~cluster,data=tmp,median)
      table(tmp2[,2]>0.8)
      table(tmp2[,2]>0.7)
      table(tmp2[,2]>0.6)
      table(tmp2[,2]>0.5)
    }
    
    matEffectiveSize=.myEffSizePropMat(adj2)
    matEffectiveSize=matEffectiveSize$effective_sample_size
    matEffectiveSize=matEffectiveSize[match(row.names(tmp),names(matEffectiveSize))]
    tmp$size=matEffectiveSize
    tmp2=NULL
    for(itrclust in unique(tmp$cluster)){
      tmp2=rbind(tmp2,data.frame(cluster=itrclust,purity=sum(tmp$size[tmp$cluster==itrclust]*tmp$purity[tmp$cluster==itrclust])/sum(tmp$size[tmp$cluster==itrclust]),stringsAsFactors = F))
    }
    table(tmp2[,2]>0.8)
    table(tmp2[,2]>0.7)
    table(tmp2[,2]>0.6)
    table(tmp2[,2]>0.5)
    summary(tmp2[,2])
  }
  
  
  c_c_aff=c_c_aff[sl_ps_list,]
  #c_c_aff=Matrix::Diagonal(x=1/qlcMatrix::rowMax(c_c_aff)) %*% c_c_aff
  #c_c_aff=Matrix::drop0(c_c_aff,tol=0.05)
  c_c_aff <- Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  
  #c_c_aff=Matrix::drop0(c_c_aff,tol=0.05)
  prop_mat_list=list()
  rowMax_vals=c()
  for(i in seq(1,ncol(prop_mat),50000)){
    tmp=c_c_aff %*% prop_mat[colnames(c_c_aff),i:min(i+49999,ncol(prop_mat))]
    tmp_max=as.numeric(qlcMatrix::rowMax(tmp))
    prop_mat_list=c(prop_mat_list,list(tmp))
    if(length(rowMax_vals)>0){
      rowMax_vals=pmax(rowMax_vals,tmp_max)
    } else {
      rowMax_vals=tmp_max
    }
  }
  
  for(i in 1:length(prop_mat_list)){
    tmp=prop_mat_list[[i]]
    #tmp=tmp %*% Matrix::Diagonal(x=1/(qlcMatrix::colMax(tmp)+0.00001))
    tmp2=Matrix::Diagonal(x = 1 / rowMax_vals) %*% tmp
    tmp2=Matrix::drop0(tmp2,tol=0.05)
    tmp2@x=rep(1,length(tmp2@x))
    tmp=tmp*tmp2
    
    #tmp=Matrix::drop0(tmp,tol=0.1)
    prop_mat_list[[i]]=tmp
  }
  
  rm(tmp,tmp2)
  
  prop_mat=prop_mat_list[[1]]
  prop_mat=prop_mat[,colSums(prop_mat)>0]
  if(length(prop_mat_list)>1){
    for(i in 2:length(prop_mat_list)){
      tmp=prop_mat_list[[i]]
      tmp=tmp[,colSums(tmp)>0]
      prop_mat=cbind(prop_mat,tmp)
    }
  }
  
  
  
  tmpLables=as.data.frame(labels)[,label_col]
  tmpLables[labels$status=="test_set"]="unknown_toexclude"
  one_hot_labels=.myOneHotFn(tmpLables)
  training_col_labels=unique(as.data.frame(labels)[training_idx,label_col])
  if(length(which(colnames(one_hot_labels)=="unknown_toexclude"))>0){
    one_hot_labels=one_hot_labels[,-which(colnames(one_hot_labels)=="unknown_toexclude")]
  }
  
  
  res=as.matrix(prop_mat %*% as.matrix(one_hot_labels))
  res=round(res,3)
  
  res2=as.data.frame(res)
  colnames(res)=paste0("inferred_",colnames(res))
  
  colnames(res2)=colnames(res)
  res=cbind(labels,res2)
  
  res_testset=res[res$status=="test_set",]
  res_training=res[res$status=="training_set",]
  
  return(list(test_labels=res_testset,training_labels=res_training,combined_labels=res))
}


.myKnnLabelTransferFn_list=function(inputPCAembeddings,meta_data,training_idx,label_col_list,n.adaptiveKernel=5,nPropIter=3){
  
  #transfers a list of covariates from the reference dataset to the query dataset
  #training_idx: the row index of training data in the inputPCAembedings
  #label_col: the column in the meta_data that specifies the cell labels (labels of test cells can be set as unknown; the labels of test sets would be excluded from the propagation step)
  
  training_data=inputPCAembeddings[training_idx,]
  training_labels=meta_data[training_idx,]
  
  test_data=inputPCAembeddings[-training_idx,]
  test_labels=meta_data[-training_idx,]
  
  nn.ranked.1 <- RANN::nn2(data = training_data,query = test_data, k = n.adaptiveKernel*12, eps = 0)
  nn.ranked.1$nn.dists=cbind(rep(0,nrow(nn.ranked.1$nn.dists)),nn.ranked.1$nn.dists)
  nn.ranked.1$nn.idx=cbind(1:nrow(nn.ranked.1$nn.idx)+nrow(training_data),nn.ranked.1$nn.idx)
  nn.ranked.1$nn.dists=nn.ranked.1$nn.dists[,1:(n.adaptiveKernel*12)]
  nn.ranked.1$nn.idx=nn.ranked.1$nn.idx[,1:(n.adaptiveKernel*12)]
  
  nn.ranked.purityIndex=RANN::nn2(data = test_data, k = n.adaptiveKernel*3, eps = 0)
  nn.ranked.purityIndex=nn.ranked.purityIndex$nn.dists[,n.adaptiveKernel*3]
  
  nn.ranked.2 <- RANN::nn2(data = training_data,query = training_data, k = n.adaptiveKernel*12, eps = 0)
  nn.training.purityIndx=nn.ranked.2$nn.dists[,n.adaptiveKernel]
  
  test_labels$status="test_set"
  training_labels$status="training_set"
  labels=rbind(training_labels,test_labels)
  rm(test_labels,training_labels,test_data,training_data)
  
  dists=nn.ranked.1$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/nn.ranked.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
  affinities2=affinities-affinitiesThr
  affinities[affinities2<=0]=0
  #removing the outlier matches
  toRM= which(scale(affinities[,2])<(-3))
  if(length(toRM)>0){
    for(i in 1:ncol(affinities)){
      affinities[toRM,i]=0
    }
  }
  nn.ranked.1$affinities=affinities
  rm(affinities)
  
  dists=nn.ranked.2$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/nn.training.purityIndx[x])^2))
  affinities=do.call("rbind",affinities)
  affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
  affinities2=affinities-affinitiesThr
  affinities[affinities2<=0]=0
  nn.ranked.2$affinities=affinities
  rm(affinities)
  
  nn.ranked.org=nn.ranked.2
  nn.ranked.org$nn.idx=rbind(nn.ranked.2$nn.idx,nn.ranked.1$nn.idx)
  nn.ranked.org$nn.dists=rbind(nn.ranked.2$nn.dists,nn.ranked.1$nn.dists)
  nn.ranked.org$affinities=rbind(nn.ranked.2$affinities,nn.ranked.1$affinities)
  rm(nn.ranked.1,nn.ranked.2)
  
  
  nn.ranked=nn.ranked.org
  affinities=nn.ranked.org$affinities
  
  affCounts=apply(affinities,2,function(x) sum(x==0))
  affinities=affinities[,-which(affCounts==nrow(affinities))]
  nn.ranked$nn.idx=nn.ranked$nn.idx[,-which(affCounts==nrow(affinities))]
  j <- as.numeric(t(nn.ranked$nn.idx))
  i <- ((1:length(j)) - 1)%/%ncol(affinities) + 1
  k=as.numeric(t(affinities))
  
  graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(inputPCAembeddings), nrow(inputPCAembeddings)))
  rm(i,j,k)
  graph <- Matrix::Diagonal(x = 1 / (rowSums(graph))) %*% graph
  
  if(class(graph)!="dgCMatrix"){
    graph=as(graph,"dgCMatrix")
  }
  
  rownames(graph) <- rownames(inputPCAembeddings)
  colnames(graph) <- rownames(inputPCAembeddings)
  
  for(i in 1:nPropIter){
    print(i)
    if(i==1){
      prop_mat=graph
    } else {
      prop_mat=prop_mat%*%graph
    }
  }
  
  res_testset=NULL
  for(ilbl in label_col_list){
    tmpLables=as.data.frame(labels)[,ilbl]
    tmpLables[labels$status=="test_set"]="unknown_toexclude"
    one_hot_labels=.myOneHotFn(tmpLables)
    training_col_labels=unique(as.data.frame(labels)[training_idx,ilbl])
    if(length(which(colnames(one_hot_labels)=="unknown_toexclude"))>0){
      one_hot_labels=one_hot_labels[,-which(colnames(one_hot_labels)=="unknown_toexclude")]
    }
    
    
    tmp_res=as.matrix(prop_mat %*% as.matrix(one_hot_labels))
    tmp_res=round(tmp_res,3)
    
    tmp_res2=as.data.frame(tmp_res)
    colnames(tmp_res)=paste0("inferred_",colnames(tmp_res))
    
    colnames(tmp_res2)=colnames(tmp_res)
    tmp_res=as.data.frame(tmp_res2)
    tmp_res$true_label=labels[,ilbl]
    tmp_res$status=labels[,"status"]
    
    tmp_res_testset=tmp_res[tmp_res$status=="test_set",]
    tmp_res_training=tmp_res[tmp_res$status=="training_set",]
    
    colnames(tmp_res_testset)=paste0(ilbl,"_",colnames(tmp_res_testset))
    colnames(tmp_res_training)=paste0(ilbl,"_",colnames(tmp_res_training))
    colnames(tmp_res)=paste0(ilbl,"_",colnames(tmp_res))
    
    if(!is.null(res_testset)){
      res_testset=cbind(res_testset,tmp_res_testset)
      res_training=cbind(res_training,tmp_res_training)
      res=cbind(res,tmp_res)
      
    } else {
      res_testset=tmp_res_testset
      res_training=tmp_res_training
      res=tmp_res
      
    }
    
  }
  
  return(list(test_labels=res_testset,training_labels=res_training,combined_labels=res))
}

.myLabelTransfer_liger=function(dataset_source,dataset_target,source_label_col,target_label_col=NULL,doubletGroups=NULL,nPCs=30,n.adaptiveKernel=5,nPropIter=3){
  
  require(rliger)
  require(Matrix)
  require(patchwork)
  require(ggplot2)
  require(ggalluvial)
  
  if(any(row.names(dataset_source)!=row.names(dataset_target))){
    stop("Row names should match between the two datasets")
  }
  
  if(class(dataset_source)!="SingleCellExperiment"){
    stop("Source dataset should be in the format of SingleCellExperiment")
  }
  
  if(class(dataset_target)!="SingleCellExperiment"){
    stop("Target dataset should be in the format of SingleCellExperiment")
  }
  
  if(is.null(target_label_col)){
    dataset_target$target_lbl="unknown"
    target_label_col="target_lbl"
  }
  
  
  if(!is.null(doubletGroups)){
    dataset_source=.extraDoubletMakerFn(inputData=dataset_source,label_col=source_label_col,sel_labels=doubletGroups)
  }
  
  fdSource=as.data.frame(rowData(dataset_source))
  fdTarget=as.data.frame(rowData(dataset_target))
  if(length(which(!colnames(fdTarget) %in% colnames(fdSource)))>1){
    fd=cbind(fdSource,fdTarget[,which(!colnames(fdTarget) %in% colnames(fdSource))])
  } else {
    fd=fdSource
  }
  
  row.names(fd)=rowData(dataset_source)$ensembl_gene_id
  
  
  pdSource=as.data.frame(colData(dataset_source))
  pdTarget=as.data.frame(colData(dataset_target))
  meta_data=plyr::rbind.fill(pdSource,pdTarget)
  meta_data$madeCluster=c(colData(dataset_source)[,source_label_col],cellType_target=colData(dataset_target)[,target_label_col])
  row.names(meta_data)=c(colnames(dataset_source),colnames(dataset_target))
  
  ligerex = rliger::createLiger(list(training = counts(dataset_source), test = counts(dataset_target))) #Can also pass in more than 2 datasets
  ligerex = rliger::normalize(ligerex)
  ligerex = rliger::selectGenes(ligerex, var.thresh = 0.1)
  ligerex = rliger::scaleNotCenter(ligerex)
  
  ligerex = rliger::optimizeALS(ligerex, k = nPCs) 
  ligerex = rliger::quantile_norm(ligerex) #SNF clustering and quantile alignment
  
  pcaData=ligerex@H.norm
  meta_data=meta_data[match(row.names(pcaData),row.names(meta_data)),]
  
  training_idx=which(row.names(meta_data) %in% colnames(dataset_source))
  
  res=.myKnnLabelTransferFn(inputPCAembeddings=pcaData,meta_data=meta_data,training_idx=training_idx,label_col="madeCluster",n.adaptiveKernel=n.adaptiveKernel,nPropIter=nPropIter)
  
  resTest=res$test_labels
  
  slCols=colnames(resTest)[which(grepl("inferred_",colnames(resTest)))]
  
  plotTestDf=NULL
  for(i in slCols){
    tmp=data.frame(label=resTest[,'madeCluster'],inferred=gsub("inferred_","",i),weight=resTest[,i],stringsAsFactors = F)
    plotTestDf=rbind(plotTestDf,tmp)
  }
  
  plotTestDf2=aggregate(weight~label+inferred,data=plotTestDf,sum)
  
  plotTrainingDf=NULL
  for(i in slCols){
    tmp=data.frame(label=res$training_labels[,'madeCluster'],inferred=gsub("inferred_","",i),weight=res$training_labels[,i],stringsAsFactors = F)
    plotTrainingDf=rbind(plotTrainingDf,tmp)
  }
  
  plotTrainingDf2=aggregate(weight~label+inferred,data=plotTrainingDf,sum)
  
  return(c(res,list(dSource=plotTrainingDf2,dTarget=plotTestDf2,source_label_col=source_label_col,liger_embeddings=ligerex@H.norm)))
}

.myLabelTransfer_seurat=function(dataset_source,dataset_target,source_label_col,target_label_col=NULL,doubletGroups=NULL,nPCs=30){
  
  require(Seurat)
  require(ggplot2)
  require(ggalluvial)
  
  if(any(row.names(dataset_source)!=row.names(dataset_target))){
    stop("Row names should match between the two datasets")
  }
  
  if(class(dataset_source)!="SingleCellExperiment"){
    stop("Source dataset should be in the format of SingleCellExperiment")
  }
  
  if(class(dataset_target)!="SingleCellExperiment"){
    stop("Target dataset should be in the format of SingleCellExperiment")
  }
  
  if(is.null(target_label_col)){
    target_label_col=rep("unknown",ncol(dataset_target))
  }
  
  if(!is.null(doubletGroups)){
    dataset_source=.extraDoubletMakerFn(inputData=dataset_source,label_col=source_label_col,sel_labels=doubletGroups)
  }
  
  ds1=.extraExport2SeuratFn(dataset_source)
  ds2=.extraExport2SeuratFn(dataset_target)
  
  ds1 <- NormalizeData(ds1, verbose = FALSE)
  ds1 <- FindVariableFeatures(ds1, selection.method = "vst", 
                                         nfeatures = 2000, verbose = FALSE)
  
  ds2 <- NormalizeData(ds2, verbose = FALSE)
  ds2 <- FindVariableFeatures(ds2, selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
  
  
  reference.list <- list(source=ds1,target=ds2)
  data.anchors <- FindTransferAnchors(reference = ds1, query = ds2, dims = 1:nPCs)
  
  predictions <- TransferData(anchorset = data.anchors, refdata = as.data.frame(colData(dataset_source))[,source_label_col], 
                              dims = 1:nPCs)
  query <- cbind(as.data.frame(colData(dataset_target)),predictions=predictions)
  
  colnames(query)=gsub("predictions\\.prediction\\.score\\.","inferred_",colnames(query))
  
  res=list(test_labels=query,training_labels=as.data.frame(colData(dataset_source)),combined_labels="")
  resTest=res$test_labels
  
  slCols=colnames(resTest)[which(grepl("inferred_",colnames(resTest)))]
  slCols=setdiff(slCols,"inferred_max")
  
  plotTestDf=NULL
  for(i in slCols){
    tmp=data.frame(label=resTest[,target_label_col],inferred=gsub("inferred_","",i),weight=resTest[,i],stringsAsFactors = F)
    plotTestDf=rbind(plotTestDf,tmp)
  }
  
  plotTestDf2=aggregate(weight~label+inferred,data=plotTestDf,sum)
  
  return(c(res,list(dTarget=plotTestDf2,source_label_col=source_label_col)))
}

.myMapToHuman_org=function(inputData,server=F){
  
  {
    #library(googleCloudStorageR)
    if(!file.exists("../data/gazestani_serverFiles/mouse_human_orthologs.rda")){
      stop("Error! file ../data/gazestani_serverFiles/mouse_human_orthologs.rda is missing")
    }
    
    load("!/data/gazestani_serverFiles/mouse_human_orthologs.rda")
  }
  
  
  #inputData=inputData[which(rowData(inputData)$ensembl_gene_id %in% mapping$mouse),]
  if(sum(rowData(inputData)$ensembl_gene_id %in% mapping$mouse,na.rm=T)==0){
    stop("ensembl_gene_ids were not found in the input data!")
  }
  
  inputData=inputData[!duplicated(rowData(inputData)$ensembl_gene_id),]
  
  mapping=mapping[which(mapping$mouse %in% rowData(inputData)$ensembl_gene_id),]
  mapping=mapping[!duplicated(mapping$mouse),]
  
  mapping=mapping[match(rowData(inputData)$ensembl_gene_id,mapping$mouse),]
  if(sum(is.na(mapping$human))>0){
    mapping$human[is.na(mapping$human)]=rowData(inputData)$ensembl_gene_id[is.na(mapping$human)]
  }
  
  if(sum(is.na(mapping$mouse))>0){
    mapping$mouse[is.na(mapping$mouse)]=rowData(inputData)$ensembl_gene_id[is.na(mapping$mouse)]
  }
  
  
  if(any(mapping$mouse!=rowData(inputData)$ensembl_gene_id)){
    stop("Error in the mapping!")
  }
  
  row.names(inputData)=mapping$human
  rowData(inputData)$mouse_ensembl_gene_id=rowData(inputData)$ensembl_gene_id
  rowData(inputData)$ensembl_gene_id=mapping$human
  
  return(inputData)
}

.myMapToHuman=function(inputData,server=F,conservativeMapping=F,oldMapping=F,MGIsymbol=F){
  
  {
    #library(googleCloudStorageR)
    if(conservativeMapping|oldMapping){
      if(!file.exists("../data/gazestani_serverFiles/ortholog_mapping_old.rda")){
        stop("Error! file ../data/gazestani_serverFiles/ortholog_mapping_old.rda is missing")
      }
      load("../data/gazestani_serverFiles/ortholog_mapping.rda")
      mapping_old=mapping
    } else if(MGIsymbol){
      if(!file.exists("../data/gazestani_serverFiles/mouse_human_orthologs_MGISymbol_03062021.rda")){
        stop("Error! file ../data/gazestani_serverFiles/mouse_human_orthologs_MGISymbol_03062021.rda is missing")
      }
      load("../data/gazestani_serverFiles/mouse_human_orthologs_MGISymbol_03062021.rda")
      mapping_old=mapping
    }
    
    if(oldMapping|MGIsymbol){
      mapping=mapping_old
    } else {
      if(!file.exists("../data/gazestani_serverFiles/ortholog_mapping.rda")){
        stop("Error! file ../data/gazestani_serverFiles/ortholog_mapping.rda is missing")
      }
      load("../data/gazestani_serverFiles/ortholog_mapping.rda")
      
      if(conservativeMapping){
        mapping=mapping[which(paste0(mapping$human,"_",mapping$mouse) %in% paste0(mapping_old$human,"_",mapping_old$mouse)),]
      }
    }
    
  }
  
  
  #inputData=inputData[which(rowData(inputData)$ensembl_gene_id %in% mapping$mouse),]
  if(sum(rowData(inputData)$ensembl_gene_id %in% mapping$mouse,na.rm=T)==0){
    stop("ensembl_gene_ids were not found in the input data!")
  }
  
  rowData(inputData)$ensembl_gene_id=as.character(rowData(inputData)$ensembl_gene_id)
  rowData(inputData)$ensembl_gene_id[is.na(rowData(inputData)$ensembl_gene_id)]=toupper(row.names(inputData))[is.na(rowData(inputData)$ensembl_gene_id)]
  if(sum(is.na(rowData(inputData)$ensembl_gene_id))>0){
    rowData(inputData)$ensembl_gene_id[is.na(rowData(inputData)$ensembl_gene_id)]=row.names(inputData)[is.na(rowData(inputData)$ensembl_gene_id)]
  }
  
  inputData=inputData[!duplicated(rowData(inputData)$ensembl_gene_id),]
  
  mapping=mapping[which(mapping$mouse %in% rowData(inputData)$ensembl_gene_id),]
  mapping=mapping[!duplicated(mapping$mouse),]
  
  mapping=mapping[match(rowData(inputData)$ensembl_gene_id,mapping$mouse),]
  if(sum(is.na(mapping$mouse))>0){
    mapping$mouse[is.na(mapping$mouse)]=rowData(inputData)$ensembl_gene_id[is.na(mapping$mouse)]
  }
  if(sum(is.na(mapping$human))>0){
    mapping$human[is.na(mapping$human)]=rowData(inputData)$ensembl_gene_id[is.na(mapping$human)]
  }
  
  
  
  
  if(any(mapping$mouse!=rowData(inputData)$ensembl_gene_id)){
    stop("Error in the mapping!")
  }
  
  row.names(inputData)=mapping$human
  rowData(inputData)$mouse_ensembl_gene_id=rowData(inputData)$ensembl_gene_id
  rowData(inputData)$ensembl_gene_id=mapping$human
  
  return(inputData)
}

.myIDMapToHuman=function(inputGeneNames,server=F,conservativeMapping=F,oldMapping=F,MGIsymbol=F){
  
  {
    #library(googleCloudStorageR)
    if(conservativeMapping|oldMapping){
      if(!file.exists("../data/gazestani_serverFiles/ortholog_mapping_old.rda")){
        stop("Error! file ../data/gazestani_serverFiles/ortholog_mapping_old.rda is missing")
      }
      load("../data/gazestani_serverFiles/ortholog_mapping.rda")
      mapping_old=mapping
    } else if(MGIsymbol){
      if(!file.exists("../data/gazestani_serverFiles/mouse_human_orthologs_MGISymbol_03062021.rda")){
        stop("Error! file ../data/gazestani_serverFiles/mouse_human_orthologs_MGISymbol_03062021.rda is missing")
      }
      load("../data/gazestani_serverFiles/mouse_human_orthologs_MGISymbol_03062021.rda")
      mapping_old=mapping
    }
    
    if(oldMapping|MGIsymbol){
      mapping=mapping_old
    } else {
      if(!file.exists("../data/gazestani_serverFiles/ortholog_mapping.rda")){
        stop("Error! file ../data/gazestani_serverFiles/ortholog_mapping.rda is missing")
      }
      load("../data/gazestani_serverFiles/ortholog_mapping.rda")
      
      if(conservativeMapping){
        mapping=mapping[which(paste0(mapping$human,"_",mapping$mouse) %in% paste0(mapping_old$human,"_",mapping_old$mouse)),]
      }
    }
    
  }
  
  
  #inputData=inputData[which(rowData(inputData)$ensembl_gene_id %in% mapping$mouse),]
  sl_gene_anno=.extraMouseGeneAnnoAdderFn(inputGeneNames = inputGeneNames,server = server)
  if(sum(sl_gene_anno$ensembl_gene_id %in% mapping$mouse,na.rm=T)==0){
    stop("ensembl_gene_ids were not found in the input data!")
  }
  
  sl_gene_anno$ensembl_gene_id=as.character(sl_gene_anno$ensembl_gene_id)
  sl_gene_anno$ensembl_gene_id[is.na(sl_gene_anno$ensembl_gene_id)]=toupper(sl_gene_anno)[is.na(sl_gene_anno$ensembl_gene_id)]
  if(sum(is.na(sl_gene_anno$ensembl_gene_id))>0){
    sl_gene_anno$ensembl_gene_id[is.na(sl_gene_anno$ensembl_gene_id)]=row.names(sl_gene_anno)[is.na(sl_gene_anno$ensembl_gene_id)]
  }
  
  sl_gene_anno=sl_gene_anno[!duplicated(sl_gene_anno$ensembl_gene_id),]
  
  mapping=mapping[which(mapping$mouse %in% sl_gene_anno$ensembl_gene_id),]
  mapping=mapping[!duplicated(mapping$mouse),]
  
  mapping=mapping[match(sl_gene_anno$ensembl_gene_id,mapping$mouse),]
  if(sum(is.na(mapping$mouse))>0){
    mapping$mouse[is.na(mapping$mouse)]=sl_gene_anno$ensembl_gene_id[is.na(mapping$mouse)]
  }
  if(sum(is.na(mapping$human))>0){
    mapping$human[is.na(mapping$human)]=sl_gene_anno$ensembl_gene_id[is.na(mapping$human)]
  }
  
  
  
  
  if(any(mapping$mouse!=sl_gene_anno$ensembl_gene_id)){
    stop("Error in the mapping!")
  }
  
  row.names(sl_gene_anno)=mapping$human
  sl_gene_anno$mouse_ensembl_gene_id=sl_gene_anno$ensembl_gene_id
  sl_gene_anno$ensembl_gene_id=mapping$human
  
  return(sl_gene_anno)
}


.myRiverPlotFn=function(inputData){
  require(ggalluvial)
  
  inputData=inputData$test_labels
  resArranged=NULL
  for(icol in colnames(inputData)[grepl("inferred",colnames(inputData))]){
    resArranged=rbind(resArranged,data.frame(label=inputData$madeCluster,inferred=icol,weight=inputData[,icol],stringsAsFactors = F))
  }
  inputData=resArranged
  inputData=aggregate(weight~label+inferred,data=inputData,sum)
  
  inputData=inputData[order(inputData$weight,decreasing = T),]
  inputData$label=factor(as.character(inputData$label),levels=unique(as.character(inputData$label)))
  inputData$inferred=factor(as.character(inputData$inferred),levels=unique(as.character(inputData$inferred)))
  
  
  pSource=ggplot(inputData,
                 aes(y = weight, axis1 = label, axis2 = inferred)) +
    geom_alluvium(aes(fill = label), width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", infer.label = TRUE) +
    scale_x_discrete(limits = c(axis1, axis2), expand = c(.05, .05))+theme_bw()+theme(panel.grid = element_blank(),legend.position = "none",axis.text.y = element_blank(),
                                                                                      axis.ticks = element_blank())+ylab("")
  return(pSource)
}


#cols = c("blue", "red"); col.min = -2.5; col.max = 2.5; dot.min = 0; dot.scale = 6; 
#idents = NULL; group.by = NULL; split.by = NULL; cluster.idents = FALSE; 
#scale = TRUE; scale.by = "radius"; scale.min = NA; scale.max = NA

.myDotPlot=function (object, assay = NULL,gene_name_col="gene_short_name", features, cols = c("blue", 
                                                                                              "red"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                     idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE, 
                     scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) 
{
  if(is.null(assay)){
    assay=DefaultAssay(object = object)
  }
  
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% rownames(x = brewer.pal.info))
  scale.func <- switch(EXPR = scale.by, size = scale_size, radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE, 
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  
  map_features=NULL
  if(!is.null(gene_name_col)){
    map_features=as.data.frame(object@assays$RNA@meta.features)
    map_features=map_features[which(toupper(map_features[,gene_name_col]) %in% toupper(features)),]
    converted_features=row.names(map_features)
  } else {
    converted_features=features
  }
  
  converted_features=gsub("_","-",converted_features)
  data.features <- FetchData(object = object, vars = converted_features, 
                             cells = cells)
  
  if(!is.null(map_features)){
    if(sum(duplicated(map_features[,gene_name_col]))>0){
      for(i in unique(map_features[duplicated(map_features[,gene_name_col]),gene_name_col])){
        map_features[which(map_features[,gene_name_col]==i),gene_name_col]=paste0(i,"_",1:sum(map_features[,gene_name_col]==i))
      }
    }
    
    for(i in 1:nrow(map_features)){
      colnames(data.features)[which(colnames(data.features) %in% row.names(map_features)[i])]=map_features[i,gene_name_col]
    }
  }
  
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  }
  else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = Seurat:::PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  if (length(x = levels(x = data.plot$id)) == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(X = as.character(x = data.plot$id), 
                         FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((", 
                                                                                          paste(sort(x = levels(x = object), decreasing = TRUE), 
                                                                                                collapse = "|"), ")_)"), replacement = "", 
                         USE.NAMES = FALSE)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                       levels = unique(x = feature.groups))
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
                                                        y = "id")) + geom_point(mapping = aes_string(size = "pct.exp", 
                                                                                                     color = color.by)) + scale.func(range = c(0, dot.scale), 
                                                                                                                                     limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
                                                                                                                                                                               axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
    labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
                                    yes = "Identity", no = "Split Identity")) + theme_cowplot()
  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(facets = ~feature.groups, scales = "free_x", 
                              space = "free_x", switch = "y") + theme(panel.spacing = unit(x = 1, 
                                                                                           units = "lines"), strip.background = element_blank())
  }
  if (split.colors) {
    plot <- plot + scale_color_identity()
  }
  else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  }
  else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  }
  return(plot)
}

.myMetaMarkerFn=function(inputClusteringRes,inputSource,inputTarget,minSimProb=0.8,n.cores=5,only.pos.markers=F,...){
  #minSimProb: cell type assignments with probability below this score are ignored
  #inputSource: the UMI counts from the source dataset
  #inputTarget: the UMI counts from the target
  
  colnames(inputClusteringRes$test_labels)=gsub("infered","inferred",colnames(inputClusteringRes$test_labels))
  lbls=colnames(inputClusteringRes$test_labels)
  lbls=lbls[!grepl("_max",lbls)]
  lbls=lbls[grepl("inferred_",lbls)]
  lblIndx=unlist(apply(inputClusteringRes$test_labels[,lbls],1,function(x) if(any(!is.na(x))){if(max(x,na.rm = T)>minSimProb){which(x==max(x))} else {NA}} else {NA}))
  lbls=gsub("inferred_","",lbls)[lblIndx]
  lbls[is.na(lbls)]="Unassigned"
  
  inputSource$cluster=inputClusteringRes$training_labels[,inputClusteringRes$source_label_col]
  inputTarget$cluster=lbls
  
  inputSource=inputSource[,which(!inputSource$cluster %in% c("","Unassigned"))]
  lbls=lbls[which(!inputTarget$cluster %in% c("","Unassigned"))]
  inputTarget=inputTarget[,which(!inputTarget$cluster %in% c("","Unassigned"))]
  
  inputSource=.extraExport2SeuratFn(inputSource)#[,sample(ncol(inputSource),3000)])
  inputTarget=.extraExport2SeuratFn(inputData = inputTarget)#[,sample(ncol(inputTarget),3000)])
  
  inputSource = NormalizeData(inputSource, normalization.method = "LogNormalize", scale.factor = 10000)
  inputTarget = NormalizeData(inputTarget, normalization.method = "LogNormalize", scale.factor = 10000)
  
  Idents(inputSource)=(as.character(inputSource$cluster))
  inputTarget$cluster=lbls
  Idents(inputTarget)=(as.character(inputTarget$cluster))
  
  source.markers = FindAllMarkers(inputSource, only.pos = only.pos.markers,verbose = F, ...)
  target.markers = FindAllMarkers(inputTarget, only.pos = only.pos.markers,verbose = F, ...)
  source.markers$cluster=as.character(source.markers$cluster)
  target.markers$cluster=as.character(target.markers$cluster)
  
  dfCombined=data.frame(gene=c(source.markers$gene,target.markers$gene),cluster=c(source.markers$cluster,target.markers$cluster),stringsAsFactors = F)
  dfCombined=dfCombined[!duplicated(paste0(dfCombined$gene,"_",dfCombined$cluster)),]
  
  myAUCfn=function(indx,inputExp,inputGeneInfo){
    tmp=.extra_auc_score(geneExp = as.numeric(inputExp@assays$RNA@data[row.names(inputExp)==inputGeneInfo$gene[indx],]),cluster = as.numeric(inputExp$cluster==as.character(inputGeneInfo$cluster[indx])))
    return(data.frame(gene=inputGeneInfo$gene[indx],cluster=inputGeneInfo$cluster[indx],AUC=tmp,stringsAsFactors = F))
  }
  
  resAUC_source=parallel::mclapply(1:nrow(dfCombined),myAUCfn,inputExp=inputSource,inputGeneInfo=dfCombined,mc.cores = n.cores)
  resAUC_source2=do.call("rbind",resAUC_source)
  source.markers=merge(source.markers,resAUC_source2,by=c("gene","cluster"))
  
  resAUC_target=parallel::mclapply(1:nrow(dfCombined),myAUCfn,inputExp=inputTarget,inputGeneInfo=dfCombined,mc.cores = n.cores)
  resAUC_target2=do.call("rbind",resAUC_target)
  target.markers=merge(target.markers,resAUC_source2,by=c("gene","cluster"))
  
  colnames(source.markers)[!colnames(source.markers) %in% c("gene","cluster")]=paste0("source_",colnames(source.markers)[!colnames(source.markers) %in% c("gene","cluster")])
  colnames(target.markers)[!colnames(target.markers) %in% c("gene","cluster")]=paste0("target_",colnames(target.markers)[!colnames(target.markers) %in% c("gene","cluster")])
  source.markers$source_zscore=qnorm(source.markers$source_p_val,lower.tail = F)
  target.markers$target_zscore=qnorm(target.markers$target_p_val,lower.tail = F)
  
  res_meta=merge(source.markers,target.markers,by=c("gene","cluster"),all=T)
  res_meta$stouffer_zscore=apply(res_meta[,c("target_zscore","source_zscore")],1,function(x) sum(pmin(x,10),na.rm = T)/sqrt(2))
  res_meta$stouffer_zscore[which(res_meta$source_avg_logFC&res_meta$target_avg_logFC<0)]=0
  
  res_meta=res_meta[order(res_meta$stouffer_zscore,decreasing = T),]
  
  
  
  resCorData=NULL
  for(i in unique(as.character(res_meta$cluster))){
    tmpData=res_meta[which(res_meta$cluster==i),]
    
    if(sum(!is.na(tmpData$target_avg_logFC))>5&sum(!is.na(tmpData$source_avg_logFC))>5){
      tmpCor=cor(tmpData$target_avg_logFC,tmpData$source_avg_logFC,use="complete.obs")
    } else {
      tmpCor=0
    }
    resCorData=rbind(resCorData,data.frame(cluster=i,correlation=tmpCor,source_num_genes=sum(!is.na(tmpData$source_avg_logFC)),target_num_genes=sum(!is.na(tmpData$target_avg_logFC)),intersect_num_genes=sum(!is.na(tmpData$source_avg_logFC)&!is.na(tmpData$target_avg_logFC)),stringsAsFactors = F))
  }
  return(list(stouffer_mata_analysis=res_meta,marker_correlations=resCorData))
}


.my2d_dotPlot=function(inputDataList,single_gene,inputPd,cluster_col,normalizeExp=T,gene_name_col="gene_short_name",dataset_name_col="ds_batch",ignore_missing_anno=F,return_plot=T,...){
  require(Seurat)
  require(ggplot2)
  require(cowplot)
  
  p=NULL
  p_data=NULL
  
  if(class(inputDataList[[1]])!="Seurat"){
    for(i in 1:length(inputDataList)){
      tmp=.extraExport2SeuratFn(inputDataList[[i]])
      tmp=Seurat::NormalizeData(tmp,verbose=F)
      inputDataList[[i]]=tmp
    }
  } else {
    if(normalizeExp){
      for(i in 1:length(inputDataList)){
        tmp=inputDataList[[i]]
        tmp=Seurat::NormalizeData(tmp,verbose=F)
        inputDataList[[i]]=tmp
      }
    }
  }
  
  for(i in 1:length(inputDataList)){
    if(sum(toupper(inputDataList[[i]]@assays$RNA@meta.features[,gene_name_col])==toupper(single_gene),na.rm = T)>0&sum(colnames(inputDataList[[i]]) %in% row.names(inputPd))>0){
      
      tmp_data=inputDataList[[i]]
      if(ignore_missing_anno){
        tmp_data=tmp_data[,which(colnames(tmp_data) %in% row.names(inputPd))]
      }
      if(ncol(tmp_data)>0&sum(toupper(inputDataList[[i]]@assays$RNA@meta.features[,gene_name_col])==toupper(single_gene),na.rm = T)>0){
        tmp_pd=inputPd[match(colnames(tmp_data),row.names(inputPd)),]
        
        if(sum(is.na(tmp_pd[,cluster_col]))>0){
          stop("Error in matching dataset with pd")
        }
        
        tmp_p=tmp_data
        Idents(tmp_p)=factor(as.character(tmp_pd[,cluster_col]))
        
        #object=tmp_p;gene_name_col=gene_name_col; features=single_gene
        tmp_p=.myDotPlot(object=tmp_p,gene_name_col=gene_name_col, features=single_gene) 
        
        if(is.null(dataset_name_col)){
          tmp_p$data$dataset=names(inputDataList)[i]
        } else {
          tmp_p$data$dataset=inputDataList[[i]]@meta.data[1,dataset_name_col]
        }
        
        if(is.null(p)){
          p=tmp_p
        }
        
        p_data=rbind(p_data,tmp_p$data)
      }
      
    }
    
  }
  
  
  if(!is.null(p_data)){
    if(nrow(p_data)>0){
      p_data$features.plot=single_gene
      p_data$avg.exp.scaled[which(is.nan(p_data$avg.exp.scaled))]=0
      tmp_avg.exp=aggregate(avg.exp~id+features.plot+dataset,data=p_data,function(x) mean(x,na.rm = T))
      tmp_pct.exp=aggregate(pct.exp~id+features.plot+dataset,data=p_data,function(x) mean(x,na.rm = T))
      tmp_avg.exp.scaled=aggregate(avg.exp.scaled~id+features.plot+dataset,data=p_data,function(x) mean(x,na.rm = T))
      
      tmp_data=merge(tmp_avg.exp,tmp_avg.exp.scaled,by=c("id","features.plot","dataset"),all=T)
      tmp_data=merge(tmp_data,tmp_pct.exp,by=c("id","features.plot","dataset"),all=T)
      
      
      tmp_feature.plot=tmp_data$features.plot
      tmp_data$features.plot=tmp_data$dataset
      tmp_data$gene=tmp_feature.plot
      p$data=tmp_data
      p=p+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5))+xlab("Dataset")
    } else {
      p=ggplot()
    }
  } else {
    tmp_data=NULL
    p=ggplot()
  }
  
  if(!return_plot){
    p=list(plot=p,p_data=tmp_data)
  }
  
  return(p)
  
}


.myProp_geneDensityPlot=function(inputDataList,single_gene,inputPd,lbl_inf_col,gene_name_col="gene_short_name",dataset_name_col="ds_batch",...){
  require(Seurat)
  require(ggplot2)
  require(cowplot)
  
  p=NULL
  p_data=NULL
  
  if(class(inputDataList[[1]])!="Seurat"){
    for(i in 1:length(inputDataList)){
      tmp=.extraExport2SeuratFn(inputDataList[[i]])
      tmp=Seurat::NormalizeData(tmp,verbose=F)
      inputDataList[[i]]=tmp
    }
  } else {
      for(i in 1:length(inputDataList)){
        tmp=inputDataList[[i]]
        tmp=Seurat::NormalizeData(tmp,verbose=F)
        inputDataList[[i]]=tmp
      }
    
  }
  
  
  for(i in 1:length(inputDataList)){
    if(sum(toupper(inputDataList[[i]]@assays$RNA@meta.features[,gene_name_col])==toupper(single_gene),na.rm = T)>0&sum(grepl(toupper("iNPH"),toupper(colnames(inputDataList[[i]]))))==0){
      
      tmp_pd=inputPd[match(colnames(inputDataList[[i]]),row.names(inputPd)),]
      
      if(sum(is.na(tmp_pd[,lbl_inf_col]))/nrow(tmp_pd)>0.2){
        stop("Error in matching dataset with pd")
      }
      
      tmp_p=inputDataList[[i]]
      if(!is.null(gene_name_col)){
        tmp_gene_ind=which(toupper(tmp_p@assays$RNA@meta.features[,gene_name_col])==toupper(single_gene))
        
      } else {
        tmp_gene_ind=which(toupper(row.names(tmp_p))==toupper(single_gene))
      }
      
      sl_cells=NULL
      if(length(tmp_gene_ind)>0){
        tmp_p=Seurat::ScaleData(tmp_p,features=row.names(tmp_p)[tmp_gene_ind],verbose=F)
        tmp_p2=tmp_p@assays$RNA@scale.data
        sl_cells=colnames(tmp_p2)[tmp_p2[1,]>0]
      }
      
      if(!is.null(sl_cells)){
        if(length(sl_cells)>0){
          p_data=rbind(p_data,data.frame(sample=sl_cells,lbl_inf_score=inputPd[sl_cells,lbl_inf_col],batch=inputPd[sl_cells,dataset_name_col],stringsAsFactors = F))
        }
      }
      
    }
    
  }
  
  
  p_plot_all=ggplot()
  p_plot_ind=ggplot()
  if(!is.null(p_data)){
    if(nrow(p_data)>0){
      p_data$lbl_inf_score[is.na(p_data$lbl_inf_score)]=0
      p_data$cdf=ecdf(p_data$lbl_inf_score)(p_data$lbl_inf_score)
      p_data2=p_data[!duplicated(p_data$cdf),]
      p_data2=p_data2[order(p_data2$lbl_inf_score,decreasing = T),]
      p_plot_all=ggplot(p_data2,aes(lbl_inf_score,cdf))+geom_path()+theme_classic()+scale_y_continuous(limits = c(0,1))+scale_x_reverse()+xlab("Propagation score")+ylab("CDF")+ggtitle("All datasets")
      
      
      p_data_list=split(p_data,p_data$batch)
      for(i in 1:length(p_data_list)){
        tmp_pd=p_data_list[[i]]
        tmp_pd$cdf=ecdf(tmp_pd$lbl_inf_score)(tmp_pd$lbl_inf_score)
        p_data2=tmp_pd[!duplicated(tmp_pd$cdf),]
        p_data2=p_data2[order(p_data2$lbl_inf_score,decreasing = T),]
        tmp_pd=ggplot(p_data2,aes(lbl_inf_score,cdf))+geom_path()+theme_classic()+scale_y_continuous(limits = c(0,1))+scale_x_reverse()+xlab("Propagation score")+ylab("CDF")+ggtitle(tmp_pd$batch[1])
        if(i==1){
          p_plot_ind=tmp_pd
        } else {
          p_plot_ind=p_plot_ind+tmp_pd
        }
      }
      
    }
  }
  
  
  return(list(p_all=p_plot_all,p_ind=p_plot_ind))
  
}

.my2d_dotPlot2=function(inputDataList,single_gene,cells.1,normalizeExp=T,gene_name_col="gene_short_name",dataset_name_col="ds_batch",...){
  require(Seurat)
  require(ggplot2)
  require(cowplot)
  
  p=NULL
  p_data=NULL
  
  if(class(inputDataList[[1]])!="Seurat"){
    for(i in 1:length(inputDataList)){
      tmp=.extraExport2SeuratFn(inputDataList[[i]])
      tmp=Seurat::NormalizeData(tmp,verbose=F)
      inputDataList[[i]]=tmp
    }
  } else {
    if(normalizeExp){
      for(i in 1:length(inputDataList)){
        tmp=inputDataList[[i]]
        tmp=Seurat::NormalizeData(tmp,verbose=F)
        inputDataList[[i]]=tmp
      }
    }
  }
  for(i in 1:length(inputDataList)){
    if(sum(toupper(inputDataList[[i]]@assays$RNA@meta.features[,gene_name_col])==toupper(single_gene),na.rm = T)>0){
      
      
      tmp_p=inputDataList[[i]]
      cluster_name=rep("other",ncol(tmp_p))
      cluster_name[colnames(tmp_p) %in% cells.1]="inferred"
      Idents(tmp_p)=factor(cluster_name)
      
      tmp_p=.myDotPlot(object=tmp_p,gene_name_col=gene_name_col, features=single_gene) 
      
      if(is.null(dataset_name_col)){
        tmp_p$data$dataset=names(inputDataList)[i]
      } else {
        tmp_p$data$dataset=inputDataList[[i]]@meta.data[1,dataset_name_col]
      }
      
      if(is.null(p)){
        p=tmp_p
      }
      
      p_data=rbind(p_data,tmp_p$data)
    }
    
  }
  if(nrow(p_data)>0){
    p_data$features.plot=single_gene
    p_data$avg.exp.scaled[which(is.nan(p_data$avg.exp.scaled))]=0
    tmp_avg.exp=aggregate(avg.exp~id+features.plot+dataset,data=p_data,function(x) mean(x,na.rm = T))
    tmp_pct.exp=aggregate(pct.exp~id+features.plot+dataset,data=p_data,function(x) mean(x,na.rm = T))
    tmp_avg.exp.scaled=aggregate(avg.exp.scaled~id+features.plot+dataset,data=p_data,function(x) mean(x,na.rm = T))
    
    tmp_data=merge(tmp_avg.exp,tmp_avg.exp.scaled,by=c("id","features.plot","dataset"),all=T)
    tmp_data=merge(tmp_data,tmp_pct.exp,by=c("id","features.plot","dataset"),all=T)
    
    
    tmp_feature.plot=tmp_data$features.plot
    tmp_data$features.plot=tmp_data$dataset
    tmp_data$gene=tmp_feature.plot
    p$data=tmp_data
    p=p+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5))+xlab("Dataset")
  } else {
    p=ggplot()
  }
  
  return(p)
  
}


.extraDoubletMakerFn=function(inputData,label_col,sel_labels,cells_per_group=30){
  print("Adding the doublet data...")
  selCells=list()
  for(i in unique(sel_labels)){
    tmp=inputData[,colData(inputData)[,label_col]==i]
    if(ncol(tmp)>0){
      tmp=tmp[,sample(ncol(tmp),min(ncol(tmp),cells_per_group))]
      selCells=c(selCells,list(new=tmp))
      names(selCells)[length(selCells)]=i
    }
  }
  
  doublet_exp=list()
  lbls=c()
  for(i in 1:(length(selCells)-1)){
    ds1=selCells[[i]]
    for(j in (i+1):length(selCells)){
      ds2=selCells[[j]]
      for(ik in 1:ncol(ds1)){
        for(jk in 1:ncol(ds2)){
          x=cbind(counts(ds1)[,ik],counts(ds2)[,jk])
          x=rowSums(x)
          x=data.frame(new=x)
          colnames(x)=paste0(names(selCells)[i],"_",names(selCells)[j])
          doublet_exp=c(doublet_exp,list(x))
          lbls=c(lbls,paste0(names(selCells)[i],"_",names(selCells)[j]))
        }
      }
    }
  }
  
  doublet_exp2=do.call("cbind",doublet_exp)
  
  colnames(doublet_exp2)=paste0(colnames(doublet_exp2),"_",LETTERS[sample(26,ncol(doublet_exp2),replace = T)])
  while(sum(duplicated(colnames(doublet_exp2)))>0){
    colnames(doublet_exp2)[duplicated(colnames(doublet_exp2))]=paste0(colnames(doublet_exp2)[duplicated(colnames(doublet_exp2))],LETTERS[sample(26,sum(duplicated(colnames(doublet_exp2))),replace = T)])
  }
  
  pd=as.data.frame(colData(inputData))
  for(i in 1:ncol(pd)){
    if(class(pd[,i])=="character"){
      pd[,i]="synthetic"
    } else{
      pd[,i]=NA
    }
  }
  
  pd=pd[1:ncol(doublet_exp2),]
  
  pd[,label_col]=lbls
  pd$is_doublet="YES"
  inputData$is_doublet="NO"
  
  doublet_exp=Matrix(as.matrix(doublet_exp2),sparse=T)
  
  fd=as.data.frame(rowData(inputData))
  if(!all(row.names(inputData)==row.names(doublet_exp))){
    stop("Error in making the doublet data!")
  }
  
  expInput=counts(inputData)
  expM=cbind(expInput,doublet_exp)
  pdM=rbind(as.data.frame(colData(inputData)),pd)
  
  res=SingleCellExperiment(assays = list(counts = expM),colData = pdM,rowData=fd)
  
  print("Doublet data is added.")
  return(res)
}

.extraKnnSimNetFn=function(inputPCAembeddings,n.adaptiveKernel=5,nPropIter=4,creatNet=TRUE,n.adaptiveKernel.coefficient=3){
  #nPropIter: number of propagation iterations
  #input: PCA res
  #CreateNet: create nn2 network? if FALSE, it is assumed that inputPCAembeddings is a nn2 network
  n.cells <- nrow(inputPCAembeddings)
  if(creatNet){
    nn.ranked <- RANN::nn2(data = inputPCAembeddings, k = n.adaptiveKernel*6, eps = 0)
    
    
    dists=nn.ranked$nn.dists
    x=apply(dists,1,function(x) sum(x==0))
    if(sum(x>10)>0){
      inputPCAembeddings=inputPCAembeddings[-which(x>10),]
      nn.ranked <- RANN::nn2(data = inputPCAembeddings, k = n.adaptiveKernel*6, eps = 0)
      print("Duplicate points in the dataset were found!")
    }
  } else {
    nn.ranked=inputPCAembeddings
  }
  
  
  
  dists=nn.ranked$nn.dists
  affinities=t(apply(dists,1,function(x) exp((-1)*(x/x[n.adaptiveKernel])^2)))
  affinitiesThr=apply(affinities,1,function(x) x[n.adaptiveKernel.coefficient*n.adaptiveKernel+1])
  affinities2=affinities-affinitiesThr
  affinities[affinities2<=0]=0
  
  affCounts=apply(affinities,2,function(x) sum(x==0))
  affinities=affinities[,-which(affCounts==nrow(affinities))]
  nn.ranked$nn.idx=nn.ranked$nn.idx[,-which(affCounts==nrow(affinities))]
  j <- as.numeric(t(nn.ranked$nn.idx))
  i <- ((1:length(j)) - 1)%/%ncol(affinities) + 1
  k=as.numeric(t(affinities))
  
  
  graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(dists), nrow(dists)))
  graph=(graph+t(graph))/2
  graph <- Matrix::Diagonal(x = 1 / (rowSums(graph))) %*% graph
  
  if(class(graph)!="dgCMatrix"){
    graph=as(graph,"dgCMatrix")
  }
  
  rownames(graph) <- rownames(inputPCAembeddings)
  colnames(graph) <- rownames(inputPCAembeddings)
  
  for(i in 1:nPropIter){
    print(i)
    if(i==1){
      prop_mat=graph
    } else {
      prop_mat=prop_mat%*%graph
    }
  }
  
  
  return(list(prop_mat=prop_mat,PCA_embeddings=inputPCAembeddings))
}

FindNeighbors_archive=function (object, reduction="pca", dims=20, k.param = 20, compute.SNN = TRUE,compute.MAN=TRUE, 
                        jaccard.indx.thr = 1/15, nn.eps = 0, verbose = TRUE) {
  
  object=.myFindNeighbors(object=object, reduction=reduction, dims=dims, k.param = k.param, compute.SNN = compute.SNN,compute.MAN=compute.MAN, 
                                    jaccard.indx.thr = jaccard.indx.thr, nn.eps = nn.eps, verbose = verbose)
  
  return(object)
}

.myFindNeighbors=function (object, reduction="pca", dims=20, k.param = 20, compute.SNN = TRUE,compute.MAN=TRUE, 
                        jaccard.indx.thr = 1/15, nn.eps = 0, verbose = TRUE) {
  
  if(class(object)!="Seurat"){
    stop("Input data is supposed to be a seurat object")
  }
  
  inputPCAscores=Seurat:::Embeddings(object = object[[reduction]])
  
  n.cells <- nrow(x = inputPCAscores)
  
  if (verbose) {
    message("Computing nearest neighbor graph")
  }
  nn.ranked_org <- RANN::nn2(data = inputPCAscores, k = k.param*3, eps = nn.eps)
  nn.ranked <- nn.ranked_org$nn.idx[,1:k.param]
  
  j <- as.numeric(t(nn.ranked))
  i <- ((1:length(j)) - 1)%/%k.param + 1
  nn.matrix <- sparseMatrix(i = i, j = j, x = 1,dims = c(nrow(inputPCAscores), nrow(inputPCAscores)))
  rownames(nn.matrix) <- rownames(inputPCAscores)
  colnames(nn.matrix) <- rownames(inputPCAscores)
  neighbor.graphs <- list(nn = nn.matrix)
  if (compute.SNN) {
    if (verbose) {
      message("Computing SNN")
    }
    snn.matrix <- Seurat:::ComputeSNN(nn_ranked = nn.ranked, prune = jaccard.indx.thr)
    rownames(x = snn.matrix) <- rownames(inputPCAscores)
    colnames(x = snn.matrix) <- rownames(inputPCAscores)
    snn.matrix <- as.Graph(x = snn.matrix)
    
    neighbor.graphs[["snn"]] <- snn.matrix
  }
  
  if(compute.MAN){
    netMAN=.extraKnnSimNetFn(inputPCAembeddings=nn.ranked_org,n.adaptiveKernel=5,nPropIter=3,creatNet=FALSE,n.adaptiveKernel.coefficient=4)
    rownames(netMAN$prop_mat) <- rownames(inputPCAscores)
    colnames(netMAN$prop_mat) <- rownames(inputPCAscores)
    neighbor.graphs[["man"]] <- netMAN$prop_mat
  }
  
  assay=DefaultAssay(object = object[[reduction]])
  graph.name <- paste0(assay, "_", names(x = neighbor.graphs))
  for (ii in 1:length(x = graph.name)) {
    tmp=neighbor.graphs[[ii]]
    tmp=Seurat::as.Graph(tmp)
    DefaultAssay(object = tmp) <- assay
    object[[graph.name[[ii]]]] <- tmp
    rm(tmp)
  }
  
  
  return(object)
}

.myMarkerBasedAnalysisFn=function(inputData,inputMarkers,all_positive=F,nPCs=30,nPropIter=4,n.adaptiveKernel=5){
  
  #inputMarkers: a data.frame with the column gene representing the geneNames that are consistent with the row.names of the inputData and the second column (if present), is the FCs
  
  inputData=inputData[!duplicated(rowData(inputData)$ensembl_gene_id),]
  inputData=inputData[!is.na(rowData(inputData)$ensembl_gene_id),]
  row.names(inputData)=rowData(inputData)$ensembl_gene_id
  
  if(class(inputData)!="Seurat"){
    inputData=.extraExport2SeuratFn(inputData = inputData)
  }
  inputData <- NormalizeData(inputData, normalization.method = "LogNormalize", scale.factor = 10000,verbose = F)
  inputData <- FindVariableFeatures(inputData, selection.method = "vst", nfeatures = 2000,verbose = F)
  inputData <- ScaleData(inputData,verbose = F)
  inputData <- RunPCA(inputData)
  similarityMat_org=.extraKnnSimNetFn(inputPCAembeddings = inputData@reductions$pca@cell.embeddings[,1:nPCs],nPropIter = nPropIter,n.adaptiveKernel = n.adaptiveKernel)
  similarityMat=similarityMat_org$prop_mat
  inputData=inputData[,colnames(inputData) %in% row.names(similarityMat_org$PCA_embeddings)]
  
  if(all_positive){
    
    res=.netAUCellFn(inputExp = inputData@assays$RNA@data,inputGmt = list(geneset=inputMarkers[,1]),numAUCgenes=min(round(0.2*nrow(inputData)),median(inputData$QC_Gene_unique_count)))
    res=similarityMat %*% matrix(res[1,],ncol=1)
    res=data.frame(sample=colnames(inputData),score=res[,1],stringsAsFactors = F)
    
  } else {
    
    inputExp=inputData@assays$RNA@data[row.names(inputData) %in% inputMarkers[,1],]
    inputExp=as.matrix(inputExp)
    
    gMean=rowMeans(inputExp)
    FCs=sweep(inputExp,1,gMean,"-")
    
    inputMarkers=inputMarkers[match(row.names(FCs),inputMarkers[,1]),]
    
    FCs_cor=apply(FCs,2,function(x) cor(x,inputMarkers[,2],use="complete.obs"))
    
    FCs_cor2=as.matrix(similarityMat %*% matrix(FCs_cor,ncol=1))
    res=data.frame(sample=row.names(FCs_cor2),score=FCs_cor2[,1],stringsAsFactors = F)
  }
  
  return(res)
}


.mycellAssignHeatmap_binary=function(input_labelTransfer_object,confidenceLevel=0.8,target_cluster_col="Cluster"){
  require(ggplot2)
  lbls=input_labelTransfer_object$test_labels[,grepl("inferred_",colnames(input_labelTransfer_object$test_labels))&!grepl("_max",colnames(input_labelTransfer_object$test_labels))]
  lblIndx=apply(lbls,1,function(x) if(sum(!is.na(x))>0){if(max(x,na.rm = T)>confidenceLevel){which(x==max(x))[1]}else{NA}}else{NA})
  lbls=gsub("inferred_","",colnames(lbls))[lblIndx]
  
  lbls=data.frame(cluster=input_labelTransfer_object$test_labels[,target_cluster_col],inferred=lbls,stringsAsFactors = F)
  lbls=as.data.frame(table(lbls$cluster,lbls$inferred,useNA='ifany'))
  colnames(lbls)=c("cluster","inferred","count")
  lbls$proportion=0
  for(i in 1:nrow(lbls)){
    lbls$proportion[i]=lbls$count[i]/sum(lbls$count[lbls$cluster==lbls$cluster[i]],na.rm = T)
  }
  p=ggplot(lbls,aes(inferred,cluster,fill=proportion))+geom_tile(color="black")+scale_fill_gradient(low="white",high="red")+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=1))
  return(p)
}

.mycellAssignHeatmap_prob=function(input_labelTransfer_object,confidenceLevel=0.8,target_cluster_col="Cluster"){
  require(ggplot2)
  lbls=input_labelTransfer_object$test_labels[,grepl("inferred_",colnames(input_labelTransfer_object$test_labels))&!grepl("_max",colnames(input_labelTransfer_object$test_labels))]
  lbls[lbls<confidenceLevel]=0
  lbls$cluster_name=input_labelTransfer_object$test_labels[,target_cluster_col]
  lbls=reshape2::melt(lbls,"cluster_name")
  lbls=aggregate(value~variable+cluster_name,data=lbls,sum)
  
  lbls_total=aggregate(value~cluster_name,data=lbls,sum)
  colnames(lbls_total)=c("cluster_name","total_count")
  lbls=merge(lbls,lbls_total,by="cluster_name")
  lbls$fraction=lbls$value/lbls$total_count
  colnames(lbls)[colnames(lbls)=="variable"]="Inferred"
  lbls$Inferred=gsub("inferred_","",lbls$Inferred)
  
  p=ggplot(lbls,aes(Inferred,cluster_name,fill=fraction))+geom_tile(color="black")+scale_fill_gradient(low="white",high="red")+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=1))
  return(p)
}

.extra_auc_score = function(geneExp, cluster) {
  if(length(table(cluster)) != 2) {
    return(0.5)
  }
  
  pred = ROCR::prediction(predictions=geneExp, cluster)
  
  
  auc.perf = ROCR::performance(pred, measure = 'auc')
  res = auc.perf@y.values
  
  return(res[[1]])
}

.extraClusteringOptimization_detailed=function(inputIndx,inputDf,inputData,reductionList,reduction_dimList,cluster_algorithmList,cluster_resolutionList,cov){
  resAnalysis=NULL
  isize=inputDf$varGeneSizes[inputIndx]
  icol=inputDf$rankCol[inputIndx]
  tmp_var=row.names(hvg_anno)[hvg_anno[,icol]<=isize]
  tmpData2=inputData
  tmpData2@assays$RNA@var.features=tmp_var
  for(ireduction in reductionList){
        
        if(ireduction=="ica"){
          tmpData2 <- RunICA(tmpData2, features = tmp_var)
        } else {
          tmpData2 <- RunPCA(tmpData2, features = tmp_var)
        }
        for(ireductionDim in reduction_dimList){
          tmpData2 <- .myFindNeighbors(tmpData2, dims = 1:ireductionDim,reduction = ireduction)
          tmpGraphName=names(tmpData2@graphs)
          
          for(igraph_name in tmpGraphName){
            for(icluster_algorithm in cluster_algorithmList){
              for(icluster_resolution in cluster_resolutionList){
                tmpData2 <- FindClusters(tmpData2, resolution = icluster_resolution,algorithm = icluster_algorithm,graph.name = igraph_name,group.singletons=T)
                tmp_clusters=data.frame(sample=names(Seurat::Idents(tmpData2)),cluster=Seurat::Idents(tmpData2),stringsAsFactors = F)
                
                tmp_clusters=merge(tmp_clusters,data.frame(sample=colnames(inputData),ref=as.character(inputData@meta.data[,clusterCol])),by="sample")
                tmp_clusters$cluster=as.character(tmp_clusters$cluster)
                ari_res=CrossClustering::ari(table(tmp_clusters$ref,tmp_clusters$cluster))
                attr(ari_res$ari,"ci")[2]
                tmpRes=data.frame(cov=paste0(cov,collapse = ","),var_gene_method=icol,var_gene_count=isize,reduction=ireduction,reduction_dim=ireductionDim,graph=igraph_name,cluster_algorithm=icluster_algorithm,cluster_resolution=icluster_resolution,total_cluster_count=length(unique(tmp_clusters$cluster)),ari=ari_res$ari[1],stringsAsFactors = F)
                for(icluster in unique(tmp_clusters$ref)){
                  tmpRes$cluster_count=length(unique(tmp_clusters$cluster[tmp_clusters$ref==icluster]))
                  tmpPurity=c()
                  for(jcluster in unique(tmp_clusters$cluster[tmp_clusters$ref==icluster])){
                    tmpJ=tmp_clusters[tmp_clusters$cluster==jcluster,]
                    tmpPurity=c(tmpPurity,round(sum(tmpJ$ref==icluster,na.rm = T)/sum(!is.na(tmpJ$ref)),3))
                  }
                  tmpRes$cluster_purity=paste(tmpPurity,collapse = ",")
                  colnames(tmpRes)[colnames(tmpRes)=="cluster_purity"]=paste0(icluster,"_cluster_purity")
                  colnames(tmpRes)[colnames(tmpRes)=="cluster_count"]=paste0(icluster,"_cluster_count")
                }
                resAnalysis=rbind(resAnalysis,tmpRes)
              }
            }
          }
          
        }
        
      }
    
  
  return(resAnalysis)
}

.extraClusteringOptimization_detailed_type2=function(inputIndx,inputDf,inputData,batch_covariate){
  resAnalysis=NULL
  
  network_nPCs=inputDf[inputIndx,"reduction_dim"]
  icluster_algorithm=inputDf[inputIndx,"cluster_algorithm"]
  icluster_resolution=inputDf[inputIndx,'cluster_resolution']
  
  tmpData2 <- Seurat::FindNeighbors(inputData[,1:network_nPCs],verbose = F)
  tmpData2 <- FindClusters(tmpData2$snn, resolution = icluster_resolution,algorithm = icluster_algorithm,group.singletons=T,verbose = F)
  res_clusters=data.frame(sample=row.names(tmpData2),cluster=tmpData2[,1],stringsAsFactors = F)
  
  if(is.null(batch_covariate)){
    tmp_dist=dist(inputData[,1:network_nPCs])
    tmp_sillhoutte=cluster::silhouette(as.numeric(as.character(res_clusters$cluster)),tmp_dist)
    tmp_sillhoutte=summary(tmp_sillhoutte)
    tmp_clusters=data.frame(cluster_name=names(tmp_sillhoutte$clus.avg.widths),cluster_sillhoute=tmp_sillhoutte$clus.avg.widths,stringsAsFactors=F)
    
    resAnalysis=data.frame(reduction_dim=network_nPCs,cluster_algorithm=icluster_algorithm,cluster_resolution=icluster_resolution,total_cluster_count=length(unique(tmp_clusters$cluster_name)),min_sillhoutte=min(tmp_clusters$cluster_sillhoute),mean_sillhoutte=mean(tmp_clusters$cluster_sillhoute),stringsAsFactors = F)
    
  } else {
    resAnalysis=list()
    for(i in unique(batch_covariate)){
      tmpData=inputData[batch_covariate==i,1:network_nPCs]
      tmp_dist=dist(tmpData)
      tmp_sillhoutte=cluster::silhouette(as.numeric(as.character(res_clusters$cluster)[batch_covariate==i]),tmp_dist)
      tmp_sillhoutte=summary(tmp_sillhoutte)
      tmp_clusters=data.frame(cluster_name=names(tmp_sillhoutte$clus.avg.widths),cluster_sillhoute=tmp_sillhoutte$clus.avg.widths,stringsAsFactors=F)
      tmp_clusters=data.frame(reduction_dim=network_nPCs,cluster_algorithm=icluster_algorithm,cluster_resolution=icluster_resolution,total_cluster_count=length(unique(tmp_clusters$cluster_name)),min_sillhoutte=min(tmp_clusters$cluster_sillhoute),mean_sillhoutte=mean(tmp_clusters$cluster_sillhoute),batch_covariate=i,stringsAsFactors = F)
      resAnalysis=c(resAnalysis,list(tmp_clusters))
    }
    resAnalysis=do.call("rbind",resAnalysis)
  }
  return(resAnalysis)
}

.myGuidedClusteringOptimizerFn=function(inputData,clusterCol,organism,vars_to_regress_out=NULL,server=F,ncores=10,varGeneSizes=c(1000,1500,2000,3000,4000),reduction=c("pca","ica"),reduction_dim=seq(10,50,5),cluster_algorithm=seq(1,3),cluster_resolution=seq(0.1,1.5,0.1)){
  #Given prior information on the cell identities, it tries to optimize the clustering
  hvg_data=.myHVGadder(inputData = inputData,organism = organism,server=server)
  
  if(class(inputData)=="SingleCellExperiment"){
    inputData=.extraExport2SeuratFn(inputData = inputData)
  }
  hvg_anno=hvg_data@assays$RNA@meta.features
  rm(hvg_data)
  
  if(!is.null(vars_to_regress_out)){
    vars_to_regress_out=expand.grid(vars_to_regress_out,vars_to_regress_out)
    vars_to_regress_out$Var1=as.character(vars_to_regress_out$Var1)
    vars_to_regress_out$Var2=as.character(vars_to_regress_out$Var2)
    Var1=vars_to_regress_out$Var1
    Var1[vars_to_regress_out$Var2>vars_to_regress_out$Var1]=vars_to_regress_out$Var2[vars_to_regress_out$Var2>vars_to_regress_out$Var1]
    Var2=vars_to_regress_out$Var2
    Var2[vars_to_regress_out$Var1<vars_to_regress_out$Var2]=vars_to_regress_out$Var1[vars_to_regress_out$Var1<vars_to_regress_out$Var2]
    vars_to_regress_out=vars_to_regress_out[!duplicated(paste0(Var1,"_",Var2)),]
    vars_to_regress_out=rbind(vars_to_regress_out,data.frame(Var1="",Var2="",stringsAsFactors = F))
    rm(Var1,Var2)
  }
  
  allVarGenes=c()
  for(i in colnames(hvg_anno)[grepl("_rank",colnames(hvg_anno))]){
    allVarGenes=c(allVarGenes,row.names(hvg_anno)[hvg_anno[,i]<=max(varGeneSizes)])
    
  }
  allVarGenes=unique(allVarGenes)
  
  inputData=Seurat::NormalizeData(inputData)
  inputData=Seurat::FindVariableFeatures(inputData)
  inputData@assays$RNA@var.features=allVarGenes
  
  resAnalysis=list()
  
  for(icov in 1:nrow(vars_to_regress_out)){
    tmp_var_to_reg_out=unique(as.character(vars_to_regress_out[icov,]))
    
    tmpData=inputData
    if(all(tmp_var_to_reg_out=="")){
      tmp_var_to_reg_out=NULL
    } else {
      if(sum(tmp_var_to_reg_out=="")>0){
        tmp_var_to_reg_out=tmp_var_to_reg_out[tmp_var_to_reg_out!=""]
      }
    }
    tmpData <- ScaleData(tmpData, features = allVarGenes,vars.to.regress = tmp_var_to_reg_out)
    
    df_hvg_data=expand.grid(varGeneSizes,colnames(hvg_anno)[grepl("rank_",colnames(hvg_anno))])
    df_hvg_data$Var1=as.numeric(df_hvg_data$Var1)
    df_hvg_data$Var2=as.character(df_hvg_data$Var2)
    colnames(df_hvg_data)=c("varGeneSizes","rankCol")
    
    tmpres=parallel::mclapply(1:nrow(df_hvg_data),.extraClusteringOptimization_detailed,inputDf=df_hvg_data,inputData=tmpData,reductionList=reduction,reduction_dimList=reduction_dim,cluster_algorithmList=cluster_algorithm,cluster_resolutionList=cluster_resolution,cov=tmp_var_to_reg_out,mc.cores = ncores)
    resAnalysis=c(resAnalysis,list(tmpres))
  }
  
  return(resAnalysis)
}

.myClusteringOptimizerFn=function(inputPCA,reduction_dim=seq(10,50,5),cluster_algorithm=seq(1,3),cluster_resolution=seq(0.1,1.5,0.1),batch_covariate=NULL,ncores=30){
  #Examines different clustering parameters, and their impacts on the clustering results
  
  df=expand.grid(cluster_algorithm,cluster_resolution,reduction_dim)
  colnames(df)=c('cluster_algorithm','cluster_resolution','reduction_dim')
  
  tmpres=parallel::mclapply(1:nrow(df),.extraClusteringOptimization_detailed_type2,inputDf=df,inputData=inputPCA,batch_covariate=batch_covariate,mc.cores = ncores)
  resAnalysis=do.call(get("rbind.fill",asNamespace('plyr')),tmpres)
  return(resAnalysis)
}

.myVlnPlot=function(inputSeurat,inputGenes,inputGroups=NULL,geneNameCol="gene_short_name",...){
  
  fd=as.data.frame(inputSeurat@assays$RNA@meta.features)
  pd=as.data.frame(inputSeurat@meta.data)
  counts=inputSeurat@assays$RNA@counts
  
  if(sum(toupper(fd[,geneNameCol]) %in% toupper(inputGenes))/length(inputGenes)<0.5){
    print("Less than 50% were found!")
  }
  
  sampleGenes=unique(toupper(c(inputGenes,sample(fd[,geneNameCol],20))))
  sampleGenes=sampleGenes[!is.na(sampleGenes)]
  sampleGenes=sampleGenes[sampleGenes!=""]
  
  
  expData=SingleCellExperiment(assays = list(counts = counts),colData = pd,rowData=fd)
  rowData(expData)[,geneNameCol]=toupper(rowData(expData)[,geneNameCol])
  expData=expData[which(toupper(rowData(expData)[,geneNameCol]) %in% sampleGenes),]
  expData=expData[!duplicated(toupper(rowData(expData)[,geneNameCol])),]
  expData2=inputSeurat@assays$RNA@data
  expData2=expData2[row.names(expData2) %in% row.names(expData),]
  if(!all(row.names(expData2)==row.names(expData))){
    stop("Error in matching!")
  }
  row.names(expData2)=toupper(rowData(expData)[,geneNameCol])
  row.names(expData)=toupper(rowData(expData)[,geneNameCol])
  row.names(rowData(expData))=toupper(rowData(expData)[,geneNameCol])
  
  expData=.extraExport2SeuratFn(expData)
  expData@assays$RNA@data=expData2
  
  clusterings=data.frame(sample=colnames(inputSeurat),cluster=as.character(Idents(inputSeurat)),stringsAsFactors = F)
  clusterings=clusterings[match(colnames(expData),clusterings$sample),]
  
  if(!all(colnames(expData)==clusterings$sample,na.rm = F)){
    stop("Error!")
  }
  Idents(expData)=factor(clusterings$cluster,levels = levels(Idents(inputSeurat)))
  inputSeurat=expData
  
  
  if(!is.null(inputGroups)){
    Idents(inputSeurat)=inputGroups
  }
  p=VlnPlot(inputSeurat, features = toupper(inputGenes),...)
  return(p)
}

#inputGroupsCol=NULL;geneNameCol="gene_short_name";reductionName="umap";slIndx=NULL;ncol=6;internalNcol=6;combine_figs=T
#inputSeurat = tmp;inputDimData = pd2;inputGenes = Tushar_markers;combine_figs = F
.myFeaturePlot=function(inputSeurat,inputDimData,inputGenes,inputGroupsCol=NULL,geneNameCol="gene_short_name",reductionName="umap",slIndx=NULL,ncol=6,internalNcol=6,combine_figs=T,order=T,...){
  
  #inputGroups: the cell attributes to be colored with
  #geneNameCol: gene name column that will be used to color cells based on the expression level
  inputGenes=toupper(as.character(inputGenes))
  if(!is.null(geneNameCol)){
    fd=as.data.frame(inputSeurat@assays$RNA@meta.features)
    if(sum(is.na(fd[,geneNameCol]))>0){
      fd[which(is.na(fd[,geneNameCol])),geneNameCol]=row.names(fd)[which(is.na(fd[,geneNameCol]))]
    }
    if(sum(duplicated(fd[,geneNameCol]))>0){
      fd[which(duplicated(fd[,geneNameCol])),geneNameCol]=row.names(fd)[which(duplicated(fd[,geneNameCol]))]
    }
    
    pd=as.data.frame(inputSeurat@meta.data)
    counts=inputSeurat@assays$RNA@counts
    
    row.names(counts)=toupper(row.names(counts))
    row.names(fd)=toupper(row.names(fd))
    fd[,geneNameCol]=toupper(fd[,geneNameCol])
    inputGenes=toupper(inputGenes)
    id_order=row.names(fd)[match(inputGenes,toupper(fd[,geneNameCol]))]
    id_order=match(c(intersect(inputGenes,row.names(counts)),setdiff(row.names(counts),inputGenes)),row.names(counts))
    counts=counts[id_order,]
    fd=fd[id_order,]
    
    
    if(sum(toupper(fd[,geneNameCol]) %in% toupper(unlist(strsplit(inputGenes,","))))/length(inputGenes)<0.5){
      print("Less than 50% were found!")
    }
    
    
    inputGenesOrg=inputGenes
    if(sum(grepl(",",inputGenes))>0){
      inputGenes=unlist(strsplit(inputGenes,","))
    }
    inputGenes=unique(inputGenes)
    
    sampleGenes=unique(toupper(c(inputGenes,sample(fd[,geneNameCol],20))))
    sampleGenes=sampleGenes[!is.na(sampleGenes)]
    sampleGenes=sampleGenes[sampleGenes!=""]
    inputGenes=inputGenes[inputGenes!=""]
    
    sampleGenes=toupper(sampleGenes)
    inputGenes=toupper(inputGenes)
    inputGenesOrg=toupper(inputGenesOrg)
    fd[,geneNameCol]=toupper(fd[,geneNameCol])
    
    slGeneIndx=which(fd[,geneNameCol] %in% sampleGenes)
    counts=counts[slGeneIndx,]
    fd=fd[slGeneIndx,]
    
    slGeneIndx=which(!duplicated(fd[,geneNameCol]))
    counts=counts[slGeneIndx,]
    fd=fd[slGeneIndx,]
    
    expData2=inputSeurat@assays$RNA@data
    expData2=expData2[match(toupper(row.names(counts)),toupper(row.names(expData2))),]
    if(!all(toupper(row.names(expData2))==toupper(row.names(counts)))){
      stop("Error in matching!")
    }
    
    if(!is.null(slIndx)){
      counts=counts[,slIndx]
      expData2=expData2[,slIndx]
      inputDimData=inputDimData[slIndx,]
      pd=pd[slIndx,]
    }
    
    row.names(expData2)=fd[,geneNameCol]
    row.names(counts)=fd[,geneNameCol]
    row.names(fd)=fd[,geneNameCol]
    
    if(sum(grepl(",",inputGenesOrg))>0){
      tmpIndx=which(grepl(",",inputGenesOrg))
      for(igene in tmpIndx){
        if(length(which(row.names(counts) %in% unlist(strsplit(inputGenesOrg[igene],","))))>0){
          tmpExp1=matrix(counts[which(row.names(counts) %in% unlist(strsplit(inputGenesOrg[igene],","))&!duplicated(row.names(counts))),],ncol=ncol(counts))
          
          tmpGene1=paste(row.names(counts)[which(row.names(counts) %in% unlist(strsplit(inputGenesOrg[igene],","))&!duplicated(row.names(counts)))],collapse = ",")
          inputGenesOrg[igene]=tmpGene1
          fd1=as.data.frame(t(apply(fd[which(row.names(fd) %in% unlist(strsplit(inputGenesOrg[igene],","))),],2,function(x) paste(x,collapse = ","))))
          row.names(fd1)=paste(row.names(counts)[which(row.names(counts) %in% unlist(strsplit(inputGenesOrg[igene],","))&!duplicated(row.names(counts)))],collapse = ",")
          
          tmpExp1=matrix(apply(tmpExp1,2,sum),nrow=1)
          tmpExp1=as(tmpExp1, "dgCMatrix")
          row.names(tmpExp1)=tmpGene1
          counts=rbind(counts,tmpExp1)
          fd=rbind(fd,fd1)
          
          tmpExp2=matrix(expData2[which(row.names(expData2) %in% unlist(strsplit(inputGenesOrg[igene],","))&!duplicated(row.names(expData2))),],ncol=ncol(counts))
          tmpGene2=paste(row.names(expData2)[which(row.names(expData2) %in% unlist(strsplit(inputGenesOrg[igene],","))&!duplicated(row.names(expData2)))],collapse = ",")
          tmpExp2=matrix(apply(tmpExp2,2,sum),nrow=1)
          tmpExp2=as(tmpExp2, "dgCMatrix")
          row.names(tmpExp2)=tmpGene2
          expData2=rbind(expData2,tmpExp2)
        }
      }
    }
    
    fd[,geneNameCol]=row.names(fd)
    row.names(expData2)=gsub(",",".",row.names(expData2))
    row.names(counts)=gsub(",",".",row.names(counts))
    row.names(fd)=gsub(",",".",row.names(fd))
    inputGenes=gsub(",",".",inputGenesOrg)
    
    row.names(expData2)=gsub(",",".",row.names(expData2))
    row.names(counts)=gsub(",",".",row.names(counts))
    row.names(fd)=gsub(",",".",row.names(fd))
    inputGenes=gsub(",",".",inputGenesOrg)
    if(sum(grepl("^[0-9]",inputGenes))>0){
      row.names(expData2)[grepl("^[0-9]",row.names(expData2))]=paste0("X",row.names(expData2)[grepl("^[0-9]",row.names(expData2))])
      row.names(counts)[grepl("^[0-9]",row.names(counts))]=paste0("X",row.names(counts)[grepl("^[0-9]",row.names(counts))])
      row.names(fd)[grepl("^[0-9]",row.names(fd))]=paste0("X",row.names(fd)[grepl("^[0-9]",row.names(fd))])
      inputGenes[grepl("^[0-9]",inputGenes)]=paste0("X",inputGenes[grepl("^[0-9]",inputGenes)])
    }
    
    expData=SingleCellExperiment(assays = list(counts = counts),colData = pd,rowData=fd)
    
    
    expData=.extraExport2SeuratFn(expData)
    expData@assays$RNA@data=expData2
    inputSeurat=expData
    if(!is.null(inputGroupsCol)){
      Idents(inputSeurat)=inputSeurat@meta.data[,inputGroupsCol]
    }
  } else {
    if(!is.null(inputGroupsCol)){
      Idents(inputSeurat)=inputSeurat@meta.data[,inputGroupsCol]
    }
  }
  
  p=""
  if(length(which(inputGenes %in% row.names(inputSeurat)))>0){
    if(!is.null(inputGroupsCol)){
      data.list <- SplitObject(inputSeurat, split.by = inputGroupsCol)
      p = lapply(1:length(data.list),.extraFeaturePlot2,expDataList=data.list,umapData=inputDimData, features=toupper(inputGenes[inputGenes!=""]),reductionName=reductionName,ncol=internalNcol,order=order,combine_figs=combine_figs,...)
      
      
    } else {
      p=.extraFeaturePlot(expData=inputSeurat,umapData=inputDimData, features=inputGenes[inputGenes!=""],reductionName=reductionName,combine = combine_figs,order=order,...)
    }
  }
 
  
  return(p)
}

.extraFeaturePlot2=function(index,expDataList,umapData,features,combine_figs=T,order=T,...){
  
  umapData=umapData[row.names(umapData) %in% colnames(expDataList[[index]]),]
  if(length(features)==1){
    p=.extraFeaturePlot(expDataList[[index]],umapData = umapData,title_prefix=names(expDataList)[index],features = features,combine = combine_figs,order=order,...)
  } else {
    p=.extraFeaturePlot(expDataList[[index]],umapData = umapData,title_prefix=names(expDataList)[index],features = features,combine = combine_figs,order=order,...)
  }
   
  return(p)
}

#dims = c(1, 2); cells = NULL; cols=c("lightgrey","#0D0887FF","#9C179EFF","#ED7953FF","#F0F921FF"); pt.size = NULL; order = TRUE; min.cutoff = NA; max.cutoff = NA; 
#reduction = NULL; split.by = NULL; shape.by = NULL; slot = "data"; 
#blend = FALSE; blend.threshold = 0.5; label = FALSE; label.size = 4; 
#repel = FALSE; ncol = NULL; coord.fixed = FALSE; by.col = TRUE; 
#sort.cell = NULL; title_prefix=NULL
#expData=inputSeurat;umapData=inputDimData; features=inputGenes[inputGenes!=""];reductionName=reductionName;combine = combine_figs

.extraFeaturePlot=function (expData,umapData, features,reductionName="umap", dims = c(1, 2), cells = NULL, cols=c("lightgrey","#0D0887FF","#9C179EFF","#ED7953FF","#F0F921FF"), pt.size = NULL, order = TRUE, min.cutoff = NA, max.cutoff = NA, 
reduction = NULL, split.by = NULL, shape.by = NULL, slot = "data", 
blend = FALSE, blend.threshold = 0.5, label = FALSE, label.size = 4, 
repel = FALSE, ncol = NULL, coord.fixed = FALSE, by.col = TRUE, 
sort.cell = NULL, combine = TRUE,title_prefix=NULL) 
{
  legend ="yes"
  if(is.null(cols)){
    cols = if (blend) {
      c("lightgrey", "#ff0000", "#00ff00")
    } else {
      c("lightgrey", "blue")
    }
    legend ="none"
  }
  require(cowplot)
  require(ggplot2)
  
  features=row.names(expData)[toupper(row.names(expData)) %in% toupper(features)]
  
  object=expData
  if (!is.null(x = sort.cell)) {
    warning("The sort.cell parameter is being deprecated. Please use the order ", 
            "parameter instead for equivalent functionality.", call. = FALSE, immediate. = TRUE)
    if (isTRUE(x = sort.cell)) {
      order <- sort.cell
    }
  }
  no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
                    axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold", size = 14, margin = margin(r = 7)))
  
  reduction=reductionName
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  if (blend && length(x = features) != 2) {
    stop("Blending feature plots only works with two features")
  }
  if (blend) {
    default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
    cols <- switch(EXPR = as.character(x = length(x = cols)), 
                   `0` = {
                     warning("No colors provided, using default colors", 
                             call. = FALSE, immediate. = TRUE)
                     default.colors
                   }, `1` = {
                     warning("Only one color provided, assuming specified is double-negative and augmenting with default colors", 
                             call. = FALSE, immediate. = TRUE)
                     c(cols, default.colors[2:3])
                   }, `2` = {
                     warning("Only two colors provided, assuming specified are for features and agumenting with '", 
                             default.colors[1], "' for double-negatives", 
                             call. = FALSE, immediate. = TRUE)
                     c(default.colors[1], cols)
                   }, `3` = cols, {
                     warning("More than three colors provided, using only first three", 
                             call. = FALSE, immediate. = TRUE)
                     cols[1:3]
                   })
  }
  if (blend && length(x = cols) != 3) {
    stop("Blending feature plots only works with three colors; first one for negative cells")
  }
  dims <- toupper(paste0(reduction,"_", dims))
  if(is.null(cells)){
    cells = colnames(x = object)
  }
  
  if(!all(colnames(object)==row.names(umapData))){
    stop("Error in aligning the two datasets")
  }
  
  if(length(features)==1){
    tmpData=object@assays$RNA@data
    data=data.frame(x=as.matrix(tmpData)[which(row.names(tmpData)==features),],stringsAsFactors = F)
    colnames(data)=features
  } else {
    data <- as.data.frame(Seurat::FetchData(object = object, vars = features, cells = cells, slot = slot))
  }
  
  if (ncol(x = data) < 1) {
    stop("None of the requested features were found: ", paste(features, 
                                                              collapse = ", "), " in slot ", slot, call. = FALSE)
  }
  features <- colnames(x = data)
  features=gsub("-",".",features)
  colnames(data)=gsub("-",".",colnames(data))
  
  
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = min(data[,feature]), no = cutoff))
  }, cutoff = min.cutoff, feature = features)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = max(data[,feature]), no = cutoff))
  }, cutoff = max.cutoff, feature = features)
  check.lengths <- unique(x = vapply(X = list(features, min.cutoff, max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols, ]$maxcolors, no = length(x = cols))
  data <- sapply(X = 1:ncol(x = data), 
                 FUN = function(index) {
                   data.feature <- as.vector(x = data[, index])
                   min.use <- Seurat:::SetQuantile(cutoff = min.cutoff[index], data.feature)
                   max.use <- Seurat:::SetQuantile(cutoff = max.cutoff[index], data.feature)
                   data.feature[which(data.feature < min.use)] <- min.use
                   data.feature[which(data.feature > max.use)] <- max.use
                   return(data.feature)
                 })
  if(class(data)[1]=="list"){
    data=as.data.frame(do.call("cbind",data))
  }
  data=as.data.frame(data)
  colnames(x = data) <- features
  
  if(sum(colnames(umapData)=="ident")>0){
    umapData=umapData[,-which(colnames(umapData)=="ident")]
  }
  data=data.frame(umapData,ident=Seurat::Idents(object = object),data)
  data=as.data.frame(data)
  
  
  if (!all(dims %in% colnames(x = data))) {
    stop("The dimensions requested were not found", call. = FALSE)
  }
  rownames(x = data) <- cells
  
  data$split <- if (is.null(x = split.by)) {
    Seurat:::RandomName()
  }else {
    switch(EXPR = split.by, ident = Idents(object = object)[cells], 
           object[[split.by, drop = TRUE]][cells])
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  plots <- vector(mode = "list", length = ifelse(test = blend, 
                                                 yes = 4, no = length(x = features) * length(x = levels(x = data$split))))
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[,dims[1]])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[,dims[2]])))
  if (blend) {
    color.matrix <- BlendMatrix(two.colors = cols[2:3], col.threshold = blend.threshold, negative.color = cols[1])
    cols <- cols[2:3]
    colors <- list(color.matrix[, 1], color.matrix[1, ], as.vector(x = color.matrix))
  }
  
  for (i in 1:length(x = levels(x = data$split))) {
    ident <- levels(x = data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident,, drop = FALSE]
    if (blend) {
      features <- features[1:2]
      no.expression <- features[colMeans(x = data.plot[, features]) == 0]
      if (length(x = no.expression) != 0) {
        stop("The following features have no value: ", paste(no.expression, collapse = ", "), call. = FALSE)
      }
      data.plot <- cbind(data.plot[, c(dims, "ident")],Seurat:::BlendExpression(data = data.plot[, features[1:2]]))
      features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
    }
    for (j in 1:length(x = features)) {
      feature <- features[j]
      if (blend) {
        cols.use <- as.numeric(x = as.character(x = data.plot[,feature])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      } else {
        cols.use <- NULL
      }
      data.single <- data.plot[, c(dims, "ident", feature,shape.by)]
      plot <- Seurat:::SingleDimPlot(data = data.single, dims = dims, col.by = feature, order = order, pt.size = pt.size, 
                                     cols = cols.use, shape.by = shape.by, label = FALSE) + 
        scale_x_continuous(limits = xlims) + scale_y_continuous(limits = ylims) + 
        theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
      if (label) {
        plot <- LabelClusters(plot = plot, id = "ident", repel = repel, size = label.size)
      }
      if (length(x = levels(x = data$split)) > 1) {
        plot <- plot + theme(panel.border = element_rect(fill = NA, 
                                                         colour = "black"))
        plot <- plot + if (i == 1) {
          if(is.null(title_prefix)){
            labs(title = feature)
          } else {
            labs(title = paste0(title_prefix,"_",feature))
          }
          
        }
        else {
          labs(title = title_prefix)
        }
        if (j == length(x = features) && !blend) {
          suppressMessages(expr = plot <- plot + scale_y_continuous(sec.axis = dup_axis(name = ident), limits = ylims) + no.right)
        }
        if (j != 1) {
          plot <- plot + theme(axis.line.y = element_blank(), 
                               axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
                               axis.title.y.left = element_blank())
        }
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(axis.line.x = element_blank(), 
                               axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
                               axis.title.x = element_blank())
        }
      } else {
        if(is.null(title_prefix)){
          plot <- plot + labs(title = feature)
        } else {
          plot <- plot + labs(title = paste0(title_prefix,"_",feature))
        }
        
      }
      if (!blend) {
        plot <- plot + guides(col=NULL)
        cols.grad <- cols
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        } else if (length(x = cols) > 1) {
          unique.feature.exp <- unique(data.plot[, feature])
          if (length(unique.feature.exp) == 1) {
            warning("All cells have the same value (", unique.feature.exp, ") of ", feature, ".")
            if (unique.feature.exp == 0) {
              cols.grad <- cols[1]
            }
            else {
              cols.grad <- cols
            }
          }
          plot <- suppressMessages(expr = plot + scale_color_gradientn(colors = cols.grad, guide = "colorbar"))
        }
      }
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      plot <- plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
      names(plots)[(length(x = features) * (i - 1)) + j]=features[j]
    }
  }
  if (blend) {
    blend.legend <- BlendMap(color.matrix = color.matrix)
    for (ii in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(x = plots, 
                                              values = list(blend.legend + scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) > 1, yes = levels(x = data$split)[ii], no = "")), 
                                                                                              expand = c(0, 0)) + labs(x = features[1], y = features[2], title = if (ii == 1) {
                                                                                                paste("Color threshold:", blend.threshold)
                                                                                              } else {NULL}) + no.right), after = 4 * ii - 1))
    }
  }
  plots <- Filter(f = Negate(f = is.null), x = plots)
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = features) == 1) {
      ncol <- 1
    }
    if (length(x = features) > 6) {
      ncol <- 3
    }
    if (length(x = features) > 9) {
      ncol <- 4
    }
  }
  
  if (combine & length(plots)>1) {
    if (by.col && !is.null(x = split.by) && !blend) {
      plots <- lapply(X = plots, FUN = function(x) {
        return(suppressMessages(expr = x + theme_cowplot() + 
                                  ggtitle("") + scale_y_continuous(sec.axis = dup_axis(name = ""),limits = ylims) + no.right))
      })
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) + 
                 1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(expr = plots[[i]] + 
                                         scale_y_continuous(sec.axis = dup_axis(name = features[[idx]]), limits = ylims) + no.right)
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots)%%length(x = features) == 1)) {
        plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) + 
          theme(plot.title = element_text(hjust = 0.5))
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) + 
            theme(plot.title = element_text(hjust = 0.5))
          idx <- idx + 1
        }
        ncol <- 1
        nrow <- nsplits
      }
      else {
        nrow <- length(x = levels(x = data$split))#split.by %iff% length(x = levels(x = data$split))
      }
      plots <- plots[c(do.call(what = rbind, args = split(x = 1:length(x = plots), 
                                                          f = ceiling(x = seq_along(along.with = 1:length(x = plots))/length(x = features)))))]
      plots <- wrap_plots(plots, ncol = ncol)
      if (!is.null(x = legend) && legend == "none") {
        plots <- plots & NoLegend()
      }
    }
    else {
      plots <- wrap_plots(plots, ncol = ncol, length(x = levels(x = data$split)))#nrow = split.by %iff% length(x = levels(x = data$split)))
    }
    if (!is.null(x = legend) && legend == "none") {
      plots <- plots & NoLegend()
    }
  } else {
    if(combine&length(plots)==1){
      plots=plots[[1]]
    }
  }
  return(plots)
}

.extraAutoPointSize=function (data, raster = NULL) 
{
  return(ifelse(test = isTRUE(x = raster), yes = 1, no = min(1583/nrow(x = data), 
                                                             1)))
}

.mySingleDimPlot=function (data, dims, col.by = NULL, cols = NULL, pt.size = NULL, 
          shape.by = NULL, alpha.by = NULL, order = NULL, label = FALSE, 
          repel = FALSE, label.size = 4, cells.highlight = NULL, cols.highlight = "#DE2D26", 
          sizes.highlight = 1, na.value = "grey50", raster = NULL) 
{
  if(is.null(pt.size)){
    pt.size=.extraAutoPointSize(data = data, raster = raster)
  }
  
  if ((nrow(x = data) > 1e+05) & !isFALSE(raster)) {
    message("Rasterizing points since number of points exceeds 100,000.", 
            "\nTo disable this behavior set `raster=FALSE`")
  }
  if(is.null(raster)){
    raster=(nrow(x = data) > 1e+05)
  }
  
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  if (!is.data.frame(x = data)) {
    data <- as.data.frame(x = data)
  }
  if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
    stop("Cannot find dimensions to plot in data")
  } else if (is.numeric(x = dims)) {
    dims <- colnames(x = data)[dims]
  }
  if (!is.null(x = cells.highlight)) {
    highlight.info <- SetHighlight(cells.highlight = cells.highlight, 
                                   cells.all = rownames(x = data), sizes.highlight = sizes.highlight %||% 
                                     pt.size, cols.highlight = cols.highlight, col.base = cols[1] %||% 
                                     "#C3C3C3", pt.size = pt.size)
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- "highlight"
    pt.size <- highlight.info$size
    cols <- highlight.info$color
  }
  if (!is.null(x = order) && !is.null(x = col.by)) {
    if (typeof(x = order) == "logical") {
      if (order) {
        data <- data[order(!is.na(x = data[, col.by]), 
                           data[, col.by]), ]
      }
    } else {
      order <- rev(x = c(order, setdiff(x = unique(x = data[, 
                                                            col.by]), y = order)))
      data[, col.by] <- factor(x = data[, col.by], levels = order)
      new.order <- order(x = data[, col.by])
      data <- data[new.order, ]
      if (length(x = pt.size) == length(x = new.order)) {
        pt.size <- pt.size[new.order]
      }
    }
  }
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find ", col.by, " in plotting data, not coloring plot")
    col.by <- NULL
  }
  else {
    col.index <- match(x = col.by, table = colnames(x = data))
    if (grepl(pattern = "^\\d", x = col.by)) {
      col.by <- paste0("x", col.by)
    }
    else if (grepl(pattern = "-", x = col.by)) {
      col.by <- gsub(pattern = "-", replacement = ".", 
                     x = col.by)
    }
    colnames(x = data)[col.index] <- col.by
  }
  if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
    warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
  }
  if (!is.null(x = alpha.by) && !alpha.by %in% colnames(x = data)) {
    warning("Cannot find alpha variable ", alpha.by, " in data, setting to NULL", 
            call. = FALSE, immediate. = TRUE)
    alpha.by <- NULL
  }
  plot <- ggplot(data = data)
  plot <- if (isTRUE(x = raster)) {
    plot + scattermore::geom_scattermore(mapping = aes_string(x = dims[1], 
                                                 y = dims[2], color = paste0("`", col.by, "`"), shape = shape.by, 
                                                 alpha = alpha.by), pointsize = pt.size)
  }
  else {
    plot + geom_point(mapping = aes_string(x = dims[1], y = dims[2], 
                                           color = paste0("`", col.by, "`"), shape = shape.by, 
                                           alpha = alpha.by), size = pt.size)
  }
  plot <- plot + guides(color = guide_legend(override.aes = list(size = 3))) + 
    labs(color = NULL, title = col.by) + CenterTitle()
  if (label && !is.null(x = col.by)) {
    plot <- LabelClusters(plot = plot, id = col.by, repel = repel, 
                          size = label.size)
  }
  if (!is.null(x = cols)) {
    if (length(x = cols) == 1 && (is.numeric(x = cols) || 
                                  cols %in% rownames(x = brewer.pal.info))) {
      scale <- scale_color_brewer(palette = cols, na.value = na.value)
    }
    else if (length(x = cols) == 1 && (cols %in% c("alphabet", 
                                                   "alphabet2", "glasbey", "polychrome", "stepped"))) {
      colors <- DiscretePalette(length(unique(data[[col.by]])), 
                                palette = cols)
      scale <- scale_color_manual(values = colors, na.value = na.value)
    }
    else {
      scale <- scale_color_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot + theme_cowplot()
  return(plot)
}

.myDimPlotFn=function (object, dimCols = c("UMAP_1", "UMAP_2"), attCol = NULL, 
                       pt.size = NULL, label = TRUE, label.size = 4, raster=T,
                       repel = FALSE, cells.highlight = NULL, cols.highlight = "#DE2D26", 
                       sizes.highlight = 1, na.value = "grey50", ncol = NULL, combine = TRUE,set_col=T,col_vals=NULL) {
  #adapted from Seurat
  require(Seurat)
  require(hues)
  require(cowplot)
  #object: data.frame that includes the coordinates (dims) and the attribute column (only one attribute!)
  #dimCols: the column name of the coordinate data
  #attCol: the attribute column name that will be used for the coloring of the cells
  
  dims=dimCols
  colName=attCol
  
  
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  
  data=object
  
  if(!is.null(attCol)){
    data[,colName[1]]=as.factor(data[,colName[1]])
  }
  
  if(set_col){
    if(is.null(col_vals)&!is.null(attCol)){
      col_vals=hues::iwanthue(length(unique(data[,colName[1]])))
    }
    
    #data = data[, c(dims, colName[1], NULL, NULL)];
    #dims = dims; col.by = colName[1]; cols = col_vals; 
    #pt.size = pt.size; shape.by = NULL; order = NULL; 
    #label = FALSE; cells.highlight = cells.highlight; 
    #cols.highlight = cols.highlight; sizes.highlight = sizes.highlight; 
    #na.value = na.value
    
    if(!is.null(attCol)){
      plot <- .mySingleDimPlot(data = data[, c(dims, colName[1], NULL, NULL)],
                               dims = dims, col.by = colName[1], cols = col_vals, 
                               pt.size = pt.size, shape.by = NULL, order = NULL, 
                               label = FALSE, cells.highlight = cells.highlight, 
                               cols.highlight = cols.highlight, sizes.highlight = sizes.highlight, 
                               na.value = na.value,raster = raster)
    } else {
      plot <- .mySingleDimPlot(data = data[, c(dims, NULL, NULL)],
                               dims = dims, col.by = NULL, cols = col_vals, 
                               pt.size = pt.size, shape.by = NULL, order = NULL, 
                               label = FALSE, cells.highlight = cells.highlight, 
                               cols.highlight = cols.highlight, sizes.highlight = sizes.highlight, 
                               na.value = na.value,raster = raster)
    }
    
  } else {
    if(!is.null(attCol)){
      plot <- .mySingleDimPlot(data = data[, c(dims, colName[1], NULL, NULL)],
                               dims = dims, col.by = colName[1], cols = NULL, 
                               pt.size = pt.size, shape.by = NULL, order = NULL, 
                               label = FALSE, cells.highlight = cells.highlight, raster=raster,
                               cols.highlight = cols.highlight, sizes.highlight = sizes.highlight, 
                               na.value = na.value)
    } else {
      plot <- .mySingleDimPlot(data = data[, c(dims,  NULL, NULL)],
                               dims = dims, col.by = NULL, cols = NULL, 
                               pt.size = pt.size, shape.by = NULL, order = NULL, 
                               label = FALSE, cells.highlight = cells.highlight, raster=raster, 
                               cols.highlight = cols.highlight, sizes.highlight = sizes.highlight, 
                               na.value = na.value)
    }
    
  }
  
  if (label) {
    plot <- LabelClusters(plot = plot, id = colName[1], repel = repel, 
                          size = label.size, split.by = NULL)
  }
  
  return(plot)
}


.myheatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      ColSideColorsSize = 1,
                      RowSideColorsSize = 1,
                      KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}

#library(Rcpp)
#sourceCpp("~/code/matmult.cpp")

.scatterPlot_summary2d=function(object,reductionCols,n=300){
  
  if(sum(colnames(object) %in% reductionCols)!=2){
    stop("Specified reduction columns couldn't be found ")
  }
  
  
  yval=object[,reductionCols]
  xval=yval[,1]
  yval=yval[,2]
  
  rescoord=expand.grid(seq(min(xval),max(xval),length.out = n),seq(min(yval),max(yval),length.out = n))
  colnames(rescoord)=reductionCols
  
  knet=RANN::nn2(data=rescoord,query = object[,reductionCols],k=1,eps=0)
  counts=as.data.frame(table(knet$nn.idx[,1]))
  counts$Var1=as.character(counts$Var1)
  rescoord$id=1:nrow(rescoord)
  rescoord=merge(rescoord,counts,by.x="id",by.y="Var1")
  rescoord=rescoord[,-which(colnames(rescoord)=="id")]
  
  return(rescoord)
}

.myEvalMarkers=function (object, cells.1, cells.2=NULL, slot = "data", features = NULL,thresh.min = 0, pseudocount.use = 1,cluster_name=NULL,cells.1.weight.col=NULL) {
  require(Matrix)
  require(Seurat)
  features <- if(is.null(features)){
    features=rownames(x = object)
  }
  
  cells.1=as.character(cells.1)
  
  if(class(cells.2)==class(1)|class(cells.2)=="factor"){cells.2=as.character(cells.2)}
  
  if(length(cells.1)==1){
    if(is.null(cluster_name)){
      cluster_name=cells.1
    }
    cells.1=colnames(object)[which(Idents(object)==cells.1)]
  }
  
  if(is.null(cells.2)){
    cells.2=setdiff(colnames(object),cells.1)
  } else if(length(cells.2)==1&class(cells.2)==class("")){
    cells.2=colnames(object)[which(Idents(object)==cells.2)]
  }
  
  if(F){
    if(cells.1.weight<1){
      cells.1.weight=cells.1.weight/(1-cells.1.weight)*(length(cells.2))/(length(cells.1))
    }
  }
  
  
  
  data <- GetAssay(object = object,assay = "RNA")
  data=GetAssayData(data,slot = slot)
  
  cell.count.1 <- rowSums(x = data[features, cells.1, drop = FALSE] > thresh.min)
  cell.count.2 <- rowSums(x = data[features, cells.2, drop = FALSE] > thresh.min)
  
  pct.1 <- round(x = cell.count.1/(length(x = cells.1)), digits = 3)
  pct.2 <- round(x = cell.count.2/length(x = cells.2),  digits = 3)
  
  if(!is.null(cells.1.weight.col)){
    weight_data=object@meta.data[,cells.1.weight.col]
    
    mean.fxn <- switch(EXPR = slot, data = function(x,weight_data) {
      return(log(x = sparseMatrixStats::rowWeightedMeans(expm1(x = x),w=weight_data) + pseudocount.use))
    }, function(x) {
      return(log(x = sparseMatrixStats::rowWeightedMeans(x,w=weight_data) + pseudocount.use))
    })
    mean.data.1 <- mean.fxn(data[features, cells.1, drop = FALSE],weight_data=weight_data[match(cells.1,colnames(object))])
    mean.data.2 <- mean.fxn(data[features, cells.2, drop = FALSE],weight_data=weight_data[match(cells.2,colnames(object))])
    
  } else {
    mean.fxn <- switch(EXPR = slot, data = function(x) {
      return(log(x = rowMeans(x = (expm1(x = x))) + pseudocount.use))
    }, function(x) {
      return(log(x = rowMeans(x = x) + pseudocount.use))
    })
    mean.data.1 <- mean.fxn(data[features, cells.1, drop = FALSE])
    mean.data.2 <- mean.fxn(data[features, cells.2, drop = FALSE])
  }
  
  
  
  
  
  
  avg_logFC <- (mean.data.1 - mean.data.2)
  
  de.results=as.data.frame(cbind(pct.1,pct.2,mean.data.1,mean.data.2,avg_logFC,cell.count.1,cell.count.2))
  
  if(!is.null(cluster_name)){
    de.results$cluster_name=cluster_name
  }
  
  return(de.results)
}


.myEvalMarkers2=function (object1,object2, slot = "data", features = NULL,thresh.min = 0, pseudocount.use = 1,cluster_name=NULL) {
  require(Matrix)
  require(Seurat)
  features <- if(is.null(features)){
    features=rownames(x = object1)
  }
  
  cells.1=colnames(object1)
  cells.2=NULL
  if(class(object2)!=class(list())){
    cells.2=colnames(object2)
  }
  
  
  
  
  data1 <- GetAssay(object = object1,assay = "RNA")
  data1=GetAssayData(data1,slot = slot)
  pct.1 <- round(x = rowSums(x = data1[features, cells.1, drop = FALSE] > thresh.min)/length(x = cells.1), digits = 3)
  
  mean.fxn <- switch(EXPR = slot, data = function(x) {
    return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use))
  }, function(x) {
    return(log(x = rowMeans(x = x) + pseudocount.use))
  })
  
  if(class(object2)==class(list(""))){
    mean.fxn2 <- switch(EXPR = slot, data = function(x) {
      return(rowMeans(x = expm1(x = x)) + pseudocount.use)
    }, function(x) {
      return(rowMeans(x = x) + pseudocount.use)
    })
    
    sizeList=c()
    pct.2List=list()
    cell.count.2List=list()
    mean.data.2List=list()
    for(iobj in 1:length(object2)){
      data2 <- GetAssay(object = object2[[iobj]],assay = "RNA")
      data2=GetAssayData(data2,slot = slot)
      data3=Matrix::drop0(data2,tol=thresh.min)
      #pct.2List <- c(pct.2List,list(rowSums(x = data3)))
      sizeList=c(sizeList,length(x = ncol(data2)))
      cell.count.2List <- c(cell.count.2List,list(rowSums(x = data3)))
      mean.data.2List <- c(mean.data.2List,list(mean.fxn2(data2)*ncol(data2)))
    }
    
    for(icell in 1:length(cell.count.2List)){
      tmpcell=cell.count.2List[[icell]]
      tmpcell=tmpcell[match(features,names(tmpcell))]
      tmpcell[is.na(tmpcell)]=0
      cell.count.2List[[icell]]=tmpcell
    }
    
  } else {
    data2 <- GetAssay(object = object2,assay = "RNA")
    data2=GetAssayData(data2,slot = slot)
    cell.count.2 <- rowSums(x = data2[features, cells.2, drop = FALSE] > thresh.min)
    pct.2 <- round(cell.count.2/length(x = cells.2),  digits = 3)
    mean.data.2 <- mean.fxn(data2[features, cells.2, drop = FALSE])
  }
  
  cell.count.1 <- rowSums(x = data1[features, cells.1, drop = FALSE] > thresh.min)
  mean.data.1 <- mean.fxn(data1[features, cells.1, drop = FALSE])
  
  avg_logFC <- (mean.data.1 - mean.data.2)
  
  de.results=as.data.frame(cbind(pct.1,pct.2,mean.data.1,mean.data.2,avg_logFC,cell.count.1,cell.count.2))
  
  if(!is.null(cluster_name)){
    de.results$cluster_name=cluster_name
  }
  
  return(de.results)
}


.myEvalMarkers_contrastMaker=function (object2, slot = "data", features,thresh.min = 0, pseudocount.use = 1,cluster_name=NULL) {
  require(Matrix)
  require(Seurat)
  
  if(is.null(features)){
    stop("Wrong features input format!")
  }
  
  
  if(class(object2)!=class(list())){
    stop("Wrong input format!")
  }
  
  
  mean.fxn2 <- switch(EXPR = slot, data = function(x) {
    return(rowSums(x = expm1(x = x)))
  }, function(x) {
    return(rowSums(x = x))
  })
  
  sizeList=c()
  cell.count.2List=list()
  mean.data.2List=list()
  for(iobj in 1:length(object2)){
    data2 <- GetAssay(object = object2[[iobj]],assay = "RNA")
    data2=GetAssayData(data2,slot = slot)
    data3=Matrix::drop0(data2,tol=thresh.min)
    data3@x=rep(1,length(data3@x))
    sizeList=c(sizeList,ncol(data2))
    cell.count.2List <- c(cell.count.2List,list(rowSums(x = data3)))
    mean.data.2List <- c(mean.data.2List,list(mean.fxn2(data2)))
  }
  
  mean.data=cell.count=matrix(0,ncol=length(cell.count.2List),nrow=length(features))
  
  for(icell in 1:length(cell.count.2List)){
    tmpcell=cell.count.2List[[icell]]
    tmpcell=tmpcell[match(features,names(tmpcell))]
    if(sum(is.na(tmpcell))>0){
      tmpcell[is.na(tmpcell)]=0
    }
    cell.count[,icell]=tmpcell
    
    tmpcell=mean.data.2List[[icell]]
    tmpcell=tmpcell[match(features,names(tmpcell))]
    if(sum(is.na(tmpcell))>0){
      tmpcell[is.na(tmpcell)]=0
    }
    mean.data[,icell]=tmpcell
    
  }
  
  row.names(cell.count)=row.names(mean.data)=features
  
  
  return(list(cell.count=cell.count,sizeList=sizeList,mean.data=mean.data))
}

.myEvalMarkers_wcontrast=function (object, cells.1, cells.2=NULL,contrast_dataList, slot = "data", features = NULL,thresh.min = 0, pseudocount.use = 1,cluster_name=NULL,secondAnno=T) {
  require(Matrix)
  require(Seurat)
  #secondAnno: if the input object is composed of two or more classes; ie, if there is a need to further split the object
  
  features <- if(is.null(features)){
    features=rownames(x = object)
  }
  
  cells.1=as.character(cells.1)
  
  if(class(cells.2)==class(1)|class(cells.2)=="factor"){cells.2=as.character(cells.2)}
  
  if(length(cells.1)==1){
    if(is.null(cluster_name)){
      cluster_name=cells.1
    }
    cells.1=colnames(object)[which(Idents(object)==cells.1)]
  }
  
  if(secondAnno){
    if(is.null(cells.2)){
      cells.2=setdiff(colnames(object),cells.1)
    } else if(length(cells.2)==1&class(cells.2)==class("")){
      cells.2=colnames(object)[which(Idents(object)==cells.2)]
    }
  } else {
    cells.2=c()
  }
  
  
  data <- GetAssay(object = object,assay = "RNA")
  data=GetAssayData(data,slot = slot)
  
  cell.count.1 <- rowSums(x = data[features, cells.1, drop = FALSE] > thresh.min)
  
  if(secondAnno){
    cell.count.2 <- rowSums(x = data[features, cells.2, drop = FALSE] > thresh.min)
  } else {
    cell.count.2=rep(0,length(cell.count.1))
  }
  
  
  pct.1 <- round(x = cell.count.1/length(x = cells.1), digits = 3)
  #pct.2 <- round(x = cell.count.2/length(x = cells.2),  digits = 3)
  
  
  
  mean.fxn <- switch(EXPR = slot, data = function(x) {
    return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use))
  }, function(x) {
    return(log(x = rowMeans(x = x) + pseudocount.use))
  })
  
  mean.fxn2 <- switch(EXPR = slot, data = function(x) {
    return(rowSums(x = expm1(x = x)))
  }, function(x) {
    return(rowSums(x = x))
  })
  
  
  mean.data.1 <- mean.fxn(data[features, cells.1, drop = FALSE])
  
  if(secondAnno){
    mean.data.2 <- mean.fxn2(data[features, cells.2, drop = FALSE])
  } else {
    mean.data.2=rep(0,length(mean.data.1))
  }
  
  
  cell.count.contrast=rowSums(contrast_dataList$cell.count)
  cell.count.contrast=cell.count.contrast[match(features,names(cell.count.contrast))]
  if(sum(is.na(cell.count.contrast))>0){
    cell.count.contrast[is.na(cell.count.contrast)]=0
  }
  cell.count.2=cell.count.2+cell.count.contrast
  
  mean.data.contrast=contrast_dataList$mean.data
  mean.data.contrast=mean.data.contrast[match(features,row.names(mean.data.contrast)),]
  if(sum(is.na(mean.data.contrast))>0){
    mean.data.contrast[is.na(mean.data.contrast)]=0
  }
  mean.data.2=cbind(mean.data.2,mean.data.contrast)
  #mean.data.2=sweep(mean.data.2,2,c(length(cells.1),contrast_dataList$sizeList),"*")
  mean.data.2=rowSums(mean.data.2)/sum(c(length(cells.2),contrast_dataList$sizeList))
  mean.data.2=log(mean.data.2+ pseudocount.use)
  
  pct.2=cell.count.2/sum(c(length(cells.2),contrast_dataList$sizeList))
  
  
  avg_logFC <- (mean.data.1 - mean.data.2)
  
  de.results=as.data.frame(cbind(pct.1,pct.2,mean.data.1,mean.data.2,avg_logFC,cell.count.1,cell.count.2))
  
  if(!is.null(cluster_name)){
    de.results$cluster_name=cluster_name
  }
  
  return(de.results)
}


.myEvalMarkers_singleCluster=function (object, cells.1=NULL, slot = "data", features = NULL,thresh.min.1 = 0,thresh.min.2 = 4, pseudocount.use = 1,cluster_name=NULL) {
  require(Matrix)
  require(Seurat)
  #secondAnno: if the input object is composed of two or more classes; ie, if there is a need to further split the object
  
  features <- if(is.null(features)){
    features=rownames(x = object)
  }
  
  if(is.null(cells.1)){
    cells.1=colnames(object)
  }
  
  cells.1=as.character(cells.1)
  
  
  if(length(cells.1)==1){
    if(is.null(cluster_name)){
      cluster_name=cells.1
    }
    cells.1=colnames(object)[which(Idents(object)==cells.1)]
  }
  
  
  
  data <- GetAssay(object = object,assay = "RNA")
  data=GetAssayData(data,slot = slot)
  
  cell.count.1 <- rowSums(x = data[features, cells.1, drop = FALSE] > thresh.min.1)
  cell.count.2 <- rowSums(x = data[features, cells.1, drop = FALSE] > thresh.min.2)
  total.cell.count.1=rep(ncol(data),length(cell.count.1))
  
  pct.1 <- round(x = cell.count.1/length(x = cells.1), digits = 3)
  pct.2 <- round(x = cell.count.2/length(x = cells.1),  digits = 3)
  
  
  
  mean.fxn <- switch(EXPR = slot, data = function(x) {
    return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use))
  }, function(x) {
    return(log(x = rowMeans(x = x) + pseudocount.use))
  })
  
  mean.data.1 <- mean.fxn(data[features, cells.1, drop = FALSE])
  
  
  de.results=as.data.frame(cbind(pct.1,pct.2,mean.data.1,cell.count.1,cell.count.2,total.cell.count.1))
  colnames(de.results)[colnames(de.results)=="pct.1"]=paste0("pct.thr",thresh.min.1+1)
  colnames(de.results)[colnames(de.results)=="pct.2"]=paste0("pct.thr",thresh.min.2+1)
  colnames(de.results)[colnames(de.results)=="cell.count.1"]=paste0("cell.count.thr",thresh.min.1+1)
  colnames(de.results)[colnames(de.results)=="cell.count.2"]=paste0("cell.count.thr",thresh.min.2+1)
  
  
  if(!is.null(cluster_name)){
    de.results$cluster_name=cluster_name
  }
  
  return(de.results)
}


.myPseudoCellfn_v2=function(inputExpData,inputPCAembeddings=NULL,resolution=10, n.adaptiveKernel=5, nPropIter=3,nPCs=30,verbose=T,return_exp_data=T,ncores=10){
  #nPropIter: number of propagation iterations
  #inputPCAembeddings: PCA res
  #inputExpData: the expression data is expected to be raw counts (Seurat object)
  
  
  #There are two slightly different implementation of the algorithm
  #algorithm v2 accepts the resolution; the first version of the algorithm does a 10-to-1 merging
  #If there is a specific/optimized PC embedding available, algorithm will make use of it; otherwise, it will do the PC analysis
  #The input data is expected to be a Seurat object. Algorithm makes use of only the raw counts
  #resolution sets the number of cells to be merged together in an ideal scenario
  #Algorithm performs propagation on the similarity network to prune it and identify highly similar nodes
  #n.adaptiveKernel sets the radius for the propagation step. Larger numbers means larger sharing of data, and hece larger smoothing
  #nPropIter sets the number of propagation rounds to identify similar cells
  
  #Note: since the merging process is heuristic in nature, some iterations of the first algorithm may not converge. Therefore, an iteration is needed to restart it.
  
  #Output is a list of three components:
  #mappingRes: Each pseudo-cell is represented in a star structure in which Fnode is the center, and the Snode is the other cells included in the pseudo-cell.
  #notMapped: singletons
  #collapsedExpData: collapsed expression data in Seurat object; cell annotations are based on the available information on the center of the pseudocell; singletons and groups are specified in the colNames
  
  
  resolution=resolution-1
  nPCs=min(nPCs,ncol(inputExpData)-1)
  
  if(class(inputExpData)!="Seurat"){
    stop("Input exp data is expected to be a Seurat object!")
  }
  
  
  if(ncol(inputExpData)<100){
    n.adaptiveKernel=4
  }
  if(ncol(inputExpData)<50){
    n.adaptiveKernel=3
  }
  if(ncol(inputExpData)<25){
    n.adaptiveKernel=2
  }
  
  if(ncol(inputExpData)<18){
    n.adaptiveKernel=1
  }
  
  if(is.null(inputPCAembeddings)){
    tmpData=inputExpData
    tmpData = NormalizeData(tmpData)
    tmpData = FindVariableFeatures(tmpData, selection.method = "vst", nfeatures = 2000)
    tmpData <- ScaleData(tmpData)
    tmpData <- RunPCA(tmpData,npcs = nPCs)
    inputPCAembeddings=tmpData@reductions$pca@cell.embeddings[,1:nPCs]
  } else {
    nPCs=min(nPCs,ncol(inputPCAembeddings))
    inputPCAembeddings=inputPCAembeddings[,1:nPCs]
  }
  
  n.cells <- nrow(inputPCAembeddings)
  
  nn.ranked.org <- RANN::nn2(data = inputPCAembeddings, k = min(ncol(inputExpData),n.adaptiveKernel*12), eps = 0)
  
  nn.ranked=nn.ranked.org
  
  dists=nn.ranked$nn.dists
  affinities=t(apply(dists,1,function(x) exp((-1)*(x/x[max(n.adaptiveKernel,3)])^2)))
  affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
  affinities2=affinities-affinitiesThr
  affinities[affinities2<=0]=0
  
  affCounts=apply(affinities,2,function(x) sum(x==0))
  affinities=affinities[,-which(affCounts==nrow(affinities))]
  nn.ranked$nn.idx=nn.ranked$nn.idx[,-which(affCounts==nrow(affinities))]
  j <- as.numeric(t(nn.ranked$nn.idx))
  i <- ((1:length(j)) - 1)%/%ncol(affinities) + 1
  k=as.numeric(t(affinities))
  
  
  graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(inputPCAembeddings), nrow(inputPCAembeddings)))
  graph=(graph+t(graph))/2
  graph <- Matrix::Diagonal(x = 1 / (rowSums(graph))) %*% graph
  
  if(class(graph)!="dgCMatrix"){
    graph=as(graph,"dgCMatrix")
  }
  
  rownames(graph) <- rownames(inputPCAembeddings)
  colnames(graph) <- rownames(inputPCAembeddings)
  
  for(i in 1:nPropIter){
    if(verbose){
      print(paste("Prop",i))
    }
    
    if(i==1){
      prop_mat=graph
    } else {
      prop_mat=prop_mat%*%graph
    }
  }
  
  maxScores=spam::diag(prop_mat)
  
  affinityNorm <- Matrix::Diagonal(x = 1 / maxScores) %*% prop_mat
  
  affinityNorm2=affinityNorm * t(affinityNorm)
  affinityNorm2=sqrt(affinityNorm2)
  
  
  df.org=list()
  for(i in 2:ncol(nn.ranked.org$nn.idx)){
    tmp=data.frame(id1=nn.ranked.org$nn.idx[,1],id2=nn.ranked.org$nn.idx[,i],value=1)
    df.org=c(df.org,list(tmp))
  }
  df.org=do.call("rbind",df.org)
  
  df.org <- sparseMatrix(i = df.org$id1, j = df.org$id2, x = df.org$value,dims = c(nrow(affinityNorm), ncol(affinityNorm)))
  
  df.org=df.org*affinityNorm2
  df.org=(df.org+t(df.org))/2
  
  df.org_counts=t(df.org)
  df.org_counts=df.org_counts@p
  df.org_counts=df.org_counts[-1]- df.org_counts[-length(df.org_counts)]
  
  
  similarityQuantile=1
  df_counts=0
  while(mean(df_counts)<resolution&(round(similarityQuantile,1)>0.31)){
    df=summary(df.org)
    similarityQuantile=similarityQuantile-0.1
    simThr=as.numeric(quantile(df$x,similarityQuantile))
    df=df[which(df$x>=simThr),]
    df = sparseMatrix(i = df$i, j = df$j, x = df$x,dims = c(nrow(affinityNorm), ncol(affinityNorm)))
    
    df_counts=t(df)
    df_counts=df_counts@p
    df_counts=df_counts[-1]- df_counts[-length(df_counts)]
  }
  
  
  .PseudocellAssigner=function(df,resolution){
    
    df_counts=t(df)
    df_counts=df_counts@p
    df_counts=df_counts[-1]- df_counts[-length(df_counts)]
    
    mapping_res=data.frame(Fnode=0,Snode=0)
    for(i in order(df_counts,decreasing = F)){
      
      if(df_counts[i]>0& sum(c(mapping_res[,1],mapping_res[,2]) %in% i)==0){
        slNodes=setdiff(c(i,which(df[i,]>0)),c(mapping_res[,1],mapping_res[,2]))
        if(length(slNodes)>1){
          if(sum(slNodes %in% mapping_res[,1])==0){
            tmp_res=data.frame(Fnode=i,Snode=setdiff(slNodes,i),stringsAsFactors = F)
            tmp_res=tmp_res[order(df_counts[tmp_res$Snode],decreasing = F),]
            tmp_res=tmp_res[1: min(resolution,nrow(tmp_res)),]
            mapping_res=rbind(mapping_res,tmp_res)
          } else if(sum(slNodes %in% unique(mapping_res[,1]))==1){
            keyNode=intersect(slNodes , unique(mapping_res[,1]))
            if(sum(mapping_res[,1] %in% keyNode)<resolution){
              tmp_res=data.frame(Fnode=i,Snode=setdiff(slNodes,i))
              tmp_res=tmp_res[order(df_counts[tmp_res$Snode],decreasing = F),]
              tmp_res=tmp_res[1: min((resolution- sum(mapping_res[,1] %in% keyNode)),nrow(tmp_res)),]
              tmp_res=unique(c(tmp_res$Fnode,tmp_res$Snode))
              tmp_res=setdiff(tmp_res,keyNode)
              mapping_res=rbind(mapping_res,data.frame(Fnode=keyNode,Snode=tmp_res,stringsAsFactors = F))
            } else {
              tmp_res=data.frame(Fnode=i,Snode=setdiff(slNodes,keyNode),stringsAsFactors = F)
              tmp_res=tmp_res[order(df_counts[tmp_res$Snode],decreasing = F),]
              tmp_res=tmp_res[1: min(resolution,nrow(tmp_res)),]
              
              mapping_res=rbind(mapping_res,tmp_res)
            }
          } else {
            keyNodes=intersect(slNodes , unique(mapping_res[,1]))
            for(j in keyNodes){
              if(sum(mapping_res[,1]==j)<resolution){
                mapping_res=rbind(mapping_res,data.frame(Fnode=j,Snode=i,stringsAsFactors = F))
                break;
              }
            }
          }
        } else {
          keyNodes=intersect(which(df[i,]>0) , unique(mapping_res[,1]))
          if(length(keyNodes)==0){
            tmpNodes=intersect(which(df[i,]>0) , unique(mapping_res[,2]))
            keyNodes=unique(mapping_res[mapping_res[,2] %in% tmpNodes,1])
          }
          
          for(j in keyNodes){
            if(sum(mapping_res[,1]==j)<resolution){
              mapping_res=rbind(mapping_res,data.frame(Fnode=j,Snode=i,stringsAsFactors = F))
              break;
            }
          }
          
        }
      } else {
        slNodes=setdiff(c(i,which(df[i,]>0)),c(mapping_res[,1],mapping_res[,2]))
        if(length(slNodes)>0){
          slNodes=slNodes[order(df_counts[slNodes],decreasing = F)]
          keyNodes=c(intersect(i,df[,1]),mapping_res[mapping_res[,2]==i,1])
          for(j in keyNodes){
            if(length(slNodes)>0){
              if(sum(mapping_res[,1]==j)<resolution){
                tmpNode=slNodes[1:min(length(slNodes),(resolution - sum(mapping_res[,1]==j)))]
                mapping_res=rbind(mapping_res,data.frame(Fnode=j,Snode=tmpNode,stringsAsFactors = F))
                slNodes=setdiff(slNodes,tmpNode)
              }
            }
            
          }
        }
        
      }
    }
    mapping_res=mapping_res[-1,]
    
    included_cells=unique(c(mapping_res[,1],mapping_res[,2]))
    length(included_cells)
    gene_counts=as.data.frame(table(mapping_res[,1]))
    gene_counts[,1]=as.numeric(as.character(gene_counts[,1]))
    for(i in 1:max(resolution-2,1)){
      tmp_counts=gene_counts[which(gene_counts[,2]==i),1]
      for(j in tmp_counts){
        tmp_genes=which(df[j,]>0)
        key_list=unique(c(intersect(mapping_res[,1],tmp_genes),mapping_res[mapping_res[,2] %in% tmp_genes,1]))
        key_list=setdiff(key_list,j)
        if(length(key_list)>0){
          tmp_counts2=gene_counts[gene_counts[,1] %in% key_list,]
          tmp_counts2=tmp_counts2[order(tmp_counts2[,2],decreasing = F),]
          
          mapping_Snode=c(j,mapping_res[which(mapping_res[,1]==j),2])
          for(ic in 1:nrow(tmp_counts2)){
            if(tmp_counts2[ic,2]<(resolution+1)){
              
              to_include=mapping_Snode[1:min((resolution+3-tmp_counts2[ic,2]),length(mapping_Snode))]
              new_mapping=mapping_res[-which(mapping_res[,1]==j),]
              if(length(to_include)>0){
                mapping_res=rbind(new_mapping,data.frame(Fnode=tmp_counts2[ic,1],Snode=to_include,stringsAsFactors = F))
                break;
              }
            }
          }
          gene_counts=as.data.frame(table(mapping_res[,1]))
          gene_counts[,1]=as.numeric(as.character(gene_counts[,1]))
        }
        
      }
      
    }
    
    new_mapping=mapping_res[which(mapping_res[,1] %in% gene_counts[gene_counts[,2]>=0.7*resolution,1]),]
    excluded_genes=setdiff(1:ncol(df),c(new_mapping[,1],new_mapping[,2]))
    
    pre_excluded_genes=""
    while(abs(length(pre_excluded_genes)-length(excluded_genes))>(max(length(pre_excluded_genes),length(excluded_genes))*0.05)){
      pre_excluded_genes=excluded_genes
      for(i in excluded_genes){
        if(sum(df[i,]>0)>0){
          slNodes=which(df[i,]>0)
          slNodes=c(new_mapping[new_mapping[,1] %in% slNodes,1],new_mapping[new_mapping[,2] %in% slNodes,1])
          if(length(slNodes)>0){
            slNodes=as.data.frame(table(slNodes))
            slNodes[,1]=as.numeric(as.character(slNodes[,1]))
            slNodes=slNodes[slNodes[,1] %in% gene_counts[gene_counts[,2]<=resolution*1.2,1],]
            if(sum(slNodes[,1] %in% gene_counts[gene_counts[,2]<=resolution,1])>0){
              slNodes=slNodes[slNodes[,1] %in% gene_counts[gene_counts[,2]<=resolution,1],]
            }
            
            if(nrow(slNodes)>0){
              slNodes=slNodes[order(slNodes[,2],decreasing = T),]
              new_mapping=rbind(new_mapping,data.frame(Fnode=slNodes[1,1],Snode=i,stringsAsFactors = F))
              gene_counts[which(gene_counts[,1]==slNodes[1,1]),2]=gene_counts[which(gene_counts[,1]==slNodes[1,1]),2]+1
            }
          }
          
        }
        
      }
      
      mapping_res=new_mapping
      new_mapping=mapping_res[which(mapping_res[,1] %in% gene_counts[gene_counts[,2]>=0.7*resolution,1]),]
      excluded_genes=setdiff(1:ncol(df),c(new_mapping[,1],new_mapping[,2]))
    }
    return(mapping_res)
  }
  
  if(ncol(inputExpData)<(resolution*1.2)){
    mapping_res=data.frame(Fnode=1,Snode=2:ncol(inputExpData))
  } else {
    mapping_res=.PseudocellAssigner(df=df,resolution=resolution)
    excluded_genes=setdiff(1:ncol(df),c(mapping_res[,1],mapping_res[,2]))
    if(length(excluded_genes)>resolution/2){
      mapping_res2=.PseudocellAssigner(df=df[excluded_genes,excluded_genes],resolution=resolution)
      if(nrow(mapping_res2)>0){
        mapping_res2$Fnode=excluded_genes[mapping_res2$Fnode]
        mapping_res2$Snode=excluded_genes[mapping_res2$Snode]
        mapping_res=rbind(mapping_res,mapping_res2)
      }
    }
  }
  
  
  
  
  
  mapping_res$center=NA
  mapping_res$center_score=NA
  for(i in unique(mapping_res[,1])){
    tmp_genes=c(i,mapping_res[mapping_res[,1]==i,2])
    tmp=df[tmp_genes,tmp_genes]
    tmp=apply(tmp,1,function(x) sum(x>0))
    tmp_center=tmp_genes[which(tmp==max(tmp))[1]]
    mapping_res$center[mapping_res[,1]==i]=tmp_center
    mapping_res$center_score[mapping_res[,1]==i]=max(tmp)
    mapping_res$module_score[mapping_res[,1]==i]=sum(tmp!=0)/(length(tmp_genes)^2-length(tmp_genes))
  }
  center_ind=which(mapping_res$Snode==mapping_res$center)
  mapping_res$Snode[center_ind]=mapping_res$Fnode[center_ind]
  mapping_res$Fnode=mapping_res$center
  mapping_res=mapping_res[,-which(colnames(mapping_res)=="center")]
  
  mapping_res$node_score=NA
  for(i in unique(mapping_res[,1])){
    tmp_genes=mapping_res[mapping_res[,1]==i,2]
    tmp=df[c(tmp_genes),c(i,tmp_genes)]
    if(class(tmp)=="numeric"){
      tmp=sum(tmp>0)
    } else {
      tmp=apply(tmp,1,function(x) sum(x>0))
    }
    
    mapping_res$node_score[mapping_res[,1]==i]=tmp
  }
  
  
  
  mappingList=mapping_res
  
  df=summary(df.org)
  df=df[!df$i %in% c(mappingList$Fnode,mappingList$Snode),]
  
  if(sum(duplicated(c(unique(mappingList$Fnode),mappingList$Snode)))>0){
    stop("error in the mapping file2!")
  }
  
  df=df[!df$i %in% c(mappingList$Fnode,mappingList$Snode),]
  df=df[!df$j %in% c(mappingList$Fnode,mappingList$Snode),]
  if(sum(duplicated(c(unique(mappingList$Fnode),mappingList$Snode)))>0){
    stop("error in the mapping file3!")
  }
  
  if(verbose){
    print(paste("Coverage:",round(length(c(unique(mappingList$Fnode),mappingList$Snode))/nrow(inputPCAembeddings),3)))
    
  }
  df=unique(df$i)
  
  collapsedData=NULL
  if(return_exp_data){
    tstMap1=aggregate(Snode~Fnode,data=mappingList,function(x) x[1:min(length(x),20)])
    tstMap=cbind(tstMap1[,1],tstMap1[,2])
    tstMap=apply(tstMap,1,function(x) list(x[!is.na(x)]))
    
    tstMap=c(tstMap,lapply(df,function(x) x))
    
    inputExpMat=inputExpData@assays$RNA@counts
    .CollapseFn=function(indexList,inputExpMat){
      if(length(indexList)==1){
        x=indexList[[1]]
      } else {
        x=indexList
      }
      if(length(unlist(x))>1){
        values=rowSums(inputExpMat[,unlist(x)])
        values=data.frame(values=values,stringsAsFactors = F)
        colnames(values)=x[[1]]
      } else {
        values=as.numeric(inputExpMat[,unlist(x)])
        values=data.frame(values=values,stringsAsFactors = F)
        colnames(values)=unlist(x)
      }
      
      return(values)
    }
    if(ncores==1){
      tst=lapply(tstMap,.CollapseFn,inputExpMat=inputExpMat)
    } else {
      tst=parallel::mclapply(tstMap,.CollapseFn,inputExpMat=inputExpMat,mc.cores =ncores)
    }
    
    tst=do.call("cbind",tst)
    
    
    
    
    tstColnames=rep("singleton_",ncol(tst))
    tstColnames[colnames(tst) %in% tstMap1[,1]]="group_"
    tstColnames=paste0(tstColnames,colnames(tst))
    colnames(tst)=tstColnames
    
    tmp5=unlist(lapply(strsplit(colnames(tst),"_"),function(x) x[2]))
    tmp5=as.numeric(tmp5)
    pdata=inputExpData@meta.data[tmp5,]
    fdata=inputExpData@assays$RNA@meta.features
    row.names(pdata)=colnames(tst)
    collapsedData=Seurat::CreateSeuratObject(counts=as.matrix(tst),
                                             project = "SeuratProject",
                                             assay = "RNA",
                                             min.cells = 0,
                                             min.features = 0,
                                             names.field = 1,
                                             names.delim = "-",
                                             meta.data = pdata)
    
    collapsedData@assays$RNA@meta.features=fdata
    
  }
  
  mappingList$Fnode=colnames(inputExpData)[mappingList$Fnode]
  mappingList$Snode=colnames(inputExpData)[mappingList$Snode]
  
  return(list(mappingRes=mappingList,notMapped=df,collapsedExpData=collapsedData))
}

.myReadGMT=function (file,ignore.format=F) {
  if (!grepl("\\.gmt$", file)[1]&(!ignore.format)) {
    stop("Pathway information must be a .gmt file")
  }
  geneSetDB = readLines(file)
  geneSetDB = strsplit(geneSetDB, "\t")
  names(geneSetDB) = sapply(geneSetDB, "[", 1)
  geneSetDB = lapply(geneSetDB, "[", -1:-2)
  geneSetDB = lapply(geneSetDB, function(x) {
    x[which(x != "")]
  })
  return(geneSetDB)
}

.mySplitObject=function(object,colName){
  res=list()
  if(class(object)=="Seurat"){
    res=Seurat::SplitObject(object,split.by = colName)
  } else{
    if(F){
      pd=colData(object)[,colName]
      res=summary(counts(object))
      res=split(res,pd[res[,2]])
      res_pd=split(as.data.frame(colData(object)),pd)
      res=lapply(1:length(res),function(x){
        y=Matrix::sparseMatrix(i = res[[x]][,1],
                               j = as.numeric(as.factor(as.character(res[[x]][,2]))),
                               x = res[[x]][,3],dims = c(nrow(object),nrow(res_pd[[x]])))
        row.names(y)=row.names(object)
        colnames(y)=row.names(res_pd[[x]])
        y=SingleCellExperiment(assays = list(counts = y),colData = res_pd[[x]],rowData=as.data.frame(rowData(object)))
        return(y)
        
      })
    }
    
    
    if(T){
      pd=colData(object)[,colName]
      for(i in unique(pd)){
        res=c(res,list(object[,which(pd==i)]))
        names(res)[length(res)]=i
      }
    }
    
  }
  return(res)
}

.mySplitObject_v2=function(object,colName,min_dataset_size,ncores=5){

  if(class(object)=="Seurat"){
    
    res=Seurat::SplitObject(object,split.by = colName)

  } else{
    pd=colData(object)[,colName]

    res=as.data.frame(summary(counts(object)))    
    res$anno=pd[res[,2]]

    res=split(res,pd[res[,2]])
    res_pd=split(as.data.frame(colData(object)), pd)

    res=parallel::mclapply(1:length(res),function(x){
      
      y=Matrix::sparseMatrix(i = res[[x]][,1],
                             j = as.numeric(as.factor(as.character(res[[x]][,2]))),
                             x = res[[x]][,3],
                             dims = c(nrow(object), nrow(res_pd[[x]])))

      #res_pd[[x]]$pseudocell_size=nrow(res_pd[[x]])
      if(ncol(y)==1){#.5*min_dataset_size){
        #y=as(matrix(rowSums(y),ncol=1),"dgCMatrix")
        
        y=SingleCellExperiment(assays = list(counts = y),colData = res_pd[[x]][1,],rowData=as.data.frame(rowData(object)))
        
      } else {
        row.names(y)=row.names(object)
        colnames(y)=row.names(res_pd[[x]])
        y=SingleCellExperiment(assays = list(counts = y),colData = res_pd[[x]],rowData=as.data.frame(rowData(object)))
        
      }
      return(y)
      
    },mc.cores=ncores)


    if(F){
      for(i in unique(pd)){
        res=c(res,list(object[,which(pd==i)]))
        names(res)[length(res)]=i
      }
    }
    
  }
  
  size_dist=unlist(lapply(res,ncol))
  res=res[which(size_dist>=min_dataset_size)]

  return(res)
}

.myLigerToExpSet=function(inputData,organism,limit=NULL,add_anno=T){
  require(Matrix)
  expList=list()
  for(i in 1:length(inputData@raw.data)){
    tmp=inputData@raw.data[[i]]
    if(!is.null(limit)){
      if(ncol(tmp)>limit){
        tmp=tmp[,sample(ncol(tmp),min(limit,ncol(tmp)))]
      }
    }
    exp=SingleCellExperiment(assays = list(counts = tmp),colData = data.frame(sampleName=colnames(tmp),stringsAsFactors = F),rowData=data.frame(gene=row.names(tmp),stringsAsFactors = F))
    expList=c(expList,list(exp))
    names(expList)[length(expList)]=names(inputData@raw.data)[i]
  }
  
  expList=.mycBindFn(expList,batchNames = names(expList))
  
  pd=as.data.frame(inputData@cell.data)
  clusters=inputData@clusters
  clusters=clusters[match(row.names(pd),names(clusters))]
  if(!all(row.names(pd)==names(clusters))){
    stop("Error!")
  }
  if(length(inputData@clusters)>0){
    pd$clusters=as.character(inputData@clusters)
  }
  
  pd=pd[match(expList$sampleName,row.names(pd)),]
  if(!all(row.names(pd)==expList$sample)){
    stop("Error!")
  }
  
  
  pd$anno_batch=expList$anno_batch
  
  row.names(pd)=as.character(expList$sampleName)
  
  data=.myExpSetCreatorFn(inputExpData=counts(expList),organism=organism,minExpCells=0,inputPdata=pd,inputFdata=as.data.frame(rowData(expList)),addExtraAnno=add_anno,server=T)
  return(data)
}

#input_seurat_object=tmp;suggest_annot=T;max_cluster_size=NA;min_cluster_size=10;ncores=5;marker_count_thr=100
.myFindAllMarkers=function(input_seurat_object,suggest_annot=T,max_cluster_size=NA,min_cluster_size=10,ncores=5,marker_count_thr=100){
  require(Matrix)
  require(plyr)
  
  myLocalMarkerSummaryfn=function(input_pct_list){
    
    res_markers=list()
    for(iclust in names(input_pct_list)){
      tmp_main=input_pct_list[[iclust]]
      tmp_main$gene=row.names(tmp_main)
      tmp_contrast=input_pct_list[names(input_pct_list)!=iclust]
      
      contrast_res=tmp_contrast[[1]]
      contrast_res$gene=row.names(contrast_res)
      contrast_res$pct.1=contrast_res$pct.1*contrast_res$total.cell.count.1
      contrast_res$mean.data.1=(exp(x = contrast_res$mean.data.1)-1)*contrast_res$total.cell.count.1
      if(length(tmp_contrast)>1){
        for(jclust in 2:length(tmp_contrast)){
          tmp=tmp_contrast[[jclust]]
          tmp$pct.1=tmp$pct.1*tmp$total.cell.count.1
          tmp$mean.data.1=(exp(x = tmp$mean.data.1)-1)*tmp$total.cell.count.1
          colnames(tmp)=paste0("new_",colnames(tmp))
          tmp$gene=row.names(tmp)
          tmp=merge(contrast_res,tmp,all=T,by="gene")
          for(j in c("pct.1","mean.data.1","cell.count.1","total.cell.count.1","new_pct.1","new_mean.data.1","new_cell.count.1","new_total.cell.count.1")){
            tmp[is.na(tmp[,j]),j]=0
          }
          tmp$pct.1=tmp$pct.1+tmp$new_pct.1
          tmp$mean.data.1=tmp$mean.data.1+tmp$new_mean.data.1
          tmp$total.cell.count.1=tmp$total.cell.count.1+tmp$new_total.cell.count.1
          tmp$cell.count.1=tmp$cell.count.1+tmp$new_cell.count.1
          tmp=tmp[,!grepl("^new_",colnames(tmp))]
          contrast_res=tmp
        }
      }
      contrast_res$pct.1=contrast_res$pct.1/contrast_res$total.cell.count.1
      contrast_res$mean.data.1=log((contrast_res$mean.data.1)/contrast_res$total.cell.count.1+1)
      colnames(contrast_res)=gsub("\\.1",".2",colnames(contrast_res))
      
      tmp_main=merge(tmp_main,contrast_res,by="gene",all=T)
      tmp_main$avg_logFC =tmp_main$mean.data.1 - tmp_main$mean.data.2
      tmp_main=tmp_main[,!grepl("^total\\.",colnames(tmp_main))]
      tmp_main$cluster=iclust
      res_markers=c(res_markers,list(tmp_main))
      names(res_markers)[length(res_markers)]=iclust
    }
    res_markers=do.call("rbind",res_markers)
    
    row.names(res_markers)=paste0(res_markers$cluster,"_",res_markers$gene)
    return(res_markers)
  }
  
  if(length(unique(as.character(Idents(input_seurat_object))))<2){
    stop("Idents of the input_seurat_object is not set")
  } else if(sum(is.na(Idents(input_seurat_object)))>0){
    stop("NA values were found in the object Ident!")
  }
  
  
  
  if(class(input_seurat_object@assays$RNA@data@x)=="dgCMatrix"){
    if(max(input_seurat_object@assays$RNA@data@x)>100){
      input_seurat_object=Seurat::NormalizeData(input_seurat_object)
    }
  }
  
  
  
  input_seurat_object$cluster=as.character(Idents(input_seurat_object))
  datalist=.mySplitObject(input_seurat_object,colName = "cluster")
  datalist_size=unlist(lapply(datalist,ncol))
  if(!is.na(min_cluster_size)){
    datalist=datalist[datalist_size>=min_cluster_size]
  }
  
  if(!is.na(max_cluster_size)){
    datalist=lapply(datalist,function(x){
      if(ncol(x)>max_cluster_size){
        x=x[,sample(ncol(x),max_cluster_size)] 
      }
      x
    })
  }
  
  if(length(datalist)==1){
    return(NULL)
  }
  
  res=list()
  for(i in seq(1,length(datalist),ncores)){
    tmp_res=parallel::mclapply(datalist[(i:min((i+ncores-1),length(datalist)))],.myEvalMarkers_singleCluster,mc.cores = ncores)
    names(tmp_res)=names(datalist)[(i:min((i+ncores-1),length(datalist)))]
    res=c(res,tmp_res)
  }
  
  for(i in 1:length(res)){
    tmp=res[[i]]
    tmp$cell.count.1=tmp$cell.count.thr1
    tmp$pct.1=tmp$pct.thr1
    tmp=tmp[,!colnames(tmp) %in% c("pct.thr1","pct.thr5","cell.count.thr1","cell.count.thr5")]
    res[[i]]=tmp
  }
  
  res_markers=myLocalMarkerSummaryfn(input_pct_list=res)
  
  fd=input_seurat_object@assays$RNA@meta.features
  row.names(fd)=gsub("\\_","-",row.names(fd))
  fd$rwname=row.names(fd)
  res_markers=merge(res_markers,fd,by.x="gene",by.y="rwname",all.x=T)
  
  
  col_order=c("cluster","gene")
  if(sum(colnames(res_markers)=="gene_short_name")>0){
    col_order=c(col_order,"gene_short_name")
  }
  col_order=c(col_order,"avg_logFC","pct.1","pct.2","mean.data.1","mean.data.2","cell.count.1","cell.count.2")
  col_order=c(col_order,setdiff(colnames(res_markers),col_order))
  col_order=col_order[col_order %in% colnames(res_markers)]
  res_markers=res_markers[,col_order]
  
  if(suggest_annot){
    res_markers=.myAnnotateFn(marker_data=res_markers,marker_count_thr=marker_count_thr,ncores=ncores)
  }
  
  res_markers=res_markers[order(res_markers$pct.2 - res_markers$pct.1, (-1)*res_markers$pct.1),]
  
  return(res_markers)
}

.myAnnotateFn=function(marker_data,marker_count_thr=100,pct.2_thr=0.5,ncores=10,gene_count_per_cluster=500){
  
  InternalORtestfn=function(marker_res,gs,marker_count_thr=100){
    
    or_res=lapply(names(gs),function(x) {
      if(length(gs[[x]])<3){
        return(NULL)
      }
      tmp=marker_res
      tmp$gs="no"
      tmp$gs[toupper(tmp$gene_short_name) %in% toupper(gs[[x]])]="yes"
      if(sum(tmp$gs=="yes")==0){
        return(NULL)
      }
      tmp_markers=tmp$gene_short_name[(1:nrow(tmp))<=marker_count_thr&tmp$gs=="yes"]
      tmp=fisher.test(table((1:nrow(tmp))<=marker_count_thr,tmp$gs))
      return(data.frame(gs=x,OR=tmp$estimate,pval=tmp$p.value,markers=paste(tmp_markers,collapse = ","),stringsAsFactors = F))
    })
    or_res=do.call("rbind",or_res)
    
    if(!is.null(or_res)){
      or_res=or_res[order(or_res$OR,decreasing=T),]
      if(sum(or_res$OR>2)==0){
        or_res=NULL
      } else {
        or_res=or_res[which(or_res$OR>2),]
        or_res=or_res[1,]
      }
    }
    return(or_res)
  }
  
  
  if(!file.exists("../data/gazestani_serverFiles/enrichr_azimuth_2021.gmt")){
    stop("Error! file ../data/gazestani_serverFiles/enrichr_azimuth_2021.gmt is missing")
  }
  
  ctype_markers=.sconline.GSEA.readGMT(file="../data/gazestani_serverFiles/enrichr_azimuth_2021.gmt",bkg_genes=NULL,min.gs.size=NULL,max.gs.size=NULL)
  
  marker_data$gene_short_name=toupper(marker_data$gene_short_name)
  marker_data_org=marker_data[order(marker_data$pct.2 - marker_data$pct.1,(-1)*marker_data$pct.1),]
  marker_data=marker_data_org[!is.na(marker_data_org$gene_short_name),]
  if(!is.null(pct.2_thr)){
    marker_data$gene_short_name[marker_data$pct.2>pct.2_thr]=paste0("gene_",1:sum(marker_data$pct.2>pct.2_thr))
  }
  
  if(F){
    for(iclust in rev(unique(as.character(marker_data$cluster)))){
      tmp_res=InternalORtestfn(marker_res = marker_data[which(marker_data$cluster==iclust),],gs=ctype_markers,marker_count_thr = marker_count_thr)
      print(paste(iclust,class(tmp_res)))
      tmp_res$cluster=iclust
      tmp_res
    }
  }
  anno_res=parallel::mclapply(unique(as.character(marker_data$cluster)),function(iclust,marker_data,ctype_markers,marker_count_thr){
    tmp_res=InternalORtestfn(marker_res = marker_data[which(marker_data$cluster==iclust),],gs=ctype_markers,marker_count_thr = marker_count_thr)
    if(is.null(tmp_res)){
      tmp_res=data.frame(cluster=iclust)
    } else {
      tmp_res$cluster=iclust
    }
    
    tmp_res
  },marker_data=marker_data,ctype_markers=ctype_markers,marker_count_thr=marker_count_thr,mc.cores = ncores)
  
  
  
  sl_ind=which(unlist(lapply(anno_res,function(x) !(class(x)[1] %in% c(class(data.frame()))))))
  if(length(sl_ind)>0){
    anno_res2=parallel::mclapply(unique(as.character(marker_data$cluster))[sl_ind],function(iclust,marker_data,ctype_markers,marker_count_thr){
      tmp_res=InternalORtestfn(marker_res = marker_data[which(marker_data$cluster==iclust),],gs=ctype_markers,marker_count_thr = marker_count_thr)
      tmp_res$cluster=iclust
      tmp_res
    },marker_data=marker_data,ctype_markers=ctype_markers,marker_count_thr=marker_count_thr,mc.cores = 1)
    anno_res[sl_ind]=anno_res2
    
  }
  library(plyr)
  
  anno_res=do.call("rbind.fill",anno_res)
  
  colnames(anno_res)[colnames(anno_res)=="gs"]="inferred_ctype_marker"
  if(length(setdiff(unique(as.character(marker_data$cluster)),anno_res$cluster))>0){
    anno_res=rbind(anno_res,
                   data.frame(cluster=setdiff(unique(as.character(marker_data$cluster)),anno_res$cluster),inferred_ctype="not identified",OR=1,pval=1,markers="",stringsAsFactors = F))
  }
  
  
  marker_data_org=merge(marker_data_org,anno_res,by="cluster",all.x=T)
  
  marker_data_org=marker_data_org[order(marker_data_org$pct.2 - marker_data_org$pct.1,(-1)*marker_data_org$pct.1),]
  return(marker_data_org)
}
.myAnnotateFn_markers=function(marker_data,marker_count_thr=100,pct.2_thr=0.5,ncores=10){
  
  InternalORtestfn=function(marker_res,gs,marker_count_thr=100){
    
    or_res=lapply(names(gs),function(x) {
      if(length(gs[[x]])<3){
        return(NULL)
      }
      tmp=marker_res
      tmp$gs="no"
      tmp$gs[toupper(tmp$gene_short_name) %in% toupper(gs[[x]])]="yes"
      if(sum(tmp$gs=="yes")==0){
        return(NULL)
      }
      tmp_markers=tmp$gene_short_name[(1:nrow(tmp))<=marker_count_thr&tmp$gs=="yes"]
      tmp=fisher.test(table((1:nrow(tmp))<=marker_count_thr,tmp$gs))
      return(data.frame(gs=x,OR=tmp$estimate,pval=tmp$p.value,markers=paste(tmp_markers,collapse = ","),stringsAsFactors = F))
    })
    or_res=do.call("rbind",or_res)
    
    if(!is.null(or_res)){
      or_res=or_res[order(or_res$OR,decreasing=T),]
      if(sum(or_res$OR>2)==0){
        or_res=NULL
      } else {
        or_res=or_res[which(or_res$OR>2),]
        or_res=or_res[1,]
      }
    }
    return(or_res)
  }
  
  
  if(!file.exists("../data/gazestani_serverFiles/enrichr_azimuth_2021.gmt")){
    stop("Error! file ../data/gazestani_serverFiles/enrichr_azimuth_2021.gmt is missing")
  }
  
  ctype_markers=.sconline.GSEA.readGMT(file="../data/gazestani_serverFiles/enrichr_azimuth_2021.gmt",bkg_genes=NULL,min.gs.size=NULL,max.gs.size=NULL)
  
  marker_data_org=marker_data[order(marker_data$pct.2 - marker_data$pct.1,(-1)*marker_data$pct.1),]
  marker_data=marker_data_org[!is.na(marker_data_org$gene_short_name),]
  if(!is.null(pct.2_thr)){
    marker_data$gene_short_name[marker_data$pct.2>pct.2_thr]=paste0("gene_",1:sum(marker_data$pct.2>pct.2_thr))
  }
  
  if(F){
    for(iclust in rev(unique(as.character(marker_data$cluster)))){
      tmp_res=InternalORtestfn(marker_res = marker_data[which(marker_data$cluster==iclust),],gs=ctype_markers,marker_count_thr = marker_count_thr)
      print(paste(iclust,class(tmp_res)))
      tmp_res$cluster=iclust
      tmp_res
    }
  }
  anno_res=parallel::mclapply(unique(as.character(marker_data$cluster)),function(iclust,marker_data,ctype_markers,marker_count_thr){
    tmp_res=InternalORtestfn(marker_res = marker_data[which(marker_data$cluster==iclust),],gs=ctype_markers,marker_count_thr = marker_count_thr)
    if(is.null(tmp_res)){
      tmp_res=data.frame(cluster=iclust)
    } else {
      tmp_res$cluster=iclust
    }
    
    tmp_res
  },marker_data=marker_data,ctype_markers=ctype_markers,marker_count_thr=marker_count_thr,mc.cores = ncores)
  
  
  
  sl_ind=which(unlist(lapply(anno_res,function(x) !(class(x)[1] %in% c(class(data.frame()))))))
  if(length(sl_ind)>0){
    anno_res2=parallel::mclapply(unique(as.character(marker_data$cluster))[sl_ind],function(iclust,marker_data,ctype_markers,marker_count_thr){
      tmp_res=InternalORtestfn(marker_res = marker_data[which(marker_data$cluster==iclust),],gs=ctype_markers,marker_count_thr = marker_count_thr)
      tmp_res$cluster=iclust
      tmp_res
    },marker_data=marker_data,ctype_markers=ctype_markers,marker_count_thr=marker_count_thr,mc.cores = 1)
    anno_res[sl_ind]=anno_res2
    
  }
  library(plyr)
  
  anno_res=do.call("rbind.fill",anno_res)
  colnames(anno_res)[colnames(anno_res)=="gs"]="inferred_ctype"
  if(length(setdiff(unique(as.character(marker_data$cluster)),anno_res$cluster))>0){
    anno_res=rbind(anno_res,
                   data.frame(cluster=setdiff(unique(as.character(marker_data$cluster)),anno_res$cluster),inferred_ctype="not identified",OR=1,pval=1,markers="",stringsAsFactors = F))
  }
  marker_data_org=merge(marker_data_org,anno_res,by="cluster",all.x=T)
  
  marker_data_org=marker_data_org[order(marker_data_org$pct.2 - marker_data_org$pct.1,(-1)*marker_data_org$pct.1),]
  return(marker_data_org)
}

.myMarkerCountSelectorFn=function(markers_org,markers2_org,marker_data_org){
  res=NULL
  for(gene_count_per_cluster in c(10,15,20,25,30,50,75,100,150,200,250,300,350,400,450,500)){
    markers=markers_org
    markers2=markers2_org
    marker_data=marker_data_org
    sl_genes=lapply(1:ncol(markers2),function(x){
      y=row.names(markers2)[order(markers2[,x],decreasing = T)]
      y[1:min(gene_count_per_cluster,sum(markers2[,x]>0))]
    })
    sl_genes=unlist(sl_genes)
    if(F){
      sl_genes=as.data.frame(table(sl_genes))
      sl_genes=as.character(sl_genes[sl_genes[,2]<3,1])
    } else {
      sl_genes=unique(sl_genes)
    }
    
    markers=markers[row.names(markers) %in% sl_genes,]
    
    tmp=marker_data
    tmp$gene_short_name=toupper(tmp$gene_short_name)
    row.names(markers)=toupper(row.names(markers))
    tmp=tmp[which(tmp$gene_short_name %in% row.names(markers)),]
    tmp=tmp[!duplicated(paste0(tmp$gene_short_name,"_",tmp$cluster)),]
    tmp=reshape2::dcast(gene_short_name~cluster,value.var = "pct.1",data=tmp)
    tmp[is.na(tmp)]=0
    row.names(tmp)=tmp[,1]
    tmp=as.matrix(tmp[,-1])
    
    
    tmp_markers=markers[match(row.names(tmp),row.names(markers)),]
    cosine_mat=matrix(0,nrow=ncol(tmp),ncol=ncol(tmp_markers))
    row.names(cosine_mat)=colnames(tmp)
    colnames(cosine_mat)=colnames(tmp_markers)
    
    {
      tmp_cosine=lsa::cosine(cbind(tmp,tmp_markers))
      tmp_cosine=tmp_cosine[1:ncol(tmp),(ncol(tmp)+1):ncol(tmp_cosine)]
      cosine_mat=tmp_cosine
    }
    
    
    tmp_cosine=(scale(cosine_mat)+t(scale(t(cosine_mat))))/sqrt(2)
    res=rbind(res,data.frame(thr=gene_count_per_cluster,scores=sum(apply(tmp_cosine,1,function(x) sum(x)>2)==1)))
  }
  res=res[order(res$score,decreasing = T),]
  return(res$thr[1])
}

.myIntegrative_oneline=function(input_data1,input_data2=NULL,batch_variable,include_human,include_mouse,do.harmony=F,HVG_gene_count=2000,nPCs=30,indScaling=T,min_dataset_size=30,create_seurat_obj=T,ncores=5){
  
  organism="unknown"
  if(include_human){
    organism="Human"
  } else if(include_mouse){
    organism="Mouse"
  }
  
  
  
  argList=list(indScaling=indScaling,commonExpressed=T,conservativeMapping=F,oldMapping=F,MGIsymbol=T,allGenesFraction=0.5,HVG_list=c(),HVG_method="vst",ncores=ncores,min_ds_size=min_dataset_size,includeHuman=include_human,includeMouse=include_mouse,newRun=T,excludeHighExp=F)
  if(!is.null(input_data2)){
    input_data1=.extra_sconline.exp_creatorFn(argList=argList,inputExpData=input_data1,batch_variable=batch_variable,organism = organism,addAnno=F)
    input_data1$data=lapply(1:length(input_data1$data),function(x){x=input_data1$data[[x]]
    x$ds_batch=paste0("data1")
    return(x)})
    input_data2=.extra_sconline.exp_creatorFn(argList=argList,inputExpData=input_data2,batch_variable=batch_variable,organism = organism,addAnno=F)
    input_data2$data=lapply(1:length(input_data2$data),function(x){x=input_data2$data[[x]]
    x$ds_batch=paste0("data2")
    return(x)})
    input_data=list(data=c(input_data1$data,input_data2$data),data_m=NULL)
  } else {
    input_data=.extra_sconline.exp_creatorFn(argList=argList,inputExpData=input_data1,batch_variable=batch_variable,organism = organism,addAnno=F)
    input_data$data=lapply(1:length(input_data$data),function(x){y=input_data$data[[x]]
    y$ds_batch=paste0("data_",x)
    return(y)})
    input_data=list(data=input_data$data,data_m=NULL)
  }
  
  
  
  input_data=.myHighVarGeneSlFn(input_data,dataorganism=organism,argList = argList,HVG_gene_count=HVG_gene_count)
  argList$HVG_list=.sconline.select_HVG_thr(data=input_data,nGene=min(HVG_gene_count,nrow(input_data$varFeatures)))[,1]
  
  #inputExpData=input_data;argList=argList;batch_variable=batch_variable;organism=organism;addAnno=F;hierarchical_mode=T;set_HVG_count_thr=HVG_gene_count
  pca_res=.sconline.pca(inputExpData=input_data,argList=argList,batch_variable=batch_variable,organism=organism,addAnno=F,hierarchical_mode=T,set_HVG_count_thr=HVG_gene_count)
  
  harmony_embeddings=NULL
  if(do.harmony){
    harmony_embeddings = harmony::HarmonyMatrix(pca_res$pca_res[,1:nPCs], pca_res$pd, batch_variable, do_pca = FALSE, verbose=FALSE)
    resUMAP=.reductionUMAPFn(harmony_embeddings,umap.method="uwot",testPCAembeddings=NULL)
  } else {
    resUMAP=.reductionUMAPFn(pca_res$pca_res[,1:nPCs],umap.method="uwot",testPCAembeddings=NULL)
  }
  resUMAP=resUMAP$embedding
  
  pd=pca_res$pd
  pd$UMAP_1=resUMAP[,1]
  pd$UMAP_2=resUMAP[,2]
  
  if(create_seurat_obj){
    inputExpData=.extraExport2SeuratFn(.mycBindFn(input_data$data))
    pd=pd[match(colnames(inputExpData),row.names(pd)),]
    umap.reduction <- Seurat::CreateDimReducObject(embeddings = as.matrix(pd[,c("UMAP_1","UMAP_2")]), 
                                                   key = "UMAP_", assay = "RNA", global = TRUE)
    
    inputExpData[["umap"]]=umap.reduction
    
    #Adding in the embedding data
    if(is.null(harmony_embeddings)){
      harmony_embeddings=pca_res$pca_res
    }
    harmony_embeddings=harmony_embeddings[match(colnames(inputExpData),row.names(harmony_embeddings)),]
    colnames(harmony_embeddings)=gsub("PC_","",colnames(harmony_embeddings))
    reduction.data <- Seurat::CreateDimReducObject(embeddings = harmony_embeddings, 
                                                   assay = "RNA",  
                                                   key = "PC_")
    inputExpData[["pca"]] <- reduction.data
    
    return(inputExpData)
    
  } else {
    return(list(pca=pca_res$pca_res,pd=pd,harmony_embeddings=harmony_embeddings))
  }
}


cat("Main Functions:\n.myRead10X()\n.myRead10X_h5()\n.myLigerToExpSet()\n.mycBindFn()\n.myExpSetCreatorFn()\n.myIntegrative_oneline()\n.myFindAllMarkers()\n.myAnnotateFn()\n.my2dPlot()\n.myPseudoCellfn_v2()\n.myLabelTransfer_harmony()\n.myLabelTransfer_liger()\n.myMapToHuman()\n.myRiverPlotFn()\n.myClusteringOptimizerFn()\n.myMarkerBasedAnalysisFn()\n.mycellAssignHeatmap()\n.myMetaMarkerFn()\n.myFindNeighbors()\n.myVlnPlot()\n.myFeaturePlot()\n.myheatmap.3()\n.myEvalMarkers()\n.myReadGMT()\n.mySplitObject()\n.myRTNgsea(two-sided GSEA implementation)\n")

