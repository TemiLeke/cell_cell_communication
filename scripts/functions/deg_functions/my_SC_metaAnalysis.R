.myMouseFileList=function(){
  humanFilePath="~/myBucket/microgliaData/Human/"
  mouseFilePath="~/myBucket/microgliaData/Mouse/"
  mouseFile=data.frame(fileName="Bhattacherjee_31519873",path=file.path(mouseFilePath,"Bhattacherjee_31519873/microglia_data.rds"),stringsAsFactors = F)
  mouseFile=rbind(mouseFile,data.frame(fileName="Sousa_30206190",path=file.path(mouseFilePath,"Sousa_30206190/data.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="Tepe_30517858",path=file.path(mouseFilePath,"Tepe_30517858/microglia_data_arranged.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="Chen_28355573",path=file.path(mouseFilePath,"Chen_28355573/microglia_data.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="Chen_30773317",path=file.path(mouseFilePath,"Chen_30773317/microglia_data.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="Zywitza_30485812",path=file.path(mouseFilePath,"Zywitza_30485812/microglia_data.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="Masuda_GSE120629",path=file.path(mouseFilePath,"Masuda/GSE120629_arranged_data.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="Masuda_GSE120744",path=file.path(mouseFilePath,"Masuda/GSE120744_arranged_data.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="Masuda_GSE120745",path=file.path(mouseFilePath,"Masuda/GSE120745_arranged_data.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="Hammond_30471926",path=file.path(mouseFilePath,"Hammond_scRNA_30471926/Hammond_data_selectedCells_inPaper.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="Keren-Shaul_28602351",path=file.path(mouseFilePath,"Keren-Shaul_scRNA_28602351/Shaul_data.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="Lau",path=file.path(mouseFilePath,"Lau/data_arranged.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="SalaFrigerio_exp1_31018141",path=file.path(mouseFilePath,"SalaFrigerio_scRNA_31018141/exp1/arranged_data.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="SalaFrigerio_exp2_31018141",path=file.path(mouseFilePath,"SalaFrigerio_scRNA_31018141/exp2/arranged_data.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="Sierksma_31951107",path=file.path(mouseFilePath,"Sierksma_31951107/arrangedData.rds"),stringsAsFactors = F))
  #mouseFile=rbind(mouseFile,data.frame(fileName="Van_Hove_App_PS1_16mo",path=file.path(mouseFilePath,"VanHove/App_PS1/16mo/data_arranged.rds"),stringsAsFactors = F))
  #mouseFile=rbind(mouseFile,data.frame(fileName="Van_Hove_App_PS1_9mo",path=file.path(mouseFilePath,"VanHove/App_PS1/9mo/data_arranged.rds"),stringsAsFactors = F))
  #mouseFile=rbind(mouseFile,data.frame(fileName="Van_Hove_CD45",path=file.path(mouseFilePath,"VanHove/CD45/data_arranged.rds"),stringsAsFactors = F))
  #mouseFile=rbind(mouseFile,data.frame(fileName="Van_Hove_Irf8",path=file.path(mouseFilePath,"VanHove/Irf8/data_arranged.rds"),stringsAsFactors = F))
  #mouseFile=rbind(mouseFile,data.frame(fileName="Van_Hove_dissociation",path=file.path(mouseFilePath,"VanHove/dissociation_induced/data_arranged.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="Wu_29024657",path=file.path(mouseFilePath,"Wu_29024657/data_microglia_GSE103976_MeA_Act-Seq_vs_Conventional_DGE.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="Ximerakis_31551601",path=file.path(mouseFilePath,"Ximerakis_scRNA_31551601/data_microglia.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="zeisel",path=file.path(mouseFilePath,"zeisel/microglia_data.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="zhou",path=file.path(mouseFilePath,"Zhou/data_microglia.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="mathys_29020624",path=file.path(mouseFilePath,"Mathys/arranged_data.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="dulken_31270459",path=file.path(mouseFilePath,"Dulken/data_microglia.rds"),stringsAsFactors = F))
  mouseFile=rbind(mouseFile,data.frame(fileName="Mancuso_chimera",path=file.path(humanFilePath,"Mancuso/Mancuso_data_mouse_GSE137444_chimera.rds"),stringsAsFactors = F))
  return(mouseFile)
}
.myHumanFileList_old=function(){
  humanFilePath="~/myBucket/microgliaData/Human/"

  humanFiles=data.frame(fileName="Hodge_31435019",path=file.path(humanFilePath,"Hodge_31435019/brainSpan_microglia_Human_smartSeq.rds"),stringsAsFactors = F)
  humanFiles=rbind(humanFiles,data.frame(fileName="Lake_VIP_29227469",path=file.path(humanFilePath,"Lake_29227469/microglia_data_VIp.rds"),stringsAsFactors = F))
  humanFiles=rbind(humanFiles,data.frame(fileName="Grubman_31768052",path=file.path(humanFilePath,"Grubman/microglia_data.rds"),stringsAsFactors = F))
  humanFiles=rbind(humanFiles,data.frame(fileName="Leng_33432193",path=file.path(humanFilePath,"Leng_33432193/microglia_data.rds"),stringsAsFactors = F))
  #humanFiles=rbind(humanFiles,data.frame(fileName="Tushar",path=file.path(humanFilePath,"Tushar/microglia_data.rds"),stringsAsFactors = F))
  #humanFiles=rbind(humanFiles,data.frame(fileName="Mancuso_human_patient",path=file.path(humanFilePath,"Mancuso/Mancuso_data_GSE137444_human_patient.rds"),stringsAsFactors = F))
  #humanFiles=rbind(humanFiles,data.frame(fileName="Mancuso_invitro_mh_mc",path=file.path(humanFilePath,"Mancuso/Mancuso_data_GSE137444_invitro_mh_mc.rds"),stringsAsFactors = F))
  humanFiles=rbind(humanFiles,data.frame(fileName="Velmeshev_31097668",path=file.path(humanFilePath,"Velmeshev_31097668/microglia_data.rds"),stringsAsFactors = F))
  humanFiles=rbind(humanFiles,data.frame(fileName="Mathys_31042697",path=file.path(humanFilePath,"Mathys_31042697/microglia_data.rds"),stringsAsFactors = F))
  humanFiles=rbind(humanFiles,data.frame(fileName="zhou_31932797",path=file.path(humanFilePath,"zhou_31932797/microglia_data.rds"),stringsAsFactors = F))
  return(humanFiles)
}

.extra_matrix_rowNorm=function(input_mat,rowValues=NULL){
  require(Matrix)
  if(is.null(rowValues)){
    mat_diag=Matrix::Diagonal(x=1/(rowSums(input_mat)+10E-10))
  } else {
    mat_diag=Matrix::Diagonal(x=rowValues)
  }
  
  row.names(mat_diag)=colnames(mat_diag)=row.names(input_mat)
  return(mat_diag %*% input_mat)
}

.extra_matrix_colNorm=function(input_mat,colValues){
  require(Matrix)
  mat_diag=Matrix::Diagonal(x=colValues)
  row.names(mat_diag)=colnames(mat_diag)=colnames(input_mat)
  return(input_mat %*% mat_diag)
}

.myCellClass_FileList=function(inputClass,inputOrganism=NULL,do.exclusion=T,second_phase=F){
  #inputClass=argList$cellClass;inputOrganism="Human";do.exclusion=T
  fileList=read.table("~/myBucket/SC_data/data_alignment_file_organized.txt",sep = "\t",header = T,stringsAsFactors = F)
  
  if(do.exclusion){
    fileList=fileList[which(fileList$Exclude=="No"),]
  }
  
  fileList=fileList[grepl("indScaling",fileList$prefix),]
  fileList$batches="No"
  fileList$batches[grepl("_batches",fileList$prop_file)]="Yes"
  fileList$batches[grepl("allcells",tolower(fileList$prop_file))]="Yes"
  fileList=fileList[order(fileList$batches,decreasing = T),]
  fileList$isLake=0
  fileList$isLake[grepl("Lake",fileList$dsName)]=1
  fileList=fileList[(!duplicated(paste0(fileList$dir_path,"_",fileList$dsName)))|fileList$isLake==1,]
  
  
  for(i in 1:nrow(fileList)){
    if(!grepl("_batch",tolower(fileList$prop_file[i]))){
      fileList$exp_file[i]=paste0(gsub("\\.rda","",fileList$exp_file[i]),"_final.rda")
    } else {
      fileList$exp_file[i]=paste0(gsub("\\.rda","",fileList$exp_file[i]),"_final_batches.rda")
    }
  }
  
  
  fileList=fileList[,c("dir_path","dsName","exp_file")]
  
  
  fileList=rbind(fileList,data.frame(dir_path="~/myBucket/SC_data/Human/brain/snRNA",dsName="Kamath_35513515",exp_file="data_arranged_updatedId_final_batches.rda",stringsAsFactors = F))
  fileList=rbind(fileList,data.frame(dir_path="~/myBucket/SC_data/Human/brain/snRNA",dsName="iNPH",exp_file="data_arranged_updatedId_final_batches.rda",stringsAsFactors = F))
  fileList=rbind(fileList,data.frame(dir_path="~/myBucket/SC_data/Human/glia/snRNAseq",dsName="Smajic_34919646",exp_file="data_arranged_updatedId_final_batches.rda",stringsAsFactors = F))
  fileList=rbind(fileList,data.frame(dir_path="~/myBucket/SC_data/Human/brain/snRNA",dsName="Gerrits_33609158",exp_file="data_arranged_updatedId_final_batches.rda",stringsAsFactors = F))
  fileList=rbind(fileList,data.frame(dir_path="~/myBucket/SC_data/Human/brain/snRNA",dsName="Hodge_31435019",exp_file="data_arranged_updatedId_final.rda",stringsAsFactors = F))
  fileList=rbind(fileList,data.frame(dir_path="~/myBucket/SC_data/Human/brain/snRNA",dsName="Bakken_34616062",exp_file="data_arranged_updatedId_final_batches.rda",stringsAsFactors = F))
  if(second_phase){
    fileList=rbind(fileList,data.frame(dir_path="~/myBucket/SC_data/Human/brain/snRNA",dsName="Allen_seaAD",exp_file="data_arranged_updatedId_final_batches.rda",stringsAsFactors = F))
    fileList=rbind(fileList,data.frame(dir_path="~/myBucket/SC_data/Human/glia/snRNAseq",dsName="Dolan_iMGL",exp_file="data_arranged_updatedId_final_batches.rda",stringsAsFactors = F))
    fileList=rbind(fileList,data.frame(dir_path="~/myBucket/SC_data/Mouse/brain/scRNA",dsName="Dhandapani_35649351",exp_file="data_arranged_updatedId_final_batches.rda",stringsAsFactors = F))
    fileList=rbind(fileList,data.frame(dir_path="~/myBucket/SC_data/Human/brain/snRNA",dsName="NPH_new_batch",exp_file="data_arranged_updatedId_final_batches.rda",stringsAsFactors = F))
    fileList=rbind(fileList,data.frame(dir_path="~/myBucket/SC_data/Human/brain/snRNA",dsName="Yang_35165441",exp_file="data_arranged_updatedId_final_batches.rda",stringsAsFactors = F))
    fileList=rbind(fileList,data.frame(dir_path="~/myBucket/SC_data/Human/brain/snRNA",dsName="Garcia_35158371_HD",exp_file="data_arranged_updatedId_final_batches.rda",stringsAsFactors = F))
    fileList=rbind(fileList,data.frame(dir_path="~/myBucket/SC_data/Human/brain/snRNA",dsName="Garcia_35158371_exvivo",exp_file="data_arranged_updatedId_final_batches.rda",stringsAsFactors = F))
    fileList=rbind(fileList,data.frame(dir_path="~/myBucket/SC_data/Human/glia/snRNAseq",dsName="Smith_34767070",exp_file="data_arranged_updatedId_final_batches.rda",stringsAsFactors = F))
    fileList=rbind(fileList,data.frame(dir_path="~/myBucket/SC_data/Human/brain/snRNA",dsName="Morabito_34239132",exp_file="data_arranged_updatedId_final_batches.rda",stringsAsFactors = F))
    fileList=rbind(fileList,data.frame(dir_path="~/myBucket/SC_data/Human/brain/snRNA",dsName="Sadik_35381189",exp_file="data_arranged_updatedId_final_batches.rda",stringsAsFactors = F))
    
  }
  
  includeMouse=T
  if(!is.null(inputOrganism)){
    if(tolower(inputOrganism)=="human"){
      includeMouse=F
    }
    fileList=fileList[grepl(tolower(inputOrganism),tolower(fileList$dir_path)),]
  }
  
  fileList$exp_file=gsub("\\.rda","\\.qs",fileList$exp_file)
  fileList$exp_file=paste0(inputClass,"_",fileList$exp_file)
  
  fileList$path=paste0(fileList$dir_path,"/",fileList$dsName,"/",fileList$exp_file)
  
  fileList$exists=F
  for(i in 1:nrow(fileList)){
    if(file.exists(fileList$path[i])){
      fileList$exists[i]=T
    }
  }
  
  fileList=fileList[fileList$exists,]
  
  
  colnames(fileList)[colnames(fileList)=="dsName"]="fileName"
  
  fileList$fileName[which(grepl("CerebellarHem",fileList$path))]="Lake_29227469_CerebellarHem"
  fileList$fileName[which(grepl("FrontalCortex",fileList$path))]="Lake_29227469_FrontalCortex"
  fileList$fileName[which(grepl("VisualCortex",fileList$path))]="Lake_29227469_VisualCortex"
  
  if(tolower(inputClass)=="mg"&includeMouse){
    
    mouseFilePath="~/myBucket/microgliaData/Mouse"
    humanFilePath="~/myBucket/microgliaData/Human"
    mouseFile=NULL
    #fileList=rbind(fileList,data.frame(dir_path="",fileName="",exp_file="",path="",exists="",stringsAsFactors = F))
    mouseFile=rbind(mouseFile,data.frame(fileName="Chen_28355573",path=file.path("~/myBucket/SC_data/organized_data/Mouse/brain/scRNA","Chen_28355573"),exp_file="microglia_data.qs",stringsAsFactors = F))
    mouseFile=rbind(mouseFile,data.frame(fileName="Zywitza_30485812",path=file.path("~/myBucket/SC_data/organized_data/Mouse/brain/scRNA","Zywitza_30485812"),exp_file="microglia_data.qs",stringsAsFactors = F))
    mouseFile=rbind(mouseFile,data.frame(fileName="Hammond_30471926",path=file.path("~/myBucket/SC_data/organized_data/Mouse/glia/scRNA","Hammond_30471926"),exp_file="Hammond_data_selectedCells_inPaper.qs",stringsAsFactors = F))
    mouseFile=rbind(mouseFile,data.frame(fileName="Keren-Shaul_28602351",path=file.path("~/myBucket/SC_data/organized_data/Mouse/glia/scRNA","Keren-Shaul_28602351"),exp_file="Shaul_data.qs",stringsAsFactors = F))
    mouseFile=rbind(mouseFile,data.frame(fileName="SalaFrigerio_exp1_31018141",path=file.path("~/myBucket/SC_data/organized_data/Mouse/glia/scRNA","SalaFrigerio_31018141/exp1"),exp_file="arranged_data.qs",stringsAsFactors = F))
    mouseFile=rbind(mouseFile,data.frame(fileName="SalaFrigerio_exp2_31018141",path=file.path("~/myBucket/SC_data/organized_data/Mouse/glia/scRNA","SalaFrigerio_31018141/exp2"),exp_file="arranged_data.qs",stringsAsFactors = F))
    mouseFile=rbind(mouseFile,data.frame(fileName="Sierksma_31951107",path=file.path("~/myBucket/SC_data/organized_data/Mouse/glia/scRNA","Sierksma_31951107"),exp_file="arrangedData.qs",stringsAsFactors = F))
    #mouseFile=rbind(mouseFile,data.frame(fileName="Van_Hove_App_PS1_16mo",path=file.path(mouseFilePath,"VanHove/App_PS1/16mo"),exp_file="data_arranged.qs",stringsAsFactors = F))
    #mouseFile=rbind(mouseFile,data.frame(fileName="Van_Hove_App_PS1_9mo",path=file.path(mouseFilePath,"VanHove/App_PS1/9mo"),exp_file="data_arranged.qs",stringsAsFactors = F))
    #mouseFile=rbind(mouseFile,data.frame(fileName="Van_Hove_CD45",path=file.path(mouseFilePath,"VanHove/CD45"),exp_file="data_arranged.qs",stringsAsFactors = F))
    #mouseFile=rbind(mouseFile,data.frame(fileName="Van_Hove_Irf8",path=file.path(mouseFilePath,"VanHove/Irf8"),exp_file="data_arranged.qs",stringsAsFactors = F))
    #mouseFile=rbind(mouseFile,data.frame(fileName="Van_Hove_dissociation",path=file.path(mouseFilePath,"VanHove/dissociation_induced"),exp_file="data_arranged.qs",stringsAsFactors = F))
    mouseFile=rbind(mouseFile,data.frame(fileName="dulken_31270459",path=file.path("~/myBucket/SC_data/organized_data/Mouse/brain/scRNA","Dulken_31270459"),exp_file="data_microglia.qs",stringsAsFactors = F))
    mouseFile=rbind(mouseFile,data.frame(fileName="Masuda_GSE120745",path=file.path("~/myBucket/SC_data/organized_data/Mouse/glia/scRNA","Masuda_30760929"),exp_file="GSE120745_arranged_data.qs",stringsAsFactors = F))
    mouseFile=rbind(mouseFile,data.frame(fileName="Masuda_GSE120629",path=file.path("~/myBucket/SC_data/organized_data/Mouse/glia/scRNA","Masuda_30760929"),exp_file="GSE120629_arranged_data.qs",stringsAsFactors = F))
    mouseFile=rbind(mouseFile,data.frame(fileName="Masuda_GSE120744",path=file.path("~/myBucket/SC_data/organized_data/Mouse/glia/scRNA","Masuda_30760929"),exp_file="GSE120744_arranged_data.qs",stringsAsFactors = F))
    
    
    
    colnames(mouseFile)=c("fileName","dir_path","exp_file")
    mouseFile$path=file.path(mouseFile$dir_path,mouseFile$exp_file)
    mouseFile$exists=T
    fileList=rbind(fileList,mouseFile)
  }
  
  return(fileList)
}

.myHumanFileList=function(){
  humanFilePath="~/data/data/MG/MG_data_arranged"
  fileList=dir(humanFilePath)
  fileList=fileList[grepl("\\.rds",fileList)]
  humanFiles=NULL
  counter=0
  for(i in fileList){
    counter=counter+1
    humanFiles=rbind(humanFiles,data.frame(fileName=paste0("batch",floor(counter/4.01),"_",gsub("\\.rds","",i)),path=file.path(humanFilePath,i),stringsAsFactors = F))
  }
  humanFiles=rbind(humanFiles,.myHumanFileList_old())
  return(humanFiles)
}

.myReadDataFn=function(argList,expGenePercentage=NULL,rankedBased=F,convert_to_gene_level=T){
  #convert_to_gene_level: if multiple ensembl gene ids were present for the same entrez id, sum them up for each sample and report one representative
  
  removeHighExp = argList$excludeHighExp
  humanFiles=.myHumanFileList()
  mouseFile=.myMouseFileList()
  
  .myIndFileReaderFn=function(dsName,dfPath,organism,convertToHuman=F,rm_artifacts=F,rankedBased){
    index=strsplit(dsName,"_")
    index=unlist(lapply(index,function(x) paste(x[-1],collapse = "_")))
    index=which(dfPath$fileName==index)
    tmp=readRDS(dfPath$path[index])
    tmp=tmp[!is.na(rowData(tmp)$ensembl_gene_id),]
    x=apply(counts(tmp),1,function(x) sum(x>0))
    tmp=tmp[which(x>5),]
    if(sum(colnames(colData(tmp))=="QC_MT.pct")>0){
        if(sum(!is.na(tmp$QC_MT.pct))>10){
          if(organism=="Human"){
            tmp=tmp[,which(tmp$QC_MT.pct<10)]
          } else {
            tmp=tmp[,which(tmp$QC_MT.pct<5)]
          }
          
        }
        
      }
      
    if(sum(colnames(rowData(tmp))=="QC_mtGenes")>0){
      tmp=tmp[which(rowData(tmp)$QC_mtGenes=="No"),]
    }
    
    tmp=tmp[,which(tmp$QC_Gene_unique_count>500)]
    expMedian=median(as.numeric(as.character(tmp$QC_Gene_unique_count)),na.rm = T)
    expMad=mad(as.numeric(as.character(tmp$QC_Gene_unique_count)),na.rm = T)
    tmp=tmp[,which(tmp$QC_Gene_unique_count<(expMedian+3*expMad))]
    
      #tmp=tmp[,which(tmp$QC_Gene_unique_count>200&tmp$QC_Gene_unique_count<6000)]
    if(sum(duplicated(rowData(tmp)$ensembl_gene_id))>0){
      o=order(apply(counts(tmp),1,sum))
      tmp=tmp[o,]
      tmp=tmp[!duplicated(rowData(tmp)$ensembl_gene_id),]
    }
    if(rm_artifacts){
      artifact_genes=c("SST","PVALB","GAD2","OLIG1","LAMP5","VIP","SLC17A7")
      expdata=counts(tmp)[which(toupper(rowData(tmp)$gene_short_name) %in% artifact_genes),]
      if(class(expdata)=="numeric"){
        expdata=matrix(expdata,nrow=1)
      } else {
        expdata=as.matrix(expdata)
      }
      
      expdata=apply(expdata,2,sum)
      
      tmp=tmp[,expdata<1]
    }
    
    
      
    row.names(tmp)[!is.na(rowData(tmp)$ensembl_gene_id)]=rowData(tmp)$ensembl_gene_id
    if(convertToHuman){
      tmp=.myMapToHuman(tmp,server=T,conservativeMapping=argList$conservativeMapping,oldMapping=argList$oldMapping,MGIsymbol=argList$MGIsymbol)
    }
    
    if(rankedBased){
      tmpExp=as.matrix(counts(tmp))
      tmpExp=apply(tmpExp,2,function(x) {x[x!=0]=rank(x[x!=0]); return(x)})
      tmpExp=Matrix(tmpExp, sparse = TRUE)  
      counts(tmp)=tmpExp
    }
      
      row.names(tmp)=rowData(tmp)$ensembl_gene_id
      row.names(rowData(tmp))=rowData(tmp)$ensembl_gene_id
    if(nrow(tmp)<10000){
      print(paste("less than 10k genes in the dataset:",index))
    }
    tmp$anno_batch=dsName
    colnames(tmp)=paste0(dsName,"_",colnames(tmp))
    if(sum(colnames(colData(tmp))=="QC_Gene_total_count")>0){
      tmp$depth_per_gene=tmp$QC_Gene_total_count/(tmp$QC_Gene_unique_count+1)
    }
    
    return(tmp)
  }
  
  .myHighExpRemoverFn=function(inputList){
    top50list=c()
    for(i in 1:length(inputList)){
      top50list=c(top50list,rowData(inputList[[i]])$ensembl_gene_id[rowData(inputList[[i]])$QC_top50_expressed=="Yes"])
    }
    
    top50list=as.data.frame(table(top50list))
    top50list=top50list[order(top50list$Freq,decreasing = T),]
    top50list=as.character(top50list$top50list)[top50list$Freq>=(length(inputList)/2)]
    for(i in 1:length(inputList)){
      tmp2=inputList[[i]]
      tmp2=tmp2[-which(rowData(tmp2)$ensembl_gene_id %in% top50list),]
      inputList[[i]]=tmp2
    }
    
    return(inputList)
  }
  
  data=list()
  if(argList$includeHuman){
    #dsName=paste0("human_",humanFiles$fileName)[1];dfPath=humanFiles;convertToHuman=F;organism="Human";rankedBased=rankedBased
    tmp=parallel::mclapply(paste0("human_",humanFiles$fileName),.myIndFileReaderFn,dfPath=humanFiles,convertToHuman=F,organism="Human",rankedBased=rankedBased,mc.cores = argList$ncores)
    if(removeHighExp){
      tmp=.myHighExpRemoverFn(inputList = tmp)
    }
    data=c(data,tmp)
  }
  
  if(argList$includeMouse){
    #dfPath=rbind(mouseFile,.myHumanFileList_old());convertToHuman=argList$includeHuman;organism="Mouse";rankedBased=rankedBased
    tmp=parallel::mclapply(paste0("mouse_",mouseFile$fileName),.myIndFileReaderFn,dfPath=rbind(mouseFile,.myHumanFileList_old()),convertToHuman=argList$includeHuman,organism="Mouse",rankedBased=rankedBased,mc.cores = argList$ncores)
    if(removeHighExp){
      tmp=.myHighExpRemoverFn(inputList = tmp)
    }
    data=c(data,tmp)
    
    }
  
  for(i in 1:length(data)){
    tmp=data[[i]]
    tmp$ds_batch=tmp$anno_batch
    data[[i]]=tmp
  }
  
  
  #dsName=paste0("mouse_",mouseFile$fileName[4]);dfPath=mouseFile;affix="mouse";convertToHuman=.includeHuman
  
  if(argList$includeMouse){
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="mouse_Sierksma_31951107"))>0){
      tmpInd=which(tmpInd=="mouse_Sierksma_31951107")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpBatch=tmpData$batch
      tmpBatch[is.na(tmpBatch)]="NA"
      
      tmpDataList=list()
      for(ik in unique(tmpBatch)){
        tmp3=tmpData[,which(tmpBatch==ik)]
        tmp3$anno_batch=paste0(tmp3$anno_batch,"_",ik)
        tmpDataList=c(tmpDataList,list(tmp3))
      }
      data=c(data,tmpDataList)
    } else {
      warning("Sierksma_31951107 mouse dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="mouse_Sousa_30206190"))>0){
      tmpInd=which(tmpInd=="mouse_Sousa_30206190")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpBatch=tmpData$anno_batch
      tmpBatch[is.na(tmpBatch)]="NA"
      
      tmpDataList=list()
      for(ik in unique(tmpBatch)){
        tmp3=tmpData[,which(tmpBatch==ik)]
        tmp3$anno_batch=paste0(tmp3$anno_batch,"_",gsub("Sousa_30206190_","",ik))
        tmpDataList=c(tmpDataList,list(tmp3))
      }
      data=c(data,tmpDataList)
    } else {
      warning("Sousa (30206190) mouse dataset was not identified!")
    }
    
    {
      tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
      if(length(which(tmpInd=="mouse_Hammond_30471926"))>0){
        tmpInd=which(tmpInd=="mouse_Hammond_30471926")
        
        tmpData=data[[tmpInd]]
        data=data[-tmpInd]
        tmpBatch=paste0(as.character(tmpData$anno_age),"_",as.character(tmpData$anno_condition))
        
        tmpDataList=list()
        for(ik in unique(tmpBatch)){
          tmp3=tmpData[,which(tmpBatch==ik)]
          tmp3$anno_batch=paste0(tmp3$anno_batch,"_",ik)
          tmpDataList=c(tmpDataList,list(tmp3))
        }
        data=c(data,tmpDataList)
      } else {
        warning("Hammond mouse dataset was not identified!")
      }
    }
    
  }
  
  
  if(argList$includeHuman){
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Mancuso_invitro_mh_mc"))>0){
      tmpInd=which(tmpInd=="human_Mancuso_invitro_mh_mc")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpBatch=tmpData$anno_condition
      
      tmpDataList=list()
      for(ik in unique(tmpBatch)){
        tmp3=tmpData[,which(tmpBatch==ik)]
        tmp3$anno_batch=paste0(tmp3$anno_batch,"_",ik)
        tmpDataList=c(tmpDataList,list(tmp3))
      }
      data=c(data,tmpDataList)
    } else {
      warning("Human Mancuso_invitro_mh_mc dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_zhou_31932797"))>0){
      tmpInd=which(tmpInd=="human_zhou_31932797")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpBatch=tmpData$libraryBatch
      
      tmpDataList=list()
      for(ik in unique(tmpBatch)){
        tmp3=tmpData[,which(tmpBatch==ik)]
        tmp3$anno_batch=paste0(tmp3$anno_batch,"_",ik)
        tmpDataList=c(tmpDataList,list(tmp3))
      }
      data=c(data,tmpDataList)
    } else {
      warning("Human zhou_31932797 dataset was not identified!")
    }
    
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Leng_33432193"))>0){
      tmpInd=which(tmpInd=="human_Leng_33432193")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpBatch=tmpData$anno_region
      
      tmpDataList=list()
      for(ik in unique(tmpBatch)){
        tmp3=tmpData[,which(tmpBatch==ik)]
        tmp3$anno_batch=paste0(tmp3$anno_batch,"_",ik)
        tmpDataList=c(tmpDataList,list(tmp3))
      }
      data=c(data,tmpDataList)
    } else {
      warning("Human Leng_33432193 dataset was not identified!")
    }
    
    if(F){
      tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
      if(length(which(tmpInd=="human_Leng_33432193"))>0){
        tmpInd=which(tmpInd=="human_Leng_33432193")
        
        tmpData=data[[tmpInd]]
        data=data[-tmpInd]
        tmpBatch=tmpData$anno_batch
        
        tmpDataList=list()
        for(ik in unique(tmpBatch)){
          tmp3=tmpData[,which(tmpBatch==ik)]
          tmp3$anno_batch=paste0("human_",ik)
          tmpDataList=c(tmpDataList,list(tmp3))
        }
        data=c(data,tmpDataList)
      } else {
        warning("Human Leng_33432193 dataset was not identified!")
      }
    }
    
    
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Tushar"))>0){
      tmpInd=which(tmpInd=="human_Tushar")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpBatch=tmpData$anno_batch
      
      
      tmpDataList=list()
      for(ik in unique(tmpBatch)){
        tmp3=tmpData[,which(tmpBatch==ik)]
        tmp3$anno_batch=paste0(tmp3$anno_batch,"_",ik)
        tmpDataList=c(tmpDataList,list(tmp3))
      }
      data=c(data,tmpDataList)
    } else {
      warning("Human Finnish dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Mancuso_human_patient"))>0){
      tmpInd=which(tmpInd=="human_Mancuso_human_patient")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpBatch=as.character(tmpData$title)
      
      tmpDataList=list()
      for(ik in unique(tmpBatch)){
        tmp3=tmpData[,which(tmpBatch==ik)]
        tmp3$anno_batch=paste0(tmp3$anno_batch,"_",ik)
        tmpDataList=c(tmpDataList,list(tmp3))
      }
      data=c(data,tmpDataList)
    } else {
      warning("Mancuso_human_patients dataset was not identified!")
    }
    
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Grubman_31768052"))>0){
      tmpInd=which(tmpInd=="human_Grubman_31768052")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpBatch=as.character(tmpData$title)
      
      tmpDataList=list()
      for(ik in unique(tmpBatch)){
        tmp3=tmpData[,which(tmpBatch==ik)]
        tmp3$anno_batch=paste0(tmp3$anno_batch,"_",ik)
        tmpDataList=c(tmpDataList,list(tmp3))
      }
      data=c(data,tmpDataList)
    } else {
      warning("Human Grubman_31768052 dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Velmeshev_31097668"))>0){
      tmpInd=which(tmpInd=="human_Velmeshev_31097668")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpBatch=as.character(tmpData$Seqbatch)
      
      tmpDataList=list()
      for(ik in unique(tmpBatch)){
        tmp3=tmpData[,which(tmpBatch==ik)]
        tmp3$anno_batch=paste0(tmp3$anno_batch,"_",ik)
        tmpDataList=c(tmpDataList,list(tmp3))
      }
      data=c(data,tmpDataList)
    } else {
      warning("Human Velmeshev_31097668 dataset was not identified!")
    }
    
  }
  
  data2=data
  if(!is.null(expGenePercentage)){
    gnames=unlist(lapply(data2,function(x) row.names(x)))
    gnames=as.data.frame(table(gnames))
    gnames=gnames[order(gnames$Freq,decreasing = T),]
    gnames=gnames[gnames$Freq>=expGenePercentage*length(data2),]
    gnames=as.character(gnames$gnames)
    for(idata in 1:length(data2)){
      tmp=data[[idata]]
      tmp=tmp[row.names(tmp) %in% gnames,]
      data2[[idata]]=tmp
    }
  }
  
  data_m=.mycBindFn(inputList = data2,batchNames = names(data2)) #44858 347459
  
  if(convert_to_gene_level){
    #library(googleCloudStorageR)
    if(!file.exists("~/serverFiles/mapping_external_gene_names.rda")){
      stop("Error! file ~/serverFiles/mapping_external_gene_names.rda is missing")
    }
    
    load("~/serverFiles/mapping_external_gene_names.rda")
    
    fd=as.data.frame(rowData(data_m))
    fd=merge(data.frame(name=row.names(fd),fd),mapping,by="ensembl_gene_id",all.x=T)
    fd$external_gene_name=toupper(fd$external_gene_name)
    
    fd=fd[!is.na(fd$external_gene_name),]
    fd=fd[fd$external_gene_name!="",]
    fd=fd[fd$external_gene_name %in% fd$external_gene_name[duplicated(fd$external_gene_name)],]
    fd=fd[order(fd$external_gene_name),]
    
    res_m=matrix(0,nrow=length(unique(fd$external_gene_name)),ncol=ncol(data_m))
    counter=0
    nameList=list()
    for(igene in unique(fd$external_gene_name)){
      counter=counter+1
      tmp=counts(data_m)[which(rowData(data_m)$ensembl_gene_id %in% fd$ensembl_gene_id[which(fd$external_gene_name==igene)]),]
      tmp_fd=rowData(data_m)[which(rowData(data_m)$ensembl_gene_id %in% fd$ensembl_gene_id[which(fd$external_gene_name==igene)]),]
      nameList=c(nameList,list(as.data.frame(tmp_fd[1,])))
      if(rankedBased){
        tmp=as.matrix(tmp)
        tmp=apply(tmp,2,function(x) max(x))
      } else {
        tmp=Matrix::colSums(tmp)
      }
      
      res_m[counter,]=tmp
    }
    
    fd_m=do.call("rbind",nameList)
    fd$name=as.character(fd$name)
    fd_c=as.data.frame(rowData(data_m))
    fd_c=fd_c[-which(row.names(fd_c) %in% fd$name),]
    exp=counts(data_m)[-which(row.names(data_m) %in% fd$name),]
    exp=rbind(exp,as(res_m, "sparseMatrix"))
    fd_c=rbind(fd_c,fd_m)
    pd=as.data.frame(colData(data_m))
    row.names(exp)=row.names(fd_c)
    data_m=SingleCellExperiment(assays = list(counts = exp),colData = pd,rowData=fd_c)
    
    for(i in 1:length(data)){
      tmp=data[[i]]
      tmp=tmp[row.names(tmp) %in% row.names(data_m),]
      data[[i]]=tmp
    }
    
  }
  
  organism=rep("human",ncol(data_m))
  organism[grepl("^mouse",data_m$anno_batch)]="mouse"
  
  colData(data_m)=colData(data_m)[,-which(grepl("organism",colnames(colData(data_m))))]
  data_m$organism=organism
  return(list(data=data,data_m=data_m))
}

#argList=.ArgList;expGenePercentage=NULL;rankedBased=F;convert_to_gene_level=F;returnDataM=F;humanFiles=NULL;mouseFile=NULL
#expGenePercentage=NULL;rankedBased=F;convert_to_gene_level=T;returnDataM=T;breakSubjects=T;humanFiles=NULL;mouseFile=NULL;second_phase=T
.myReadDataFn_cellClass=function(argList,humanFiles=NULL,mouseFile=NULL,expGenePercentage=NULL,rankedBased=F,convert_to_gene_level=T,returnDataM=T,breakSubjects=T,second_phase=T){
  #convert_to_gene_level: if multiple ensembl gene ids were present for the same entrez id, sum them up for each sample and report one representative
  
  if(sum(names(argList)=="Leng200")==0){
    argList$Leng200=F
  }
  
  if(sum(names(argList)=="Jakel200")==0){
    argList$Jakel200=F
  }
  
  removeHighExp = argList$excludeHighExp
  if(is.null(humanFiles)){
    humanFiles=.myCellClass_FileList(inputClass=argList$cellClass,inputOrganism="Human",do.exclusion = T,second_phase=second_phase)
  }
  
  if(is.null(mouseFile)){
    mouseFile=.myCellClass_FileList(inputClass=argList$cellClass,inputOrganism="Mouse",do.exclusion = T,second_phase=second_phase)
  }
  
  .myIndFileReaderFn=function(dsName,dfPath,organism,convertToHuman=F,rm_artifacts=F,rankedBased,Leng200=F,Jakel200=F){
    require(qs)
    index=strsplit(dsName,"_")
    index=unlist(lapply(index,function(x) paste(x[-1],collapse = "_")))
    index=which(dfPath$fileName==index)
    tmplist=qread(dfPath$path[index])
    if(class(tmplist)!=class(list())){
      tmplist=list(tmplist)
    }
    
    for(ids in 1:length(tmplist)){
      tmp=tmplist[[ids]]
      tmp=tmp[!is.na(rowData(tmp)$ensembl_gene_id),]
      tmp$QC_totalCell_count=ncol(tmp)
      
      detectionRate=c()
      if(sum(colnames(rowData(tmp))=="QC_detectionRate")>0){
        detectionRate=rowData(tmp)$QC_detectionRate
        if(ncol(tmp)>10){
          ind=which(detectionRate>max(0.0001,5/ncol(tmp)))
          tmp=tmp[ind,]
          detectionRate=detectionRate[ind]
        }
        
      } else {
        
        for(i in seq(1,nrow(tmp),3000)){
          detectionRate=c(detectionRate,apply(as.matrix(counts(tmp)[i:min(i+3000-1,nrow(tmp)),]),1,function(x) sum(x>0)))
        }
        if(length(detectionRate)!=nrow(tmp)){
          stop("Error in counting genes expressed")
        }
        if(ncol(tmp)>5){
          ind=which(detectionRate>5)
          tmp=tmp[ind,]
          detectionRate=detectionRate[ind]
        }
        
      }
      
      if(T){
        if(sum(colnames(colData(tmp))=="QC_MT.pct")>0){
          if(sum(!is.na(tmp$QC_MT.pct))>10){
            if(organism=="Human"){
              tmp=tmp[,which(tmp$QC_MT.pct<10)]
            } else {
              tmp=tmp[,which(tmp$QC_MT.pct<5)]
            }
            
          }
          
        }
      }
      
      
      if(sum(colnames(rowData(tmp))=="QC_mtGenes")>0){
        ind=which(rowData(tmp)$QC_mtGenes=="No")
        tmp=tmp[ind,]
        detectionRate=detectionRate[ind]
      }
      
      if(ncol(tmp)>10){
        if(dsName!="human_Mathys_31042697"){
          if(Leng200&dsName=="human_Leng_33432193"){
            tmp=tmp[,which(tmp$QC_Gene_unique_count>200)]
            expMedian=median(as.numeric(as.character(tmp$QC_Gene_unique_count)),na.rm = T)
            expMad=mad(as.numeric(as.character(tmp$QC_Gene_unique_count)),na.rm = T)
            tmp=tmp[,which(tmp$QC_Gene_unique_count<=(expMedian+3*expMad))]
          } else if(Jakel200&dsName=="human_Jakel_30747918") {
            tmp=tmp[,which(tmp$QC_Gene_unique_count>200)]
            expMedian=median(as.numeric(as.character(tmp$QC_Gene_unique_count)),na.rm = T)
            expMad=mad(as.numeric(as.character(tmp$QC_Gene_unique_count)),na.rm = T)
            tmp=tmp[,which(tmp$QC_Gene_unique_count<=max(expMedian+3*expMad,1200))]
          } else if(dsName=="mouse_Keren-Shaul_28602351") {
            tmp=tmp[,which(tmp$QC_Gene_unique_count>200)]
            expMedian=median(as.numeric(as.character(tmp$QC_Gene_unique_count)),na.rm = T)
            expMad=mad(as.numeric(as.character(tmp$QC_Gene_unique_count)),na.rm = T)
            tmp=tmp[,which(tmp$QC_Gene_unique_count<=max(expMedian+3*expMad,1200))]
          } else {
            tmp=tmp[,which(tmp$QC_Gene_unique_count>500)]
            expMedian=median(as.numeric(as.character(tmp$QC_Gene_unique_count)),na.rm = T)
            expMad=mad(as.numeric(as.character(tmp$QC_Gene_unique_count)),na.rm = T)
            tmp=tmp[,which(tmp$QC_Gene_unique_count<=max(expMedian+3*expMad,1200))]
          }
        } else {
          tmp=tmp[,which(tmp$QC_Gene_unique_count>200)]
          expMedian=median(as.numeric(as.character(tmp$QC_Gene_unique_count)),na.rm = T)
          expMad=mad(as.numeric(as.character(tmp$QC_Gene_unique_count)),na.rm = T)
          tmp=tmp[,which(tmp$QC_Gene_unique_count<=max(expMedian+3*expMad,1200))]
        }
        
      }
      
      
      #tmp=tmp[,which(tmp$QC_Gene_unique_count>200&tmp$QC_Gene_unique_count<6000)]
      if(sum(duplicated(rowData(tmp)$ensembl_gene_id))>0){
        o=order(detectionRate,decreasing = T)
        tmp=tmp[o,]
        ind=which(!duplicated(rowData(tmp)$ensembl_gene_id))
        tmp=tmp[ind,]
        detectionRate=detectionRate[ind]
      }
      if(rm_artifacts){
        artifact_genes=c("SST","PVALB","GAD2","OLIG1","LAMP5","VIP","SLC17A7")
        expdata=counts(tmp)[which(toupper(rowData(tmp)$gene_short_name) %in% artifact_genes),]
        if(class(expdata)=="numeric"){
          expdata=matrix(expdata,nrow=1)
        } else {
          expdata=as.matrix(expdata)
        }
        
        expdata=apply(expdata,2,sum)
        
        tmp=tmp[,expdata<1]
      }
      
      row.names(tmp)[!is.na(rowData(tmp)$ensembl_gene_id)]=rowData(tmp)$ensembl_gene_id
      if(convertToHuman){
        tmp=.myMapToHuman(tmp,server=T,conservativeMapping=argList$conservativeMapping,oldMapping=argList$oldMapping,MGIsymbol=argList$MGIsymbol)
      }
      
      if(rankedBased){
        tmpExp=as.matrix(counts(tmp))
        tmpExp=apply(tmpExp,2,function(x) {x[x!=0]=rank(x[x!=0]); return(x)})
        tmpExp=Matrix(tmpExp, sparse = TRUE)  
        counts(tmp)=tmpExp
      }
      
      row.names(tmp)=rowData(tmp)$ensembl_gene_id
      row.names(rowData(tmp))=rowData(tmp)$ensembl_gene_id
      if(nrow(tmp)<10000){
        print(paste("less than 10k genes in the dataset:",dsName,"(genes:",nrow(tmp),")"))
      }
      if(sum(colnames(colData(tmp))=="anno_batch")>0){
        tmp$anno_batch2=tmp$anno_batch
      }
      
      tmp$anno_batch=dsName
      
      colnames(tmp)=paste0(dsName,"_",colnames(tmp))
      if(sum(colnames(colData(tmp))=="QC_Gene_total_count")>0){
        tmp$depth_per_gene=tmp$QC_Gene_total_count/(tmp$QC_Gene_unique_count+1)
      }
      
      
      tmp$QC_pass_basicQC=ncol(tmp)
      tmplist[[ids]]=tmp
    }
    
    return(tmplist)
  }
  
  
  .myHighExpRemoverFn=function(inputList){
    top50list=c()
    for(i in 1:length(inputList)){
      top50list=c(top50list,rowData(inputList[[i]])$ensembl_gene_id[rowData(inputList[[i]])$QC_top50_expressed=="Yes"])
    }
    
    top50list=as.data.frame(table(top50list))
    top50list=top50list[order(top50list$Freq,decreasing = T),]
    top50list=as.character(top50list$top50list)[top50list$Freq>=(length(inputList)/2)]
    for(i in 1:length(inputList)){
      tmp2=inputList[[i]]
      tmp2=tmp2[-which(rowData(tmp2)$ensembl_gene_id %in% top50list),]
      inputList[[i]]=tmp2
    }
    
    return(inputList)
  }
  #QC_totalCell_count
  #QC_pass_basicQC
  data=list()
  if(argList$includeHuman){
    #dsName=paste0("human_",humanFiles$fileName)[1];dfPath=humanFiles;convertToHuman=F;organism="Human";rankedBased=rankedBased;Leng200=F;Jakel200=F
    tmp=parallel::mclapply(paste0("human_",humanFiles$fileName),.myIndFileReaderFn,dfPath=humanFiles,convertToHuman=F,organism="Human",Leng200=argList$Leng200,rankedBased=rankedBased,mc.cores = argList$ncores)
    if(removeHighExp){
      tmp=.myHighExpRemoverFn(inputList = tmp)
    }
    data=c(data,tmp)
  }
  
  if(argList$includeMouse){
    #dsName=paste0("mouse_",mouseFile$fileName)[4];dfPath=mouseFile;convertToHuman=argList$includeHuman;organism="Mouse";rankedBased=rankedBased
    
    tmp=parallel::mclapply(paste0("mouse_",mouseFile$fileName),.myIndFileReaderFn,dfPath=mouseFile,convertToHuman=argList$includeHuman,organism="Mouse",rankedBased=rankedBased,mc.cores = argList$ncores)
    if(removeHighExp){
      tmp=.myHighExpRemoverFn(inputList = tmp)
    }
   
    data=c(data,tmp)
  }
  
  res_arranged=list()
  for(ids in 1:length(data)){
    tmp=data[[ids]]
    res_arranged=c(res_arranged,tmp)
  }
  data=res_arranged
  rm(res_arranged)
  if(F){
    res=NULL
    for(i in 1:length(data)){
      tmp=data[[i]]
      if(sum(colnames(colData(tmp))=="QC_MT.pct")>0){
        tmp_val=tmp$QC_MT.pct
        tmp_val2=tmp$QC_Gene_total_count
        tmp_val3=tmp$QC_Gene_unique_count
        tmp_val=tmp_val[!is.na(tmp_val)]
        tmp_val2=tmp_val2[!is.na(tmp_val2)]
        tmp_val3=tmp_val3[!is.na(tmp_val3)]
        res=rbind(res,data.frame(dataset=tmp$anno_batch[1],med_gene=median(tmp_val3),med_UMI=median(tmp_val2),min=min(tmp_val),qt_1=quantile(tmp_val,0.25),median=median(tmp_val,0.5),mean=mean(tmp_val),qt_3=quantile(tmp_val,0.75),max=max(tmp_val),stringsAsFactors = F,RNAseq_method=tmp$anno_RNAseq_type[1],RNAseq_method2=tmp$anno_RNAseq_type2[1]))
      }
      
    }
  }
  
  x_name=unlist(lapply(data,function(x) x$anno_batch[1]))
  x_total=unlist(lapply(data,function(x) x$QC_totalCell_count[1]))
  x_included=unlist(lapply(data,ncol))
  
  summary_statistic=data.frame(name=x_name,total=x_total,included=x_included,fraction=x_included/x_total,stringsAsFactors = F)
  
  
  #dsName=paste0("mouse_",mouseFile$fileName[4]);dfPath=mouseFile;affix="mouse";convertToHuman=.includeHuman
  
  
  for(i in 1:length(data)){
    tmp=data[[i]]
    tmp$ds_batch=tmp$anno_batch
    data[[i]]=tmp
  }
  
  if(argList$includeMouse&breakSubjects){
    
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="mouse_Dhandapani_35649351"))>0){
      tmpInd=which(tmpInd=="mouse_Dhandapani_35649351")
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData$anno_batch=paste0("mouse_Dhandapani_35649351_",tmpData$anno_batch2)
      tmpData=.mySplitObject(tmpData,"anno_batch")
      
      data=c(data,tmpData)
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="mouse_Habib_32341542"))>0){
      tmpInd=which(tmpInd=="mouse_Habib_32341542")
      
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData$anno_batch=paste0(tmpData$anno_batch,"_Mouse",tmpData$anno_batch2)
      tmpData=.mySplitObject(tmpData,"anno_batch")
      
      data=c(data,tmpData)
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="mouse_XPO"))>0){
      tmpInd=which(tmpInd=="mouse_XPO")
      
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      if(sum(colnames(colData(tmpData))=="dataset")>0){
        tmpData$anno_batch=tmpData$dataset
      }
      
      tmpData=.mySplitObject(tmpData,"anno_batch")
      
      data=c(data,tmpData)
    }
    
    tmpInd=as.character(unlist(lapply(data,function(x) as.character(x$anno_batch)[1])))
    if(length(which(tmpInd=="mouse_pfc_dev"))>0){
      tmpInd=which(tmpInd=="mouse_pfc_dev")
      
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData$anno_batch=tmpData$age
      tmpData=.mySplitObject(tmpData,"anno_batch")
      
      data=c(data,tmpData)
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="mouse_Sierksma_31951107"))>0){
      tmpInd=which(tmpInd=="mouse_Sierksma_31951107")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpBatch=tmpData$batch
      tmpBatch[is.na(tmpBatch)]="NA"
      
      tmpDataList=list()
      for(ik in unique(tmpBatch)){
        tmp3=tmpData[,which(tmpBatch==ik)]
        tmp3$anno_batch=paste0(tmp3$anno_batch,"_",ik)
        tmpDataList=c(tmpDataList,list(tmp3))
      }
      data=c(data,tmpDataList)
    } else {
      warning("Sierksma_31951107 mouse dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="mouse_Yao_34616066_BICCN_sn"))>0){
      tmpInd=which(tmpInd=="mouse_Yao_34616066_BICCN_sn")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData=.mySplitObject(tmpData,"anno_batch2")
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0(tmp3$anno_batch,"_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
      
    } else {
      warning("mouse BICCN dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="mouse_Yao_34616066_BICCN_sc10Xv2"))>0){
      tmpInd=which(tmpInd=="mouse_Yao_34616066_BICCN_sc10Xv2")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData=.mySplitObject(tmpData,"anno_batch2")
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0(tmp3$anno_batch,"_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
      
    } else {
      warning("mouse Yao_34616066_BICCN_sc10Xv2 dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="mouse_Yao_34616066_BICCN_sc10Xv3"))>0){
      tmpInd=which(tmpInd=="mouse_Yao_34616066_BICCN_sc10Xv3")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData=.mySplitObject(tmpData,"anno_batch2")
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0(tmp3$anno_batch,"_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
      
    } else {
      warning("mouse Yao_34616066_BICCN_sc10Xv3 dataset was not identified!")
    }
    
    if(F){
      #Excluded datasets
      tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
      if(length(which(tmpInd=="mouse_Sousa_30206190"))>0){
        tmpInd=which(tmpInd=="mouse_Sousa_30206190")
        
        tmpData=data[[tmpInd]]
        data=data[-tmpInd]
        tmpBatch=tmpData$anno_batch
        tmpBatch[is.na(tmpBatch)]="NA"
        
        tmpDataList=list()
        for(ik in unique(tmpBatch)){
          tmp3=tmpData[,which(tmpBatch==ik)]
          tmp3$anno_batch=paste0(tmp3$anno_batch,"_",gsub("Sousa_30206190_","",ik))
          tmpDataList=c(tmpDataList,list(tmp3))
        }
        data=c(data,tmpDataList)
      } else {
        warning("Sousa (30206190) mouse dataset was not identified!")
      }
      
      tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
      if(length(which(tmpInd=="mouse_Ding_32341560"))>0){
        tmpInd=which(tmpInd=="mouse_Ding_32341560")
        
        tmpData=data[[tmpInd]]
        data=data[-tmpInd]
        tmpData=.mySplitObject(tmpData,'anno_batch2')
        
        for(ik in 1:length(tmpData)){
          tmp3=tmpData[[ik]]
          tmp3$anno_batch=paste0(tmp3$anno_batch,"_",names(tmpData)[ik])
          data=c(data,tmp3)
        }
        
      } else {
        warning("Ding (32341560) mouse dataset was not identified!")
      }
      
      
      tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
      if(length(which(tmpInd=="mouse_Tepe_30517858"))>0){
        tmpInd=which(tmpInd=="mouse_Tepe_30517858")
        
        tmpData=data[[tmpInd]]
        data=data[-tmpInd]
        tmpTitle=as.character(tmpData$title)
        tmpTitle=unlist(lapply(strsplit(tmpTitle,"\\:"),function(x) x[1]))
        tmpData$title=tmpTitle
        tmpData=.mySplitObject(tmpData,'title')
        
        for(ik in 1:length(tmpData)){
          tmp3=tmpData[[ik]]
          tmp3$anno_batch=paste0(tmp3$anno_batch,"_",names(tmpData)[ik])
          data=c(data,tmp3)
        }
      } else {
        warning("mouse Tepe_30517858 dataset was not identified!")
      }
    }
    
    
   
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="mouse_Zhou_31932797"))>0){
      tmpInd=which(tmpInd=="mouse_Zhou_31932797")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData=.mySplitObject(tmpData,'anno_batch2')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0(tmp3$anno_batch,"_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
    } else {
      warning("mouse_Zhou_31932797 dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="mouse_dulken_31270459"))>0){
      tmpInd=which(tmpInd=="mouse_dulken_31270459")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData=.mySplitObject(tmpData,'sampleId')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0(tmp3$anno_batch,"_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
    } else {
      warning("mouse_Dulken dataset was not identified!")
    }
    
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="mouse_Zywitza_30485812"))>0){
      tmpInd=which(tmpInd=="mouse_Zywitza_30485812")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      
      tmpData=.mySplitObject(tmpData,'anno_batch2')
      names(tmpData)=paste0("mouse_",names(tmpData))
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=names(tmpData)[ik]
        data=c(data,tmp3)
      }
    } else {
      warning("mouse Zywitza_30485812 dataset was not identified!")
    }
    
    {
      tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
      if(length(which(tmpInd=="mouse_Hammond_30471926"))>0){
        tmpInd=which(tmpInd=="mouse_Hammond_30471926")
        
        tmpData=data[[tmpInd]]
        data=data[-tmpInd]
        #tmpBatch=paste0(as.character(tmpData$anno_age),"_",gsub(" ","_",as.character(tmpData$anno_condition)))
        tmpBatch=gsub(" ","_",as.character(tmpData$title))
        
        tmpDataList=list()
        for(ik in unique(tmpBatch)){
          tmp3=tmpData[,which(tmpBatch==ik)]
          tmp3$anno_batch=paste0(tmp3$anno_batch,"_",ik)
          tmpDataList=c(tmpDataList,list(tmp3))
        }
        data=c(data,tmpDataList)
      } else {
        warning("Hammond mouse dataset was not identified!")
      }
    }
    
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="mouse_Bhattacherjee_31519873"))>0){
      tmpInd=which(tmpInd=="mouse_Bhattacherjee_31519873")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData$anno_condition=tmpData$treatment
      tmpData$anno_condition[is.na(tmpData$anno_condition)]="WT"
      #tmpData$Sample=paste0(tmpData$anno_batch,'_',tmpData$anno_condition)
      tmpData=.mySplitObject(tmpData,'anno_condition')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0("mouse_Bhattacherjee_31519873_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
    } else {
      warning("Mouse Bhattacherjee_31519873 dataset was not identified!")
    }
    
  }
  
  
  if(argList$includeHuman&breakSubjects){
    
    if(F){
      #excluded datasets
      tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
      if(length(which(tmpInd=="human_Mancuso_invitro_mh_mc"))>0){
        tmpInd=which(tmpInd=="human_Mancuso_invitro_mh_mc")
        
        tmpData=data[[tmpInd]]
        data=data[-tmpInd]
        tmpBatch=tmpData$anno_condition
        
        tmpDataList=list()
        for(ik in unique(tmpBatch)){
          tmp3=tmpData[,which(tmpBatch==ik)]
          tmp3$anno_batch=paste0(tmp3$anno_batch,"_",ik)
          tmpDataList=c(tmpDataList,list(tmp3))
        }
        data=c(data,tmpDataList)
      } else {
        warning("Human Mancuso_invitro_mh_mc dataset was not identified!")
      }
      
      tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
      if(length(which(tmpInd=="human_Mancuso_human_patient"))>0){
        tmpInd=which(tmpInd=="human_Mancuso_human_patient")
        
        tmpData=data[[tmpInd]]
        data=data[-tmpInd]
        tmpBatch=as.character(tmpData$title)
        
        tmpDataList=list()
        for(ik in unique(tmpBatch)){
          tmp3=tmpData[,which(tmpBatch==ik)]
          tmp3$anno_batch=paste0(tmp3$anno_batch,"_",ik)
          tmpDataList=c(tmpDataList,list(tmp3))
        }
        data=c(data,tmpDataList)
      } else {
        warning("Mancuso_human_patients dataset was not identified!")
      }
      
      
      tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
      if(length(which(tmpInd=="human_Polioudakis_31303374"))>0){
        tmpInd=which(tmpInd=="human_Polioudakis_31303374")
        
        tmpData=data[[tmpInd]]
        data=data[-tmpInd]
        tmpData=.mySplitObject(tmpData,'Donor')
        
        for(ik in 1:length(tmpData)){
          tmp3=tmpData[[ik]]
          tmp3$anno_batch=paste0("human_Polioudakis_31303374_",names(tmpData)[ik])
          data=c(data,tmp3)
        }
      } else {
        warning("Human Polioudakis_31303374 dataset was not identified!")
      }
      
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Allen_seaAD"))>0){
      tmpInd=which(tmpInd=="human_Allen_seaAD")
      
      tmpData=data[tmpInd]
      data=data[-tmpInd]
      
      for(ibatch in 1:length(tmpData)){
        tmp=tmpData[[ibatch]]
        tmp$anno_batch=paste0("human_Allen_seaAD_",tmp$anno_batch2)
        data=c(data,list(tmp))
      }
     
      
    } else {
      warning("Human Allen_seaAD dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Dolan_iMGL"))>0){
      tmpInd=which(tmpInd=="human_Dolan_iMGL")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      
      tmpData=.mySplitObject(tmpData,'anno_batch2')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0("human_Dolan_iMGL_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
      
    } else {
      warning("human_Dolan_iMGL dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Garcia_35158371_HD"))>0){
      tmpInd=which(tmpInd=="human_Garcia_35158371_HD")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      
      tmpData=.mySplitObject(tmpData,'anno_batch2')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0("human_Garcia_35158371_HD_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
      
    } else {
      warning("human_Garcia_35158371_HD dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Garcia_35158371_exvivo" ))>0){
      tmpInd=which(tmpInd=="human_Garcia_35158371_exvivo" )
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      
      tmpData=.mySplitObject(tmpData,'anno_batch2')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0("human_Garcia_35158371_exvivo_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
      
    } else {
      warning("human_Garcia_35158371_exvivo dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Smith_34767070"))>0){
      tmpInd=which(tmpInd=="human_Smith_34767070")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      
      tmpData=.mySplitObject(tmpData,'anno_batch2')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0("human_Smith_34767070_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
      
    } else {
      warning("human_Garcia_35158371_exvivo dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Morabito_34239132"))>0){
      tmpInd=which(tmpInd=="human_Morabito_34239132")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      
      tmpData=.mySplitObject(tmpData,'SampleID')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0("human_Morabito_34239132_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
      
    } else {
      warning("human_Morabito_34239132 dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_NPH_new_batch"))>0){
      tmpInd=which(tmpInd=="human_NPH_new_batch")
      
      tmpData=data[tmpInd]
      data=data[-tmpInd]
      
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0("human_NPH_new_batch_",tmp3$anno_batch2)
        data=c(data,tmp3)
      }
      
    } else {
      warning("human_NPH_new_batch dataset was not identified!")
    }  
                    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Zhou_31932797"))>0){
      tmpInd=which(tmpInd=="human_Zhou_31932797")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      
      tmpData=.mySplitObject(tmpData,'anno_batch2')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0("human_Zhou_31932797_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
      
    } else {
      warning("Human Zhou_31932797 dataset was not identified!")
    }
    
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Leng_33432193"))>0){
      tmpInd=which(tmpInd=="human_Leng_33432193")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData=.mySplitObject(tmpData,'anno_batch2')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0("human_Leng_33432193_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
    } else {
      warning("Human Leng_33432193 dataset was not identified!")
    }
    
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Gerrits_33609158"))>0){
      tmpInd=which(tmpInd=="human_Gerrits_33609158")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData=.mySplitObject(tmpData,'anno_batch2')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0("human_Gerrits_33609158_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
    } else {
      warning("Human Gerrits_33609158 dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Smajic_34919646"))>0){
      tmpInd=which(tmpInd=="human_Smajic_34919646")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData=.mySplitObject(tmpData,'anno_batch2')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0("human_Smajic_34919646_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
    } else {
      warning("Human Smajic_34919646 dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_OteroGarcia_35882228"))>0){
      tmpInd=which(tmpInd=="human_OteroGarcia_35882228")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData=.mySplitObject(tmpData,'anno_batch2')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0("human_OteroGarcia_35882228_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
    } else {
      warning("Human OteroGarcia_35882228 dataset was not identified!")
    }
    
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_iNPH"))>0){
      tmpInd=which(tmpInd=="human_iNPH")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData=.mySplitObject(tmpData,'anno_batch2')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0("human_Tushar_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
    } else {
      warning("Human Finnish dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Kamath_35513515"))>0){
      tmpInd=which(tmpInd=="human_Kamath_35513515")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData=.mySplitObject(tmpData,'anno_batch2')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0("human_Kamath_35513515_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
    } else {
      warning("Human PD dataset was not identified!")
    }
    
    
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Grubman_31768052"))>0){
      tmpInd=which(tmpInd=="human_Grubman_31768052")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData$geo_accession=as.character(tmpData$geo_accession)
      tmpData=.mySplitObject(tmpData,'geo_accession')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0("human_Grubman_31768052_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
      
    } else {
      warning("Human Grubman_31768052 dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Velmeshev_31097668"))>0){
      tmpInd=which(tmpInd=="human_Velmeshev_31097668")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData=.mySplitObject(tmpData,'anno_batch2')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0("human_Velmeshev_31097668_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
    } else {
      warning("Human Velmeshev_31097668 dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Mathys_31042697"))>0){
      tmpInd=which(tmpInd=="human_Mathys_31042697")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData=.mySplitObject(tmpData,'anno_batch2')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0("human_Mathys_31042697_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
    } else {
      warning("Human Mathys_31042697 dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Lau_32989152"))>0){
      tmpInd=which(tmpInd=="human_Lau_32989152")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData=.mySplitObject(tmpData,'anno_batch2')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0("human_Lau_32989152_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
    } else {
      warning("Human Lau_32989152 dataset was not identified!")
    }
    
    tmpInd=unlist(lapply(data,function(x) x$anno_batch[1]))
    if(length(which(tmpInd=="human_Jakel_30747918"))>0){
      tmpInd=which(tmpInd=="human_Jakel_30747918")
      
      tmpData=data[[tmpInd]]
      data=data[-tmpInd]
      tmpData=.mySplitObject(tmpData,'Sample')
      
      for(ik in 1:length(tmpData)){
        tmp3=tmpData[[ik]]
        tmp3$anno_batch=paste0("human_Jakel_30747918_",names(tmpData)[ik])
        data=c(data,tmp3)
      }
    } else {
      warning("Human Jakel_30747918 dataset was not identified!")
    }
    
    
  }
  
  
  ds_sizes=unlist(lapply(data,ncol))
  total_size=sum(ds_sizes)
  total_batch_count=length(data)
  data=data[which(ds_sizes>=argList$min_ds_size)]
  ds_sizes=unlist(lapply(data,ncol))
  ds_gene_sizes=unlist(lapply(data,nrow))
  data=data[which(ds_gene_sizes>=argList$min_ds_genes)]
  
  ds_batch=unlist(lapply(data,function(x) x$ds_batch[1]))
  retained_dataset_size=unlist(lapply(data,ncol))
  df=data.frame(name=ds_batch,retained_samples=retained_dataset_size,stringsAsFactors = F)
  df=aggregate(retained_samples~name,data=df,sum)
  df$name=as.character(df$name)
  dim(df)
  
  summary_statistic=merge(summary_statistic,df,by="name",all=T)
  if(sum(is.na(summary_statistic$retained_samples))>0){
    summary_statistic$retained_samples[is.na(summary_statistic$retained_samples)]=0
  }
  if(sum(colnames(summary_statistic)=="fraction")>0){
    summary_statistic=summary_statistic[,-which(colnames(summary_statistic)=="fraction")]
  }
  
  if(sum(colnames(summary_statistic)=="included")>0){
    colnames(summary_statistic)[which(colnames(summary_statistic)=="included")]="passed_basicQC"
  }
  
  
  retained_size=sum(ds_sizes)
  retained_batch_count=length(data)
  print("***********************")
  print(paste0("Retained ",sum(summary_statistic$retained_samples)," samples out of ",sum(summary_statistic$total)," (",round(sum(summary_statistic$retained_samples)/sum(summary_statistic$total),2)*100,"%)"))
  print("***********************")
  
  data2=data
  if(!is.null(expGenePercentage)){
    gnames=unlist(lapply(data2,function(x) row.names(x)))
    gnames=as.data.frame(table(gnames))
    gnames=gnames[order(gnames$Freq,decreasing = T),]
    gnames=gnames[gnames$Freq>=expGenePercentage*length(data2),]
    gnames=as.character(gnames$gnames)
    for(idata in 1:length(data2)){
      tmp=data[[idata]]
      tmp=tmp[row.names(tmp) %in% gnames,]
      data2[[idata]]=tmp
    }
  }
  
  data_m=NULL
  if(returnDataM){
    data_m=.mycBindFn(inputList = data2,batchNames = names(data2)) #44858 347459
    
    if(convert_to_gene_level){
      
      
      load("~/serverFiles/mapping_external_gene_names.rda")
      
      fd=as.data.frame(rowData(data_m))
      fd=merge(data.frame(name=row.names(fd),fd),mapping,by="ensembl_gene_id",all.x=T)
      fd$external_gene_name=toupper(fd$external_gene_name)
      
      fd=fd[!is.na(fd$external_gene_name),]
      fd=fd[fd$external_gene_name!="",]
      fd=fd[fd$external_gene_name %in% fd$external_gene_name[duplicated(fd$external_gene_name)],]
      fd=fd[order(fd$external_gene_name),]
      
      res_m=matrix(0,nrow=length(unique(fd$external_gene_name)),ncol=ncol(data_m))
      counter=0
      nameList=list()
      for(igene in unique(fd$external_gene_name)){
        counter=counter+1
        tmp=counts(data_m)[which(rowData(data_m)$ensembl_gene_id %in% fd$ensembl_gene_id[which(fd$external_gene_name==igene)]),]
        tmp_fd=rowData(data_m)[which(rowData(data_m)$ensembl_gene_id %in% fd$ensembl_gene_id[which(fd$external_gene_name==igene)]),]
        nameList=c(nameList,list(as.data.frame(tmp_fd[1,])))
        if(rankedBased){
          tmp=as.matrix(tmp)
          tmp=apply(tmp,2,function(x) max(x))
        } else {
          tmp=Matrix::colSums(tmp)
        }
        
        res_m[counter,]=tmp
      }
      
      fd_m=do.call("rbind",nameList)
      fd$name=as.character(fd$name)
      fd_c=as.data.frame(rowData(data_m))
      fd_c=fd_c[-which(row.names(fd_c) %in% fd$name),]
      exp=counts(data_m)[-which(row.names(data_m) %in% fd$name),]
      exp=rbind(exp,as(res_m, "sparseMatrix"))
      fd_c=rbind(fd_c,fd_m)
      pd=as.data.frame(colData(data_m))
      row.names(exp)=row.names(fd_c)
      data_m=SingleCellExperiment(assays = list(counts = exp),colData = pd,rowData=fd_c)
      
      for(i in 1:length(data)){
        tmp=data[[i]]
        tmp=tmp[row.names(tmp) %in% row.names(data_m),]
        data[[i]]=tmp
      }
      
      
    }
    
    organism=rep("human",ncol(data_m))
    organism[grepl("^mouse",data_m$anno_batch)]="mouse"
    
    colData(data_m)=colData(data_m)[,-which(grepl("organism",colnames(colData(data_m))))]
    data_m$organism=organism
  }
  
  return(list(data=data,data_m=data_m,summary_statistic=summary_statistic))
}

#argList=.ArgList;ncores=3
.mycellClassMerger=function(argList,ncores=3,second_phase=F){
  #argList: any of the arglists used for the integrative analysis
  
  expGenePercentage=NULL
  
  .myIndFileReaderFn=function(dsName,dfPath,organism,convertToHuman=F,rm_artifacts=F,Leng200=F){
    require(qs)
    index=strsplit(dsName,"_")
    index=unlist(lapply(index,function(x) paste(x[-1],collapse = "_")))
    index=which(dfPath$fileName==index)
    tmp=qread(dfPath$path[index])
    tmp=tmp[!is.na(rowData(tmp)$ensembl_gene_id),]
    tmp$QC_totalCell_count=ncol(tmp)
    
    detectionRate=c()
    if(sum(colnames(rowData(tmp))=="QC_detectionRate")>0){
      detectionRate=rowData(tmp)$QC_detectionRate
      
      
    } else {
      
      for(i in seq(1,nrow(tmp),3000)){
        detectionRate=c(detectionRate,apply(as.matrix(counts(tmp)[i:min(i+3000-1,nrow(tmp)),]),1,function(x) sum(x>0)))
      }
      if(length(detectionRate)!=nrow(tmp)){
        stop("Error in counting genes expressed")
      }
    }
    
    if(T){
      if(sum(colnames(colData(tmp))=="QC_MT.pct")>0){
        if(sum(!is.na(tmp$QC_MT.pct))>10){
          if(organism=="Human"){
            tmp=tmp[,which(tmp$QC_MT.pct<10)]
          } else {
            tmp=tmp[,which(tmp$QC_MT.pct<5)]
          }
          
        }
        
      }
    }
    
    if(ncol(tmp)>10){
      if(dsName!="human_Mathys_31042697"){
        if(Leng200&dsName=="human_Leng_33432193"){
          tmp=tmp[,which(tmp$QC_Gene_unique_count>200)]
          expMedian=median(as.numeric(as.character(tmp$QC_Gene_unique_count)),na.rm = T)
          expMad=mad(as.numeric(as.character(tmp$QC_Gene_unique_count)),na.rm = T)
          tmp=tmp[,which(tmp$QC_Gene_unique_count<=(expMedian+3*expMad))]
        } else {
          tmp=tmp[,which(tmp$QC_Gene_unique_count>500)]
          expMedian=median(as.numeric(as.character(tmp$QC_Gene_unique_count)),na.rm = T)
          expMad=mad(as.numeric(as.character(tmp$QC_Gene_unique_count)),na.rm = T)
          tmp=tmp[,which(tmp$QC_Gene_unique_count<=(expMedian+3*expMad))]
        }
      } else {
        tmp=tmp[,which(tmp$QC_Gene_unique_count>200)]
        expMedian=median(as.numeric(as.character(tmp$QC_Gene_unique_count)),na.rm = T)
        expMad=mad(as.numeric(as.character(tmp$QC_Gene_unique_count)),na.rm = T)
        tmp=tmp[,which(tmp$QC_Gene_unique_count<=(expMedian+3*expMad))]
      }
      
    }
    
    
    #tmp=tmp[,which(tmp$QC_Gene_unique_count>200&tmp$QC_Gene_unique_count<6000)]
    if(sum(duplicated(rowData(tmp)$ensembl_gene_id))>0){
      o=order(detectionRate,decreasing = T)
      tmp=tmp[o,]
      ind=which(!duplicated(rowData(tmp)$ensembl_gene_id))
      tmp=tmp[ind,]
      detectionRate=detectionRate[ind]
    }
    
    
    
    row.names(tmp)[!is.na(rowData(tmp)$ensembl_gene_id)]=rowData(tmp)$ensembl_gene_id
    if(convertToHuman){
      tmp=.myMapToHuman(tmp,server=T,conservativeMapping=argList$conservativeMapping,oldMapping=argList$oldMapping,MGIsymbol=argList$MGIsymbol)
    }
    
    
    row.names(tmp)=rowData(tmp)$ensembl_gene_id
    row.names(rowData(tmp))=rowData(tmp)$ensembl_gene_id
    if(nrow(tmp)<10000){
      print(paste("less than 10k genes in the dataset:",dsName,"(genes:",nrow(tmp),")"))
    }
    if(sum(colnames(colData(tmp))=="anno_batch")>0){
      tmp$anno_batch2=tmp$anno_batch
    }
    
    tmp$anno_batch=dsName
    
    colnames(tmp)=paste0(dsName,"_",colnames(tmp))
    if(sum(colnames(colData(tmp))=="QC_Gene_total_count")>0){
      tmp$depth_per_gene=tmp$QC_Gene_total_count/(tmp$QC_Gene_unique_count+1)
    }
    
    
    tmp$QC_pass_basicQC=ncol(tmp)
    return(tmp)
  }
  
  .myHighExpRemoverFn=function(inputList){
    top50list=c()
    for(i in 1:length(inputList)){
      top50list=c(top50list,rowData(inputList[[i]])$ensembl_gene_id[rowData(inputList[[i]])$QC_top50_expressed=="Yes"])
    }
    
    top50list=as.data.frame(table(top50list))
    top50list=top50list[order(top50list$Freq,decreasing = T),]
    top50list=as.character(top50list$top50list)[top50list$Freq>=(length(inputList)/2)]
    for(i in 1:length(inputList)){
      tmp2=inputList[[i]]
      tmp2=tmp2[-which(rowData(tmp2)$ensembl_gene_id %in% top50list),]
      inputList[[i]]=tmp2
    }
    
    return(inputList)
  }
  
  if(sum(names(argList)=="Leng200")==0){
    argList$Leng200=F
  }
  
  removeHighExp = argList$excludeHighExp
  cellClass_data=list()
  for(cellClass in c("Oligo","ExN","InN","Astro","MG","Endo","OPC")){
    argList$cellClass=cellClass
    humanFiles=.myCellClass_FileList(inputClass=cellClass,inputOrganism="Human",do.exclusion = T,second_phase = second_phase)
    mouseFile=.myCellClass_FileList(inputClass=cellClass,inputOrganism="Mouse",do.exclusion = T,second_phase = second_phase)
    
    data=list()
    #dsName=paste0("human_",humanFiles$fileName)[1];dfPath=humanFiles;convertToHuman=F;organism="Human";Leng200=argList$Leng200
    tmp=parallel::mclapply(paste0("human_",humanFiles$fileName),.myIndFileReaderFn,dfPath=humanFiles,convertToHuman=F,organism="Human",Leng200=argList$Leng200,mc.cores = ncores)
    if(removeHighExp){
      tmp=.myHighExpRemoverFn(inputList = tmp)
    }
    data=c(data,tmp)
    
    tmp=parallel::mclapply(paste0("mouse_",setdiff(mouseFile$fileName,"Ximerakis_scRNA_31551601")),.myIndFileReaderFn,dfPath=mouseFile,convertToHuman=argList$includeHuman,organism="Mouse",mc.cores = ncores)
    if(removeHighExp){
      tmp=.myHighExpRemoverFn(inputList = tmp)
    }
    
    data=c(data,tmp)
    
    for(i in 1:length(data)){
      tmp=data[[i]]
      tmp$ds_batch=tmp$anno_batch
      data[[i]]=tmp
    }
    
    cellClass_data=c(cellClass_data,list(data))
    names(cellClass_data)[length(cellClass_data)]=cellClass
    
  }
  
  datasetList=c()
  for(iclass in 1:length(cellClass_data)){
    tmp=cellClass_data[[iclass]]
    names(tmp)=unlist(lapply(tmp,function(x) x$ds_batch[1]))
    datasetList=c(datasetList,names(tmp))
    cellClass_data[[iclass]]=tmp
  }
  
  if(!dir.exists("~/myBucket/SC_data/datasets")){
    dir.create("~/myBucket/SC_data/datasets",recursive = T)
  }
  
  
  for(ids in unique(datasetList)){
    if(!file.exists(paste0("~/myBucket/SC_data/datasets/",ids,".qs"))){
      tmp=list()
      for(i in 1:length(cellClass_data)){
        if(sum(names(cellClass_data[[i]])==ids)>0){
          tmp2=cellClass_data[[i]][[ids]]
          tmp=c(tmp,list(tmp2))
        }
      }
      tmp=.mycBindFn(tmp)
      qsave(tmp,file=paste0("~/myBucket/SC_data/datasets/",ids,".qs"))
    }
    
    #tmp2=.mycBindFn(tmp[1:4])
    #tmp3=.mycBindFn(tmp[6:7])
    #res_dataset=c(res_dataset,list(tmp))
    #names(res_dataset)[length(res_dataset)]=ids
  }
  
  
  #QC_totalCell_count
  #QC_pass_basicQC
  
  
  
  
  return(list(data=data,data_m=data_m,summary_statistic=summary_statistic))
}


.myExpMatrixCreatorFn=function(cellClass,second_phase=F){
  #convert_to_gene_level: if multiple ensembl gene ids were present for the same entrez id, sum them up for each sample and report one representative
  
  humanFiles=.myCellClass_FileList(inputClass=cellClass,inputOrganism="Human",do.exclusion = T,second_phase = second_phase)
  
  
  mouseFile=.myCellClass_FileList(inputClass=cellClass,inputOrganism="Mouse",do.exclusion = T,second_phase = second_phase)
  
  
  for(i in 1:nrow(humanFiles)){
    tmp_path=unlist(lapply(strsplit(humanFiles$path[i],"/"),function(x) paste(x[-length(x)],collapse = "/")))
    tmp_fileName=gsub("\\.qs$","",unlist(lapply(strsplit(humanFiles$path[i],"/"),function(x) paste(x[length(x)],collapse = "/"))))
    if(!dir.exists(paste0(tmp_path,"/expMatrix_txt"))){
      dir.create(paste0(tmp_path,"/expMatrix_txt"),recursive = T)
    }
    
    expData=qread(humanFiles$path[i])
    expData_count=summary(counts(expData))
    pd=as.data.frame(colData(expData))
    fd=as.data.frame(rowData(expData))
    write.table(expData_count,file=paste0(tmp_path,"/expMatrix_txt/",tmp_fileName,"_count.txt"),sep="\t",row.names = F,quote = F)
    write.table(pd,file=paste0(tmp_path,"/expMatrix_txt/",tmp_fileName,"_pd.txt"),sep="\t",row.names = F,quote = F)
    write.table(fd,file=paste0(tmp_path,"/expMatrix_txt/",tmp_fileName,"_fd.txt"),sep="\t",row.names = F,quote = F)
  }
  
  for(i in 1:nrow(mouseFile)){
    tmp_path=unlist(lapply(strsplit(mouseFile$path[i],"/"),function(x) paste(x[-length(x)],collapse = "/")))
    tmp_fileName=gsub("\\.qs$","",unlist(lapply(strsplit(mouseFile$path[i],"/"),function(x) paste(x[length(x)],collapse = "/"))))
    if(!dir.exists(paste0(tmp_path,"/expMatrix_txt"))){
      dir.create(paste0(tmp_path,"/expMatrix_txt"),recursive = T)
    }
    
    expData=qread(mouseFile$path[i])
    expData_count=summary(counts(expData))
    pd=as.data.frame(colData(expData))
    fd=as.data.frame(rowData(expData))
    write.table(expData_count,file=paste0(tmp_path,"/expMatrix_txt/",tmp_fileName,"_count.txt"),sep="\t",row.names = F,quote = F)
    write.table(pd,file=paste0(tmp_path,"/expMatrix_txt/",tmp_fileName,"_pd.txt"),sep="\t",row.names = F,quote = F)
    write.table(fd,file=paste0(tmp_path,"/expMatrix_txt/",tmp_fileName,"_fd.txt"),sep="\t",row.names = F,quote = F)
  }
  
  
  return("Done")
}


.myReadData_spliterFn=function(inputData,removeHighExp=F){
  
  
  if(ncol(rowData(inputData))==0){
    rowData(inputData)=data.frame(gene=row.names(inputData),gene_short_name=row.names(inputData))
  }
  
  if(sum(colnames(rowData(inputData))=="gene_short_name")==0){
    rowData(inputData)$gene_short_name=row.names(inputData)
  }
  
  rowData(inputData)$gene_short_name=toupper(rowData(inputData)$gene_short_name)
  
  .myHighExpRemoverFn=function(inputList){
    top50list=c()
    for(i in 1:length(inputList)){
      top50list=c(top50list,rowData(inputList[[i]])$ensembl_gene_id[rowData(inputList[[i]])$QC_top50_expressed=="Yes"])
    }
    
    top50list=as.data.frame(table(top50list))
    top50list=top50list[order(top50list$Freq,decreasing = T),]
    top50list=as.character(top50list$top50list)[top50list$Freq>=(length(inputList)/2)]
    for(i in 1:length(inputList)){
      tmp2=inputList[[i]]
      tmp2=tmp2[-which(rowData(tmp2)$ensembl_gene_id %in% top50list),]
      inputList[[i]]=tmp2
    }
    
    return(inputList)
  }
  
  
  tmpBatch=inputData$anno_batch
  tmpBatch[is.na(tmpBatch)]="NA"
  
  data=list()
  for(ik in unique(tmpBatch)){
    tmp3=inputData[,which(tmpBatch==ik)]
    data=c(data,list(tmp3))
    names(data)[length(data)]=ik
  }
  
  if(removeHighExp){
    data=.myHighExpRemoverFn(inputList = data)
    data_m=.mycBindFn(inputList = data,batchNames = names(data))
  } else {
    data_m=inputData
  }
  
  
  
  return(list(data=data,data_m=data_m))
}

.myHighVarGeneSlFn=function(data,argList,dataorganism,batch_variable="anno_batch",HVG_gene_count="automatic"){
  #min sample size in each dataset for that dataset to be included in the analysis
  
  
  minDatasetSize=argList$min_ds_size
  
  mySeuratFn_vst=function(inputData,dataorganism,HVG_gene_count){
    #vstGenes=.extraHVGvstFn(counts(inputData))
    if(HVG_gene_count=="automatic"){
      if(tolower(dataorganism)==tolower("Human")){
        HVG_gene_count=2500
      } else {
        HVG_gene_count=2000
      }
    }
    inputData=.extraExport2SeuratFn(inputData = inputData)
    if(ncol(inputData)>10){
      inputData<- Seurat::NormalizeData(inputData, verbose = FALSE)
      if(length(argList$input_highly_var_genes)==0){
        argList$input_highly_var_genes=NULL
      }
      if(is.null(argList$input_highly_var_genes)){
        inputData <- Seurat::FindVariableFeatures(inputData, selection.method = "vst", nfeatures = HVG_gene_count, verbose = FALSE)
      } else {
        inputData@assays$RNA@var.features=unique(as.character(argList$input_highly_var_genes))
      }
    }
    
    return(inputData)
  }
  
  myVSTFn_vst=function(inputData){
    if(ncol(inputData)>10){
      vstGenes=.extraHVGvstFn(counts(inputData))
      vstGenes=vstGenes[order(vstGenes$vst.variance.standardized,decreasing = T),]
      vstGenes=data.frame(gene=row.names(vstGenes),rank=1:nrow(vstGenes),stringsAsFactors = F)
    } else {
      vstGenes=data.frame(gene=row.names(inputData),rank=1:nrow(inputData),stringsAsFactors = F)
    }
    
    return(vstGenes)
  }
  
  mySeuratFn_sctransform=function(inputData,HVG_gene_count){
    inputData=.extraExport2SeuratFn(inputData = inputData)
    inputData <- SCTransform(inputData, vars.to.regress = NULL, verbose = FALSE,variable.features.n = HVG_gene_count)
    return(inputData)
  }
  
  
  varFeatures=c()
  allGenes=c()
  
  if(sum(colnames(colData(data$data[[1]]))=="ds_batch")==0){
    tmpBatchNames_complete=unlist(lapply(data$data,function(x) colData(x)[,batch_variable][1]))
    tmpBatchSizes=unlist(lapply(data$data,ncol))
    tmpBatchNames=strsplit(tmpBatchNames_complete,"_")
    tmpBatchNames=unlist(lapply(tmpBatchNames,function(x) if(length(x)>2){if(length(x)>3){paste(x[-((length(x)-1):length(x))],collapse = "_")} else {if(length(x)>1) {paste(x[-(length(x))],collapse = "_")} else {x}}} else {paste(x,collapse = "_")}))
    
    if(sum(grepl("human_Mancuso",tmpBatchNames))>0){
      tmpBatchNames[grepl("human_Mancuso",tmpBatchNames)]="human_Mancuso"
    }
  } else {
    tmpBatchNames_complete=unlist(lapply(data$data,function(x) colData(x)[,batch_variable][1]))
    tmpBatchNames=unlist(lapply(data$data,function(x) colData(x)$ds_batch[1]))
    tmpBatchSizes=unlist(lapply(data$data,ncol))
  }
  
  if(is.null(argList$input_highly_var_genes)){
    if(argList$HVG_method=="sctransform"){
      stop("needs to be revised!")
      data$data=parallel::mclapply(data$data,mySeuratFn_sctransform,HVG_gene_count=HVG_gene_count,mc.cores = argList$ncores) #22 datasets
      gc()
      #save(data,data_m,file="myBucket/microgliaData/combined_data_arranged.rda")
      
      for(i in 1:length(data$data)){
        varFeatures=c(varFeatures,data$data[[i]]@assays$SCT@var.features)
        allGenes=c(allGenes,row.names(data$data[[i]]))
      }
      varFeatures=table(as.character(varFeatures))
      varFeatures=as.data.frame(varFeatures)
    } else {
      vstScores=parallel::mclapply(data$data,myVSTFn_vst,mc.cores = argList$ncores)
      #data$data=parallel::mclapply(data$data[1],mySeuratFn_vst,dataorganism=dataorganism,mc.cores = argList$ncores)
      gc()
      
      varScores=NULL
      varFeatures=NULL
      allGenes=NULL
      for(i in 1:length(data$data)){
        if(ncol(data$data[[i]])>minDatasetSize){
          tmp=vstScores[[i]]
          tmp=tmp[order(tmp$rank,decreasing = F),]
          tmp=tmp[1:HVG_gene_count,]
          varScores=rbind(varScores,data.frame(dataset=as.character(tmpBatchNames_complete)[i],dataSource=as.character(tmpBatchNames)[i],varGenes=tmp$gene,varRank=((max(tmp$rank)-tmp$rank)/1000),weight=1/sum(tmpBatchNames==tmpBatchNames[i]&tmpBatchSizes>minDatasetSize),stringsAsFactors = F))
        }
        
      }
      varScores$weightedRank=varScores$varRank*varScores$weight
      varScores=aggregate(weightedRank~varGenes,data=varScores,sum)
      varScores=varScores[order(varScores$weightedRank,decreasing = F),]
      #save(data,data_m,file="myBucket/microgliaData/combined_data_arranged.rda")
      
      hvg_count_thr=HVG_gene_count
      if(hvg_count_thr=="automatic"){
        if(tolower(dataorganism)=="human"){
          hvg_count_thr=2500
        } else {
          hvg_count_thr=2000
        }
      }
      
      for(i in 1:length(data$data)){
        if(ncol(data$data[[i]])>minDatasetSize){
          varFeatures=rbind(varFeatures,data.frame(dataset=as.character(tmpBatchNames_complete)[i],dataSource=as.character(tmpBatchNames)[i],varGenes=vstScores[[i]][vstScores[[i]][,2]<=hvg_count_thr,1],weight=1/sum(sum(tmpBatchNames==tmpBatchNames[i]&tmpBatchSizes>minDatasetSize)),stringsAsFactors = F))
          allGenes=rbind(allGenes,data.frame(dataset=as.character(tmpBatchNames_complete)[i],dataSource=as.character(tmpBatchNames)[i],Genes=row.names(data$data[[i]]),weight=1/sum(sum(tmpBatchNames==tmpBatchNames[i]&tmpBatchSizes>minDatasetSize)),stringsAsFactors = F))
        }
        
      }
      
      varFeatures=aggregate(weight~varGenes,data=varFeatures,sum)
      colnames(varFeatures)=c("Var1","Freq")
      
      allGenes=aggregate(weight~Genes,data=allGenes,sum)
      colnames(allGenes)=c("allGenes","Freq")
      if(argList$includeHuman&argList$includeMouse){
        consMapping=F
        oldMapping=F
        if(sum(names(argList)=="conservativeMapping")>0){
          if(argList$conservativeMapping){
            
            load("~/serverFiles/ortholog_mapping.rda")
            mapping_old=mapping
            consMapping=T
          } else if(argList$oldMapping){
            
            load("~/serverFiles/ortholog_mapping.rda")
            mapping_old=mapping
            oldMapping=T
            
          } else if(argList$MGIsymbol){
            
            load("~/serverFiles/mouse_human_orthologs_MGISymbol_03062021.rda")
            mapping_old=mapping
            oldMapping=T
          }
          
        }
        
        
        load("~/serverFiles/mouse_human_orthologs.rda")
        if(consMapping){
          mapping=mapping[which(paste0(mapping$human,"_",mapping$mouse) %in% paste0(mapping_old$human,"_",mapping_old$mouse)),]
        } else if(oldMapping){
          mapping=mapping_old
        }
        
        
        allGenes=allGenes[which(as.character(allGenes$allGenes) %in% mapping$human),]
        
      }
    }
  } else {
    #data$data=parallel::mclapply(data$data,mySeuratFn_vst,dataorganism=dataorganism,mc.cores = argList$ncores) #22 datasets
    #if(sum(names(data)=="data_rearranged")>0){
    #  data$data_rearranged=parallel::mclapply(data$data_rearranged,mySeuratFn_vst,dataorganism=dataorganism,mc.cores = argList$ncores)
    #}
    allGenes=NULL
    for(i in 1:length(data$data)){
      if(ncol(data$data[[i]])>minDatasetSize){
        allGenes=rbind(allGenes,data.frame(dataset=as.character(tmpBatchNames_complete)[i],dataSource=as.character(tmpBatchNames)[i],Genes=row.names(data$data[[i]]),weight=1/sum(sum(tmpBatchNames==tmpBatchNames[i]&tmpBatchSizes>minDatasetSize)),stringsAsFactors = F))
      }
    }
    
    allGenes=aggregate(weight~Genes,data=allGenes,sum)
    colnames(allGenes)=c("allGenes","Freq")
    if(argList$includeHuman&argList$includeMouse){
      consMapping=F
      oldMapping=F
      if(sum(names(argList)=="conservativeMapping")>0){
        if(argList$conservativeMapping){
          
          load("~/serverFiles/ortholog_mapping.rda")
          mapping_old=mapping
          consMapping=T
        } else if(argList$oldMapping){
          
          load("~/serverFiles/ortholog_mapping.rda")
          mapping_old=mapping
          oldMapping=T
          
        } else if(argList$MGIsymbol){
          
          load("~/serverFiles/mouse_human_orthologs_MGISymbol_03062021.rda")
          mapping_old=mapping
          oldMapping=T
        }
        
      }
      
      
      load("~/serverFiles/mouse_human_orthologs.rda")
      if(consMapping){
        mapping=mapping[which(paste0(mapping$human,"_",mapping$mouse) %in% paste0(mapping_old$human,"_",mapping_old$mouse)),]
      } else if(oldMapping){
        mapping=mapping_old
      }
      
      
      allGenes=allGenes[which(as.character(allGenes$allGenes) %in% mapping$human),]
      
    }
  }
  
  
  
  allGenes$allGenes=as.character(allGenes$allGenes)
  allGenes=allGenes[which(allGenes$Freq>=(length(unique(tmpBatchNames[which(tmpBatchSizes>minDatasetSize)]))*as.numeric(as.character(argList$allGenesFraction)))),]
  
  if(is.null(argList$input_highly_var_genes)){
    varScores=varScores[as.character(varScores$varGenes) %in% as.character(allGenes$allGenes),]
    
    #varFeatures=varFeatures[order(varFeatures$Freq,decreasing = T),]
    #.varFeatures=varFeatures
    
    #varFeatures=.varFeatures
    #argList$HVG_count=9
    #varFeatures=varFeatures[which(varFeatures$Freq>(argList$HVG_count-1)),]
    varFeatures$Var1=as.character(varFeatures$Var1)
    {
      if(!is.null(argList$external_DE_path)){
        load(argList$external_DE_path)
        varFeatures=data.frame(Var1=intersect(external_de_genes,allGenes$allGenes),Freq=100)
      } else {
        varFeatures=varFeatures[which(varFeatures$Var1 %in% allGenes$allGenes),]
      }
      
    }
  } else {
    varFeatures=data.frame(Var1=intersect(argList$input_highly_var_genes,allGenes$allGenes),Freq=100)
  }
  
  #length(varFeatures)
  #sum(varFeatures %in% mapping[,1])
  #if(argList$onlyRankBased){
  #  varFeatures=as.character(varScores$varGenes[which(varScores$weightedRank<=argList$varScore.thr)])
  #} else {
  #  varFeatures=intersect(varFeatures,as.character(varScores$varGenes)[which(varScores$weightedRank<=argList$varScore.thr)])
  #}
  if(!is.null(data$data_m)){
    data$data_m=.extraExport2SeuratFn(data$data_m)
    data$data_m=NormalizeData(data$data_m)
    data$data_m@assays$RNA@var.features=intersect(varFeatures$Var1[which(varFeatures$Freq>(argList$HVG_list-1))],row.names(data$data_m))
    if(sum(colnames(data$data_m@meta.data)=="QC_Gene_total_count")>0){
      data$data_m$depth_per_gene=data$data_m$QC_Gene_total_count/(data$data_m$QC_Gene_unique_count+1)
    }
  }
  
  colnames(varFeatures)=c("Gene","Freq")
  
  return(c(data,list(varFeatures=varFeatures,allGenes=allGenes)))
}

.sconline.select_HVG_thr=function(data,nGene){
  tmp=data$varFeatures
  tmp=tmp[order(tmp$Freq,decreasing = T),]
  res=NULL
  for(icount in nGene){
    res=rbind(res,data.frame(HVG_threshold=round(tmp[icount,2]+1,3),number_of_genes=icount,stringsAsFactors = F))
  }
  return(res)
}

.myFilePathMakerFn=function(filename,argList=NULL,pdf=F,pdf2=F,conservativeMapping=NULL,oldMapping=NULL,MGIsymbol=NULL,expData=F,varGenes=F,HVG_method=NULL,HVG_list=NULL,do.split.prop=NULL,commonExpressed=NULL,covariates=NULL,indScaling=NULL,saveDir=NULL,saveDirGlobal=NULL,excludeHighExp=NULL,UMI_cor_thr=NULL,onlyRankBased=NULL,varScore.thr=NULL,allGenesFraction=NULL,uniformZscore=NULL,uniformImportant=F,dist_zscore_gamma=NULL,dist_zscore_nbinom=NULL,regularize=NULL,geoMean=NULL,returnDir=F,makePlotDir=F,DE_nPropIter=NULL,pseudocell_size=NULL,qsFormat=F,prop.n.neighbors=NULL,propImportant=F,pseudoImportant=T,pseudocell_selection_method=NULL,include.singletons=NULL,colNormalize=NULL,...){
  #conservativeMapping: use a conservative method to map the orthologs
  #oldMapping: use the previous mapping file to map the orthologs
  
  if(is.null(HVG_method)){
    HVG_method=argList$HVG_method
  }
  
  if(is.null(HVG_list)){
    HVG_list=argList$HVG_list
  }
  
  if(length(HVG_list)>1){
    stop("Only one HVG_list value should be specified!")
  }
  
  if(is.null(uniformZscore)){
    uniformZscore=argList$uniformZscore
  }
  
  
  if(is.null(commonExpressed)){
    commonExpressed=argList$commonExpressed
  }
  
  if(is.null(covariates)){
    covariates=argList$covariates
  }
  
  if(is.null(excludeHighExp)){
    excludeHighExp=argList$excludeHighExp
  }
  
  if(is.null(colNormalize)){
    colNormalize=argList$colNormalize
  }
  
  if(is.null(include.singletons)){
    include.singletons=argList$include.singletons
  }
  
  if(is.null(pseudocell_selection_method)){
    pseudocell_selection_method=argList$pseudocell_selection_method
  }
  
  if(is.null(prop.n.neighbors)){
    if(sum(names(argList)=="prop.n.neighbors")>0){
      prop.n.neighbors=argList$prop.n.neighbors
    } else {
      prop.n.neighbors=4
    }
    
  }
  
  if(is.null(UMI_cor_thr)){
    UMI_cor_thr=argList$UMI_cor_thr
  }
  
  if(is.null(conservativeMapping)){
    conservativeMapping=argList$conservativeMapping
  }
  
  if(is.null(oldMapping)){
    oldMapping=argList$oldMapping
  }
  
  if(is.null(MGIsymbol)){
    MGIsymbol=argList$MGIsymbol
  }
  
  if(is.null(allGenesFraction)){
    allGenesFraction=argList$allGenesFraction
  }
  
  if(is.null(dist_zscore_gamma)){
    dist_zscore_gamma=argList$dist_zscore_gamma
  }
  
  if(is.null(dist_zscore_nbinom)){
    dist_zscore_nbinom=argList$dist_zscore_nbinom
  }
  
  if(is.null(regularize)){
    regularize=argList$regularize
  }
  
  if(is.null(geoMean)){
    geoMean=argList$geoMean
  }
  
  if(is.null(DE_nPropIter)){
    DE_nPropIter=argList$DE_nPropIter
  }
  
  if(is.null(pseudocell_size)){
    pseudocell_size=argList$pseudocell_size
  }
  
  
  if(!is.null(covariates)){
    if(all(covariates=='3')){
      covariates=c('QC_Gene_total_count','QC_top50_pct','QC_Gene_unique_count')
    } else if(all(covariates=='2')) {
      covariates=c('QC_top50_pct','QC_Gene_unique_count')
    } else if(all(covariates=='1')){
      covariates=c('QC_Gene_unique_count')
    } else if(all(covariates=="")){
      covariates=NULL
    }
  }
  
  
  if(is.null(indScaling)){
    indScaling=argList$indScaling
  }
  
  if(is.null(saveDir)){
    saveDir=argList$saveDir
  }
  
  if(is.null(saveDirGlobal)){
    saveDirGlobal=argList$saveDirGlobal
  }
  
  if(is.null(do.split.prop)){
    do.split.prop=argList$do.split.prop
  }
  
  if(is.null(onlyRankBased)){
    onlyRankBased=argList$onlyRankBased
  }
  
  if(is.null(varScore.thr)){
    varScore.thr=argList$varScore.thr
  }
  
  if(expData){
    filename=file.path(saveDirGlobal,filename)
    
    
    filename=paste0(filename,"_UMIcor_",UMI_cor_thr)
    
    if(conservativeMapping){
      filename=paste0(filename,"_consOrth")
    }
    
    if(oldMapping){
      filename=paste0(filename,"_oldMapping")
    }
    
    if(MGIsymbol){
      filename=paste0(filename,"_MGIsymbol")
    }
    
    if(!qsFormat){
      filename=paste0(filename,".rda")
    } else{
      filename=paste0(filename,".qs")
    }
    
  } else if(varGenes){
    filename=file.path(saveDirGlobal,paste0(filename,"_",HVG_method,"_HVG",HVG_list))
    if(commonExpressed){
      filename=paste0(filename,"_commonExp")
    } else {
      filename=paste0(filename,"_noCommonExp")
    }
    if(!is.null(excludeHighExp)){
      if(excludeHighExp){
        filename=paste0(filename,"_exHighVar")
      }
    }
    
    
    if(conservativeMapping){
      filename=paste0(filename,"_consOrth")
    }
    
    if(oldMapping){
      filename=paste0(filename,"_oldMapping")
    }
    
    if(MGIsymbol){
      filename=paste0(filename,"_MGIsymbol")
    }
    
    if(!is.null(onlyRankBased)){
      if(onlyRankBased){
        filename=paste0(filename,"_onlyRankBased")
      }
    }
    
    
    filename=paste0(filename,"_varRankThr",varScore.thr)
    
    filename=paste0(filename,"_UMIcor_",UMI_cor_thr)
    
    if(allGenesFraction!=0.5){
      filename=paste0(filename,"_allGenesFraction",allGenesFraction)
    }
    
    if(!qsFormat){
      filename=paste0(filename,".rda")
    } else {
      filename=paste0(filename,".qs")
    }
    
  } else {
    filename_org=filename
    if(!is.null(uniformZscore)){
      if(uniformZscore&uniformImportant){
        filename=paste0(filename,"_uniformZscore")
      }
    }
    
    filename=paste0(filename,"_",HVG_method,"_HVG",HVG_list)
    
    if(!is.null(dist_zscore_gamma)){
      if(dist_zscore_gamma&uniformImportant){
        filename=paste0(filename,"_gamma")
      }
    }
    
    if(!is.null(dist_zscore_nbinom)){
      if(dist_zscore_nbinom&uniformImportant){
        filename=paste0(filename,"_nbinom")
      }
    }
    
    if(!is.null(prop.n.neighbors)){
      if(prop.n.neighbors!=4&propImportant){
        filename=paste0(filename,"_propN",prop.n.neighbors)
      }
    }
    
    if(!is.null(do.split.prop)){
      if(!do.split.prop&propImportant){
        filename=paste0(filename,"_noSplitPorp")
      }
    }
    
    if(!is.null(regularize)){
      if(regularize&uniformImportant){
        filename=paste0(filename,"_regularized")
        if(geoMean&(dist_zscore_nbinom|dist_zscore_gamma)){
          filename=paste0(filename,"_geoMean")
        }
      }
    }
    
    if(!is.null(DE_nPropIter)){
      if(DE_nPropIter!=1&uniformImportant){
        filename=paste0(filename,"_propIter",DE_nPropIter)
      }
    }
    
    if(!is.null(covariates)){
      filename=paste0(filename,"_cov",paste(covariates,collapse = "_"))
    }
    
    if(indScaling){
      filename=paste0(filename,"_indScaling")
    }
    
    if(commonExpressed){
      filename=paste0(filename,"_commonExp")
    } else {
      filename=paste0(filename,"_noCommonExp")
    }
    
    if(!is.null(excludeHighExp)){
      if(excludeHighExp){
        filename=paste0(filename,"_exHighVar")
      }
    }
    
    if(!is.null(onlyRankBased)){
      if(onlyRankBased){
        filename=paste0(filename,"_onlyRankBased")
      }
    }
    
    
    if(conservativeMapping){
      filename=paste0(filename,"_consOrth")
    }
    
    if(oldMapping){
      filename=paste0(filename,"_oldMapping")
    }
    
    if(MGIsymbol){
      filename=paste0(filename,"_MGIsymbol")
    }
    
    filename=paste0(filename,"_varRankThr",varScore.thr)
    
    filename=paste0(filename,"_UMIcor_",UMI_cor_thr)
    if(allGenesFraction!=0.5){
      filename=paste0(filename,"_allGenesFraction",allGenesFraction)
    }
    
    if(pseudoImportant){
      if(!is.null(pseudocell_size)){
        filename=paste0(filename,"_pseudocell",pseudocell_size)
        if(!is.null(pseudocell_selection_method)){
          if(pseudocell_selection_method!="kmeans"){
            filename=paste0(filename,"_densitypeak")
          }
        }
      }
      
    }
    
    if(propImportant){
      if(!is.null(include.singletons)){
        if(!include.singletons){
          filename=paste0(filename,"_NoSingleton")
        }
      }
      
      
      if(!is.null(colNormalize)){
        if(!colNormalize){
          filename=paste0(filename,"_NoColNorm")
        }
      }
      
    }
    
    
    if(returnDir){
      filename=gsub(paste0("^",filename_org,"_"),"",filename)
      filename=paste0(filename,"/")
    } else if(makePlotDir){
      filename=gsub(paste0("^",filename_org,"_"),"",filename)
      filename=paste0(filename,"/",filename_org)
    }
    filename=file.path(saveDir,filename)
    
    if(!returnDir){
      if(pdf|makePlotDir){
        if(!pdf2){
          filename=paste0(filename,".png")
        } else {
          filename=paste0(filename,".pdf")
        }
        
      } else if(pdf2){
        filename=paste0(filename,".pdf")
      } else {
        if(!qsFormat){
          filename=paste0(filename,".rda")
        } else {
          filename=paste0(filename,".qs")
        }
      }
    }
  }
  
  
  return(filename)
}

.myRunSettings=function(){
  .dfSettings=data.frame(commonExpressed=T,nPCs=30,HVG_count=5,covariates="depth_per_gene",indScaling=F,stringsAsFactors = F)
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=5,covariates="",indScaling=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=F,nPCs=30,HVG_count=4,covariates="",indScaling=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=T,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=F,nPCs=30,HVG_count=4,covariates="depth_per_gene",indScaling=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="depth_per_gene",indScaling=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="depth_per_gene",indScaling=T,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates=3,indScaling=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=F,nPCs=30,HVG_count=5,covariates="",indScaling=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=20,HVG_count=5,covariates="",indScaling=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=F,nPCs=25,HVG_count=5,covariates="",indScaling=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=F,nPCs=35,HVG_count=4,covariates="",indScaling=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=35,HVG_count=4,covariates="",indScaling=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=35,HVG_count=4,covariates="depth_per_gene",indScaling=F,stringsAsFactors = F))
return(.dfSettings)
}

.myRunSettingsComplete_mouse=function(){
  .dfSettings=NULL
  
  #.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=5,covariates="depth_per_gene",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
#.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=5,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
#.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=F,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
#.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=T,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
#.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=F,nPCs=30,HVG_count=4,covariates="depth_per_gene",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="depth_per_gene",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
#.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="depth_per_gene",indScaling=T,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
#.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates=3,indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates=2,indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates=1,indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
#.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=F,nPCs=30,HVG_count=4,covariates=2,indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
#.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=F,nPCs=30,HVG_count=4,covariates=1,indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=T,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=5,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
#.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=20,HVG_count=5,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
#.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=F,nPCs=25,HVG_count=5,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
#.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=F,nPCs=35,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=35,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
#.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="0,1,2,3,7,8,11,14,16,18,19,20,21,23,32,33,35",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,do.liger=F,varScore.thr=10000,allGenesFraction=0.5,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=F,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates="",indScaling=F,removeHighExp=T,slMicrogliaClusters="",breakHammond=F,UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates="",indScaling=F,removeHighExp=T,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=5000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))

.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=5500,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.3,onlyRankBased=T,varScore.thr=5000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.3,onlyRankBased=T,varScore.thr=5500,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=5000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=5500,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=T,varScore.thr=5000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=T,varScore.thr=5500,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=F,UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=5500,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.6667,do.liger=F,stringsAsFactors = F))
#.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.75,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=6000,allGenesFraction=0.6667,do.liger=F,stringsAsFactors = F))
#.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=6000,allGenesFraction=0.75,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=8000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=5000,allGenesFraction=0.75,do.liger=F,stringsAsFactors = F))

.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=8,covariates="depth_per_gene",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=5000,allGenesFraction=0.7,do.liger=T,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=8,covariates=1,indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=5000,allGenesFraction=0.7,do.liger=T,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=8,covariates=2,indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=5000,allGenesFraction=0.7,do.liger=T,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=8,covariates=3,indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=5000,allGenesFraction=0.7,do.liger=T,stringsAsFactors = F))
.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=8,covariates="depth_per_gene",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=5000,allGenesFraction=0.7,do.liger=F,stringsAsFactors = F))
#.dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=35,HVG_count=4,covariates="depth_per_gene",indScaling=F,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,stringsAsFactors = F))
return(.dfSettings)
}

.myRunSettingsComplete_human=function(){
  .dfSettings=NULL
  
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates="depth_per_gene",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates=2,indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates=1,indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates="",indScaling=F,removeHighExp=T,slMicrogliaClusters="",UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=35,HVG_count=3,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=2,covariates="",indScaling=F,removeHighExp=T,slMicrogliaClusters="",UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=5000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=5500,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.3,onlyRankBased=T,varScore.thr=5000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.3,onlyRankBased=T,varScore.thr=5500,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=5000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=5500,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.2,onlyRankBased=T,varScore.thr=5000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.2,onlyRankBased=T,varScore.thr=5500,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.6667,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=6000,allGenesFraction=0.6667,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=30,HVG_count=3,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=8000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=20,HVG_count=3,covariates="",indScaling=F,removeHighExp=T,slMicrogliaClusters="",UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=20,HVG_count=4,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=25,HVG_count=3,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.2,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=20,HVG_count=3,covariates="",indScaling=F,removeHighExp=F,slMicrogliaClusters="",UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  .dfSettings=rbind(.dfSettings,data.frame(commonExpressed=T,nPCs=20,HVG_count=2,covariates="",indScaling=F,removeHighExp=T,slMicrogliaClusters="",UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=10000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  return(.dfSettings)
}

.mySaveDirMaker=function(prefix,nPCs,cov){
  res=prefix
  
  if(!is.null(cov)){
    res=paste0(res,"-cov-",paste0(cov,collapse = "-"))
  }
  res=paste0(res,"-nPCs",nPCs)
  
  return(res)
}

.mySaveFn=function(file,...){
  
  object=list(...)[[1]]
  objectName=lapply(substitute(list(...))[-1], deparse)[[1]]
  repCycle=T
  repCount=0
  
  assign(objectName,object)
  while(repCycle&repCount<4){
    repCycle=F
    repCount=repCount+1
    repCycle=tryCatch({save(list=objectName,file=file);F}, error=function(e) {return(T)})
    if(repCycle){
      Sys.sleep(30)
    }
  }
  
  if(repCycle){
    stop(file)
  }
  return("done")
}


.myArgFn=function(runIndx,nPCs=NULL,saveDir=NULL,prop.n.neighbors = 4,ncores=7,conservativeMapping=F,oldMapping=F,MGIsymbol=F,includeHuman,includeMouse,DE_supportingFractionThr,DE_n.adaptiveKernel,DE_nPropIter,uniformZscore=F,Leng200=F,dist_zscore_gamma,dist_zscore_norm,dist_zscore_nbinom,regularize,geoMean=F,prefix=NULL,newRun=F,inputDf=NULL,internal_pseudocell_count,external_DE_path=NULL,do.split.prop=T,external_DE_name=NULL){
  
  if(is.null(inputDf)){
    if(includeMouse&(!includeHuman)){
      .dfSetting=.myRunSettingsComplete_mouse()
    } else if((!includeMouse)&(includeHuman)){
      .dfSetting=.myRunSettingsComplete_human()
    } else {
      .dfSetting=.myRunSettingsComplete_mouse()
    }
  } else {
    .dfSetting=inputDf
  }
  
  
  .commonExpressed=.dfSetting$commonExpressed[runIndx]
  .includeHuman=includeHuman
  .includeMouse=includeMouse
  if(is.null(nPCs)){
    .nPCs=.dfSetting$nPCs[runIndx]
  } else {
    .nPCs=nPCs
  }
  
  .HVG_list=.dfSetting$HVG_count[runIndx] #4
  .HVG_method="vst" #"vst" and "sctransform"
  .covariates=.dfSetting$covariates[runIndx]
  .removeHighExp=.dfSetting$removeHighExp[runIndx]
  .UMI_cor_thr=.dfSetting$UMI_cor_thr[runIndx]
  
  
  if(!is.null(.covariates)){
    if(all(.covariates=='3')){
      .covariates=c('QC_Gene_total_count','QC_top50_pct','QC_Gene_unique_count')
    } else if(all(.covariates=='2')) {
      .covariates=c('QC_top50_pct','QC_Gene_unique_count')
    } else if(all(.covariates=='1')){
      .covariates=c('QC_Gene_unique_count')
    } else if(all(.covariates=="")){
      .covariates=NULL
    }
  }
  
  
  .exCellStates=c("CAM","Mnc","BAM","T cells","Neutrophils","cDC","B cells","NK cells","Proliferation","T/NKT cells","cDC2","Monocyte/Mdc","ILC","ydT cells","pDC","cDC1","migDC","MCs","Non-cl. monocytes")
  
  #"depth_per_gene"#'QC_Gene_unique_count'#c('QC_Gene_total_count','QC_top50_pct',)#, "organism","depth_per_gene" 'QC_Gene_total_count','QC_top50_pct','QC_Gene_unique_count'
  .indScaling=.dfSetting$indScaling[runIndx]
  
  if(is.null(prefix)){
    prefix=""
  }
  
  if(sum(c(.includeHuman,.includeMouse))==2){
    prefix=paste0(prefix,"human-mouse")
  } else if(.includeMouse){
    prefix=paste0(prefix,"mouse")
  } else if(.includeHuman){
    prefix=paste0(prefix,"human")
  }
  
  if(.dfSetting$do.liger[runIndx]){
    prefix=paste0(prefix,"_liger")
  }
  
  if(!is.null(external_DE_name)){
    prefix=paste0(prefix,"_",external_DE_name)
  }
  
  if(Leng200){
    prefix=paste0(prefix,"_Leng200")
  }
  
  if(is.null(saveDir)){
    .saveDir=.mySaveDirMaker(paste0("~/data/data/tmpBucket/results/",prefix),nPCs = .nPCs,cov=.covariates) 
  } else {
    .saveDir=.mySaveDirMaker(file.path(saveDir,prefix),nPCs = .nPCs,cov=.covariates) 
  }
  
  .saveDirGlobal=paste0("~/data/data/tmpBucket/results/",prefix,"-Global")
  
  
  argList=list(commonExpressed=.dfSetting$commonExpressed[runIndx],do.split.prop=do.split.prop,prop.n.neighbors = prop.n.neighbors,includeHuman=.includeHuman,includeMouse=.includeMouse,conservativeMapping=conservativeMapping,oldMapping=oldMapping,Leng200=Leng200,MGIsymbol=MGIsymbol,nPCs=.dfSetting$nPCs[runIndx],HVG_list=.dfSetting$HVG_count[runIndx],HVG_method="vst",covariates=.covariates,indScaling=.dfSetting$indScaling[runIndx],saveDir=.saveDir,saveDirGlobal=.saveDirGlobal,ncores=ncores,excludeHighExp=.removeHighExp,UMI_cor_thr=.UMI_cor_thr,onlyRankBased=.dfSetting$onlyRankBased[runIndx],varScore.thr=.dfSetting$varScore.thr[runIndx],do.liger=.dfSetting$do.liger[runIndx],allGenesFraction=.dfSetting$allGenesFraction[runIndx],uniformZscore=uniformZscore,DE_supportingFractionThr=DE_supportingFractionThr,DE_n.adaptiveKernel=DE_n.adaptiveKernel,DE_nPropIter=DE_nPropIter,dist_zscore_gamma=dist_zscore_gamma,dist_zscore_norm=dist_zscore_norm,dist_zscore_nbinom=dist_zscore_nbinom,regularize=regularize,geoMean=geoMean,prefix=prefix,internal_pseudocell_count=internal_pseudocell_count,external_DE_path=external_DE_path,external_DE_name=external_DE_name)
  
  .extraCreateDirFn(argList)
  
  .extraNewRun(argList=argList,newRun=newRun)
  
  return(argList)
  
}

.myGeneMeanSdFn=function(inputExpList,argList){
  myGeneMeanSdFn=function(tmpData){
    
    tmpCount=tmpData@assays$RNA@counts
    gm_count=apply(tmpCount,1,mean)
    gsd_count=apply(tmpCount,1,sd)
    geo_mean_count=.extra_gmean(tmpCount)
    
    
    tmpdata=tmpData@assays$RNA@data
    gm_lognorm=apply(tmpdata,1,mean)
    gsd_lognorm=apply(tmpdata,1,sd)
    geo_mean_lognorm=.extra_gmean(tmpdata)
    
    return(data.frame(gene=row.names(tmpData),meanExp_count=gm_count,sdExp_count=gsd_count,geo_mean_count=geo_mean_count,meanExp_lognorm=gm_lognorm,sdExp_lognorm=gsd_lognorm,geo_mean_lognorm=geo_mean_lognorm,stringsAsFactors = F))
  }
  resGeneMeanSd=parallel::mclapply(inputExpList,myGeneMeanSdFn,mc.cores = argList$ncores)
  resGeneMeanSd=data.table::rbindlist(resGeneMeanSd,use.names = T,idcol = T)
  colnames(resGeneMeanSd)[colnames(resGeneMeanSd)==".id"]="anno_batch"
  resGeneMeanSd=as.data.frame(resGeneMeanSd)
  return(resGeneMeanSd)
}

.myConcensusDEFn_step1=function(argList,expData=NULL,batch_variable="anno_batch"){
  
  load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  if(is.null(expData)){
    load(.myFilePathMakerFn("exp_merged",argList=argList,expData=T))
    expData=SplitObject(tmp, split.by = batch_variable)
  } else {
    names(expData)=unlist(lapply(expData,function(x) as.character(x[1,batch_variable])))
  }
  
  if(!file.exists(.myFilePathMakerFn("kmeans_clustering",argList=argList))){
    load(file=.myFilePathMakerFn("kmeans_res_clusters",argList=argList))
    
    load(.myFilePathMakerFn("pca_centroids",argList=argList))
    
    if(sum(colnames(pd)=="cluster_id")>0){
      
      pd$cluster_id_data=pd[,"cluster_id"]
      pd=pd[,-which(colnames(pd)=="cluster_id")]
    }
    
    res_clusters=merge(res_clusters,data.frame(sample=row.names(pd),pd,stringsAsFactors = F),by="sample")
    
    
    
    cluster_weights=as.data.frame(table(res_clusters$cluster_id))
    cluster_weights$weight=cluster_weights$Freq/sum(cluster_weights$Freq)
    clust_w_batch=list()
    for(i in unique(res_clusters[,batch_variable])){
      tmp_cl_weight=res_clusters[which(res_clusters[,batch_variable]==i),]
      tmp_cl_weight=as.data.frame(table(tmp_cl_weight$cluster_id))
      tmp_cl_weight$weight=tmp_cl_weight$Freq/sum(tmp_cl_weight$Freq)
      tmp_w=cluster_weights
      tmp_w=merge(tmp_w,tmp_cl_weight,by="Var1")
      tmp_w$weight.x=tmp_w$weight.x/sum(tmp_w$weight.x)
      
      tmp_w$weight=tmp_w$weight.x/tmp_w$weight.y
      tmp_w$Freq.new=tmp_w$Freq.y*tmp_w$weight/sum(tmp_w$Freq.y*tmp_w$weight)*sum(tmp_w$Freq.y)
      tmp_w$weight=tmp_w$Freq.new/tmp_w$Freq.y
      tmp_w[,batch_variable]=i
      tmp_w=tmp_w[,c("Var1",batch_variable,"weight")]
      colnames(tmp_w)[1]="cluster_id"
      colnames(tmp_w)[3]="cluster_weight"
      clust_w_batch=c(clust_w_batch,list(tmp_w))
    }
    clust_w_batch=do.call("rbind",clust_w_batch)
    
    if(length(which(clust_w_batch$cluster_weight>quantile(clust_w_batch$cluster_weight,0.95)))>0){
      clust_w_batch$cluster_weight[which(clust_w_batch$cluster_weight>quantile(clust_w_batch$cluster_weight,0.95))]=quantile(clust_w_batch$cluster_weight,0.95)
    }
    
    if(length(which(clust_w_batch$cluster_weight<quantile(clust_w_batch$cluster_weight,0.05)))>0){
      clust_w_batch$cluster_weight[which(clust_w_batch$cluster_weight<quantile(clust_w_batch$cluster_weight,0.05))]=quantile(clust_w_batch$cluster_weight,0.05)
    }
    
    res_clusters=merge(res_clusters,clust_w_batch,by=c("cluster_id",batch_variable))
    gc()
    save(res_clusters,file=.myFilePathMakerFn("kmeans_clustering",argList=argList))
  }
  
  
  if(!file.exists(.myFilePathMakerFn("geneExp_meanSd",argList=argList))){
    
    resGeneMeanSd=.myGeneMeanSdFn(inputExpList=expData,argList=argList)
    save(resGeneMeanSd,file=.myFilePathMakerFn("geneExp_meanSd",argList=argList))
    gc()
  }
  
  if(!file.exists(.myFilePathMakerFn("geneExp_meanSd_weighted",argList=argList))){
    
    myWeightedFn=function(inputExpData,weights){
      tmpCount=inputExpData@assays$RNA@counts
      weights=weights[match(colnames(tmpCount),weights$sample),]
      res_sd_count=matrixStats::rowWeightedSds(tmpCount, w = weights$cluster_weight)
      res_m_count=apply(tmpCount,1,function(x) weighted.mean(x,w=weights$cluster_weight))
      res_geo_m_count=.extra_gmean(tmpCount, eps = 1,weights=weights$cluster_weight)
      rm(tmpCount)
      tmpCount=inputExpData@assays$RNA@data
      weights=weights[match(colnames(tmpCount),weights$sample),]
      res_sd_lognorm=matrixStats::rowWeightedSds(tmpCount, w = weights$cluster_weight)
      res_m_lognorm=apply(tmpCount,1,function(x) weighted.mean(x,w=weights$cluster_weight))
      res_geo_m_lognorm=.extra_gmean(tmpCount, eps = 1,weights=weights$cluster_weight)
      
      res=data.frame(gene=row.names(inputExpData),w_sd_count=res_sd_count,w_geomean_count=res_geo_m_count,w_mean_count=res_m_count,w_sd_lognorm=res_sd_lognorm,w_geomean_lognorm=res_geo_m_lognorm,w_mean_lognorm=res_m_lognorm)
      return(res)
    }
    
    load(.myFilePathMakerFn("kmeans_clustering",argList=argList))
    res_w=parallel::mclapply(expData,myWeightedFn,weights=res_clusters,mc.cores = 3)
    res_w=do.call("rbind",res_w)
    
    ##for gamma, counts need to be 10k normalized in WeightedFn
    #res_w_1$gamma_rate = res_w_1$count_w_mean/(res_w_1$gamma_w_sd^2)
    #res_w_1$gamma_shape = (res_w_1$count_w_mean^2)/(res_w_1$gamma_w_sd^2)
    #res_w_1$gamma_rate[res_w_1$gamma_w_sd==0]=0
    #res_w_1$gamma_shape[res_w_1$gamma_w_sd==0]=0
    
    tmpbatch=unlist(lapply(strsplit(row.names(res_w),"\\."),function(x) x[1]))
    res_w$anno_batch=tmpbatch
    rm(tmpbatch)
    
    resGeneMeanSd=res_w
    save(resGeneMeanSd,file=.myFilePathMakerFn("geneExp_meanSd_weighted",argList=argList))
    
  }
  
  gc()
  return("Done")
}

.myCentroidGeneMeanSdFn=function(){
  
  load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  #load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
  
  if(!file.exists(.myFilePathMakerFn("centroids_geneExp_meanSd_weighted",argList=argList))|!file.exists(.myFilePathMakerFn("centroids_geneExp",argList=argList))){
    load(.myFilePathMakerFn("exp_merged",argList=argList,expData=T))
    expData=SplitObject(tmp, split.by = "anno_batch")
    
  }
  
  
  
  if(!file.exists(.myFilePathMakerFn("kmeans_clustering",argList=argList))){
    load(file=.myFilePathMakerFn("kmeans_res_clusters",argList=argList))
    
    load(.myFilePathMakerFn("pca_centroids",argList=argList))
    
    
    
    res_clusters=merge(res_clusters,data.frame(sample=row.names(pd),pd,stringsAsFactors = F),by="sample")
    
    
    
    cluster_weights=as.data.frame(table(res_clusters$cluster_id))
    cluster_weights$weight=cluster_weights$Freq/sum(cluster_weights$Freq)
    clust_w_batch=list()
    for(i in unique(res_clusters$anno_batch)){
      tmp_cl_weight=res_clusters[which(res_clusters$anno_batch==i),]
      tmp_cl_weight=as.data.frame(table(tmp_cl_weight$cluster_id))
      tmp_cl_weight$weight=tmp_cl_weight$Freq/sum(tmp_cl_weight$Freq)
      tmp_w=cluster_weights
      tmp_w=merge(tmp_w,tmp_cl_weight,by="Var1")
      tmp_w$weight.x=tmp_w$weight.x/sum(tmp_w$weight.x)
      
      tmp_w$weight=tmp_w$weight.x/tmp_w$weight.y
      tmp_w$Freq.new=tmp_w$Freq.y*tmp_w$weight/sum(tmp_w$Freq.y*tmp_w$weight)*sum(tmp_w$Freq.y)
      tmp_w$weight=tmp_w$Freq.new/tmp_w$Freq.y
      tmp_w$anno_batch=i
      tmp_w=tmp_w[,c("Var1","anno_batch","weight")]
      colnames(tmp_w)[1]="cluster_id"
      colnames(tmp_w)[3]="cluster_weight"
      clust_w_batch=c(clust_w_batch,list(tmp_w))
    }
    clust_w_batch=do.call("rbind",clust_w_batch)
    
    if(length(which(clust_w_batch$cluster_weight>quantile(clust_w_batch$cluster_weight,0.95)))>0){
      clust_w_batch$cluster_weight[which(clust_w_batch$cluster_weight>quantile(clust_w_batch$cluster_weight,0.95))]=quantile(clust_w_batch$cluster_weight,0.95)
    }
    
    if(length(which(clust_w_batch$cluster_weight<quantile(clust_w_batch$cluster_weight,0.05)))>0){
      clust_w_batch$cluster_weight[which(clust_w_batch$cluster_weight<quantile(clust_w_batch$cluster_weight,0.05))]=quantile(clust_w_batch$cluster_weight,0.05)
    }
    
    res_clusters=merge(res_clusters,clust_w_batch,by=c("cluster_id","anno_batch"))
    
    save(res_clusters,file=.myFilePathMakerFn("kmeans_clustering",argList=argList))
  } else {
    load(.myFilePathMakerFn("kmeans_clustering",argList=argList))
    load(.myFilePathMakerFn("pca_centroids",argList=argList))
  }
  
  if(!file.exists(.myFilePathMakerFn("centroids_geneExp",argList=argList))){
    myCentroidExpFn=function(ids,expData,res_clusters){
      tmpData=expData[[ids]]@assays$RNA@counts
      tmpData=apply(tmpData,2,function(x) x/sum(x)*10000)
      tmpData=t(tmpData)
      tmpClust=res_clusters[match(row.names(tmpData),res_clusters$sample),]
      tmpData=split(as.data.frame(tmpData),tmpClust$cluster_id)
      sumFn=function(x){
        x=apply(x,2,sum)
        x=x*10000/sum(x)
        return(x)
      }
      tmpData=lapply(tmpData,sumFn)
      tmpData2=t(do.call("rbind",tmpData))
      return(tmpData2)
    }
    centroidExpData=parallel::mclapply(1:length(expData),myCentroidExpFn,expData=expData,res_clusters=res_clusters,mc.cores = 4)
    names(centroidExpData)=names(expData)
    save(centroidExpData,file=.myFilePathMakerFn("centroids_geneExp",argList=argList))
  } else {
    load(.myFilePathMakerFn("centroids_geneExp",argList=argList))
  }
  
  if(!file.exists(.myFilePathMakerFn("centroids_geneExp_meanSd_weighted",argList=argList))){
    
    myGammaWeightedFn=function(inputExpData,weights){
      inputExpData=as.matrix(inputExpData@assays$RNA@counts)
      inputExpData=apply(inputExpData,2,function(x) x/sum(x)*10000)
      weights=weights[match(colnames(inputExpData),weights$sample),]
      res_sd=matrixStats::rowWeightedSds(inputExpData, w = weights$cluster_weight)
      res_m=apply(inputExpData,1,function(x) weighted.mean(x,w=weights$cluster_weight))
      res_geo_m=.extra_gmean(inputExpData, eps = 1,weights=weights$cluster_weight)
      res=data.frame(gene=row.names(inputExpData),gamma_w_sd=res_sd,count_w_geomean=res_geo_m,count_w_mean=res_m)
      return(res)
    }
    
    
    res_w=parallel::mclapply(expData,myGammaWeightedFn,weights=res_clusters,mc.cores = 3)
    res_w_1=do.call("rbind",res_w)
    
    res_w_1$gamma_rate = res_w_1$count_w_mean/(res_w_1$gamma_w_sd^2)
    res_w_1$gamma_shape = (res_w_1$count_w_mean^2)/(res_w_1$gamma_w_sd^2)
    res_w_1$gamma_rate[res_w_1$gamma_w_sd==0]=0
    res_w_1$gamma_shape[res_w_1$gamma_w_sd==0]=0
    
    
    {
      myWeightedFn=function(inputExpData,weights){
        weights=weights[match(colnames(inputExpData),weights$sample),]
        res_sd=matrixStats::rowWeightedSds(as.matrix(inputExpData@assays$RNA@data), w = weights$cluster_weight)
        res_m=apply(as.matrix(inputExpData@assays$RNA@data),1,function(x) weighted.mean(x,w=weights$cluster_weight))
        res=data.frame(gene=row.names(inputExpData),normalDist_w_mean=res_m,normalDist_w_sd=res_sd)
        return(res)
      }
      res_w=parallel::mclapply(expData,myWeightedFn,weights=res_clusters,mc.cores = 3)
      res_w_2=do.call("rbind",res_w)
      
    }
    
    res_w=res_w_1[match(row.names(res_w_2),row.names(res_w_1)),]
    if(sum(is.na(res_w$gene))>0){
      stop("Error in matching names in step1")
    } else {
      res_w=cbind(res_w,res_w_2[,-1])
    }
    tmpbatch=unlist(lapply(strsplit(row.names(res_w),"\\."),function(x) x[1]))
    res_w$anno_batch=tmpbatch
    rm(tmpbatch)
    save(res_w,file=.myFilePathMakerFn("centroids_geneExp_meanSd_weighted",argList=argList))
    
  }
  
  
  if(!file.exists(.myFilePathMakerFn("centroids_geneExp_meanSd",argList=argList))){
    myGeneMeanSdFn=function(tmpData){
      
      gm=apply(tmpData,1,mean)
      gsd=apply(tmpData,1,sd)
      gm_norm=apply(tmpData,1,function(x) mean(log(x+1)))
      gsd_norm=apply(tmpData,1,function(x) sd(log(x+1)))
      
      geo_mean=.extra_gmean(tmpData)
      
      return(data.frame(gene=row.names(tmpData),meanExp_count=gm,sdExp_count=gsd,geo_mean_count=geo_mean,meanExp_lognorm=gm_norm,sdExp_lognorm=gsd_norm,stringsAsFactors = F))
    }
    
    resGeneMeanSd=parallel::mclapply(centroidExpData,myGeneMeanSdFn,mc.cores = argList$ncores)
    resGeneMeanSd=do.call("rbind",resGeneMeanSd)
    
    tmpBatchNames=strsplit(row.names(resGeneMeanSd),"\\.")
    tmpBatchNames=unlist(lapply(tmpBatchNames,function(x) x[1]))
    resGeneMeanSd$anno_batch=tmpBatchNames
    
    save(resGeneMeanSd,file=.myFilePathMakerFn("centroids_geneExp_meanSd",argList=argList))
    
  }
  return("Done")
}

.myParRegFn=function (object,model_pars, genes_count_mean,genes,doGeneLogTransform, bw_adjust= 3){
  
  #doGeneLogTransform: should the program log10 transforms the genes_count_mean
  object=as.data.frame(object)
  genes <- as.character(object[,genes])
  
  if(doGeneLogTransform){
    genes_log_mean=log10(object[,genes_count_mean])
  } else {
    genes_log_mean=object[,genes_count_mean]
  }
  
  if(length(which(object[,genes_count_mean]==0))>0){
    genes_log_mean[which(object[,genes_count_mean]==0)]=0
  }
  
  genes_log_mean_org=genes_log_mean
  object_org=object
  if(length(model_pars)==1){
    outliers <- sctransform:::is_outlier(object[,model_pars],genes_log_mean)
  } else {
    outliers <- apply(object[,model_pars], 2, function(y) sctransform:::is_outlier(y,genes_log_mean))
  }
  
  if(length(ncol(outliers))>0){
    outliers <- apply(outliers, 1, any)
  }
  
  if (sum(outliers) > 0) {
    object=object[!outliers, ]
    genes_log_mean=genes_log_mean[!outliers]
  }
  
  outliers=(genes_log_mean==0)
  if (sum(outliers) > 0) {
    object=object[!outliers, ]
    genes_log_mean=genes_log_mean[!outliers]
  }
  
  if(length(model_pars)==1){
    outliers <- abs(object[,model_pars])==Inf
  } else {
    outliers <- apply(object[,model_pars], 2, function(y) abs(y)==Inf)
  }
  
  if(length(ncol(outliers))>0){
    outliers <- apply(outliers, 1, any)
  }
  
  if (sum(outliers) > 0) {
    object=object[!outliers, ]
    genes_log_mean=genes_log_mean[!outliers]
  }
  
  bw <- bw.SJ(genes_log_mean) * bw_adjust
  x_points <- genes_log_mean_org
  o <- order(x_points)
  model_pars_fit <- matrix(NA_real_, nrow=length(genes), ncol=length(model_pars), 
                           dimnames = list(genes, model_pars))
  
  if(length(model_pars)==1){
    model_pars_fit[o, 1] <- ksmooth(x = genes_log_mean, 
                                    y = object[,model_pars], x.points = x_points, bandwidth = bw, kernel = "normal")$y
  } else {
    for (i in model_pars) {
      model_pars_fit[o, i] <- ksmooth(x = genes_log_mean, 
                                      y = object[,i], x.points = x_points, bandwidth = bw, kernel = "normal")$y
    }
  }
  
  
  model_pars_fit=as.data.frame(model_pars_fit)
  model_pars_fit$count_log_mean=genes_log_mean_org
  tmpInd=which(is.na(model_pars_fit$theta))
  if(length(tmpInd)>0){
    for(ik in tmpInd){
      tmp=model_pars_fit[-tmpInd,]
      tmpInd2=which(abs(tmp$count_log_mean-model_pars_fit$count_log_mean[ik])==min(abs(tmp$count_log_mean-model_pars_fit$count_log_mean[ik])))
      if(abs(log10(tmp$count_log_mean[tmpInd2]/model_pars_fit$count_log_mean[ik]))<0.04){
        for(ij in 1:ncol(model_pars_fit)){
          model_pars_fit[ik,ij]=tmp[tmpInd2,ij]
        }
      }
    }
  }
  
  if(length(which(genes_log_mean_org==0))>0){
    for(i in 1:ncol(model_pars_fit)){
      model_pars_fit[which(genes_log_mean_org==0),i]=0
    }
  }
  
  return(model_pars_fit)
}

.myConcensusDEFn_step2_detail_prop_org_v1=function(data,centroidPCAdata,argList,n.adaptiveKernel=20,nPropIter=3){
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=3
  batchPCAdata=as.matrix(data$pcaData)[,1:argList$nPCs]
  dsName=data$dsName
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  
  nn.ranked.1 <- RANN::nn2(data = batchPCAdata,query = centroidPCAdata, k = min(2000,nrow(batchPCAdata)), eps = 0)
  nn.ranked.1$nn.dists=cbind(rep(0,nrow(nn.ranked.1$nn.dists)),nn.ranked.1$nn.dists)
  nn.ranked.1$nn.idx=cbind(1:nrow(nn.ranked.1$nn.idx)+nrow(batchPCAdata),nn.ranked.1$nn.idx)
  nn.ranked.1$nn.dists=nn.ranked.1$nn.dists[,1:min(2000,nrow(batchPCAdata))]
  nn.ranked.1$nn.idx=nn.ranked.1$nn.idx[,1:min(2000,nrow(batchPCAdata))]
  
  nn.centroid.purityIndex=RANN::nn2(data = centroidPCAdata, k = n.adaptiveKernel, eps = 0)
  nn.centroid.purityIndex=nn.centroid.purityIndex$nn.dists[,2]
  
  nn.dataset <- RANN::nn2(data = batchPCAdata,query = batchPCAdata, k = min(2000,nrow(batchPCAdata)), eps = 0)
  nn.dataset.purityIndx=nn.dataset$nn.dists[,n.adaptiveKernel]
  
  
  dists=nn.ranked.1$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
  affinities2=affinities-pmax(0.001,affinitiesThr)
  affinities[affinities2<=0]=0
  #removing the outlier matches
  #toRM= which(scale(affinities[,2])<(-3))
  toRM= which(scale(dists[,2])>(3))
  if(length(toRM)>0){
    for(i in 2:ncol(affinities)){
      affinities[toRM,i]=0
    }
  }
  
  nn.ranked.1$affinities=affinities
  rm(affinities)
  
  dists=nn.dataset$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/nn.dataset.purityIndx[x])^2))
  affinities=do.call("rbind",affinities)
  affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
  affinities2=affinities-pmax(0.001,affinitiesThr)
  affinities[affinities2<=0]=0
  nn.dataset$affinities=affinities
  rm(affinities)
  
  nn.ranked.org=nn.dataset
  nn.ranked.org$nn.idx=rbind(nn.dataset$nn.idx,nn.ranked.1$nn.idx)
  nn.ranked.org$nn.dists=rbind(nn.dataset$nn.dists,nn.ranked.1$nn.dists)
  nn.ranked.org$affinities=rbind(nn.dataset$affinities,nn.ranked.1$affinities)
  rm(nn.ranked.1,nn.dataset)
  
  
  nn.ranked=nn.ranked.org
  affinities=nn.ranked.org$affinities
  
  affCounts=apply(affinities,2,function(x) sum(x==0))
  affinities=affinities[,-which(affCounts==nrow(affinities))]
  nn.ranked$nn.idx=nn.ranked$nn.idx[,-which(affCounts==nrow(affinities))]
  j <- as.numeric(t(nn.ranked$nn.idx))
  i <- ((1:length(j)) - 1)%/%ncol(affinities) + 1
  k=as.numeric(t(affinities))
  
  graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(batchPCAdata)+nrow(centroidPCAdata), nrow(batchPCAdata)+nrow(centroidPCAdata)))
  rm(i,j,k)
  graph <- .extra_matrix_rowNorm(graph)#Matrix::Diagonal(x = 1 / (rowSums(graph))) %*% graph
  
  if(class(graph)!="dgCMatrix"){
    graph=as(graph,"dgCMatrix")
  }
  
  rownames(graph) <- c(row.names(batchPCAdata),row.names(centroidPCAdata))
  colnames(graph) <- c(row.names(batchPCAdata),row.names(centroidPCAdata))
  
  for(i in 1:nPropIter){
    #print(i)
    if(i==1){
      prop_mat=graph
    } else {
      prop_mat=prop_mat%*%graph
    }
  }
  
  
  return(list(data=data,prop_mat=prop_mat))
}

.myConcensusDEFn_step2_detail_prop_org_v2=function(data,centroidPCAdata,argList,n.adaptiveKernel=20,nPropIter=3,exCentroids=NULL){
  
  batchPCAdata=as.matrix(data$pcaData)[,1:argList$nPCs]
  dsName=data$dsName
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  if(nrow(batchPCAdata)<(3*n.adaptiveKernel+1)){
    n.adaptiveKernel=2
  }
  
  
  nn.ranked.1 <- RANN::nn2(data = batchPCAdata,query = centroidPCAdata, k = min(2000,nrow(batchPCAdata)), eps = 0)
  nn.ranked.1$nn.dists=cbind(rep(0,nrow(nn.ranked.1$nn.dists)),nn.ranked.1$nn.dists)
  nn.ranked.1$nn.idx=cbind(1:nrow(nn.ranked.1$nn.idx)+nrow(batchPCAdata),nn.ranked.1$nn.idx)
  nn.ranked.1$nn.dists=nn.ranked.1$nn.dists[,1:min(2000,nrow(batchPCAdata))]
  nn.ranked.1$nn.idx=nn.ranked.1$nn.idx[,1:min(2000,nrow(batchPCAdata))]
  
  load(.myFilePathMakerFn("kmeans_res_clusters",argList=argList))
  cluster_ids=row.names(centroidPCAdata)
  res_clusters=res_clusters[as.character(res_clusters$sample) %in% colnames(data$logNormData),]
  res_clusters=res_clusters[which(res_clusters$cluster %in% row.names(centroidPCAdata)),]
  res_clusters=rbind(res_clusters,data.frame(cluster_id=unique(cluster_ids),sample=unique(cluster_ids),stringsAsFactors = F))
  res_clusters$weight=1
  res_clusters=reshape2::dcast(cluster_id~sample,data=res_clusters,value.var = "weight")
  res_clusters[is.na(res_clusters)]=0
  row.names(res_clusters)=res_clusters[,1]
  res_clusters=res_clusters[,-1]
  
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)+nrow(centroidPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata),row.names(centroidPCAdata))
  for(i in seq_len(nrow(centroidPCAdata))) {
    adj[i,nn.ranked.1$nn.idx[i,]] = nn.ranked.1$nn.dists[i,]
  }
  
  adj=adj[match(row.names(res_clusters),row.names(adj)),match(colnames(res_clusters),colnames(adj))]
  if(!all(row.names(adj)==row.names(res_clusters))){
    stop("Error!")
  }
  if(!all(colnames(adj)==colnames(res_clusters))){
    stop("Error!")
  }
  tst=adj*res_clusters
  tst=apply(tst,1,function(x) median(x[which(x!=0)]))
  tst=tst[match(row.names(centroidPCAdata),names(tst))]
  
  nn.centroid.purityIndex=RANN::nn2(data = centroidPCAdata, k = 2, eps = 0)
  nn.centroid.purityIndex=nn.centroid.purityIndex$nn.dists[,2]
  
  tst=(tst+nn.centroid.purityIndex)/2
  tst[is.na(tst)]=nn.centroid.purityIndex[is.na(tst)]
  #nn.centroid.purityIndex=tst
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  nn.dataset <- RANN::nn2(batchPCAdata, k = min(2000,nrow(batchPCAdata)), eps = 0)
  nn.dataset.purityIndx=nn.dataset$nn.dists[,n.adaptiveKernel]
  
  
  dists=nn.ranked.1$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  if(ncol(affinities)>(3*n.adaptiveKernel+1)){
    affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
    affinities2=affinities-affinitiesThr
    affinities[affinities2<=0]=0
  }
  
  ##removing the outlier matches
  #toRM= which(scale(affinities[,2])<(-3))
  #if(length(toRM)>0){
  #  for(i in 1:ncol(affinities)){
  #    affinities[toRM,i]=0
  #  }
  #}
  
  nn.ranked.1$affinities=affinities
  rm(affinities)
  
  dists=nn.dataset$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/nn.dataset.purityIndx[x])^2))
  affinities=do.call("rbind",affinities)
  if(ncol(affinities)>(3*n.adaptiveKernel+1)){
    affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
    affinities2=affinities-affinitiesThr
    affinities[affinities2<=0]=0
  }
  nn.dataset$affinities=affinities
  rm(affinities)
  
  nn.ranked.org=nn.dataset
  nn.ranked.org$nn.idx=rbind(nn.dataset$nn.idx,nn.ranked.1$nn.idx)
  nn.ranked.org$nn.dists=rbind(nn.dataset$nn.dists,nn.ranked.1$nn.dists)
  nn.ranked.org$affinities=rbind(nn.dataset$affinities,nn.ranked.1$affinities)
  rm(nn.ranked.1,nn.dataset)
  
  nn.ranked=nn.ranked.org
  affinities=nn.ranked$affinities
  
  affCounts=apply(affinities,2,function(x) sum(x==0))
  if(length(which(affCounts==nrow(affinities)))>0){
    affinities=affinities[,-which(affCounts==nrow(affinities))]
    nn.ranked$nn.idx=nn.ranked$nn.idx[,-which(affCounts==nrow(affinities))]
  }
  
  
  j <- as.numeric(t(nn.ranked$nn.idx))
  i <- ((1:length(j)) - 1)%/%ncol(affinities) + 1
  k=as.numeric(t(affinities))
  
  graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(batchPCAdata)+nrow(centroidPCAdata), nrow(batchPCAdata)+nrow(centroidPCAdata)))
  rm(i,j,k)
  graph <- .extra_matrix_rowNorm(graph)#Matrix::Diagonal(x = 1 / (rowSums(graph))) %*% graph
  
  if(class(graph)!="dgCMatrix"){
    graph=as(graph,"dgCMatrix")
  }
  
  rownames(graph) <- c(row.names(batchPCAdata),row.names(centroidPCAdata))
  colnames(graph) <- c(row.names(batchPCAdata),row.names(centroidPCAdata))
  
  checkProp=T
  if(nPropIter<1){
    nPropIter=1
  }
  for(i in 1:nPropIter){
    if(checkProp){
      if(i==1){
        prop_mat=graph
      } else {
        prop_mat2=prop_mat%*%graph
        
        prop_check=unlist(lapply((ncol(data$logNormData)+1):nrow(prop_mat2),function(x) length(intersect(nn.ranked.org$nn.idx[x,1:50],order(prop_mat2[x,],decreasing = T)[1:50]))))
        if(median(prop_check)>=40&quantile(prop_check,0.1)>24){
          prop_mat=prop_mat2
        } else {
          print(paste("iteration:",i-1))
          checkProp=F
        }
      }
    }
    
  }
  
  
  res_clusters=res_clusters[match(colnames(prop_mat)[(ncol(data$logNormData)+1):nrow(prop_mat)],row.names(res_clusters)),]
  centroid_weights=1/rowSums(res_clusters)
  
  #prop_weights=1-centroid_weights
  #centroid_weights2=(1-(diag(as.matrix(graph))[(ncol(data$logNormData)+1):nrow(prop_mat)]))
  #centroid_weights1=1/(1+exp((12-2/(centroid_weights))/2))
  #centroid_weights=centroid_weights2*centroid_weights1
  #centroid_weightsN=centroid_weights/(prop_weights+0.01)
  #prop_mat2=sweep(prop_mat[(ncol(data$logNormData)+1):nrow(prop_mat),],1,centroid_weightsN,"*")
  #prop_mat2=rbind(prop_mat[1:ncol(data$logNormData),],prop_mat2)
  #prop_mat2[cbind((ncol(data$logNormData)+1):nrow(prop_mat),(ncol(data$logNormData)+1):nrow(prop_mat))]=(1- centroid_weights)
  
  centroid_weightsN=1/(1+exp((12-2/(centroid_weights))/2))
  prop_mat2=sweep(prop_mat[(ncol(data$logNormData)+1):nrow(prop_mat),],1,centroid_weightsN,"*")
  centroid_weightsN=1-rowSums(prop_mat2[,1:ncol(data$logNormData)])
  prop_mat2[cbind(1:nrow(prop_mat2),ncol(data$logNormData)+1:nrow(prop_mat2))]=centroid_weightsN
  prop_mat2=rbind(prop_mat[1:ncol(data$logNormData),],prop_mat2)
  
  final_check=rowSums(prop_mat2)
  if(sum(prop_mat2>1.001)>0){
    stop(paste0("Error in the construction of the propagation matrix for ",dsName," dataset!"))
  }
  
  return(list(data=data,prop_mat=prop_mat2))
}

.myConcensusDEFn_step2_detail_prop_DNN=function(data,centroidPCAdata,argList,n.adaptiveKernel=20,nPropIter=3,exCentroids=NULL,adjPurityIndx=F){
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=3;exCentroids=NULL;adjPurityIndx=F
  batchPCAdata=as.matrix(data$pcaData)[,1:argList$nPCs]
  dsName=data$dsName
  #row.names(batchPCAdata)=paste0(dsName,"_",row.names(batchPCAdata))
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  if(nrow(batchPCAdata)<(3*n.adaptiveKernel+1)){
    n.adaptiveKernel=2
  }
  
  
  nn.ranked.1 <- RANN::nn2(data = batchPCAdata,query = centroidPCAdata, k = min(2000,nrow(batchPCAdata)), eps = 0)
  nn.ranked.1$nn.dists=cbind(rep(0,nrow(nn.ranked.1$nn.dists)),nn.ranked.1$nn.dists)
  nn.ranked.1$nn.idx=cbind(1:nrow(nn.ranked.1$nn.idx)+nrow(batchPCAdata),nn.ranked.1$nn.idx)
  nn.ranked.1$nn.dists=nn.ranked.1$nn.dists[,1:min(2000,nrow(batchPCAdata))]
  nn.ranked.1$nn.idx=nn.ranked.1$nn.idx[,1:min(2000,nrow(batchPCAdata))]
  
  load(.myFilePathMakerFn("kmeans_res_clusters",argList=argList))
  cluster_ids=row.names(centroidPCAdata)
  res_clusters=res_clusters[as.character(res_clusters$sample) %in% colnames(data$logNormData),]
  res_clusters=res_clusters[which(res_clusters$cluster %in% row.names(centroidPCAdata)),]
  res_clusters=rbind(res_clusters,data.frame(cluster_id=unique(cluster_ids),sample=unique(cluster_ids),stringsAsFactors = F))
  res_clusters$weight=1
  res_clusters=reshape2::dcast(cluster_id~sample,data=res_clusters,value.var = "weight")
  res_clusters[is.na(res_clusters)]=0
  row.names(res_clusters)=res_clusters[,1]
  res_clusters=res_clusters[,-1]
  
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)+nrow(centroidPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata),row.names(centroidPCAdata))
  for(i in seq_len(nrow(centroidPCAdata))) {
    adj[i,nn.ranked.1$nn.idx[i,]] = nn.ranked.1$nn.dists[i,]
  }
  
  adj=adj[match(row.names(res_clusters),row.names(adj)),match(colnames(res_clusters),colnames(adj))]
  if(!all(row.names(adj)==row.names(res_clusters))){
    stop("Error!")
  }
  if(!all(colnames(adj)==colnames(res_clusters))){
    stop("Error!")
  }
  tst=adj*res_clusters
  tst_zero=apply(tst,1,function(x) max(x)==0)
  tst=apply(tst,1,function(x) median(x[which(x!=0)]))
  tst[tst_zero]=0
  tst=tst[match(row.names(centroidPCAdata),names(tst))]
  
  nn.centroid.purityIndex=RANN::nn2(data = centroidPCAdata, k = 2, eps = 0)
  nn.centroid.purityIndex=nn.centroid.purityIndex$nn.dists[,2]
  
  tst=(tst+nn.centroid.purityIndex)/2
  #tst[is.na(tst)]=nn.centroid.purityIndex[is.na(tst)]
  if(adjPurityIndx){
    nn.centroid.purityIndex=tst
  }
  
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  nn.dataset <- RANN::nn2(batchPCAdata, k = min(2000,nrow(batchPCAdata)), eps = 0)
  nn.dataset.purityIndx=RANN::nn2(data = centroidPCAdata,query = batchPCAdata, k = 1, eps = 0)
  nn.dataset.purityIndx=as.numeric(nn.dataset.purityIndx$nn.dist[,1])
  
  
  dists=nn.ranked.1$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  if(ncol(affinities)>(3*n.adaptiveKernel+1)){
    affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
    affinities2=affinities-affinitiesThr
    affinities[affinities2<=0]=0
  }
  
  ##removing the outlier matches
  #toRM= which(scale(affinities[,2])<(-3))
  #if(length(toRM)>0){
  #  for(i in 1:ncol(affinities)){
  #    affinities[toRM,i]=0
  #  }
  #}
  
  nn.ranked.1$affinities=affinities
  rm(affinities)
  
  dists=nn.dataset$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/(nn.dataset.purityIndx[x]+0.01))^2))
  affinities=do.call("rbind",affinities)
  if(ncol(affinities)>(3*n.adaptiveKernel+1)){
    affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
    affinities2=affinities-affinitiesThr
    affinities[affinities2<=0]=0
  }
  nn.dataset$affinities=affinities
  rm(affinities)
  
  nn.ranked.org=nn.dataset
  nn.ranked.org$nn.idx=rbind(nn.dataset$nn.idx,nn.ranked.1$nn.idx)
  nn.ranked.org$nn.dists=rbind(nn.dataset$nn.dists,nn.ranked.1$nn.dists)
  nn.ranked.org$affinities=rbind(nn.dataset$affinities,nn.ranked.1$affinities)
  rm(nn.ranked.1,nn.dataset)
  
  nn.ranked=nn.ranked.org
  affinities=nn.ranked$affinities
  
  affCounts=apply(affinities,2,function(x) sum(x==0))
  if(length(which(affCounts==nrow(affinities)))>0){
    affinities=affinities[,-which(affCounts==nrow(affinities))]
    nn.ranked$nn.idx=nn.ranked$nn.idx[,-which(affCounts==nrow(affinities))]
  }
  
  
  j <- as.numeric(t(nn.ranked$nn.idx))
  i <- ((1:length(j)) - 1)%/%ncol(affinities) + 1
  k=as.numeric(t(affinities))
  
  graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(batchPCAdata)+nrow(centroidPCAdata), nrow(batchPCAdata)+nrow(centroidPCAdata)))
  rm(i,j,k)
  graph <- .extra_matrix_rowNorm(graph)#Matrix::Diagonal(x = 1 / (rowSums(graph))) %*% graph
  
  if(class(graph)!="dgCMatrix"){
    graph=as(graph,"dgCMatrix")
  }
  
  rownames(graph) <- c(row.names(batchPCAdata),row.names(centroidPCAdata))
  colnames(graph) <- c(row.names(batchPCAdata),row.names(centroidPCAdata))
  
  checkProp=T
  if(nPropIter<1){
    nPropIter=1
  }
  
  #for(nPropIter in c(1,2,4,8,16,32)){
    init_weights=rowSums(graph[(nrow(batchPCAdata)+1):ncol(graph),1:nrow(batchPCAdata)])
    for(i in 1:nPropIter){
      if(i==1){
        prop_mat=graph
        prop_mat3=prop_mat
      } else {
        prop_mat2=prop_mat%*%(graph)
        prop_mat3=prop_mat3+prop_mat2
        prop_mat=prop_mat2
      }
    }
    
    prop_mat=prop_mat3/nPropIter
   
  #}
  #save(prop_routine,prop_DNN,file="~/prop_methods.rda")
  
  
  #save(prop_mat,file="~/prop_mat4.rda")
  
  prop_mat=prop_mat[row.names(prop_mat) %in% row.names(res_clusters),]
  res_clusters=res_clusters[match(row.names(prop_mat),row.names(res_clusters)),]
  centroid_weights=1/rowSums(res_clusters)
  
  centroid_weightsN=1/(1+exp((12-2/(centroid_weights))/2))
  prop_mat2=sweep(prop_mat,1,centroid_weightsN,"*")
  
  final_check=rowSums(prop_mat2)
  if(sum(final_check>1.001)>0){
    stop(paste0("Error in the construction of the propagation matrix for ",dsName," dataset!"))
  }
  
  return(list(data=data,prop_mat=prop_mat2))
}

.myConcensusDEFn_step2_detail_prop=function(data,centroidPCAdata,argList,n.adaptiveKernel=20,nPropIter=3,exCentroids=NULL,runIndx=6){
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=3;exCentroids=NULL
  batchPCAdata=as.matrix(data$pcaData)[,1:argList$nPCs]
  dsName=data$dsName
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  if(nrow(batchPCAdata)<(3*n.adaptiveKernel+1)){
    n.adaptiveKernel=2
  }
  
  
  nn.ranked.1 <- RANN::nn2(data = batchPCAdata,query = centroidPCAdata, k = min(2000,nrow(batchPCAdata)), eps = 0)
  nn.ranked.1$nn.dists=cbind(rep(0,nrow(nn.ranked.1$nn.dists)),nn.ranked.1$nn.dists)
  nn.ranked.1$nn.idx=cbind(1:nrow(nn.ranked.1$nn.idx)+nrow(batchPCAdata),nn.ranked.1$nn.idx)
  nn.ranked.1$nn.dists=nn.ranked.1$nn.dists[,1:min(2000,nrow(batchPCAdata))]
  nn.ranked.1$nn.idx=nn.ranked.1$nn.idx[,1:min(2000,nrow(batchPCAdata))]
  
  load(.myFilePathMakerFn("kmeans_res_clusters",argList=argList))
  cluster_ids=row.names(centroidPCAdata)
  res_clusters=res_clusters[as.character(res_clusters$sample) %in% colnames(data$logNormData),]
  res_clusters=res_clusters[which(res_clusters$cluster %in% row.names(centroidPCAdata)),]
  res_clusters=rbind(res_clusters,data.frame(cluster_id=unique(cluster_ids),sample=unique(cluster_ids),stringsAsFactors = F))
  res_clusters$weight=1
  res_clusters=reshape2::dcast(cluster_id~sample,data=res_clusters,value.var = "weight")
  res_clusters[is.na(res_clusters)]=0
  row.names(res_clusters)=res_clusters[,1]
  res_clusters=res_clusters[,-1]
  
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)+nrow(centroidPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata),row.names(centroidPCAdata))
  for(i in seq_len(nrow(centroidPCAdata))) {
    adj[i,nn.ranked.1$nn.idx[i,]] = nn.ranked.1$nn.dists[i,]
  }
  
  adj=adj[match(row.names(res_clusters),row.names(adj)),match(colnames(res_clusters),colnames(adj))]
  if(!all(row.names(adj)==row.names(res_clusters))){
    stop("Error!")
  }
  if(!all(colnames(adj)==colnames(res_clusters))){
    stop("Error!")
  }
  tst=adj*res_clusters
  tst=apply(tst,1,function(x) median(x[which(x!=0)]))
  tst=tst[match(row.names(centroidPCAdata),names(tst))]
  
  nn.centroid.purityIndex=RANN::nn2(data = centroidPCAdata, k = 2, eps = 0)
  nn.centroid.purityIndex=nn.centroid.purityIndex$nn.dists[,2]
  
  tst=(tst+nn.centroid.purityIndex)/2
  tst[is.na(tst)]=nn.centroid.purityIndex[is.na(tst)]
  #nn.centroid.purityIndex=tst
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  pd=pd[match(row.names(batchPCAdata),row.names(pd)),]
  
  nn.dataset <- RANN::nn2(batchPCAdata, k = min(2000,nrow(batchPCAdata)), eps = 0)
  nn.dataset.purityIndx=RANN::nn2(data = centroidPCAdata,query = batchPCAdata, k = 1, eps = 0)
  nn.dataset.purityIndx=as.numeric(nn.dataset.purityIndx$nn.dist[,1])
  
  
  dists=nn.ranked.1$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  if(ncol(affinities)>(3*n.adaptiveKernel+1)){
    affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
    affinities2=affinities-affinitiesThr
    affinities[affinities2<=0]=0
  }
  
  ##removing the outlier matches
  #toRM= which(scale(affinities[,2])<(-3))
  #if(length(toRM)>0){
  #  for(i in 1:ncol(affinities)){
  #    affinities[toRM,i]=0
  #  }
  #}
  
  nn.ranked.1$affinities=affinities
  rm(affinities)
  
  dists=nn.dataset$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/(nn.dataset.purityIndx[x]+0.01))^2))
  affinities=do.call("rbind",affinities)
  if(ncol(affinities)>(3*n.adaptiveKernel+1)){
    affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
    affinities2=affinities-affinitiesThr
    affinities[affinities2<=0]=0
  }
  nn.dataset$affinities=affinities
  rm(affinities)
  
  nn.ranked.org=nn.dataset
  nn.ranked.org$nn.idx=rbind(nn.dataset$nn.idx,nn.ranked.1$nn.idx)
  nn.ranked.org$nn.dists=rbind(nn.dataset$nn.dists,nn.ranked.1$nn.dists)
  nn.ranked.org$affinities=rbind(nn.dataset$affinities,nn.ranked.1$affinities)
  rm(nn.ranked.1,nn.dataset)
  
  nn.ranked=nn.ranked.org
  affinities=nn.ranked$affinities
  
  affCounts=apply(affinities,2,function(x) sum(x==0))
  if(length(which(affCounts==nrow(affinities)))>0){
    affinities=affinities[,-which(affCounts==nrow(affinities))]
    nn.ranked$nn.idx=nn.ranked$nn.idx[,-which(affCounts==nrow(affinities))]
  }
  
  
  j <- as.numeric(t(nn.ranked$nn.idx))
  i <- ((1:length(j)) - 1)%/%ncol(affinities) + 1
  k=as.numeric(t(affinities))
  
  graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(batchPCAdata)+nrow(centroidPCAdata), nrow(batchPCAdata)+nrow(centroidPCAdata)))
  rm(i,j,k)
  graph <- .extra_matrix_rowNorm(graph)#Matrix::Diagonal(x = 1 / (rowSums(graph))) %*% graph
  
  if(class(graph)!="dgCMatrix"){
    graph=as(graph,"dgCMatrix")
  }
  
  rownames(graph) <- c(row.names(batchPCAdata),row.names(centroidPCAdata))
  colnames(graph) <- c(row.names(batchPCAdata),row.names(centroidPCAdata))
  
  checkProp=T
  if(nPropIter<1){
    nPropIter=1
  }
  
  #nPropIter=4
  if(runIndx==1){
    for(i in 1:nPropIter){
      if(checkProp){
        if(i==1){
          prop_mat=graph
        } else {
          prop_mat2=prop_mat%*%graph
          
          prop_check=unlist(lapply((ncol(data$logNormData)+1):nrow(prop_mat2),function(x) length(intersect(nn.ranked.org$nn.idx[x,1:50],order(prop_mat2[x,],decreasing = T)[1:50]))))
          if(median(prop_check)>=40&quantile(prop_check,0.1)>24){
            prop_mat=prop_mat2
          } else {
            print(paste("iteration:",i-1))
            checkProp=F
          }
        }
      }
      
    }
    prop_mat=prop_mat[(ncol(data$logNormData)+1):nrow(prop_mat),1:ncol(data$logNormData)]
    
    sl_thr=0
    sl_prop_mat=NULL
    for(i in c(seq(0.01,0.15,0.01),seq(0.17,0.5,0.03))){
      if(sl_thr==0){
        prop_mat2=t(apply(prop_mat,1,function(x) {x[which(x/max(x)<i)]=0;x}))
        check_score=median(apply(prop_mat2,2,function(x) sum(x>0)))
        if(check_score>1.5&check_score<2.5){
          sl_thr=i
          sl_prop_mat=prop_mat2
        }
      }
    }
    if(sl_thr==0){
      sl_thr=i
      sl_prop_mat=t(apply(prop_mat,1,function(x) {x[which(x/max(x)<sl_thr)]=0;x}))
    }
    
    sl_prop_mat=t(apply(sl_prop_mat,1,function(x) x/sum(x)))
    sl_prop_mat=sl_prop_mat*rowSums(prop_mat)
    row.names(sl_prop_mat)=row.names(prop_mat)
    colnames(sl_prop_mat)=colnames(prop_mat)
    prop_mat=sl_prop_mat
    rm(sl_prop_mat)
    
  }
  
  if(runIndx==8){
    for(i in 1:nPropIter){
      if(checkProp){
        if(i==1){
          prop_mat=graph
        } else {
          prop_mat2=prop_mat%*%graph
        }
      }
    }
    prop_mat=prop_mat[(ncol(data$logNormData)+1):nrow(prop_mat),1:ncol(data$logNormData)]
    
    sl_thr=0
    sl_prop_mat=NULL
    for(i in c(seq(0.01,0.15,0.01),seq(0.17,0.5,0.03))){
      if(sl_thr==0){
        prop_mat2=t(apply(prop_mat,1,function(x) {x[which(x/max(x)<i)]=0;x}))
        check_score=median(apply(prop_mat2,2,function(x) sum(x>0)))
        if(check_score>1.5&check_score<2.5){
          sl_thr=i
          sl_prop_mat=prop_mat2
        }
      }
    }
    if(sl_thr==0){
      sl_thr=i
      sl_prop_mat=t(apply(prop_mat,1,function(x) {x[which(x/max(x)<sl_thr)]=0;x}))
    }
    
    sl_prop_mat=t(apply(sl_prop_mat,1,function(x) x/sum(x)))
    sl_prop_mat=sl_prop_mat*rowSums(prop_mat)
    row.names(sl_prop_mat)=row.names(prop_mat)
    colnames(sl_prop_mat)=colnames(prop_mat)
    prop_mat=sl_prop_mat
    rm(sl_prop_mat)
    
  }
  
  if(runIndx==8.1){
    for(i in 1:nPropIter){
      if(checkProp){
        if(i==1){
          prop_mat=graph
        } else {
          prop_mat2=prop_mat%*%graph
        }
      }
    }
    prop_mat=prop_mat[(ncol(data$logNormData)+1):nrow(prop_mat),1:ncol(data$logNormData)]
    
    sl_thr=0
    sl_prop_mat=NULL
    for(i in c(seq(0.01,0.15,0.01),seq(0.17,0.5,0.03))){
      if(sl_thr==0){
        prop_mat2=t(apply(prop_mat,1,function(x) {x[which(x/max(x)<i)]=0;x}))
        check_score=median(apply(prop_mat2,2,function(x) sum(x>0)))
        if(check_score>2.5&check_score<3.5){
          sl_thr=i
          sl_prop_mat=prop_mat2
        }
      }
    }
    if(sl_thr==0){
      sl_thr=i
      sl_prop_mat=t(apply(prop_mat,1,function(x) {x[which(x/max(x)<sl_thr)]=0;x}))
    }
    
    sl_prop_mat=t(apply(sl_prop_mat,1,function(x) x/sum(x)))
    sl_prop_mat=sl_prop_mat*rowSums(prop_mat)
    row.names(sl_prop_mat)=row.names(prop_mat)
    colnames(sl_prop_mat)=colnames(prop_mat)
    prop_mat=sl_prop_mat
    rm(sl_prop_mat)
    
  }
  
  if(runIndx==8.2){
    for(i in 1:nPropIter){
      if(checkProp){
        if(i==1){
          prop_mat=graph
        } else {
          prop_mat2=prop_mat%*%graph
        }
      }
    }
    prop_mat=prop_mat[(ncol(data$logNormData)+1):nrow(prop_mat),1:ncol(data$logNormData)]
    
    sl_thr=0
    sl_prop_mat=NULL
    for(i in c(seq(0.01,0.15,0.01),seq(0.17,0.5,0.03))){
      if(sl_thr==0){
        prop_mat2=t(apply(prop_mat,1,function(x) {x[which(x/max(x)<i)]=0;x}))
        check_score=median(apply(prop_mat2,2,function(x) sum(x>0)))
        if(check_score>3.5&check_score<4.5){
          sl_thr=i
          sl_prop_mat=prop_mat2
        }
      }
    }
    if(sl_thr==0){
      sl_thr=i
      sl_prop_mat=t(apply(prop_mat,1,function(x) {x[which(x/max(x)<sl_thr)]=0;x}))
    }
    
    sl_prop_mat=t(apply(sl_prop_mat,1,function(x) x/sum(x)))
    sl_prop_mat=sl_prop_mat*rowSums(prop_mat)
    row.names(sl_prop_mat)=row.names(prop_mat)
    colnames(sl_prop_mat)=colnames(prop_mat)
    prop_mat=sl_prop_mat
    rm(sl_prop_mat)
    
  }
  
  if(runIndx==9){
    if(nPropIter>1){
      for(i in nPropIter:2){
        if(i==nPropIter){
          prop_mat=graph
        } else {
          prop_mat2=graph %*% prop_mat
          
          prop_mat=(0.5*graph+0.5*prop_mat2)
        }
      }
      prop_mat=(0.5)*graph%*%prop_mat+(0.5)*graph
      prop_mat=prop_mat[(nrow(batchPCAdata)+1):ncol(graph),1:nrow(batchPCAdata)]
    } else {
      prop_mat=graph[(nrow(batchPCAdata)+1):ncol(graph),1:nrow(batchPCAdata)]
    }
    
    sl_thr=0
    sl_prop_mat=NULL
    for(i in c(seq(0.01,0.15,0.01),seq(0.17,0.5,0.03))){
      if(sl_thr==0){
        prop_mat2=t(apply(prop_mat,1,function(x) {x[which(x/max(x)<i)]=0;x}))
        check_score=median(apply(prop_mat2,2,function(x) sum(x>0)))
        if(check_score>1.5&check_score<2.5){
          sl_thr=i
          sl_prop_mat=prop_mat2
        }
      }
    }
    if(sl_thr==0){
      sl_thr=i
      sl_prop_mat=t(apply(prop_mat,1,function(x) {x[which(x/max(x)<sl_thr)]=0;x}))
    }
    
    sl_prop_mat=t(apply(sl_prop_mat,1,function(x) x/sum(x)))
    sl_prop_mat=sl_prop_mat*rowSums(prop_mat)
    row.names(sl_prop_mat)=row.names(prop_mat)
    colnames(sl_prop_mat)=colnames(prop_mat)
    prop_mat=sl_prop_mat
    rm(sl_prop_mat)
    
  }
  
  if(runIndx==7){
    for(i in 1:nPropIter){
      if(checkProp){
        if(i==1){
          prop_mat=graph
        } else {
          prop_mat2=prop_mat%*%graph
          
          prop_check=unlist(lapply((ncol(data$logNormData)+1):nrow(prop_mat2),function(x) length(intersect(nn.ranked.org$nn.idx[x,1:50],order(prop_mat2[x,],decreasing = T)[1:50]))))
          if(median(prop_check)>=40&quantile(prop_check,0.1)>24){
            prop_mat=prop_mat2
          } else {
            print(paste("iteration:",i-1))
            checkProp=F
          }
        }
      }
      
    }
    prop_mat=prop_mat[(ncol(data$logNormData)+1):nrow(prop_mat),1:ncol(data$logNormData)]
    
    sl_thr=0
    sl_prop_mat=NULL
    for(i in c(seq(0.01,0.15,0.01),seq(0.17,0.5,0.03))){
      if(sl_thr==0){
        prop_mat2=t(apply(prop_mat,1,function(x) {x[which(x/max(x)<i)]=0;x}))
        check_score=median(apply(prop_mat2,2,function(x) sum(x>0)))
        if(check_score>1.5&check_score<2.5){
          sl_thr=i
          sl_prop_mat=prop_mat2
        }
      }
    }
    if(sl_thr==0){
      sl_thr=i
      sl_prop_mat=t(apply(prop_mat,1,function(x) {x[which(x/max(x)<sl_thr)]=0;x}))
    }
    
    sl_prop_mat=t(apply(sl_prop_mat,1,function(x) x/sum(x)))
    sl_prop_mat=sl_prop_mat*rowSums(prop_mat)
    row.names(sl_prop_mat)=row.names(prop_mat)
    colnames(sl_prop_mat)=colnames(prop_mat)
    #prop_mat=sl_prop_mat
    rm(sl_prop_mat)
    
  }
  
  if(runIndx==2){
    for(i in 1:nPropIter){
      if(i==1){
        prop_mat=graph[(nrow(batchPCAdata)+1):ncol(graph),1:nrow(batchPCAdata)]
      } else {
        prop_mat2=prop_mat%*%graph[1:nrow(batchPCAdata),1:nrow(batchPCAdata)]
        
        prop_mat=(prop_mat+prop_mat2)/2
      }
    }
    
  }
  
  if(runIndx==3){
    {
      if(nPropIter>1){
        for(i in nPropIter:2){
          if(i==nPropIter){
            prop_mat=graph[1:nrow(batchPCAdata),1:nrow(batchPCAdata)]
          } else {
            prop_mat2=graph[1:nrow(batchPCAdata),1:nrow(batchPCAdata)] %*% prop_mat
            
            prop_mat=(0.5*graph[1:nrow(batchPCAdata),1:nrow(batchPCAdata)]+0.5*prop_mat2)
          }
        }
        prop_mat=(0.5)*graph[(nrow(batchPCAdata)+1):ncol(graph),1:nrow(batchPCAdata)]%*%prop_mat+(0.5)*graph[(nrow(batchPCAdata)+1):ncol(graph),1:nrow(batchPCAdata)]
      } else {
        prop_mat=graph[(nrow(batchPCAdata)+1):ncol(graph),1:nrow(batchPCAdata)]
      }
      
    }
  }
  
  if(runIndx==4){
    if(nPropIter>1){
      for(i in nPropIter:2){
        if(i==nPropIter){
          prop_mat=graph[1:nrow(batchPCAdata),1:nrow(batchPCAdata)]
        } else {
          prop_mat2=graph[1:nrow(batchPCAdata),1:nrow(batchPCAdata)] %*% prop_mat
        
          prop_mat=(0.5*graph[1:nrow(batchPCAdata),1:nrow(batchPCAdata)]+0.5*prop_mat2)
        }
      }
      prop_mat=(0.5)*graph[(nrow(batchPCAdata)+1):ncol(graph),1:nrow(batchPCAdata)]%*%prop_mat+(0.5)*graph[(nrow(batchPCAdata)+1):ncol(graph),1:nrow(batchPCAdata)]
    } else {
      prop_mat=graph[(nrow(batchPCAdata)+1):ncol(graph),1:nrow(batchPCAdata)]
    }
  
  
  }
  
  #17
  if(runIndx==5){
    if(nPropIter>1){
      for(i in nPropIter:2){
        if(i==nPropIter){
          prop_mat=graph
        } else {
          prop_mat2=graph %*% prop_mat
          
          prop_mat=(0.5*graph+0.5*prop_mat2)
        }
      }
      prop_mat=(0.5)*graph%*%prop_mat+(0.5)*graph
      prop_mat=prop_mat[(nrow(batchPCAdata)+1):ncol(graph),1:nrow(batchPCAdata)]
    } else {
      prop_mat=graph[(nrow(batchPCAdata)+1):ncol(graph),1:nrow(batchPCAdata)]
    }
    
  }
  
  #18
  if(runIndx==6){
    if(nPropIter>1){
      for(i in nPropIter:1){
        if(i==nPropIter){
          prop_mat=graph
        } else {
          prop_mat2=graph %*% prop_mat
          
          prop_mat=(prop_mat2+prop_mat)/2#(0.5*graph+0.5*prop_mat2)
        }
      }
      #prop_mat=(0.5)*graph%*%prop_mat+(0.5)*graph
      prop_mat=prop_mat[(nrow(batchPCAdata)+1):ncol(graph),1:nrow(batchPCAdata)]
      
    } else {
      prop_mat=graph[(nrow(batchPCAdata)+1):ncol(graph),1:nrow(batchPCAdata)]
    }
    
    
    
  }
  
  
  res_clusters=res_clusters[match(row.names(prop_mat),row.names(res_clusters)),]
  centroid_weights=1/rowSums(res_clusters)
  
  centroid_weightsN=1/(1+exp((12-2/(centroid_weights))/2))
  prop_mat2=sweep(prop_mat,1,centroid_weightsN,"*")
  
  final_check=rowSums(prop_mat2)
  if(sum(final_check>1.001)>0){
    stop(paste0("Error in the construction of the propagation matrix for ",dsName," dataset!"))
  }
  
  return(list(data=data,prop_mat=prop_mat2))
}

.myConcensusDEFn_step2_detail_newprop=function(data,centroidPCAdata,argList,n.adaptiveKernel=20,nPropIter=3,exCentroids=NULL,runIndx=1,n.neighbors=2){
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=3;exCentroids=NULL
  batchPCAdata=as.matrix(data$pcaData)[,1:argList$nPCs]
  dsName=data$dsName
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  if(nrow(batchPCAdata)<(3*n.adaptiveKernel+1)){
    n.adaptiveKernel=2
  }
  
  
  nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
    if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0)>0){
      stop("error!")
    }
    adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = nn.ranked.1$nn.dists[,i]
  }
  
  
  nn.centroid.purityIndex=RANN::nn2(data = centroidPCAdata, k = 2, eps = 0)
  nn.centroid.purityIndex=nn.centroid.purityIndex$nn.dists[,2]
  
 
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  pd=pd[match(row.names(batchPCAdata),row.names(pd)),]
  nn.dataset <- RANN::nn2(batchPCAdata, k = min(2000,nrow(batchPCAdata)), eps = 0)
  nn.dataset.purityIndx=RANN::nn2(data = centroidPCAdata,query = batchPCAdata, k = 1, eps = 0)
  nn.dataset.purityIndx=as.numeric(nn.dataset.purityIndx$nn.dist[,1])
  
  
  
  affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinities[which(adj==0)]=0
  centroid.affinities=affinities
  
  dists=nn.dataset$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/(nn.dataset.purityIndx[x]+0.01))^2))
  affinities=do.call("rbind",affinities)
  if(ncol(affinities)>(3*n.adaptiveKernel+1)){
    affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
    affinities2=affinities-affinitiesThr
    affinities[affinities2<=0]=0
  }
  nn.dataset$affinities=affinities
  rm(affinities)
  
  j <- as.numeric(t(nn.dataset$nn.idx))
  i <- ((1:length(j)) - 1)%/%ncol(nn.dataset$affinities) + 1
  k=as.numeric(t(nn.dataset$affinities))
  
  graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(batchPCAdata), nrow(batchPCAdata)))
  rm(i,j,k)
  graph <- .extra_matrix_rowNorm(graph)#Matrix::Diagonal(x = 1 / (rowSums(graph))) %*% graph
  
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  if(class(graph)!="dgCMatrix"){
    graph=as(graph,"dgCMatrix")
  }
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(graph) <- row.names(batchPCAdata)
  colnames(graph) <- row.names(batchPCAdata)
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  checkProp=T
  if(nPropIter<1){
    nPropIter=1
  }
  
  #nPropIter=4
  {
    for(i in nPropIter:1){
      if(i==nPropIter){
        prop_mat=centroid.affinities
      } else {
        prop_mat= prop_mat %*% graph
          
      }
    }
    
    
  }
  
  centroid_weights=1/rowSums(adj)
  
  centroid_weightsN=1/(1+exp((12-2/(centroid_weights))/2))
  prop_mat2=sweep(prop_mat,1,centroid_weightsN,"*")
  
  final_check=rowSums(prop_mat2)
  if(sum(final_check>1.001)>0){
    stop(paste0("Error in the construction of the propagation matrix for ",dsName," dataset!"))
  }
  
  return(list(data=data,prop_mat=prop_mat))
}

.myConcensusDEFn_step2_detail_newprop2=function(data,centroidPCAdata,argList,n.adaptiveKernel=20,nPropIter=3,exCentroids=NULL,runIndx=1,n.neighbors=2){
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=3;exCentroids=NULL;n.neighbors=2
  batchPCAdata=as.matrix(data$pcaData)[,1:argList$nPCs]
  dsName=data$dsName
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  if(nrow(batchPCAdata)<(3*n.adaptiveKernel+1)){
    n.adaptiveKernel=2
  }
  
  
  nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
    if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0)>0){
      stop("error!")
    }
    adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = nn.ranked.1$nn.dists[,i]
  }
  
  
  nn.centroid.purityIndex=RANN::nn2(data = centroidPCAdata, k = 2, eps = 0)
  nn.centroid.purityIndex=nn.centroid.purityIndex$nn.dists[,2]
  
  
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  pd=pd[match(row.names(batchPCAdata),row.names(pd)),]
  
  nn.dataset <- RANN::nn2(batchPCAdata, k = min(2000,nrow(batchPCAdata)), eps = 0)
  nn.dataset.purityIndx=RANN::nn2(data = centroidPCAdata,query = batchPCAdata, k = 1, eps = 0)
  nn.dataset.purityIndx=as.numeric(nn.dataset.purityIndx$nn.dist[,1])
  
  
  
  affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinities[which(adj==0)]=0
  centroid.affinities=affinities
  
  dists=nn.dataset$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/(nn.dataset.purityIndx[x]+0.01))^2))
  affinities=do.call("rbind",affinities)
  if(ncol(affinities)>(3*n.adaptiveKernel+1)){
    affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
    affinities2=affinities-affinitiesThr
    affinities[affinities2<=0]=0
  }
  nn.dataset$affinities=affinities
  rm(affinities)
  
  j <- as.numeric(t(nn.dataset$nn.idx))
  i <- ((1:length(j)) - 1)%/%ncol(nn.dataset$affinities) + 1
  k=as.numeric(t(nn.dataset$affinities))
  
  graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(batchPCAdata), nrow(batchPCAdata)))
  rm(i,j,k)
  graph <- .extra_matrix_rowNorm(graph)#Matrix::Diagonal(x = 1 / (rowSums(graph))) %*% graph
  
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  if(class(graph)!="dgCMatrix"){
    graph=as(graph,"dgCMatrix")
  }
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(graph) <- row.names(batchPCAdata)
  colnames(graph) <- row.names(batchPCAdata)
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  checkProp=T
  if(nPropIter<1){
    nPropIter=1
  }
  
  #nPropIter=4
  {
    for(i in nPropIter:1){
      if(i==nPropIter){
        prop_mat=centroid.affinities
      } else {
        prop_mat= prop_mat %*% graph
        
      }
    }
    
    
  }
  
  prop_mat=prop_mat*adj
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  centroid_weights=1/rowSums(adj)
  
  centroid_weightsN=1/(1+exp((12-2/(centroid_weights))/2))
  prop_mat2=sweep(prop_mat,1,centroid_weightsN,"*")
  
  final_check=rowSums(prop_mat2)
  if(sum(final_check>1.001)>0){
    stop(paste0("Error in the construction of the propagation matrix for ",dsName," dataset!"))
  }
  
  return(list(data=data,prop_mat=prop_mat2))
}

.myConcensusDEFn_step2_detail_newprop3=function(data,centroidPCAdata,argList,n.adaptiveKernel=20,nPropIter=1,exCentroids=NULL,runIndx=1,n.neighbors=4){
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  batchPCAdata=as.matrix(data$pcaData)[,1:argList$nPCs]
  dsName=data$dsName
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  if(nrow(batchPCAdata)<(3*n.adaptiveKernel+1)){
    n.adaptiveKernel=2
  }
  
  
  nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
    if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0)>0){
      stop("error!")
    }
    adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = nn.ranked.1$nn.dists[,i]
  }
  
  
  nn.centroid.purityIndex=RANN::nn2(data = centroidPCAdata, k = 2, eps = 0)
  nn.centroid.purityIndex=nn.centroid.purityIndex$nn.dists[,2]
  
  
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  pd=pd[match(row.names(batchPCAdata),row.names(pd)),]
  
  nn.dataset <- RANN::nn2(batchPCAdata, k = min(2000,nrow(batchPCAdata)), eps = 0)
  nn.dataset.purityIndx=RANN::nn2(data = centroidPCAdata,query = batchPCAdata, k = 1, eps = 0)
  nn.dataset.purityIndx=as.numeric(nn.dataset.purityIndx$nn.dist[,1])
  
  
  
  affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinities[which(adj==0)]=0
  centroid.affinities=affinities
  
  dists=nn.dataset$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/(nn.dataset.purityIndx[x]+0.01))^2))
  affinities=do.call("rbind",affinities)
  if(ncol(affinities)>(3*n.adaptiveKernel+1)){
    affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
    affinities2=affinities-affinitiesThr
    affinities[affinities2<=0]=0
  }
  nn.dataset$affinities=affinities
  rm(affinities)
  
  j <- as.numeric(t(nn.dataset$nn.idx))
  i <- ((1:length(j)) - 1)%/%ncol(nn.dataset$affinities) + 1
  k=as.numeric(t(nn.dataset$affinities))
  
  graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(batchPCAdata), nrow(batchPCAdata)))
  
  rm(i,j,k)
  graph <- .extra_matrix_rowNorm(graph)#Matrix::Diagonal(x = 1 / (rowSums(graph))) %*% graph
  
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  if(class(graph)!="dgCMatrix"){
    graph=as(graph,"dgCMatrix")
  }
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(graph) <- row.names(batchPCAdata)
  colnames(graph) <- row.names(batchPCAdata)
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  checkProp=T
  if(nPropIter<1){
    nPropIter=1
  }
  
  #nPropIter=4
  {
    for(i in nPropIter:1){
      if(i==nPropIter){
        prop_mat=centroid.affinities
      } else {
        prop_mat= prop_mat %*% graph
        
      }
    }
    
    
  }
  
  adj[adj>0]=1
  prop_mat=prop_mat*adj
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  c_c_aff=prop_mat %*% t(prop_mat)
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff
  
  prop_mat=c_c_aff %*% prop_mat
  
  
  return(list(data=data,prop_mat=prop_mat))
}

.myConcensusDEFn_step2_detail_newprop3_final=function(data,centroidPCAdata,argList,exCentroids=NULL,runIndx=1,n.neighbors=4,...){
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  batchPCAdata=as.matrix(data$pcaData)[,1:argList$nPCs]
  dsName=data$dsName
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  
  nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
    if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0)>0){
      stop("error!")
    }
    adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = nn.ranked.1$nn.dists[,i]
  }
  
  
  
  
  #######################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_adj=adj
    tst_adj[tst_adj>0]=1
    
    tst_adj=tst_adj %*% as.matrix(tst)
    
    .tmp_step1=tst_adj
  }
  
  #######################
  
  nn.centroid.purityIndex=apply(adj,1,function(x){
    if(sum(x>0,na.rm = T)>0){
      x=x[which(x>0)]
      x=median(x)
    } else {
      x=0
    }
    
    return(x)
  })
  
  
  affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinities[which(adj==0)]=0
  centroid.affinities=affinities
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  prop_mat=centroid.affinities
  
  
  #################################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_affinities=prop_mat %*% as.matrix(tst)
    
    apply(tst_affinities,2,mean)
    apply(tst_affinities,2,function(x) sum(x>0.5))
    
    eff_size=.myEffSizePropMat(prop_mat)
    
    x=apply(tst_affinities,1,function(x) which(x==max(x))[1])
    x=data.frame(size=eff_size$effective_sample_size,g=x)
    aggregate(size~g,data=x,summary)
    
  }
  
  
  
  ###############################
  
  
  c_c_aff=prop_mat %*% t(prop_mat)
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff
  
  prop_mat=c_c_aff %*% prop_mat
  
  prop_mat2=t(apply(prop_mat,1,function(x) x/max(x)))
  prop_mat[which(prop_mat2<0.1)]=0
  
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  
  return(list(data=data,prop_mat=prop_mat,pseudocell_sim_mat=c_c_aff))
}

.myConcensusDEFn_step2_detail_newprop3_final_v2_split_prop=function(data,centroidPCAdata,argList,exCentroids=NULL,n.neighbors=4,...){
  
  
  if(F){
    pca_centroid[155,]=(pca_centroid[165,]+pca_centroid[145,])/2
  }
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  
  batchPCAdata=as.matrix(data$pcaData)[,1:argList$nPCs]
  dsName=data$dsName
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  
  nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
    if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0)>0){
      stop("error!")
    }
    adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = nn.ranked.1$nn.dists[,i]
  }
  
  
  
  
  #######################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_adj=adj
    tst_adj[tst_adj>0]=1
    
    tst_adj=tst_adj %*% as.matrix(tst)
    
    .tmp_step1=tst_adj
  }
  
  #######################
  
  if(F){
    #old way: it doesn't work when pseuodcells are too sparse!
    nn.centroid.purityIndex=apply(adj,1,function(x){
      
      if(sum(x>0,na.rm = T)>0){
        x=x[which(x>0)]
        x=median(x)
      } else {
        x=0
      }
      
      return(x)
    })
  } else {
    tmp_df=data.frame(pseudocells=nn.ranked.1$nn.idx[,1],dist=nn.ranked.1$nn.dists[,1],stringsAsFactors = F)
    tmp_df=aggregate(dist~pseudocells,data=tmp_df,function(x) median(x,na.rm = T))
    nn.centroid.purityIndex=tmp_df$dist
    names(nn.centroid.purityIndex)=as.character(tmp_df$pseudocells)
    
    
    nn.centroid.purityIndex=nn.centroid.purityIndex[match(row.names(centroidPCAdata),names(nn.centroid.purityIndex))]
    if(sum(is.na(nn.centroid.purityIndex))>0){
      nn.centroid.purityIndex[is.na(nn.centroid.purityIndex)]=0
    }
    tmp_df=quantile(nn.centroid.purityIndex,0.95)
    nn.centroid.purityIndex[nn.centroid.purityIndex>tmp_df]=tmp_df
    rm(tmp_df)
  }
  
  
  
  affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinities[which(adj==0)]=0
  centroid.affinities=affinities
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  prop_mat=centroid.affinities
  
  
  #################################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_affinities=prop_mat %*% as.matrix(tst)
    
    #apply(tst_affinities,2,mean)
    apply(tst_affinities,2,function(x) sum(x>0.5))
    
    
    tst_pd=pd[match(colnames(prop_mat),row.names(pd)),]
    aggregate(prop_mat["8",]~tst_pd$anno_orig_cellState,sum,data=NULL)
    
    eff_size=.myEffSizePropMat(prop_mat)
    
    x=apply(tst_affinities,1,function(x) which(x==max(x))[1])
    x=data.frame(size=eff_size$effective_sample_size,g=x)
    aggregate(size~g,data=x,summary)
    
  }
  
  
  
  ###############################
  
  
  c_c_aff=prop_mat %*% t(prop_mat)
  c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
  diag(c_c_aff)=1
  c_c_aff=pmin(c_c_aff,1)
  c_c_aff=sqrt(c_c_aff*t(c_c_aff))
  #c_c_aff[which(c_c_aff<0.1)]=0
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff
  
  prop_mat=c_c_aff %*% prop_mat
  
  prop_mat2=t(apply(prop_mat,1,function(x) x/max(x)))
  prop_mat[which(prop_mat2<0.1)]=0
  
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  matWeights=.myEffSizePropMat(prop_mat)
  
  matEffectiveSize=matWeights$effective_sample_size
  matWeights=matWeights$centroid_weights
  matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  
  
  return(list(data=data,prop_mat=prop_mat,pseudocell_sim_mat=c_c_aff,matWeights=matWeights,matEffectiveSize=matEffectiveSize))
}

.myConcensusDEFn_step2_detail_newprop3_final_v2=function(dataArranged,centroidPCAdata,argList,exCentroids=NULL,n.neighbors=4,batchPCAdata=NULL,colNormalize=F,...){
  
  if(is.null(batchPCAdata)){
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    batchPCAdata=harmony_embeddings
  }
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  #batchPCAdata=harmony_embeddings
  batchPCAdata=as.matrix(batchPCAdata)[,1:argList$nPCs]
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  
  nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
    if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0)>0){
      stop("error!")
    }
    adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = (nn.ranked.1$nn.dists[,i]+0.01)
  }
  
  
  
  
  #######################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_adj=adj
    tst_adj[tst_adj>0]=1
    
    tst_adj=tst_adj %*% as.matrix(tst)
    
    .tmp_step1=tst_adj
  }
  
  #######################
  
  if(F){
    #old way: it doesn't work when pseuodcells are too sparse!
    nn.centroid.purityIndex=apply(adj,1,function(x){
      
      if(sum(x>0,na.rm = T)>0){
        x=x[which(x>0)]
        x=median(x)
      } else {
        x=0
      }
      
      return(x)
    })
  } else {
    tmp_df=data.frame(pseudocells=row.names(centroidPCAdata)[nn.ranked.1$nn.idx[,1]],dist=nn.ranked.1$nn.dists[,1],stringsAsFactors = F)
    tmp_df=aggregate(dist~pseudocells,data=tmp_df,function(x) median(x,na.rm = T))
    nn.centroid.purityIndex=tmp_df$dist
    names(nn.centroid.purityIndex)=as.character(tmp_df$pseudocells)
    
    
    nn.centroid.purityIndex=nn.centroid.purityIndex[match(row.names(centroidPCAdata),names(nn.centroid.purityIndex))]
    if(sum(is.na(nn.centroid.purityIndex))>0){
      nn.centroid.purityIndex[is.na(nn.centroid.purityIndex)]=0
    }
    tmp_df=quantile(nn.centroid.purityIndex,0.95)
    nn.centroid.purityIndex[nn.centroid.purityIndex>tmp_df]=tmp_df
    rm(tmp_df)
  }
  
  
  
  affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinities[which(adj==0)]=0
  centroid.affinities=affinities
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  
  if(F){
    dataArranged2=dataArranged
    for(i in 1:length(dataArranged2)){
      tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged2[[i]]$pcaData),colnames(prop_mat))],data=dataArranged2[[i]],pseudocell_sim_mat=prop_mat)
      dataArranged2[[i]]=tmp
    }
    p=.sconline.anno2pseudocell_tmp(res_prop=dataArranged2,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T)
    ggsave(plot=p$plot,file="~/myBucket/torm3.pdf",width=12,height = 12)
    
  }
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  prop_mat=centroid.affinities
  
  
  #################################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_affinities=prop_mat %*% as.matrix(tst)
    
    #apply(tst_affinities,2,mean)
    apply(tst_affinities,2,function(x) sum(x>0.5))
    
    
    tst_pd=pd[match(colnames(prop_mat),row.names(pd)),]
    aggregate(prop_mat["8",]~tst_pd$anno_orig_cellState,sum,data=NULL)
    
    eff_size=.myEffSizePropMat(prop_mat)
    
    x=apply(tst_affinities,1,function(x) which(x==max(x))[1])
    x=data.frame(size=eff_size$effective_sample_size,g=x)
    aggregate(size~g,data=x,summary)
    
  }
  
  
  
  ###############################
  
  c_c_aff=t(prop_mat)
  c_c_aff=prop_mat %*% c_c_aff
  c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
  diag(c_c_aff)=1
  c_c_aff@x=pmin(c_c_aff@x,1)
  c_c_aff=sqrt(c_c_aff*t(c_c_aff))
  #c_c_aff[which(c_c_aff<0.1)]=0
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff
  
  prop_mat=c_c_aff %*% prop_mat
  
  
  if(F){
    prop_mat2=t(apply(prop_mat,1,function(x) x/max(x)))
    prop_mat[which(prop_mat2<0.1)]=0
    
  } else {
    rowMax_vals=c()
    for(i in seq(1,nrow(prop_mat),1000)){
      tmp=qlcMatrix::rowMax(prop_mat[i:min(i+999,nrow(prop_mat)),])
      rowMax_vals=c(rowMax_vals,tmp)
    }
    rowMax_vals=unlist(lapply(rowMax_vals,as.numeric))
    
    prop_mat2=.extra_matrix_rowNorm(input_mat = prop_mat,rowValues = 1 / rowMax_vals)# Matrix::Diagonal(x = 1 / rowMax_vals) %*% prop_mat
    prop_mat2=Matrix::drop0(prop_mat2,tol=0.1)
    prop_mat2@x=rep(1,length(prop_mat2@x))
    
    if(colNormalize){
      cat("Column-wise normalization")
      tst=c()
      for(i in seq(1,ncol(prop_mat),10000)){
        tst1=qlcMatrix::colMax(prop_mat[,i:min(i+9999,ncol(prop_mat))])
        tst=c(tst,tst1)
      }
      tst=unlist(lapply(tst,function(x) unlist(as.numeric(x))))
      
      prop_mat3 <- .extra_matrix_colNorm(prop_mat,colValues = 1 / (tst+0.000000000001))# prop_mat %*% Matrix::Diagonal(x = 1 / (tst+0.000000000001))
      prop_mat3=drop0(prop_mat3,tol=0.1)
      
      prop_mat4=prop_mat3+prop_mat2
      prop_mat4@x=1
      prop_mat=prop_mat * prop_mat4
    } else {
      prop_mat=prop_mat * prop_mat2
    }
    
    
    rm(prop_mat2)
  }
  
  
  
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  matWeights=.myEffSizePropMat(prop_mat)
  
  matEffectiveSize=matWeights$effective_sample_size
  matWeights=matWeights$centroid_weights
  matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  
  
  if(F){
    cat("    Removing redundant pseudocells")
    tst=prop_mat
    #tst[tst>0]=1
    tst=sqrt(tst %*% t(tst))
    tst=sweep(tst,1,diag(tst),"/")
    cosine_sim=(tst+t(tst))/2
    #cosine_sim=proxyC::simil(prop_mat,margin = 1,method = "cosine")
    
    dataArranged2=dataArranged
    for(i in 1:length(dataArranged2)){
      tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged2[[i]]$pcaData),colnames(prop_mat))],data=dataArranged2[[i]],pseudocell_sim_mat=prop_mat)
      dataArranged2[[i]]=tmp
    }
    
    #res_prop=dataArranged2;annoCol="anno_orig_cellState";collapse_datasets=T;return_plot=T
    p=.sconline.anno2pseudocell_tmp(res_prop=dataArranged2,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T)
    
    
    summary(apply(cosine_sim,1,function(x) sum(x>0.9)))
    which(apply(cosine_sim,1,function(x) sum(x>0.9))>10)
    
    tst=cosine_sim[row.names(p$results),row.names(p$results)]
    summary(apply(tst,1,function(x) sum(x>0.9)))
    which(apply(tst,1,function(x) sum(x>0.9))>10)
    tst=p$results[which(tst[13,]>0.9),]
    tst2=apply(tst[,-c(1,ncol(tst))],2,max)
    tst=tst[,c(1,which(tst2>0.1),ncol(tst))]
    
    if(sum(upper.tri(cosine_sim)&cosine_sim>0.95)>0){
      ident_list=which(upper.tri(cosine_sim)&cosine_sim>0.95,arr.ind = T)
      ident_list=igraph::components(igraph::graph_from_data_frame(as.data.frame(ident_list)))
      ident_list=data.frame(group=ident_list$membership,pseudocell=names(ident_list$membership),stringsAsFactors = F)
      ident_list=split(ident_list$pseudocell,ident_list$group)
      torm_list=c()
      for(i in 1:length(ident_list)){
        x=prop_mat[ident_list[[i]],]
        x=apply(x,2,sum)
        prop_mat[ident_list[[i]][1],]=x/(sum(x)+0.000000000001)
        torm_list=c(torm_list,ident_list[[i]][2:length(ident_list[[i]])])
      }
      prop_mat=prop_mat[-which(row.names(prop_mat) %in% torm_list),]
    }
    
  }
  
  for(i in 1:length(dataArranged)){
    tmp_prop_mat=prop_mat[,match(row.names(dataArranged[[i]]$pcaData),colnames(prop_mat))]
    #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
    tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
    tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
    tmp_weights[tmp_effsize<4]=0
    tmp=list(prop_mat=tmp_prop_mat,data=dataArranged[[i]],pseudocell_sim_mat=c_c_aff,matWeights=tmp_weights,matEffectiveSize=tmp_effsize)
    dataArranged[[i]]=tmp
  }
  
  return(dataArranged)
}

.myConcensusDEFn_step2_detail_newprop3_final_v3=function(dataArranged,centroidPCAdata,argList,exCentroids=NULL,n.neighbors=4,batchPCAdata=NULL,colNormalize=F,...){
  
  if(is.null(batchPCAdata)){
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    batchPCAdata=harmony_embeddings
  }
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  #batchPCAdata=harmony_embeddings
  batchPCAdata=as.matrix(batchPCAdata)[,1:argList$nPCs]
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  
  nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
    if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0)>0){
      stop("error!")
    }
    adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = (nn.ranked.1$nn.dists[,i]+0.01)
  }
  
  
  
  
  #######################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_adj=adj
    tst_adj[tst_adj>0]=1
    
    tst_adj=tst_adj %*% as.matrix(tst)
    
    .tmp_step1=tst_adj
  }
  
  #######################
  
  if(F){
    #old way: it doesn't work when pseuodcells are too sparse!
    nn.centroid.purityIndex=apply(adj,1,function(x){
      
      if(sum(x>0,na.rm = T)>0){
        x=x[which(x>0)]
        x=median(x)
      } else {
        x=0
      }
      
      return(x)
    })
  } else {
    tmp_df=data.frame(pseudocells=row.names(centroidPCAdata)[nn.ranked.1$nn.idx[,1]],dist=nn.ranked.1$nn.dists[,1],stringsAsFactors = F)
    tmp_df=aggregate(dist~pseudocells,data=tmp_df,function(x) median(x,na.rm = T))
    nn.centroid.purityIndex=tmp_df$dist
    names(nn.centroid.purityIndex)=as.character(tmp_df$pseudocells)
    
    
    nn.centroid.purityIndex=nn.centroid.purityIndex[match(row.names(centroidPCAdata),names(nn.centroid.purityIndex))]
    if(sum(is.na(nn.centroid.purityIndex))>0){
      nn.centroid.purityIndex[is.na(nn.centroid.purityIndex)]=0
    }
    tmp_df=quantile(nn.centroid.purityIndex,0.95)
    nn.centroid.purityIndex[nn.centroid.purityIndex>tmp_df]=tmp_df
    rm(tmp_df)
  }
  
  
  
  affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinities[which(adj==0)]=0
  centroid.affinities=affinities
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  
  if(F){
    dataArranged2=dataArranged
    for(i in 1:length(dataArranged2)){
      tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged2[[i]]$pcaData),colnames(prop_mat))],data=dataArranged2[[i]],pseudocell_sim_mat=prop_mat)
      dataArranged2[[i]]=tmp
    }
    p=.sconline.anno2pseudocell_tmp(res_prop=dataArranged2,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T)
    ggsave(plot=p$plot,file="~/myBucket/torm3.pdf",width=12,height = 12)
    
  }
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  prop_mat=centroid.affinities
  
  
  #################################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_affinities=prop_mat %*% as.matrix(tst)
    
    #apply(tst_affinities,2,mean)
    apply(tst_affinities,2,function(x) sum(x>0.5))
    
    
    tst_pd=pd[match(colnames(prop_mat),row.names(pd)),]
    aggregate(prop_mat["8",]~tst_pd$anno_orig_cellState,sum,data=NULL)
    
    eff_size=.myEffSizePropMat(prop_mat)
    
    x=apply(tst_affinities,1,function(x) which(x==max(x))[1])
    x=data.frame(size=eff_size$effective_sample_size,g=x)
    aggregate(size~g,data=x,summary)
    
  }
  
  
  
  ###############################
  
  c_c_aff=t(prop_mat)
  c_c_aff=prop_mat %*% c_c_aff
  c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
  diag(c_c_aff)=1
  c_c_aff=pmin(c_c_aff,1)
  c_c_aff=sqrt(c_c_aff*t(c_c_aff))
  #c_c_aff[which(c_c_aff<0.1)]=0
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff
  
  prop_mat=c_c_aff %*% prop_mat
  
  
  if(F){
    prop_mat2=t(apply(prop_mat,1,function(x) x/max(x)))
    prop_mat[which(prop_mat2<0.1)]=0
    
  } else {
    cat("Column-wise normalization")
    tst=c()
    for(i in seq(1,ncol(prop_mat),10000)){
      tst1=qlcMatrix::colMax(prop_mat[,i:min(i+9999,ncol(prop_mat))])
      tst=c(tst,tst1)
    }
    tst=unlist(lapply(tst,function(x) unlist(as.numeric(x))))
    
    prop_mat <- .extra_matrix_colNorm(input_mat = prop_mat,colValues = 1 / (tst+0.000000000001))# prop_mat %*% Matrix::Diagonal(x = 1 / (tst+0.000000000001))
    prop_mat=drop0(prop_mat,tol=0.1)
    
  }
  
  
  
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  matWeights=.myEffSizePropMat(prop_mat)
  
  matEffectiveSize=matWeights$effective_sample_size
  matWeights=matWeights$centroid_weights
  matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  
  
  if(F){
    cat("    Removing redundant pseudocells")
    tst=prop_mat
    #tst[tst>0]=1
    tst=sqrt(tst %*% t(tst))
    tst=sweep(tst,1,diag(tst),"/")
    cosine_sim=(tst+t(tst))/2
    #cosine_sim=proxyC::simil(prop_mat,margin = 1,method = "cosine")
    
    dataArranged2=dataArranged
    for(i in 1:length(dataArranged2)){
      tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged2[[i]]$pcaData),colnames(prop_mat))],data=dataArranged2[[i]],pseudocell_sim_mat=prop_mat)
      dataArranged2[[i]]=tmp
    }
    
    #res_prop=dataArranged2;annoCol="anno_orig_cellState";collapse_datasets=T;return_plot=T
    p=.sconline.anno2pseudocell_tmp(res_prop=dataArranged2,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T)
    
    
    summary(apply(cosine_sim,1,function(x) sum(x>0.9)))
    which(apply(cosine_sim,1,function(x) sum(x>0.9))>10)
    
    tst=cosine_sim[row.names(p$results),row.names(p$results)]
    summary(apply(tst,1,function(x) sum(x>0.9)))
    which(apply(tst,1,function(x) sum(x>0.9))>10)
    tst=p$results[which(tst[13,]>0.9),]
    tst2=apply(tst[,-c(1,ncol(tst))],2,max)
    tst=tst[,c(1,which(tst2>0.1),ncol(tst))]
    
    if(sum(upper.tri(cosine_sim)&cosine_sim>0.95)>0){
      ident_list=which(upper.tri(cosine_sim)&cosine_sim>0.95,arr.ind = T)
      ident_list=igraph::components(igraph::graph_from_data_frame(as.data.frame(ident_list)))
      ident_list=data.frame(group=ident_list$membership,pseudocell=names(ident_list$membership),stringsAsFactors = F)
      ident_list=split(ident_list$pseudocell,ident_list$group)
      torm_list=c()
      for(i in 1:length(ident_list)){
        x=prop_mat[ident_list[[i]],]
        x=apply(x,2,sum)
        prop_mat[ident_list[[i]][1],]=x/(sum(x)+0.000000000001)
        torm_list=c(torm_list,ident_list[[i]][2:length(ident_list[[i]])])
      }
      prop_mat=prop_mat[-which(row.names(prop_mat) %in% torm_list),]
    }
    
  }
  
  for(i in 1:length(dataArranged)){
    tmp_prop_mat=prop_mat[,match(row.names(dataArranged[[i]]$pcaData),colnames(prop_mat))]
    #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
    tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
    tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
    tmp_weights[tmp_effsize<4]=0
    tmp=list(prop_mat=tmp_prop_mat,data=dataArranged[[i]],pseudocell_sim_mat=c_c_aff,matWeights=tmp_weights,matEffectiveSize=tmp_effsize)
    dataArranged[[i]]=tmp
  }
  
  return(dataArranged)
}

.myConcensusDEFn_step2_detail_newprop3_final_v4=function(dataArranged,centroidPCAdata,argList,exCentroids=NULL,n.neighbors=4,batchPCAdata=NULL,colNormalize=F,...){
  
  if(is.null(batchPCAdata)){
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    batchPCAdata=harmony_embeddings
  }
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  #batchPCAdata=harmony_embeddings
  batchPCAdata=as.matrix(batchPCAdata)[,1:argList$nPCs]
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  
  nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
    if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0)>0){
      stop("error!")
    }
    adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = (nn.ranked.1$nn.dists[,i]+0.01)
  }
  
  
  
  
  #######################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_adj=adj
    tst_adj[tst_adj>0]=1
    
    tst_adj=tst_adj %*% as.matrix(tst)
    
    .tmp_step1=tst_adj
  }
  
  #######################
  
  if(F){
    #old way: it doesn't work when pseuodcells are too sparse!
    nn.centroid.purityIndex=apply(adj,1,function(x){
      
      if(sum(x>0,na.rm = T)>0){
        x=x[which(x>0)]
        x=median(x)
      } else {
        x=0
      }
      
      return(x)
    })
  } else {
    tmp_df=data.frame(pseudocells=row.names(centroidPCAdata)[nn.ranked.1$nn.idx[,1]],dist=nn.ranked.1$nn.dists[,1],stringsAsFactors = F)
    tmp_df=aggregate(dist~pseudocells,data=tmp_df,function(x) median(x,na.rm = T))
    nn.centroid.purityIndex=tmp_df$dist
    names(nn.centroid.purityIndex)=as.character(tmp_df$pseudocells)
    
    
    nn.centroid.purityIndex=nn.centroid.purityIndex[match(row.names(centroidPCAdata),names(nn.centroid.purityIndex))]
    if(sum(is.na(nn.centroid.purityIndex))>0){
      nn.centroid.purityIndex[is.na(nn.centroid.purityIndex)]=0
    }
    tmp_df=quantile(nn.centroid.purityIndex,0.95)
    nn.centroid.purityIndex[nn.centroid.purityIndex>tmp_df]=tmp_df
    rm(tmp_df)
  }
  
  
  
  affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinities[which(adj==0)]=0
  centroid.affinities=affinities
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  
  if(F){
    dataArranged2=dataArranged
    for(i in 1:length(dataArranged2)){
      tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged2[[i]]$pcaData),colnames(prop_mat))],data=dataArranged2[[i]],pseudocell_sim_mat=prop_mat)
      dataArranged2[[i]]=tmp
    }
    p=.sconline.anno2pseudocell_tmp(res_prop=dataArranged2,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T)
    ggsave(plot=p$plot,file="~/myBucket/torm3.pdf",width=12,height = 12)
    
  }
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  prop_mat=centroid.affinities
  
  tst=c()
  for(i in seq(1,ncol(prop_mat),10000)){
    tst1=qlcMatrix::colMax(prop_mat[,i:min(i+9999,ncol(prop_mat))])
    tst=c(tst,tst1)
  }
  tst=unlist(lapply(tst,function(x) unlist(as.numeric(x))))
  
  prop_mat <-.extra_matrix_colNorm(input_mat = prop_mat,colValues = 1 / (tst+0.000000000001))# prop_mat %*% Matrix::Diagonal(x = 1 / (tst+0.000000000001))
  prop_mat=Matrix::drop0(prop_mat,0.1)
  prop_mat = .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  
  #################################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_affinities=prop_mat %*% as.matrix(tst)
    
    #apply(tst_affinities,2,mean)
    apply(tst_affinities,2,function(x) sum(x>0.5))
    
    
    tst_pd=pd[match(colnames(prop_mat),row.names(pd)),]
    aggregate(prop_mat["8",]~tst_pd$anno_orig_cellState,sum,data=NULL)
    
    eff_size=.myEffSizePropMat(prop_mat)
    
    x=apply(tst_affinities,1,function(x) which(x==max(x))[1])
    x=data.frame(size=eff_size$effective_sample_size,g=x)
    aggregate(size~g,data=x,summary)
    
  }
  
  
  
  ###############################
  
  c_c_aff=t(prop_mat)
  c_c_aff=prop_mat %*% c_c_aff
  c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
  diag(c_c_aff)=1
  c_c_aff@x=pmin(c_c_aff@x,1)
  c_c_aff=sqrt(c_c_aff*t(c_c_aff))
  #c_c_aff[which(c_c_aff<0.1)]=0
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff
  
  prop_mat=c_c_aff %*% prop_mat
  
  
  if(F){
    prop_mat2=t(apply(prop_mat,1,function(x) x/max(x)))
    prop_mat[which(prop_mat2<0.1)]=0
    
  } else {
    rowMax_vals=c()
    for(i in seq(1,nrow(prop_mat),1000)){
      tmp=qlcMatrix::rowMax(prop_mat[i:min(i+999,nrow(prop_mat)),])
      rowMax_vals=c(rowMax_vals,tmp)
    }
    rowMax_vals=unlist(lapply(rowMax_vals,as.numeric))
    
    prop_mat2=.extra_matrix_rowNorm(input_mat = prop_mat,rowValues = 1 / rowMax_vals)# Matrix::Diagonal(x = 1 / rowMax_vals) %*% prop_mat
    prop_mat2=Matrix::drop0(prop_mat2,tol=0.1)
    prop_mat2@x=rep(1,length(prop_mat2@x))
    prop_mat=prop_mat * prop_mat2
    
    
    rm(prop_mat2)
  }
  
  
  
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  matWeights=.myEffSizePropMat(prop_mat)
  
  matEffectiveSize=matWeights$effective_sample_size
  matWeights=matWeights$centroid_weights
  matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  
  
  if(F){
    cat("    Removing redundant pseudocells")
    tst=prop_mat
    #tst[tst>0]=1
    tst=sqrt(tst %*% t(tst))
    tst=sweep(tst,1,diag(tst),"/")
    cosine_sim=(tst+t(tst))/2
    #cosine_sim=proxyC::simil(prop_mat,margin = 1,method = "cosine")
    
    dataArranged2=dataArranged
    for(i in 1:length(dataArranged2)){
      tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged2[[i]]$pcaData),colnames(prop_mat))],data=dataArranged2[[i]],pseudocell_sim_mat=prop_mat)
      dataArranged2[[i]]=tmp
    }
    
    #res_prop=dataArranged2;annoCol="anno_orig_cellState";collapse_datasets=T;return_plot=T
    p=.sconline.anno2pseudocell_tmp(res_prop=dataArranged2,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T)
    
    
    summary(apply(cosine_sim,1,function(x) sum(x>0.9)))
    which(apply(cosine_sim,1,function(x) sum(x>0.9))>10)
    
    tst=cosine_sim[row.names(p$results),row.names(p$results)]
    summary(apply(tst,1,function(x) sum(x>0.9)))
    which(apply(tst,1,function(x) sum(x>0.9))>10)
    tst=p$results[which(tst[13,]>0.9),]
    tst2=apply(tst[,-c(1,ncol(tst))],2,max)
    tst=tst[,c(1,which(tst2>0.1),ncol(tst))]
    
    if(sum(upper.tri(cosine_sim)&cosine_sim>0.95)>0){
      ident_list=which(upper.tri(cosine_sim)&cosine_sim>0.95,arr.ind = T)
      ident_list=igraph::components(igraph::graph_from_data_frame(as.data.frame(ident_list)))
      ident_list=data.frame(group=ident_list$membership,pseudocell=names(ident_list$membership),stringsAsFactors = F)
      ident_list=split(ident_list$pseudocell,ident_list$group)
      torm_list=c()
      for(i in 1:length(ident_list)){
        x=prop_mat[ident_list[[i]],]
        x=apply(x,2,sum)
        prop_mat[ident_list[[i]][1],]=x/(sum(x)+0.000000000001)
        torm_list=c(torm_list,ident_list[[i]][2:length(ident_list[[i]])])
      }
      prop_mat=prop_mat[-which(row.names(prop_mat) %in% torm_list),]
    }
    
  }
  
  for(i in 1:length(dataArranged)){
    tmp_prop_mat=prop_mat[,match(row.names(dataArranged[[i]]$pcaData),colnames(prop_mat))]
    #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
    tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
    tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
    tmp_weights[tmp_effsize<4]=0
    tmp=list(prop_mat=tmp_prop_mat,data=dataArranged[[i]],pseudocell_sim_mat=c_c_aff,matWeights=tmp_weights,matEffectiveSize=tmp_effsize)
    dataArranged[[i]]=tmp
  }
  
  return(dataArranged)
}

.myConcensusDEFn_step2_detail_newprop3_final_v5=function(dataArranged,centroidPCAdata,argList,exCentroids=NULL,n.neighbors=4,batchPCAdata=NULL,colNormalize=F,...){
  
  if(is.null(batchPCAdata)){
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    batchPCAdata=harmony_embeddings
  }
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  #batchPCAdata=harmony_embeddings
  batchPCAdata=as.matrix(batchPCAdata)[,1:argList$nPCs]
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  
  nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
    if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0)>0){
      stop("error!")
    }
    adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = (nn.ranked.1$nn.dists[,i]+0.01)
  }
  
  
  
  
  #######################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_adj=adj
    tst_adj[tst_adj>0]=1
    
    tst_adj=tst_adj %*% as.matrix(tst)
    
    .tmp_step1=tst_adj
  }
  
  #######################
  
  if(F){
    #old way: it doesn't work when pseuodcells are too sparse!
    nn.centroid.purityIndex=apply(adj,1,function(x){
      
      if(sum(x>0,na.rm = T)>0){
        x=x[which(x>0)]
        x=median(x)
      } else {
        x=0
      }
      
      return(x)
    })
  } else {
    tmp_df=data.frame(pseudocells=row.names(centroidPCAdata)[nn.ranked.1$nn.idx[,1]],dist=nn.ranked.1$nn.dists[,1],stringsAsFactors = F)
    tmp_df=aggregate(dist~pseudocells,data=tmp_df,function(x) median(x,na.rm = T))
    nn.centroid.purityIndex=tmp_df$dist
    names(nn.centroid.purityIndex)=as.character(tmp_df$pseudocells)
    
    
    nn.centroid.purityIndex=nn.centroid.purityIndex[match(row.names(centroidPCAdata),names(nn.centroid.purityIndex))]
    if(sum(is.na(nn.centroid.purityIndex))>0){
      nn.centroid.purityIndex[is.na(nn.centroid.purityIndex)]=0
    }
    tmp_df=quantile(nn.centroid.purityIndex,0.95)
    nn.centroid.purityIndex[nn.centroid.purityIndex>tmp_df]=tmp_df
    rm(tmp_df)
  }
  
  
  
  affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinities[which(adj==0)]=0
  
  tst=c()
  for(i in seq(1,ncol(affinities),10000)){
    tst1=qlcMatrix::colMax(affinities[,i:min(i+9999,ncol(affinities))])
    tst=c(tst,tst1)
  }
  tst=unlist(lapply(tst,function(x) unlist(as.numeric(x))))
  
  affinities = .extra_matrix_colNorm(input_mat = affinities,colValues = 1 / (tst+0.000000000001)) #affinities %*% Matrix::Diagonal(x = 1 / (tst+0.000000000001))
  
  centroid.affinities=affinities
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  
  if(F){
    dataArranged2=dataArranged
    for(i in 1:length(dataArranged2)){
      tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged2[[i]]$pcaData),colnames(prop_mat))],data=dataArranged2[[i]],pseudocell_sim_mat=prop_mat)
      dataArranged2[[i]]=tmp
    }
    p=.sconline.anno2pseudocell_tmp(res_prop=dataArranged2,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T)
    ggsave(plot=p$plot,file="~/myBucket/torm3.pdf",width=12,height = 12)
    
  }
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  tst=unlist(lapply(1:nrow(centroid.affinities),function(x) sum(centroid.affinities[x,]>0)))
  
  prop_mat=centroid.affinities[tst>1,]
  
  
  #################################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_affinities=prop_mat %*% as.matrix(tst)
    
    #apply(tst_affinities,2,mean)
    apply(tst_affinities,2,function(x) sum(x>0.5))
    
    
    tst_pd=pd[match(colnames(prop_mat),row.names(pd)),]
    aggregate(prop_mat["8",]~tst_pd$anno_orig_cellState,sum,data=NULL)
    
    eff_size=.myEffSizePropMat(prop_mat)
    
    x=apply(tst_affinities,1,function(x) which(x==max(x))[1])
    x=data.frame(size=eff_size$effective_sample_size,g=x)
    aggregate(size~g,data=x,summary)
    
  }
  
  
  
  ###############################
  
  c_c_aff=t(prop_mat)
  c_c_aff=prop_mat %*% c_c_aff
  c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
  diag(c_c_aff)=1
  c_c_aff@x=pmin(c_c_aff@x,1)
  c_c_aff=sqrt(c_c_aff*t(c_c_aff))
  #c_c_aff[which(c_c_aff<0.1)]=0
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff
  
  prop_mat=c_c_aff %*% prop_mat
  
  
  if(F){
    prop_mat2=t(apply(prop_mat,1,function(x) x/max(x)))
    prop_mat[which(prop_mat2<0.1)]=0
    
  } else {
    rowMax_vals=c()
    for(i in seq(1,nrow(prop_mat),1000)){
      tmp=qlcMatrix::rowMax(prop_mat[i:min(i+999,nrow(prop_mat)),])
      rowMax_vals=c(rowMax_vals,tmp)
    }
    rowMax_vals=unlist(lapply(rowMax_vals,as.numeric))
    
    prop_mat2=.extra_matrix_rowNorm(input_mat = prop_mat,rowValues = 1 / rowMax_vals)# Matrix::Diagonal(x = 1 / rowMax_vals) %*% prop_mat
    prop_mat2=Matrix::drop0(prop_mat2,tol=0.1)
    prop_mat2@x=rep(1,length(prop_mat2@x))
    prop_mat=prop_mat * prop_mat2
    
    
    rm(prop_mat2)
  }
  
  
  
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  matWeights=.myEffSizePropMat(prop_mat)
  
  matEffectiveSize=matWeights$effective_sample_size
  matWeights=matWeights$centroid_weights
  matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  
  
  if(F){
    cat("    Removing redundant pseudocells")
    tst=prop_mat
    #tst[tst>0]=1
    tst=sqrt(tst %*% t(tst))
    tst=sweep(tst,1,diag(tst),"/")
    cosine_sim=(tst+t(tst))/2
    #cosine_sim=proxyC::simil(prop_mat,margin = 1,method = "cosine")
    
    dataArranged2=dataArranged
    for(i in 1:length(dataArranged2)){
      tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged2[[i]]$pcaData),colnames(prop_mat))],data=dataArranged2[[i]],pseudocell_sim_mat=prop_mat)
      dataArranged2[[i]]=tmp
    }
    
    #res_prop=dataArranged2;annoCol="anno_orig_cellState";collapse_datasets=T;return_plot=T
    p=.sconline.anno2pseudocell_tmp(res_prop=dataArranged2,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T)
    
    
    summary(apply(cosine_sim,1,function(x) sum(x>0.9)))
    which(apply(cosine_sim,1,function(x) sum(x>0.9))>10)
    
    tst=cosine_sim[row.names(p$results),row.names(p$results)]
    summary(apply(tst,1,function(x) sum(x>0.9)))
    which(apply(tst,1,function(x) sum(x>0.9))>10)
    tst=p$results[which(tst[13,]>0.9),]
    tst2=apply(tst[,-c(1,ncol(tst))],2,max)
    tst=tst[,c(1,which(tst2>0.1),ncol(tst))]
    
    if(sum(upper.tri(cosine_sim)&cosine_sim>0.95)>0){
      ident_list=which(upper.tri(cosine_sim)&cosine_sim>0.95,arr.ind = T)
      ident_list=igraph::components(igraph::graph_from_data_frame(as.data.frame(ident_list)))
      ident_list=data.frame(group=ident_list$membership,pseudocell=names(ident_list$membership),stringsAsFactors = F)
      ident_list=split(ident_list$pseudocell,ident_list$group)
      torm_list=c()
      for(i in 1:length(ident_list)){
        x=prop_mat[ident_list[[i]],]
        x=apply(x,2,sum)
        prop_mat[ident_list[[i]][1],]=x/(sum(x)+0.000000000001)
        torm_list=c(torm_list,ident_list[[i]][2:length(ident_list[[i]])])
      }
      prop_mat=prop_mat[-which(row.names(prop_mat) %in% torm_list),]
    }
    
  }
  
  for(i in 1:length(dataArranged)){
    tmp_prop_mat=prop_mat[,match(row.names(dataArranged[[i]]$pcaData),colnames(prop_mat))]
    #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
    tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
    tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
    tmp_weights[tmp_effsize<4]=0
    tmp=list(prop_mat=tmp_prop_mat,data=dataArranged[[i]],pseudocell_sim_mat=c_c_aff,matWeights=tmp_weights,matEffectiveSize=tmp_effsize)
    dataArranged[[i]]=tmp
  }
  
  return(dataArranged)
}

.myConcensusDEFn_step2_detail_newprop3_final_v6=function(dataArranged,centroidPCAdata,argList,exCentroids=NULL,n.neighbors=4,batchPCAdata=NULL,colNormalize=F,...){
  
  if(is.null(batchPCAdata)){
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    batchPCAdata=harmony_embeddings
  }
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  #batchPCAdata=harmony_embeddings
  batchPCAdata=as.matrix(batchPCAdata)[,1:argList$nPCs]
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  
  nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
    if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0)>0){
      stop("error!")
    }
    adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = (nn.ranked.1$nn.dists[,i]+0.01)
  }
  
  
  
  
  #######################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_adj=adj
    tst_adj[tst_adj>0]=1
    
    tst_adj=tst_adj %*% as.matrix(tst)
    
    .tmp_step1=tst_adj
  }
  
  #######################
  
  if(F){
    #old way: it doesn't work when pseuodcells are too sparse!
    nn.centroid.purityIndex=apply(adj,1,function(x){
      
      if(sum(x>0,na.rm = T)>0){
        x=x[which(x>0)]
        x=median(x)
      } else {
        x=0
      }
      
      return(x)
    })
  } else {
    tmp_df=data.frame(pseudocells=row.names(centroidPCAdata)[nn.ranked.1$nn.idx[,1]],dist=nn.ranked.1$nn.dists[,1],stringsAsFactors = F)
    tmp_df=aggregate(dist~pseudocells,data=tmp_df,function(x) median(x,na.rm = T))
    nn.centroid.purityIndex=tmp_df$dist
    names(nn.centroid.purityIndex)=as.character(tmp_df$pseudocells)
    
    
    nn.centroid.purityIndex=nn.centroid.purityIndex[match(row.names(centroidPCAdata),names(nn.centroid.purityIndex))]
    if(sum(is.na(nn.centroid.purityIndex))>0){
      nn.centroid.purityIndex[is.na(nn.centroid.purityIndex)]=0
    }
    tmp_df=quantile(nn.centroid.purityIndex,0.95)
    nn.centroid.purityIndex[nn.centroid.purityIndex>tmp_df]=tmp_df
    rm(tmp_df)
  }
  
  
  
  affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinities[which(adj==0)]=0
  
  tst=c()
  for(i in seq(1,ncol(affinities),10000)){
    tst1=qlcMatrix::colMax(affinities[,i:min(i+9999,ncol(affinities))])
    tst=c(tst,tst1)
  }
  tst=unlist(lapply(tst,function(x) unlist(as.numeric(x))))
  
  affinities = .extra_matrix_colNorm(input_mat = affinities,colValues = 1 / (tst+0.000000000001)) #affinities %*% Matrix::Diagonal(x = 1 / (tst+0.000000000001))
  affinities=Matrix::drop0(affinities,tol=0.1)
  
  centroid.affinities=affinities
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  
  if(F){
    dataArranged2=dataArranged
    for(i in 1:length(dataArranged2)){
      tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged2[[i]]$pcaData),colnames(prop_mat))],data=dataArranged2[[i]],pseudocell_sim_mat=prop_mat)
      dataArranged2[[i]]=tmp
    }
    p=.sconline.anno2pseudocell_tmp(res_prop=dataArranged2,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T)
    ggsave(plot=p$plot,file="~/myBucket/torm3.pdf",width=12,height = 12)
    
  }
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  tst=unlist(lapply(1:nrow(centroid.affinities),function(x) sum(centroid.affinities[x,]>0)))
  
  prop_mat=centroid.affinities[tst>1,]
  
  
  #################################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_affinities=prop_mat %*% as.matrix(tst)
    
    #apply(tst_affinities,2,mean)
    apply(tst_affinities,2,function(x) sum(x>0.5))
    
    
    tst_pd=pd[match(colnames(prop_mat),row.names(pd)),]
    aggregate(prop_mat["8",]~tst_pd$anno_orig_cellState,sum,data=NULL)
    
    eff_size=.myEffSizePropMat(prop_mat)
    
    x=apply(tst_affinities,1,function(x) which(x==max(x))[1])
    x=data.frame(size=eff_size$effective_sample_size,g=x)
    aggregate(size~g,data=x,summary)
    
  }
  
  
  
  ###############################
  
  c_c_aff=t(prop_mat)
  c_c_aff=prop_mat %*% c_c_aff
  c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
  diag(c_c_aff)=1
  c_c_aff@x=pmin(c_c_aff@x,1)
  c_c_aff=sqrt(c_c_aff*t(c_c_aff))
  #c_c_aff[which(c_c_aff<0.1)]=0
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff
  
  prop_mat=c_c_aff %*% prop_mat
  
  
  if(F){
    prop_mat2=t(apply(prop_mat,1,function(x) x/max(x)))
    prop_mat[which(prop_mat2<0.1)]=0
    
  } else {
    rowMax_vals=c()
    for(i in seq(1,nrow(prop_mat),1000)){
      tmp=qlcMatrix::rowMax(prop_mat[i:min(i+999,nrow(prop_mat)),])
      rowMax_vals=c(rowMax_vals,tmp)
    }
    rowMax_vals=unlist(lapply(rowMax_vals,as.numeric))
    
    prop_mat2=.extra_matrix_rowNorm(input_mat = prop_mat,rowValues = 1/rowMax_vals)# Matrix::Diagonal(x = 1 / rowMax_vals) %*% prop_mat
    prop_mat2=Matrix::drop0(prop_mat2,tol=0.1)
    prop_mat2@x=rep(1,length(prop_mat2@x))
    prop_mat=prop_mat * prop_mat2
    
    
    rm(prop_mat2)
  }
  
  
  
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  matWeights=.myEffSizePropMat(prop_mat)
  
  matEffectiveSize=matWeights$effective_sample_size
  matWeights=matWeights$centroid_weights
  matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  
  
  if(F){
    cat("    Removing redundant pseudocells")
    tst=prop_mat
    #tst[tst>0]=1
    tst=sqrt(tst %*% t(tst))
    tst=sweep(tst,1,diag(tst),"/")
    cosine_sim=(tst+t(tst))/2
    #cosine_sim=proxyC::simil(prop_mat,margin = 1,method = "cosine")
    
    dataArranged2=dataArranged
    for(i in 1:length(dataArranged2)){
      tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged2[[i]]$pcaData),colnames(prop_mat))],data=dataArranged2[[i]],pseudocell_sim_mat=prop_mat)
      dataArranged2[[i]]=tmp
    }
    
    #res_prop=dataArranged2;annoCol="anno_orig_cellState";collapse_datasets=T;return_plot=T
    p=.sconline.anno2pseudocell_tmp(res_prop=dataArranged2,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T)
    
    
    summary(apply(cosine_sim,1,function(x) sum(x>0.9)))
    which(apply(cosine_sim,1,function(x) sum(x>0.9))>10)
    
    tst=cosine_sim[row.names(p$results),row.names(p$results)]
    summary(apply(tst,1,function(x) sum(x>0.9)))
    which(apply(tst,1,function(x) sum(x>0.9))>10)
    tst=p$results[which(tst[13,]>0.9),]
    tst2=apply(tst[,-c(1,ncol(tst))],2,max)
    tst=tst[,c(1,which(tst2>0.1),ncol(tst))]
    
    if(sum(upper.tri(cosine_sim)&cosine_sim>0.95)>0){
      ident_list=which(upper.tri(cosine_sim)&cosine_sim>0.95,arr.ind = T)
      ident_list=igraph::components(igraph::graph_from_data_frame(as.data.frame(ident_list)))
      ident_list=data.frame(group=ident_list$membership,pseudocell=names(ident_list$membership),stringsAsFactors = F)
      ident_list=split(ident_list$pseudocell,ident_list$group)
      torm_list=c()
      for(i in 1:length(ident_list)){
        x=prop_mat[ident_list[[i]],]
        x=apply(x,2,sum)
        prop_mat[ident_list[[i]][1],]=x/(sum(x)+0.000000000001)
        torm_list=c(torm_list,ident_list[[i]][2:length(ident_list[[i]])])
      }
      prop_mat=prop_mat[-which(row.names(prop_mat) %in% torm_list),]
    }
    
  }
  
  for(i in 1:length(dataArranged)){
    tmp_prop_mat=prop_mat[,match(row.names(dataArranged[[i]]$pcaData),colnames(prop_mat))]
    #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
    tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
    tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
    tmp_weights[tmp_effsize<4]=0
    tmp=list(prop_mat=tmp_prop_mat,data=dataArranged[[i]],pseudocell_sim_mat=c_c_aff,matWeights=tmp_weights,matEffectiveSize=tmp_effsize)
    dataArranged[[i]]=tmp
  }
  
  return(dataArranged)
}

.myConcensusDEFn_step2_detail_newprop3_final_v7=function(dataArranged,centroidPCAdata,argList,exCentroids=NULL,n.neighbors=4,batchPCAdata=NULL,colNormalize=F,...){
  
  if(is.null(batchPCAdata)){
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    batchPCAdata=harmony_embeddings
  }
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  #batchPCAdata=harmony_embeddings
  batchPCAdata=as.matrix(batchPCAdata)[,1:argList$nPCs]
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  
  nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
    if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0)>0){
      stop("error!")
    }
    adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = (nn.ranked.1$nn.dists[,i]+0.01)
  }
  
  
  
  
  #######################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_adj=adj
    tst_adj[tst_adj>0]=1
    
    tst_adj=tst_adj %*% as.matrix(tst)
    
    .tmp_step1=tst_adj
  }
  
  #######################
  
  if(F){
    #old way: it doesn't work when pseuodcells are too sparse!
    nn.centroid.purityIndex=apply(adj,1,function(x){
      
      if(sum(x>0,na.rm = T)>0){
        x=x[which(x>0)]
        x=median(x)
      } else {
        x=0
      }
      
      return(x)
    })
  } else {
    tmp_df=data.frame(pseudocells=row.names(centroidPCAdata)[nn.ranked.1$nn.idx[,1]],dist=nn.ranked.1$nn.dists[,1],stringsAsFactors = F)
    tmp_df=aggregate(dist~pseudocells,data=tmp_df,function(x) median(x,na.rm = T))
    nn.centroid.purityIndex=tmp_df$dist
    names(nn.centroid.purityIndex)=as.character(tmp_df$pseudocells)
    
    
    nn.centroid.purityIndex=nn.centroid.purityIndex[match(row.names(centroidPCAdata),names(nn.centroid.purityIndex))]
    if(sum(is.na(nn.centroid.purityIndex))>0){
      nn.centroid.purityIndex[is.na(nn.centroid.purityIndex)]=0
    }
    tmp_df=quantile(nn.centroid.purityIndex,0.95)
    nn.centroid.purityIndex[nn.centroid.purityIndex>tmp_df]=tmp_df
    rm(tmp_df)
  }
  
  
  
  affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinities[which(adj==0)]=0
  centroid.affinities=affinities
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  
  if(F){
    dataArranged2=dataArranged
    for(i in 1:length(dataArranged2)){
      tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged2[[i]]$pcaData),colnames(prop_mat))],data=dataArranged2[[i]],pseudocell_sim_mat=prop_mat)
      dataArranged2[[i]]=tmp
    }
    p=.sconline.anno2pseudocell_tmp(res_prop=dataArranged2,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T)
    ggsave(plot=p$plot,file="~/myBucket/torm3.pdf",width=12,height = 12)
    
  }
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  prop_mat=centroid.affinities
  
  
  #################################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_affinities=prop_mat %*% as.matrix(tst)
    
    #apply(tst_affinities,2,mean)
    apply(tst_affinities,2,function(x) sum(x>0.5))
    
    
    tst_pd=pd[match(colnames(prop_mat),row.names(pd)),]
    aggregate(prop_mat["8",]~tst_pd$anno_orig_cellState,sum,data=NULL)
    
    eff_size=.myEffSizePropMat(prop_mat)
    
    x=apply(tst_affinities,1,function(x) which(x==max(x))[1])
    x=data.frame(size=eff_size$effective_sample_size,g=x)
    aggregate(size~g,data=x,summary)
    
  }
  
  
  
  ###############################
  
  c_c_aff=t(prop_mat)
  c_c_aff=prop_mat %*% c_c_aff
  c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
  diag(c_c_aff)=1
  c_c_aff=pmin(c_c_aff,1)
  c_c_aff=sqrt(c_c_aff*t(c_c_aff))
  #c_c_aff[which(c_c_aff<0.1)]=0
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff
  
  prop_mat=c_c_aff %*% prop_mat
  
  
  if(F){
    prop_mat2=t(apply(prop_mat,1,function(x) x/max(x)))
    prop_mat[which(prop_mat2<0.1)]=0
    
  } else {
    cat("Column-wise normalization")
    tst=c()
    for(i in seq(1,ncol(prop_mat),10000)){
      tst1=qlcMatrix::colMax(prop_mat[,i:min(i+9999,ncol(prop_mat))])
      tst=c(tst,tst1)
    }
    tst=unlist(lapply(tst,function(x) unlist(as.numeric(x))))
    
    prop_mat <- .extra_matrix_colNorm(input_mat = prop_mat,colValues = 1 / (colSums(tst)+0.000000000001))#prop_mat %*% Matrix::Diagonal(x = 1 / (colSums(tst)+0.000000000001))
    prop_mat=drop0(prop_mat,tol=0.1)
    
  }
  
  
  
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  matWeights=.myEffSizePropMat(prop_mat)
  
  matEffectiveSize=matWeights$effective_sample_size
  matWeights=matWeights$centroid_weights
  matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  
  
  if(F){
    cat("    Removing redundant pseudocells")
    tst=prop_mat
    #tst[tst>0]=1
    tst=sqrt(tst %*% t(tst))
    tst=sweep(tst,1,diag(tst),"/")
    cosine_sim=(tst+t(tst))/2
    #cosine_sim=proxyC::simil(prop_mat,margin = 1,method = "cosine")
    
    dataArranged2=dataArranged
    for(i in 1:length(dataArranged2)){
      tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged2[[i]]$pcaData),colnames(prop_mat))],data=dataArranged2[[i]],pseudocell_sim_mat=prop_mat)
      dataArranged2[[i]]=tmp
    }
    
    #res_prop=dataArranged2;annoCol="anno_orig_cellState";collapse_datasets=T;return_plot=T
    p=.sconline.anno2pseudocell_tmp(res_prop=dataArranged2,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T)
    
    
    summary(apply(cosine_sim,1,function(x) sum(x>0.9)))
    which(apply(cosine_sim,1,function(x) sum(x>0.9))>10)
    
    tst=cosine_sim[row.names(p$results),row.names(p$results)]
    summary(apply(tst,1,function(x) sum(x>0.9)))
    which(apply(tst,1,function(x) sum(x>0.9))>10)
    tst=p$results[which(tst[13,]>0.9),]
    tst2=apply(tst[,-c(1,ncol(tst))],2,max)
    tst=tst[,c(1,which(tst2>0.1),ncol(tst))]
    
    if(sum(upper.tri(cosine_sim)&cosine_sim>0.95)>0){
      ident_list=which(upper.tri(cosine_sim)&cosine_sim>0.95,arr.ind = T)
      ident_list=igraph::components(igraph::graph_from_data_frame(as.data.frame(ident_list)))
      ident_list=data.frame(group=ident_list$membership,pseudocell=names(ident_list$membership),stringsAsFactors = F)
      ident_list=split(ident_list$pseudocell,ident_list$group)
      torm_list=c()
      for(i in 1:length(ident_list)){
        x=prop_mat[ident_list[[i]],]
        x=apply(x,2,sum)
        prop_mat[ident_list[[i]][1],]=x/(sum(x)+0.000000000001)
        torm_list=c(torm_list,ident_list[[i]][2:length(ident_list[[i]])])
      }
      prop_mat=prop_mat[-which(row.names(prop_mat) %in% torm_list),]
    }
    
  }
  
  for(i in 1:length(dataArranged)){
    tmp_prop_mat=prop_mat[,match(row.names(dataArranged[[i]]$pcaData),colnames(prop_mat))]
    #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
    tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
    tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
    tmp_weights[tmp_effsize<4]=0
    tmp=list(prop_mat=tmp_prop_mat,data=dataArranged[[i]],pseudocell_sim_mat=c_c_aff,matWeights=tmp_weights,matEffectiveSize=tmp_effsize)
    dataArranged[[i]]=tmp
  }
  
  return(dataArranged)
}

.myConcensusDEFn_step2_detail_newprop3_final_v8=function(dataArranged,centroidPCAdata,argList,exCentroids=NULL,n.neighbors=4,batchPCAdata=NULL,colNormalize=F,...){
  
  if(is.null(batchPCAdata)){
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    batchPCAdata=harmony_embeddings
  }
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  #batchPCAdata=harmony_embeddings
  batchPCAdata=as.matrix(batchPCAdata)[,1:argList$nPCs]
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  
  nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
    if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0)>0){
      stop("error!")
    }
    adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = (nn.ranked.1$nn.dists[,i]+0.01)
  }
  
  
  
  
  #######################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_adj=adj
    tst_adj[tst_adj>0]=1
    
    tst_adj=tst_adj %*% as.matrix(tst)
    
    .tmp_step1=tst_adj
  }
  
  #######################
  
  if(F){
    #old way: it doesn't work when pseuodcells are too sparse!
    nn.centroid.purityIndex=apply(adj,1,function(x){
      
      if(sum(x>0,na.rm = T)>0){
        x=x[which(x>0)]
        x=median(x)
      } else {
        x=0
      }
      
      return(x)
    })
  } else {
    tmp_df=data.frame(pseudocells=row.names(centroidPCAdata)[nn.ranked.1$nn.idx[,1]],dist=nn.ranked.1$nn.dists[,1],stringsAsFactors = F)
    tmp_df=aggregate(dist~pseudocells,data=tmp_df,function(x) median(x,na.rm = T))
    nn.centroid.purityIndex=tmp_df$dist
    names(nn.centroid.purityIndex)=as.character(tmp_df$pseudocells)
    
    
    nn.centroid.purityIndex=nn.centroid.purityIndex[match(row.names(centroidPCAdata),names(nn.centroid.purityIndex))]
    if(sum(is.na(nn.centroid.purityIndex))>0){
      nn.centroid.purityIndex[is.na(nn.centroid.purityIndex)]=0
    }
    tmp_df=quantile(nn.centroid.purityIndex,0.95)
    nn.centroid.purityIndex[nn.centroid.purityIndex>tmp_df]=tmp_df
    rm(tmp_df)
  }
  
  
  
  affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinities[which(adj==0)]=0
  centroid.affinities=affinities
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  
  if(F){
    dataArranged2=dataArranged
    for(i in 1:length(dataArranged2)){
      tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged2[[i]]$pcaData),colnames(prop_mat))],data=dataArranged2[[i]],pseudocell_sim_mat=prop_mat)
      dataArranged2[[i]]=tmp
    }
    p=.sconline.anno2pseudocell_tmp(res_prop=dataArranged2,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T)
    ggsave(plot=p$plot,file="~/myBucket/torm3.pdf",width=12,height = 12)
    
  }
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  prop_mat=centroid.affinities
  
  
  #################################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_affinities=prop_mat %*% as.matrix(tst)
    
    #apply(tst_affinities,2,mean)
    apply(tst_affinities,2,function(x) sum(x>0.5))
    
    
    tst_pd=pd[match(colnames(prop_mat),row.names(pd)),]
    aggregate(prop_mat["8",]~tst_pd$anno_orig_cellState,sum,data=NULL)
    
    eff_size=.myEffSizePropMat(prop_mat)
    
    x=apply(tst_affinities,1,function(x) which(x==max(x))[1])
    x=data.frame(size=eff_size$effective_sample_size,g=x)
    aggregate(size~g,data=x,summary)
    
  }
  
  
  
  ###############################
  
  c_c_aff=t(prop_mat)
  c_c_aff=prop_mat %*% c_c_aff
  c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
  diag(c_c_aff)=1
  c_c_aff@x=pmin(c_c_aff@x,1)
  c_c_aff=sqrt(c_c_aff*t(c_c_aff))
  #c_c_aff[which(c_c_aff<0.1)]=0
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff
  
  prop_mat=c_c_aff %*% prop_mat
  
  
  if(F){
    prop_mat2=t(apply(prop_mat,1,function(x) x/max(x)))
    prop_mat[which(prop_mat2<0.1)]=0
    
  } else {
    rowMax_vals=c()
    for(i in seq(1,nrow(prop_mat),1000)){
      tmp=qlcMatrix::rowMax(prop_mat[i:min(i+999,nrow(prop_mat)),])
      rowMax_vals=c(rowMax_vals,tmp)
    }
    rowMax_vals=unlist(lapply(rowMax_vals,as.numeric))
    
    prop_mat2=.extra_matrix_rowNorm(input_mat = prop_mat,rowValues = 1/rowMax_vals)# Matrix::Diagonal(x = 1 / rowMax_vals) %*% prop_mat
    prop_mat2=Matrix::drop0(prop_mat2,tol=0.1)
    prop_mat2@x=rep(1,length(prop_mat2@x))
    
    if(colNormalize){
      cat("Column-wise normalization")
      tst=c()
      for(i in seq(1,ncol(prop_mat),10000)){
        tst1=qlcMatrix::colMax(prop_mat[,i:min(i+9999,ncol(prop_mat))])
        tst=c(tst,tst1)
      }
      tst=unlist(lapply(tst,function(x) unlist(as.numeric(x))))
      
      prop_mat3 <- .extra_matrix_colNorm(input_mat = prop_mat,colValues = 1/(tst+0.0000000000001))# prop_mat %*% Matrix::Diagonal(x = 1 / (tst+0.000000000001))
      prop_mat3=drop0(prop_mat3,tol=0.1)
      
      prop_mat4=prop_mat3+prop_mat2
      prop_mat4@x=1
      prop_mat=prop_mat * prop_mat4
    } else {
      prop_mat=prop_mat * prop_mat2
    }
    
    
    rm(prop_mat2)
  }
  
  
  
  
  
  adj_all=Seurat::FindNeighbors(batchPCAdata,compute.SNN=T)
  adj_all=adj_all[["snn"]]
  adj_all=adj_all[,colnames(prop_mat)]
  adj_all=adj_all %*% t(prop_mat)
  prop_mat=t(adj_all)
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  matWeights=.myEffSizePropMat(prop_mat)
  
  matEffectiveSize=matWeights$effective_sample_size
  matWeights=matWeights$centroid_weights
  matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  
  
  if(F){
    cat("    Removing redundant pseudocells")
    tst=prop_mat
    #tst[tst>0]=1
    tst=sqrt(tst %*% t(tst))
    tst=sweep(tst,1,diag(tst),"/")
    cosine_sim=(tst+t(tst))/2
    #cosine_sim=proxyC::simil(prop_mat,margin = 1,method = "cosine")
    
    dataArranged2=dataArranged
    for(i in 1:length(dataArranged2)){
      tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged2[[i]]$pcaData),colnames(prop_mat))],data=dataArranged2[[i]],pseudocell_sim_mat=prop_mat)
      dataArranged2[[i]]=tmp
    }
    
    #res_prop=dataArranged2;annoCol="anno_orig_cellState";collapse_datasets=T;return_plot=T
    p=.sconline.anno2pseudocell_tmp(res_prop=dataArranged2,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T)
    
    
    summary(apply(cosine_sim,1,function(x) sum(x>0.9)))
    which(apply(cosine_sim,1,function(x) sum(x>0.9))>10)
    
    tst=cosine_sim[row.names(p$results),row.names(p$results)]
    summary(apply(tst,1,function(x) sum(x>0.9)))
    which(apply(tst,1,function(x) sum(x>0.9))>10)
    tst=p$results[which(tst[13,]>0.9),]
    tst2=apply(tst[,-c(1,ncol(tst))],2,max)
    tst=tst[,c(1,which(tst2>0.1),ncol(tst))]
    
    if(sum(upper.tri(cosine_sim)&cosine_sim>0.95)>0){
      ident_list=which(upper.tri(cosine_sim)&cosine_sim>0.95,arr.ind = T)
      ident_list=igraph::components(igraph::graph_from_data_frame(as.data.frame(ident_list)))
      ident_list=data.frame(group=ident_list$membership,pseudocell=names(ident_list$membership),stringsAsFactors = F)
      ident_list=split(ident_list$pseudocell,ident_list$group)
      torm_list=c()
      for(i in 1:length(ident_list)){
        x=prop_mat[ident_list[[i]],]
        x=apply(x,2,sum)
        prop_mat[ident_list[[i]][1],]=x/(sum(x)+0.000000000001)
        torm_list=c(torm_list,ident_list[[i]][2:length(ident_list[[i]])])
      }
      prop_mat=prop_mat[-which(row.names(prop_mat) %in% torm_list),]
    }
    
  }
  
  for(i in 1:length(dataArranged)){
    tmp_prop_mat=prop_mat[,match(row.names(dataArranged[[i]]$pcaData),colnames(prop_mat))]
    #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
    tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
    tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
    tmp_weights[tmp_effsize<4]=0
    tmp=list(prop_mat=tmp_prop_mat,data=dataArranged[[i]],pseudocell_sim_mat=c_c_aff,matWeights=tmp_weights,matEffectiveSize=tmp_effsize)
    dataArranged[[i]]=tmp
  }
  
  return(dataArranged)
}

.myConcensusDEFn_step2_detail_newprop3_final_v9=function(dataArranged,centroidPCAdata,argList,exCentroids=NULL,n.neighbors=4,batchPCAdata=NULL,colNormalize=F,...){
  
  if(is.null(batchPCAdata)){
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    batchPCAdata=harmony_embeddings
  }
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  #batchPCAdata=harmony_embeddings
  batchPCAdata=as.matrix(batchPCAdata)[,1:argList$nPCs]
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  
  nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
    if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0)>0){
      stop("error!")
    }
    adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = (nn.ranked.1$nn.dists[,i]+0.01)
  }
  
  
  
  
  #######################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_adj=adj
    tst_adj[tst_adj>0]=1
    
    tst_adj=tst_adj %*% as.matrix(tst)
    
    .tmp_step1=tst_adj
  }
  
  #######################
  
  if(F){
    #old way: it doesn't work when pseuodcells are too sparse!
    nn.centroid.purityIndex=apply(adj,1,function(x){
      
      if(sum(x>0,na.rm = T)>0){
        x=x[which(x>0)]
        x=median(x)
      } else {
        x=0
      }
      
      return(x)
    })
  } else {
    tmp_df=data.frame(pseudocells=row.names(centroidPCAdata)[nn.ranked.1$nn.idx[,1]],dist=nn.ranked.1$nn.dists[,1],stringsAsFactors = F)
    tmp_df=aggregate(dist~pseudocells,data=tmp_df,function(x) median(x,na.rm = T))
    nn.centroid.purityIndex=tmp_df$dist
    names(nn.centroid.purityIndex)=as.character(tmp_df$pseudocells)
    
    
    nn.centroid.purityIndex=nn.centroid.purityIndex[match(row.names(centroidPCAdata),names(nn.centroid.purityIndex))]
    if(sum(is.na(nn.centroid.purityIndex))>0){
      nn.centroid.purityIndex[is.na(nn.centroid.purityIndex)]=0
    }
    tmp_df=quantile(nn.centroid.purityIndex,0.95)
    nn.centroid.purityIndex[nn.centroid.purityIndex>tmp_df]=tmp_df
    rm(tmp_df)
  }
  
  
  
  affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinities[which(adj==0)]=0
  centroid.affinities=affinities
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  
  if(F){
    dataArranged2=dataArranged
    for(i in 1:length(dataArranged2)){
      tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged2[[i]]$pcaData),colnames(prop_mat))],data=dataArranged2[[i]],pseudocell_sim_mat=prop_mat)
      dataArranged2[[i]]=tmp
    }
    p=.sconline.anno2pseudocell_tmp(res_prop=dataArranged2,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T)
    ggsave(plot=p$plot,file="~/myBucket/torm3.pdf",width=12,height = 12)
    
  }
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  prop_mat=centroid.affinities
  
  
  #################################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_affinities=prop_mat %*% as.matrix(tst)
    
    #apply(tst_affinities,2,mean)
    apply(tst_affinities,2,function(x) sum(x>0.5))
    
    
    tst_pd=pd[match(colnames(prop_mat),row.names(pd)),]
    aggregate(prop_mat["8",]~tst_pd$anno_orig_cellState,sum,data=NULL)
    
    eff_size=.myEffSizePropMat(prop_mat)
    
    x=apply(tst_affinities,1,function(x) which(x==max(x))[1])
    x=data.frame(size=eff_size$effective_sample_size,g=x)
    aggregate(size~g,data=x,summary)
    
  }
  
  
  
  ###############################
  
  c_c_aff=t(prop_mat)
  c_c_aff=prop_mat %*% c_c_aff
  c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
  diag(c_c_aff)=1
  c_c_aff@x=pmin(c_c_aff@x,1)
  c_c_aff=sqrt(c_c_aff*t(c_c_aff))
  #c_c_aff[which(c_c_aff<0.1)]=0
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff
  
  prop_mat=c_c_aff %*% prop_mat
  
  
  if(F){
    prop_mat2=t(apply(prop_mat,1,function(x) x/max(x)))
    prop_mat[which(prop_mat2<0.1)]=0
    
  } else {
    rowMax_vals=c()
    for(i in seq(1,nrow(prop_mat),1000)){
      tmp=qlcMatrix::rowMax(prop_mat[i:min(i+999,nrow(prop_mat)),])
      rowMax_vals=c(rowMax_vals,tmp)
    }
    rowMax_vals=unlist(lapply(rowMax_vals,as.numeric))
    
    prop_mat2=.extra_matrix_rowNorm(input_mat = prop_mat,rowValues = 1/rowMax_vals)#Matrix::Diagonal(x = 1 / rowMax_vals) %*% prop_mat
    prop_mat2=Matrix::drop0(prop_mat2,tol=0.1)
    prop_mat2@x=rep(1,length(prop_mat2@x))
    
    prop_mat=prop_mat * prop_mat2
    
    rm(prop_mat2)
  }
  
  
  
  
  
  adj_all=Seurat::FindNeighbors(batchPCAdata,compute.SNN=T)
  adj_all=adj_all[["snn"]]
  adj_all=adj_all[,colnames(prop_mat)]
  prop_mat <- .extra_matrix_colNorm(input_mat = prop_mat,colValues = 1/(colSums(prop_mat)+0.000000000001))# prop_mat %*% Matrix::Diagonal(x = 1 / (colSums(prop_mat)+0.000000000001))
  adj_all=adj_all %*% t(prop_mat)
  prop_mat=t(adj_all)
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  matWeights=.myEffSizePropMat(prop_mat)
  
  matEffectiveSize=matWeights$effective_sample_size
  matWeights=matWeights$centroid_weights
  matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  
  
  if(F){
    cat("    Removing redundant pseudocells")
    tst=prop_mat
    #tst[tst>0]=1
    tst=sqrt(tst %*% t(tst))
    tst=sweep(tst,1,diag(tst),"/")
    cosine_sim=(tst+t(tst))/2
    #cosine_sim=proxyC::simil(prop_mat,margin = 1,method = "cosine")
    
    dataArranged2=dataArranged
    for(i in 1:length(dataArranged2)){
      tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged2[[i]]$pcaData),colnames(prop_mat))],data=dataArranged2[[i]],pseudocell_sim_mat=prop_mat)
      dataArranged2[[i]]=tmp
    }
    
    #res_prop=dataArranged2;annoCol="anno_orig_cellState";collapse_datasets=T;return_plot=T
    p=.sconline.anno2pseudocell_tmp(res_prop=dataArranged2,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T)
    
    
    summary(apply(cosine_sim,1,function(x) sum(x>0.9)))
    which(apply(cosine_sim,1,function(x) sum(x>0.9))>10)
    
    tst=cosine_sim[row.names(p$results),row.names(p$results)]
    summary(apply(tst,1,function(x) sum(x>0.9)))
    which(apply(tst,1,function(x) sum(x>0.9))>10)
    tst=p$results[which(tst[13,]>0.9),]
    tst2=apply(tst[,-c(1,ncol(tst))],2,max)
    tst=tst[,c(1,which(tst2>0.1),ncol(tst))]
    
    if(sum(upper.tri(cosine_sim)&cosine_sim>0.95)>0){
      ident_list=which(upper.tri(cosine_sim)&cosine_sim>0.95,arr.ind = T)
      ident_list=igraph::components(igraph::graph_from_data_frame(as.data.frame(ident_list)))
      ident_list=data.frame(group=ident_list$membership,pseudocell=names(ident_list$membership),stringsAsFactors = F)
      ident_list=split(ident_list$pseudocell,ident_list$group)
      torm_list=c()
      for(i in 1:length(ident_list)){
        x=prop_mat[ident_list[[i]],]
        x=apply(x,2,sum)
        prop_mat[ident_list[[i]][1],]=x/(sum(x)+0.000000000001)
        torm_list=c(torm_list,ident_list[[i]][2:length(ident_list[[i]])])
      }
      prop_mat=prop_mat[-which(row.names(prop_mat) %in% torm_list),]
    }
    
  }
  
  for(i in 1:length(dataArranged)){
    tmp_prop_mat=prop_mat[,match(row.names(dataArranged[[i]]$pcaData),colnames(prop_mat))]
    #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
    tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
    tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
    tmp_weights[tmp_effsize<4]=0
    tmp=list(prop_mat=tmp_prop_mat,data=dataArranged[[i]],pseudocell_sim_mat=c_c_aff,matWeights=tmp_weights,matEffectiveSize=tmp_effsize)
    dataArranged[[i]]=tmp
  }
  
  return(dataArranged)
}

.myConcensusDEFn_step2_detail_newprop3_final_v10=function(dataArranged,centroidPCAdata,argList,exCentroids=NULL,n.neighbors=4,batchPCAdata=NULL,n.trees=50,NNmethod="annoy",L2Norm=T,returnPropMat=F,...){
  
  myL2normFn=function(inputMat){
    prop_mat2=rowSums(inputMat^2)
    prop_mat2=sqrt(prop_mat2)
    res=.extra_matrix_rowNorm(input_mat = inputMat,rowValues = 1/(prop_mat2+0.00000001))# Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
    return(res)
  }
  
  if(is.null(batchPCAdata)){
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    batchPCAdata=harmony_embeddings
  }
  
  if(sum(argList$singleton.method %in% c("snn","fast"))==0){
    stop("unrecognized singleton.method")
  }
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  #batchPCAdata=harmony_embeddings
  batchPCAdata=as.matrix(batchPCAdata)[,1:argList$nPCs]
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  if(NNmethod=="annoy"){
    idx=Seurat:::AnnoyBuildIndex(data = centroidPCAdata, metric = "euclidean", 
                                 n.trees = n.trees)
    nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = batchPCAdata,k=n.neighbors,include.distance = T,search.k = -1)
  } else {
    nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  }
  
  adj= nn.ranked.1$nn.idx
  j <- as.numeric(t(adj))
  i <- ((1:length(j)) - 1)%/%n.neighbors + 1
  x=as.numeric(t(nn.ranked.1$nn.dists+0.01))
  adj = t(sparseMatrix(i = i, j = j, x = x, dims = c(nrow(batchPCAdata),nrow(centroidPCAdata))))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  
  adj=adj[rowSums(adj)>0.01,]
  
  adj2=adj
  adj2@x=rep(1,length(adj2@x))
  adj=adj[rowSums(adj2)>4,]
  
  tmp_df=data.frame(pseudocells=row.names(centroidPCAdata)[nn.ranked.1$nn.idx[,1]],dist=nn.ranked.1$nn.dists[,1],stringsAsFactors = F)
  tmp_df=aggregate(dist~pseudocells,data=tmp_df,function(x) median(x,na.rm = T))
  nn.centroid.purityIndex=tmp_df$dist
  names(nn.centroid.purityIndex)=as.character(tmp_df$pseudocells)
  
  nn.centroid.purityIndex=nn.centroid.purityIndex[match(row.names(centroidPCAdata),names(nn.centroid.purityIndex))]
  if(sum(is.na(nn.centroid.purityIndex))>0){
    nn.centroid.purityIndex[is.na(nn.centroid.purityIndex)]=0
  }
  tmp_df=quantile(nn.centroid.purityIndex,0.95)
  nn.centroid.purityIndex[nn.centroid.purityIndex>tmp_df]=tmp_df
  rm(tmp_df)
  nn.centroid.purityIndex=nn.centroid.purityIndex[row.names(adj)]
  
  affinities=.extra_matrix_rowNorm(input_mat = adj,rowValues = 1/(nn.centroid.purityIndex+0.000001))#Matrix::Diagonal(x=1/(nn.centroid.purityIndex+0.000001)) %*% adj
  affinities@x=-1*affinities@x^2
  affinities@x=exp(affinities@x)
  
  #adj@x=rep(1,length(adj@x))
  #affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  #affinities=do.call("rbind",affinities)
  #affinities=affinities*adj
  
  centroid.affinities <- affinities #%*% Matrix::Diagonal(x = 1 / (colSums(affinities)+0.000000000001))
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  
  if(F){
    dataArranged2=dataArranged
    for(i in 1:length(dataArranged2)){
      tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged2[[i]]$pcaData),colnames(prop_mat))],data=dataArranged2[[i]],pseudocell_sim_mat=prop_mat)
      dataArranged2[[i]]=tmp
    }
    p=.sconline.anno2pseudocell_tmp(res_prop=dataArranged2,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T)
    ggsave(plot=p$plot,file="~/myBucket/torm3.pdf",width=12,height = 12)
    
  }
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(centroid.affinities) <- row.names(adj)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  prop_mat=centroid.affinities
  
  
  #################################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_affinities=prop_mat %*% as.matrix(tst)
    
    #apply(tst_affinities,2,mean)
    apply(tst_affinities,2,function(x) sum(x>0.5))
    
    
    tst_pd=pd[match(colnames(prop_mat),row.names(pd)),]
    aggregate(prop_mat["8",]~tst_pd$anno_orig_cellState,sum,data=NULL)
    
    eff_size=.myEffSizePropMat(prop_mat)
    
    x=apply(tst_affinities,1,function(x) which(x==max(x))[1])
    x=data.frame(size=eff_size$effective_sample_size,g=x)
    aggregate(size~g,data=x,summary)
    
  }
  
  
  
  ###############################
  
  if(L2Norm){
    prop_mat2=myL2normFn(inputMat=prop_mat)
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
  } else {
    prop_mat2=prop_mat
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
    c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
    diag(c_c_aff)=1
    c_c_aff@x=pmin(c_c_aff@x,1)
    c_c_aff=sqrt(c_c_aff*t(c_c_aff))
  }
  
  
  
  #c_c_aff[which(c_c_aff<0.1)]=0
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff
  
  
  prop_mat_list=list()
  rowMax_vals=c()
  for(i in seq(1,ncol(prop_mat),50000)){
    tmp=c_c_aff %*% prop_mat[,i:min(i+49999,ncol(prop_mat))]
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
    tmp2=.extra_matrix_rowNorm(input_mat = tmp,rowValues = 1/rowMax_vals)# Matrix::Diagonal(x = 1 / rowMax_vals) %*% tmp
    tmp2=Matrix::drop0(tmp2,tol=0.1)
    tmp2@x=rep(1,length(tmp2@x))
    
    tmp=tmp*tmp2
    tmp=Matrix::drop0(tmp)
    prop_mat_list[[i]]=tmp
  }
  
  #prop_mat=c_c_aff %*% prop_mat
  
  if(F){
    rowMax_vals=c()
    for(i in seq(1,nrow(prop_mat),1000)){
      tmp=qlcMatrix::rowMax(prop_mat[i:min(i+999,nrow(prop_mat)),])
      rowMax_vals=c(rowMax_vals,tmp)
    }
    rowMax_vals=unlist(lapply(rowMax_vals,as.numeric))
    
    prop_mat2=.extra_matrix_rowNorm(input_mat = prop_mat,rowValues = 1/rowMax_vals)#Matrix::Diagonal(x = 1 / rowMax_vals) %*% prop_mat
    prop_mat2=Matrix::drop0(prop_mat2,tol=0.1)
    prop_mat2@x=rep(1,length(prop_mat2@x))
    
    prop_mat=prop_mat * prop_mat2
    rm(prop_mat2)
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
    rm(tmp)
  }
  
  rm(prop_mat_list)
  gc()
  
  if(F){
    tmp_pd=pd[match(colnames(prop_mat),row.names(pd)),]
    anno=prop_mat %*% as.matrix(.myOneHotFn(tmp_pd$anno_orig_cellState))
    anno_group=apply(anno,1,function(x) which(x==max(x))[1])
    anno_group=colnames(anno)[anno_group]
    anno[which(anno_group=="Macrophages"),]
    
    anno[which(anno_group=="Microglia"),]
  }
  
  
  matEffectiveSize=.myEffSizePropMat(prop_mat)
  matEffectiveSize=matEffectiveSize$effective_sample_size
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  prop_mat=prop_mat[matEffectiveSize>=argList$min_cluster_size,]
  prop_mat=prop_mat[,colSums(prop_mat)>0]
  
  connected_cells=colnames(prop_mat)
  if(argList$singleton.method=="snn"){
    cat("Handling singletons by SNN\n")
    nn.ranked <- Seurat:::NNHelper(query = batchPCAdata[c(connected_cells,setdiff(row.names(batchPCAdata),connected_cells)),], data = batchPCAdata[connected_cells,], k = 20, 
                                   method = "annoy", n.trees = 50, searchtype = "standard", 
                                   eps = 0, metric = "euclidean", cache.index = F, 
                                   index = NULL)
    nn.ranked = Indices(nn.ranked)
    
    snn.matrix = Seurat:::ComputeSNN(nn_ranked = nn.ranked, prune = 1/15)
    rownames(snn.matrix) = c(connected_cells,setdiff(row.names(batchPCAdata),connected_cells))
    colnames(snn.matrix) = c(connected_cells,setdiff(row.names(batchPCAdata),connected_cells))
    
    tst=snn.matrix
    diag(tst)=0
    tst=tst[setdiff(row.names(batchPCAdata),connected_cells),]
    tst=as.numeric(qlcMatrix::rowMax(tst))
    tst=c(rep(1,length(connected_cells)),tst)
    snn.matrix=.extra_matrix_rowNorm(input_mat = snn.matrix,rowValues = 1/(tst+0.000001))#Matrix::Diagonal(x=1/(tst+0.000001))%*% snn.matrix
    diag(snn.matrix)=1
    adj_all=snn.matrix[,colnames(prop_mat)]
  } else if(argList$singleton.method=="fast") {
    cat("Handling singletons by fast method\n")
    if(length(setdiff(row.names(batchPCAdata),connected_cells))>0){
      nn.ranked <- Seurat:::NNHelper(query = batchPCAdata[setdiff(row.names(batchPCAdata),connected_cells),], data = batchPCAdata[connected_cells,], k = 1, 
                                     method = "annoy", n.trees = 50, searchtype = "standard", 
                                     eps = 0, metric = "euclidean", cache.index = F, 
                                     index = NULL)
      nn.ranked = Indices(nn.ranked)
      nn.ranked=connected_cells[nn.ranked[,1]]
      nn.ranked=data.frame(sample=c(connected_cells,setdiff(row.names(batchPCAdata),connected_cells)),map=c(connected_cells,nn.ranked),stringsAsFactors = F)
    } else {
      nn.ranked=data.frame(sample=colnames(prop_mat),map=colnames(prop_mat),stringsAsFactors = F)
    }
    
    nn.ranked=merge(nn.ranked,data.frame(map=colnames(prop_mat),id=1:ncol(prop_mat),stringsAsFactors = F),by="map")
    #nn.ranked=nn.ranked[match(colnames(prop_mat),nn.ranked$sample),]
    adj_all=Matrix::sparseMatrix(i=1:nrow(nn.ranked),j=nn.ranked$id,x=rep(1,nrow(nn.ranked)))
    row.names(adj_all)=nn.ranked$sample
    colnames(adj_all)=colnames(prop_mat)
    rm(nn.ranked)
  } else {
    stop("unrecognized singleton.method")
  }
  
  
  ###############
  #New
  #Excluded the col normalization step
  #prop_mat <- prop_mat %*% Matrix::Diagonal(x = 1 / (colSums(prop_mat)+0.000000000001))
  adj_all=adj_all %*% t(prop_mat)
  #adj_all = Matrix::Diagonal(x = 1 / (rowSums(adj_all)+0.000000000001)) %*% adj_all
  prop_mat=t(adj_all)
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  ##############
  if(L2Norm){
    prop_mat2=myL2normFn(inputMat = prop_mat)
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
  } else {
    prop_mat2=prop_mat
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
    c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
    diag(c_c_aff)=1
    c_c_aff@x=pmin(c_c_aff@x,1)
    c_c_aff=sqrt(c_c_aff*t(c_c_aff))
  }
  
  
  
  ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
  ps_merging=cutree(ps_merging,h=0.15)
  while(length(unique(ps_merging))<length(ps_merging)){
    ps_merging=data.frame(pseudocell=row.names(prop_mat),id=as.character(ps_merging),stringsAsFactors = F)
    ps_mid=ps_merging[!duplicated(ps_merging$id),]
    colnames(ps_mid)[1]="mid"
    ps_merging=merge(ps_merging,ps_mid,by="id")
    ps_merging=ps_merging[match(row.names(prop_mat),ps_merging$pseudocell),]
    ps_merging=as.matrix(.myOneHotFn(as.character(ps_merging$mid)))
    ps_merging=as(ps_merging,"dgCMatrix")
    ps_merging=t(ps_merging)
    colnames(ps_merging)=row.names(prop_mat)
    prop_mat=ps_merging %*% prop_mat
    rm(ps_merging)
    prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
    
    
    if(L2Norm){
      prop_mat2=myL2normFn(inputMat = prop_mat)
      c_c_aff=t(prop_mat2)
      c_c_aff=prop_mat2 %*% c_c_aff
    } else {
      prop_mat2=prop_mat
      c_c_aff=t(prop_mat2)
      c_c_aff=prop_mat2 %*% c_c_aff
      c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
      diag(c_c_aff)=1
      c_c_aff@x=pmin(c_c_aff@x,1)
      c_c_aff=sqrt(c_c_aff*t(c_c_aff))
    }
    
    ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
    ps_merging=cutree(ps_merging,h=0.15)
  }
  
  if(F){
    tmp_pd=pd[match(colnames(prop_mat),row.names(pd)),]
    anno=prop_mat %*% as.matrix(.myOneHotFn(tmp_pd$anno_orig_cellState))
    anno_group=apply(anno,1,function(x) which(x==max(x)))
    anno_group=colnames(anno)[anno_group]
    
    
    res_check=NULL
    for(i in ps_merging[duplicated(ps_merging)]){
      ii=which(ps_merging==i)
      for(j in 1:(length(ii)-1)){
        for(k in 2:length(ii)){
          Fscore=anno[ii[j],anno_group[ii[j]]]/anno[ii[j],anno_group[ii[k]]]
          Sscore=anno[ii[k],anno_group[ii[k]]]/anno[ii[k],anno_group[ii[j]]]
          res_check=rbind(res_check,data.frame(Fnode=anno_group[ii[j]],Snode=anno_group[ii[k]],Find=ii[j],Sind=ii[k],Fscore=Fscore,Sscore=Sscore,stringsAsFactors = F))
        }
      }
    }
    table(res_check$Fnode==res_check$Snode)
    length(unique(ps_merging))
    res_check=res_check[res_check$Fnode!=res_check$Snode,]
    table(res_check$Fscore>1.5|res_check$Sscore>1.5)
    res_check[res_check$Fscore>1.5|res_check$Sscore>1.5,]
  }
  
  ##################
  if(returnPropMat){
    return(prop_mat)
  }
  
  matWeights=.myEffSizePropMat(prop_mat)
  
  matEffectiveSize=matWeights$effective_sample_size
  matWeights=matWeights$centroid_weights
  matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  
  
  for(i in 1:length(dataArranged)){
    tmp_prop_mat=prop_mat[,match(row.names(dataArranged[[i]]$pcaData),colnames(prop_mat))]
    #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
    tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
    tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
    tmp_weights[tmp_effsize<4]=0
    tmp=list(prop_mat=tmp_prop_mat,data=dataArranged[[i]],pseudocell_sim_mat=c_c_aff,matWeights=tmp_weights,matEffectiveSize=tmp_effsize)
    dataArranged[[i]]=tmp
  }
  
  return(dataArranged)
}

.myConcensusDEFn_step2_detail_newprop3_final_v11=function(dataArranged,centroidPCAdata,argList,exCentroids=NULL,n.neighbors=4,batchPCAdata=NULL,n.trees=50,NNmethod="annoy",L2Norm=T,mergePseudocells=T,returnPropMat=F,cluster_count=NULL){
  
  myL2normFn=function(inputMat){
    prop_mat2=rowSums(inputMat^2)
    prop_mat2=sqrt(prop_mat2)
    res=.extra_matrix_rowNorm(input_mat = inputMat,rowValues = 1/(prop_mat2+0.00000001))#Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
    return(res)
  }
  
  if(is.null(batchPCAdata)){
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    batchPCAdata=harmony_embeddings
  }
  
  if(sum(argList$singleton.method %in% c("snn","fast"))==0){
    stop("unrecognized singleton.method")
  }
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  #batchPCAdata=harmony_embeddings
  batchPCAdata=as.matrix(batchPCAdata)[,1:argList$nPCs]
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  if(NNmethod=="annoy"){
    idx=Seurat:::AnnoyBuildIndex(data = centroidPCAdata, metric = "euclidean", 
                                 n.trees = n.trees)
    nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = batchPCAdata,k=n.neighbors,include.distance = T,search.k = -1)
  } else {
    nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  }
  
  adj= nn.ranked.1$nn.idx
  j <- as.numeric(t(adj))
  i <- ((1:length(j)) - 1)%/%n.neighbors + 1
  x=as.numeric(t(nn.ranked.1$nn.dists+0.01))
  adj = t(sparseMatrix(i = i, j = j, x = x, dims = c(nrow(batchPCAdata),nrow(centroidPCAdata))))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  
  adj=adj[rowSums(adj)>0.01,]
  
  adj2=adj
  adj2@x=rep(1,length(adj2@x))
  adj=adj[rowSums(adj2)>4,]
  
  tmp_df=data.frame(pseudocells=row.names(centroidPCAdata)[nn.ranked.1$nn.idx[,1]],dist=nn.ranked.1$nn.dists[,1],stringsAsFactors = F)
  tmp_df=aggregate(dist~pseudocells,data=tmp_df,function(x) median(x,na.rm = T))
  nn.centroid.purityIndex=tmp_df$dist
  names(nn.centroid.purityIndex)=as.character(tmp_df$pseudocells)
  
  nn.centroid.purityIndex=nn.centroid.purityIndex[match(row.names(centroidPCAdata),names(nn.centroid.purityIndex))]
  if(sum(is.na(nn.centroid.purityIndex))>0){
    nn.centroid.purityIndex[is.na(nn.centroid.purityIndex)]=0
  }
  tmp_df=quantile(nn.centroid.purityIndex,0.95)
  nn.centroid.purityIndex[nn.centroid.purityIndex>tmp_df]=tmp_df
  rm(tmp_df)
  nn.centroid.purityIndex=nn.centroid.purityIndex[row.names(adj)]
  
  affinities=.extra_matrix_rowNorm(input_mat = adj,rowValues = 1/(nn.centroid.purityIndex+0.000001))#Matrix::Diagonal(x=1/(nn.centroid.purityIndex+0.000001)) %*% adj
  affinities@x=-1*affinities@x^2
  affinities@x=exp(affinities@x)
  
  #adj@x=rep(1,length(adj@x))
  #affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  #affinities=do.call("rbind",affinities)
  #affinities=affinities*adj
  
  centroid.affinities <- affinities #%*% Matrix::Diagonal(x = 1 / (colSums(affinities)+0.000000000001))
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  
  if(F){
    dataArranged2=dataArranged
    for(i in 1:length(dataArranged2)){
      tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged2[[i]]$pcaData),colnames(prop_mat))],data=dataArranged2[[i]],pseudocell_sim_mat=prop_mat)
      dataArranged2[[i]]=tmp
    }
    p=.sconline.anno2pseudocell_tmp(res_prop=dataArranged2,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T)
    ggsave(plot=p$plot,file="~/myBucket/torm3.pdf",width=12,height = 12)
    
  }
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(centroid.affinities) <- row.names(adj)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  prop_mat=centroid.affinities[rowSums(centroid.affinities>0)>3,]
  
  
  #################################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_affinities=prop_mat %*% as.matrix(tst)
    
    #apply(tst_affinities,2,mean)
    apply(tst_affinities,2,function(x) sum(x>0.5))
    
    
    tst_pd=pd[match(colnames(prop_mat),row.names(pd)),]
    aggregate(prop_mat["8",]~tst_pd$anno_orig_cellState,sum,data=NULL)
    
    eff_size=.myEffSizePropMat(prop_mat)
    
    x=apply(tst_affinities,1,function(x) which(x==max(x))[1])
    x=data.frame(size=eff_size$effective_sample_size,g=x)
    aggregate(size~g,data=x,summary)
    
  }
  
  
  
  ###############################
  
  if(L2Norm){
    prop_mat2=myL2normFn(inputMat=prop_mat)
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
  } else {
    prop_mat2=prop_mat
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
    c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
    diag(c_c_aff)=1
    c_c_aff@x=pmin(c_c_aff@x,1)
    c_c_aff=sqrt(c_c_aff*t(c_c_aff))
  }
  
  
  
  #c_c_aff[which(c_c_aff<0.1)]=0
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff %*% c_c_aff
  
  
  prop_mat_list=list()
  rowMax_vals=c()
  for(i in seq(1,ncol(prop_mat),50000)){
    tmp=c_c_aff %*% prop_mat[,i:min(i+49999,ncol(prop_mat))]
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
    tmp2=.extra_matrix_rowNorm(input_mat = tmp,rowValues = 1/rowMax_vals)# Matrix::Diagonal(x = 1 / rowMax_vals) %*% tmp
    tmp2=Matrix::drop0(tmp2,tol=0.1)
    tmp2@x=rep(1,length(tmp2@x))
    
    tmp=tmp*tmp2
    tmp=Matrix::drop0(tmp)
    prop_mat_list[[i]]=tmp
  }
  
  #prop_mat=c_c_aff %*% prop_mat
  
  if(F){
    rowMax_vals=c()
    for(i in seq(1,nrow(prop_mat),1000)){
      tmp=qlcMatrix::rowMax(prop_mat[i:min(i+999,nrow(prop_mat)),])
      rowMax_vals=c(rowMax_vals,tmp)
    }
    rowMax_vals=unlist(lapply(rowMax_vals,as.numeric))
    
    prop_mat2=.extra_matrix_rowNorm(input_mat = prop_mat,rowValues = 1/rowMax_vals) #Matrix::Diagonal(x = 1 / rowMax_vals) %*% prop_mat
    prop_mat2=Matrix::drop0(prop_mat2,tol=0.1)
    prop_mat2@x=rep(1,length(prop_mat2@x))
    
    prop_mat=prop_mat * prop_mat2
    rm(prop_mat2)
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
    rm(tmp)
  }
  
  rm(prop_mat_list)
  gc()
  
  if(F){
    tmp_pd=pd[match(colnames(prop_mat),row.names(pd)),]
    anno=prop_mat %*% as.matrix(.myOneHotFn(tmp_pd$anno_orig_cellState))
    anno_group=apply(anno,1,function(x) which(x==max(x))[1])
    anno_group=colnames(anno)[anno_group]
    anno[which(anno_group=="Macrophages"),]
    
    anno[which(anno_group=="Microglia"),]
  }
  
  
  matEffectiveSize=.myEffSizePropMat(prop_mat)
  matEffectiveSize=matEffectiveSize$effective_sample_size
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  prop_mat=prop_mat[matEffectiveSize>=argList$min_cluster_size,]
  prop_mat=prop_mat[,colSums(prop_mat)>0]
  
  ##############
  if(mergePseudocells){
    if(L2Norm){
      prop_mat2=myL2normFn(inputMat = prop_mat)
      c_c_aff=t(prop_mat2)
      c_c_aff=prop_mat2 %*% c_c_aff
    } else {
      prop_mat2=prop_mat
      c_c_aff=t(prop_mat2)
      c_c_aff=prop_mat2 %*% c_c_aff
      c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
      diag(c_c_aff)=1
      c_c_aff@x=pmin(c_c_aff@x,1)
      c_c_aff=sqrt(c_c_aff*t(c_c_aff))
    }
    
    
    ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
    ps_merging=cutree(ps_merging,h=0.15)
    while(length(unique(ps_merging))<length(ps_merging)){
      ps_merging=data.frame(pseudocell=row.names(prop_mat),id=as.character(ps_merging),stringsAsFactors = F)
      ps_mid=ps_merging[!duplicated(ps_merging$id),]
      colnames(ps_mid)[1]="mid"
      ps_merging=merge(ps_merging,ps_mid,by="id")
      ps_merging=ps_merging[match(row.names(prop_mat),ps_merging$pseudocell),]
      ps_merging=as.matrix(.myOneHotFn(as.character(ps_merging$mid)))
      ps_merging=as(ps_merging,"dgCMatrix")
      ps_merging=t(ps_merging)
      colnames(ps_merging)=row.names(prop_mat)
      prop_mat=ps_merging %*% prop_mat
      rm(ps_merging)
      prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
      
      
      if(L2Norm){
        prop_mat2=myL2normFn(inputMat = prop_mat)
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
      } else {
        prop_mat2=prop_mat
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
        c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
        diag(c_c_aff)=1
        c_c_aff@x=pmin(c_c_aff@x,1)
        c_c_aff=sqrt(c_c_aff*t(c_c_aff))
      }
      
      ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
      ps_merging=cutree(ps_merging,h=0.15)
    }
    
    
    prop_mat2=prop_mat
    prop_mat2@x=rep(1,length(prop_mat2@x))
    prop_mat2=myL2normFn(prop_mat2)
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
    ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
    if(is.null(cluster_count)){
      ps_merging=cutree(ps_merging,h=0.15)
    } else {
      ps_merging=cutree(ps_merging,k=cluster_count)
    }
    
    while(length(unique(ps_merging))<length(ps_merging)){
      ps_merging=data.frame(pseudocell=row.names(prop_mat),id=as.character(ps_merging),stringsAsFactors = F)
      ps_mid=ps_merging[!duplicated(ps_merging$id),]
      colnames(ps_mid)[1]="mid"
      ps_merging=merge(ps_merging,ps_mid,by="id")
      ps_merging=ps_merging[match(row.names(prop_mat),ps_merging$pseudocell),]
      ps_merging=as.matrix(.myOneHotFn(as.character(ps_merging$mid)))
      ps_merging=as(ps_merging,"dgCMatrix")
      ps_merging=t(ps_merging)
      colnames(ps_merging)=row.names(prop_mat)
      prop_mat=ps_merging %*% prop_mat
      rm(ps_merging)
      prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
      
      
      if(L2Norm){
        prop_mat2=prop_mat
        prop_mat2@x=rep(1,length(prop_mat2@x))
        prop_mat2=myL2normFn(inputMat = prop_mat2)
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
      } else {
        prop_mat2=prop_mat
        prop_mat2@x=rep(1,length(prop_mat2@x))
        prop_mat2=prop_mat2
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
        c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
        diag(c_c_aff)=1
        c_c_aff@x=pmin(c_c_aff@x,1)
        c_c_aff=sqrt(c_c_aff*t(c_c_aff))
      }
      
      ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
      ps_merging=cutree(ps_merging,h=0.15)
    }
    
    if(F){
      tmp_pd=pd[match(colnames(prop_mat),row.names(pd)),]
      anno=prop_mat %*% as.matrix(.myOneHotFn(tmp_pd$anno_orig_cellState))
      anno_group=apply(anno,1,function(x) which(x==max(x)))
      anno_group=colnames(anno)[anno_group]
      
      matWeights=.myEffSizePropMat(prop_mat)
      matEffectiveSize=matWeights$effective_sample_size
      
      res_check=NULL
      for(i in ps_merging[duplicated(ps_merging)]){
        ii=which(ps_merging==i)
        for(j in 1:(length(ii)-1)){
          for(k in 2:length(ii)){
            Fscore=anno[ii[j],anno_group[ii[j]]]/anno[ii[j],anno_group[ii[k]]]
            Sscore=anno[ii[k],anno_group[ii[k]]]/anno[ii[k],anno_group[ii[j]]]
            res_check=rbind(res_check,data.frame(Fnode=anno_group[ii[j]],Snode=anno_group[ii[k]],effSize1=matEffectiveSize[ii[j]],effsize2=matEffectiveSize[ii[k]],Find=ii[j],Sind=ii[k],Fscore=Fscore,Sscore=Sscore,stringsAsFactors = F))
          }
        }
      }
      table(res_check$Fnode==res_check$Snode)
      length(unique(ps_merging))
      res_check=res_check[res_check$Fnode!=res_check$Snode,]
      table(res_check$Fscore>1.5|res_check$Sscore>1.5)
      res_check[res_check$Fscore>1.5|res_check$Sscore>1.5,]
    }
  }
  
  cat(paste("Number of retined pseudocells:",nrow(prop_mat),"out of",nrow(centroidPCAdata),"\n"))
  
  ##################
  
  
  
  connected_cells=colnames(prop_mat)
  if(argList$singleton.method=="snn"){
    cat("Handling singletons by SNN\n")
    nn.ranked <- Seurat:::NNHelper(query = batchPCAdata[c(connected_cells,setdiff(row.names(batchPCAdata),connected_cells)),], data = batchPCAdata[connected_cells,], k = 20, 
                                   method = "annoy", n.trees = 50, searchtype = "standard", 
                                   eps = 0, metric = "euclidean", cache.index = F, 
                                   index = NULL)
    nn.ranked = Indices(nn.ranked)
    
    snn.matrix = Seurat:::ComputeSNN(nn_ranked = nn.ranked, prune = 1/15)
    rownames(snn.matrix) = c(connected_cells,setdiff(row.names(batchPCAdata),connected_cells))
    colnames(snn.matrix) = c(connected_cells,setdiff(row.names(batchPCAdata),connected_cells))
    
    tst=snn.matrix
    diag(tst)=0
    tst=tst[setdiff(row.names(batchPCAdata),connected_cells),]
    tst=as.numeric(qlcMatrix::rowMax(tst))
    tst=c(rep(1,length(connected_cells)),tst)
    snn.matrix=.extra_matrix_rowNorm(input_mat = snn.matrix,rowValues = 1/(tst+0.000001))#Matrix::Diagonal(x=1/(tst+0.000001))%*% snn.matrix
    diag(snn.matrix)=1
    adj_all=snn.matrix[,colnames(prop_mat)]
  } else if(argList$singleton.method=="fast") {
    cat("Handling singletons by fast method\n")
    if(length(setdiff(row.names(batchPCAdata),connected_cells))>0){
      nn.ranked <- Seurat:::NNHelper(query = batchPCAdata[setdiff(row.names(batchPCAdata),connected_cells),], data = batchPCAdata[connected_cells,], k = 1, 
                                     method = "annoy", n.trees = 50, searchtype = "standard", 
                                     eps = 0, metric = "euclidean", cache.index = F, 
                                     index = NULL)
      nn.ranked = Indices(nn.ranked)
      nn.ranked=connected_cells[nn.ranked[,1]]
      nn.ranked=data.frame(sample=c(connected_cells,setdiff(row.names(batchPCAdata),connected_cells)),map=c(connected_cells,nn.ranked),stringsAsFactors = F)
    } else {
      nn.ranked=data.frame(sample=colnames(prop_mat),map=colnames(prop_mat),stringsAsFactors = F)
    }
    
    nn.ranked=merge(nn.ranked,data.frame(map=colnames(prop_mat),id=1:ncol(prop_mat),stringsAsFactors = F),by="map")
    #nn.ranked=nn.ranked[match(colnames(prop_mat),nn.ranked$sample),]
    adj_all=Matrix::sparseMatrix(i=1:nrow(nn.ranked),j=nn.ranked$id,x=rep(1,nrow(nn.ranked)))
    row.names(adj_all)=nn.ranked$sample
    colnames(adj_all)=colnames(prop_mat)
    rm(nn.ranked)
  } else {
    stop("unrecognized singleton.method")
  }
  
  
  ###############
  #New
  #Excluded the col normalization step
  #prop_mat <- prop_mat %*% Matrix::Diagonal(x = 1 / (colSums(prop_mat)+0.000000000001))
  adj_all=adj_all %*% t(prop_mat)
  #adj_all = Matrix::Diagonal(x = 1 / (rowSums(adj_all)+0.000000000001)) %*% adj_all
  prop_mat=t(adj_all)
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  if(returnPropMat){
    return(prop_mat)
  }
  
  
  matWeights=.myEffSizePropMat(prop_mat)
  
  matEffectiveSize=matWeights$effective_sample_size
  matWeights=matWeights$centroid_weights
  matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  
  
  for(i in 1:length(dataArranged)){
    tmp_prop_mat=prop_mat[,match(row.names(dataArranged[[i]]$pcaData),colnames(prop_mat))]
    #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
    tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
    tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
    tmp_weights[tmp_effsize<4]=0
    tmp=list(prop_mat=tmp_prop_mat,data=dataArranged[[i]],pseudocell_sim_mat=c_c_aff,matWeights=tmp_weights,matEffectiveSize=tmp_effsize)
    dataArranged[[i]]=tmp
  }
  
  return(dataArranged)
}

.myConcensusDEFn_step2_detail_newprop3_final_v11subset=function(dataArranged,centroidPCAdata,argList,exCentroids=NULL,n.neighbors=4,batchPCAdata=NULL,n.trees=50,NNmethod="annoy",L2Norm=T,mergePseudocells=T,returnPropMat=F,cluster_count=NULL){
  
  myL2normFn=function(inputMat){
    prop_mat2=rowSums(inputMat^2)
    prop_mat2=sqrt(prop_mat2)
    res=.extra_matrix_rowNorm(input_mat = inputMat,rowValues = 1/(prop_mat2+0.00000001))#Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
    return(res)
  }
  
  if(is.null(batchPCAdata)){
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    batchPCAdata=harmony_embeddings
  }
  
  if(sum(argList$singleton.method %in% c("snn","fast"))==0){
    stop("unrecognized singleton.method")
  }
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  #batchPCAdata=harmony_embeddings
  batchPCAdata=as.matrix(batchPCAdata)[,1:argList$nPCs]
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  if(NNmethod=="annoy"){
    idx=Seurat:::AnnoyBuildIndex(data = centroidPCAdata, metric = "euclidean", 
                                 n.trees = n.trees)
    nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = batchPCAdata,k=n.neighbors,include.distance = T,search.k = -1)
  } else {
    nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  }
  
  adj= nn.ranked.1$nn.idx
  j <- as.numeric(t(adj))
  i <- ((1:length(j)) - 1)%/%n.neighbors + 1
  x=as.numeric(t(nn.ranked.1$nn.dists+0.01))
  adj = t(sparseMatrix(i = i, j = j, x = x, dims = c(nrow(batchPCAdata),nrow(centroidPCAdata))))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  
  adj=adj[rowSums(adj)>0.01,]
  
  adj2=adj
  adj2@x=rep(1,length(adj2@x))
  adj=adj[rowSums(adj2)>4,]
  
  tmp_df=data.frame(pseudocells=row.names(centroidPCAdata)[nn.ranked.1$nn.idx[,1]],dist=nn.ranked.1$nn.dists[,1],stringsAsFactors = F)
  tmp_df=aggregate(dist~pseudocells,data=tmp_df,function(x) median(x,na.rm = T))
  nn.centroid.purityIndex=tmp_df$dist
  names(nn.centroid.purityIndex)=as.character(tmp_df$pseudocells)
  
  nn.centroid.purityIndex=nn.centroid.purityIndex[match(row.names(centroidPCAdata),names(nn.centroid.purityIndex))]
  if(sum(is.na(nn.centroid.purityIndex))>0){
    nn.centroid.purityIndex[is.na(nn.centroid.purityIndex)]=0
  }
  tmp_df=quantile(nn.centroid.purityIndex,0.95)
  nn.centroid.purityIndex[nn.centroid.purityIndex>tmp_df]=tmp_df
  rm(tmp_df)
  nn.centroid.purityIndex=nn.centroid.purityIndex[row.names(adj)]
  
  affinities=.extra_matrix_rowNorm(input_mat = adj,rowValues = 1/(nn.centroid.purityIndex+0.000001))#Matrix::Diagonal(x=1/(nn.centroid.purityIndex+0.000001)) %*% adj
  affinities@x=-1*affinities@x^2
  affinities@x=exp(affinities@x)
  
  centroid.affinities <- affinities #%*% Matrix::Diagonal(x = 1 / (colSums(affinities)+0.000000000001))
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(centroid.affinities) <- row.names(adj)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  prop_mat=centroid.affinities[rowSums(centroid.affinities>0)>3,]
  
  
  
  if(L2Norm){
    prop_mat2=myL2normFn(inputMat=prop_mat)
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
  } else {
    prop_mat2=prop_mat
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
    c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
    diag(c_c_aff)=1
    c_c_aff@x=pmin(c_c_aff@x,1)
    c_c_aff=sqrt(c_c_aff*t(c_c_aff))
  }
  
  
  
  #c_c_aff[which(c_c_aff<0.1)]=0
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff %*% c_c_aff
  
  
  prop_mat_list=list()
  rowMax_vals=c()
  for(i in seq(1,ncol(prop_mat),50000)){
    tmp=c_c_aff %*% prop_mat[,i:min(i+49999,ncol(prop_mat))]
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
    tmp2=.extra_matrix_rowNorm(input_mat = tmp,rowValues = 1/rowMax_vals)#Matrix::Diagonal(x = 1 / rowMax_vals) %*% tmp
    tmp2=Matrix::drop0(tmp2,tol=0.1)
    tmp2@x=rep(1,length(tmp2@x))
    
    tmp=tmp*tmp2
    tmp=Matrix::drop0(tmp)
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
    rm(tmp)
  }
  
  rm(prop_mat_list)
  gc()
  
  matEffectiveSize=.myEffSizePropMat(prop_mat)
  matEffectiveSize=matEffectiveSize$effective_sample_size
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  prop_mat=prop_mat[matEffectiveSize>=argList$min_cluster_size,]
  prop_mat=prop_mat[,colSums(prop_mat)>0]
  
  ##############
  if(mergePseudocells){
    if(L2Norm){
      prop_mat2=myL2normFn(inputMat = prop_mat)
      c_c_aff=t(prop_mat2)
      c_c_aff=prop_mat2 %*% c_c_aff
    } else {
      prop_mat2=prop_mat
      c_c_aff=t(prop_mat2)
      c_c_aff=prop_mat2 %*% c_c_aff
      c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
      diag(c_c_aff)=1
      c_c_aff@x=pmin(c_c_aff@x,1)
      c_c_aff=sqrt(c_c_aff*t(c_c_aff))
    }
    
    
    ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
    ps_merging=cutree(ps_merging,h=0.15)
    while(length(unique(ps_merging))<length(ps_merging)){
      ps_merging=data.frame(pseudocell=row.names(prop_mat),id=as.character(ps_merging),stringsAsFactors = F)
      ps_mid=ps_merging[!duplicated(ps_merging$id),]
      colnames(ps_mid)[1]="mid"
      ps_merging=merge(ps_merging,ps_mid,by="id")
      ps_merging=ps_merging[match(row.names(prop_mat),ps_merging$pseudocell),]
      ps_merging=as.matrix(.myOneHotFn(as.character(ps_merging$mid)))
      ps_merging=as(ps_merging,"dgCMatrix")
      ps_merging=t(ps_merging)
      colnames(ps_merging)=row.names(prop_mat)
      prop_mat=ps_merging %*% prop_mat
      rm(ps_merging)
      prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
      
      
      if(L2Norm){
        prop_mat2=myL2normFn(inputMat = prop_mat)
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
      } else {
        prop_mat2=prop_mat
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
        c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
        diag(c_c_aff)=1
        c_c_aff@x=pmin(c_c_aff@x,1)
        c_c_aff=sqrt(c_c_aff*t(c_c_aff))
      }
      
      ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
      ps_merging=cutree(ps_merging,h=0.15)
    }
    
    
    prop_mat2=prop_mat
    prop_mat2@x=rep(1,length(prop_mat2@x))
    prop_mat2=myL2normFn(prop_mat2)
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
    ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
    ps_merging=cutree(ps_merging,h=0.15)
    
    
    while(length(unique(ps_merging))<length(ps_merging)){
      ps_merging=data.frame(pseudocell=row.names(prop_mat),id=as.character(ps_merging),stringsAsFactors = F)
      ps_mid=ps_merging[!duplicated(ps_merging$id),]
      colnames(ps_mid)[1]="mid"
      ps_merging=merge(ps_merging,ps_mid,by="id")
      ps_merging=ps_merging[match(row.names(prop_mat),ps_merging$pseudocell),]
      ps_merging=as.matrix(.myOneHotFn(as.character(ps_merging$mid)))
      ps_merging=as(ps_merging,"dgCMatrix")
      ps_merging=t(ps_merging)
      colnames(ps_merging)=row.names(prop_mat)
      prop_mat=ps_merging %*% prop_mat
      rm(ps_merging)
      prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
      
      
      if(L2Norm){
        prop_mat2=prop_mat
        prop_mat2@x=rep(1,length(prop_mat2@x))
        prop_mat2=myL2normFn(inputMat = prop_mat2)
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
      } else {
        prop_mat2=prop_mat
        prop_mat2@x=rep(1,length(prop_mat2@x))
        prop_mat2=prop_mat2
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
        c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
        diag(c_c_aff)=1
        c_c_aff@x=pmin(c_c_aff@x,1)
        c_c_aff=sqrt(c_c_aff*t(c_c_aff))
      }
      
      ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
      ps_merging=cutree(ps_merging,h=0.15)
    }
    
    if(!is.null(cluster_count)&nrow(prop_mat)>cluster_count){
      prop_mat2=prop_mat
      prop_mat2@x=rep(1,length(prop_mat2@x))
      prop_mat2=myL2normFn(prop_mat2)
      c_c_aff=t(prop_mat2)
      c_c_aff=prop_mat2 %*% c_c_aff
      ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
      ps_merging=cutree(ps_merging,k=cluster_count)
      if(length(unique(ps_merging))<length(ps_merging)){
        ps_merging=data.frame(pseudocell=row.names(prop_mat),id=as.character(ps_merging),stringsAsFactors = F)
        ps_mid=ps_merging[!duplicated(ps_merging$id),]
        colnames(ps_mid)[1]="mid"
        ps_merging=merge(ps_merging,ps_mid,by="id")
        ps_merging=ps_merging[match(row.names(prop_mat),ps_merging$pseudocell),]
        ps_merging=as.matrix(.myOneHotFn(as.character(ps_merging$mid)))
        ps_merging=as(ps_merging,"dgCMatrix")
        ps_merging=t(ps_merging)
        colnames(ps_merging)=row.names(prop_mat)
        prop_mat=ps_merging %*% prop_mat
        rm(ps_merging)
        prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
      }
    }
  }
  
  cat(paste("Number of retined pseudocells:",nrow(prop_mat),"out of",nrow(centroidPCAdata),"\n"))
  
  ##################
  
  
  
  connected_cells=colnames(prop_mat)
  if(argList$singleton.method=="snn"){
    cat("Handling singletons by SNN\n")
    nn.ranked <- Seurat:::NNHelper(query = batchPCAdata[c(connected_cells,setdiff(row.names(batchPCAdata),connected_cells)),], data = batchPCAdata[connected_cells,], k = 20, 
                                   method = "annoy", n.trees = 50, searchtype = "standard", 
                                   eps = 0, metric = "euclidean", cache.index = F, 
                                   index = NULL)
    nn.ranked = Indices(nn.ranked)
    
    snn.matrix = Seurat:::ComputeSNN(nn_ranked = nn.ranked, prune = 1/15)
    rownames(snn.matrix) = c(connected_cells,setdiff(row.names(batchPCAdata),connected_cells))
    colnames(snn.matrix) = c(connected_cells,setdiff(row.names(batchPCAdata),connected_cells))
    
    tst=snn.matrix
    diag(tst)=0
    tst=tst[setdiff(row.names(batchPCAdata),connected_cells),]
    tst=as.numeric(qlcMatrix::rowMax(tst))
    tst=c(rep(1,length(connected_cells)),tst)
    snn.matrix=.extra_matrix_rowNorm(input_mat = snn.matrix,rowValues = 1/(tst+0.000001))#Matrix::Diagonal(x=1/(tst+0.000001))%*% snn.matrix
    diag(snn.matrix)=1
    adj_all=snn.matrix[,colnames(prop_mat)]
  } else if(argList$singleton.method=="fast") {
    cat("Handling singletons by fast method\n")
    if(length(setdiff(row.names(batchPCAdata),connected_cells))>0){
      nn.ranked <- Seurat:::NNHelper(query = batchPCAdata[setdiff(row.names(batchPCAdata),connected_cells),], data = batchPCAdata[connected_cells,], k = 1, 
                                     method = "annoy", n.trees = 50, searchtype = "standard", 
                                     eps = 0, metric = "euclidean", cache.index = F, 
                                     index = NULL)
      nn.ranked = Indices(nn.ranked)
      nn.ranked=connected_cells[nn.ranked[,1]]
      nn.ranked=data.frame(sample=c(connected_cells,setdiff(row.names(batchPCAdata),connected_cells)),map=c(connected_cells,nn.ranked),stringsAsFactors = F)
    } else {
      nn.ranked=data.frame(sample=colnames(prop_mat),map=colnames(prop_mat),stringsAsFactors = F)
    }
    
    nn.ranked=merge(nn.ranked,data.frame(map=colnames(prop_mat),id=1:ncol(prop_mat),stringsAsFactors = F),by="map")
    #nn.ranked=nn.ranked[match(colnames(prop_mat),nn.ranked$sample),]
    adj_all=Matrix::sparseMatrix(i=1:nrow(nn.ranked),j=nn.ranked$id,x=rep(1,nrow(nn.ranked)))
    row.names(adj_all)=nn.ranked$sample
    colnames(adj_all)=colnames(prop_mat)
    rm(nn.ranked)
  } else {
    stop("unrecognized singleton.method")
  }
  
  
  ###############
  #New
  #Excluded the col normalization step
  #prop_mat <- prop_mat %*% Matrix::Diagonal(x = 1 / (colSums(prop_mat)+0.000000000001))
  adj_all=adj_all %*% t(prop_mat)
  #adj_all = Matrix::Diagonal(x = 1 / (rowSums(adj_all)+0.000000000001)) %*% adj_all
  prop_mat=t(adj_all)
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  return(prop_mat)
  
}

.myConcensusDEFn_step2_detail_newprop3_final_v12=function(dataArranged,centroidPCAdata,argList,exCentroids=NULL,n.neighbors=4,batchPCAdata=NULL,n.trees=50,NNmethod="annoy",L2Norm=T,mergePseudocells=T,hierarchical_refinement=T,batch_variable="anno_batch",colNormalize=F,merging_strength=0.15,...){
  
  myL2normFn=function(inputMat){
    prop_mat2=rowSums(inputMat^2)
    prop_mat2=sqrt(prop_mat2)
    res=.extra_matrix_rowNorm(input_mat = inputMat,rowValues = 1/(prop_mat2+0.00000001))#Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
    return(res)
  }
  
  if(is.null(batchPCAdata)){
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    batchPCAdata=harmony_embeddings
  }
  
  if(sum(argList$singleton.method %in% c("snn","fast"))==0){
    stop("unrecognized singleton.method")
  }
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  #batchPCAdata=harmony_embeddings
  batchPCAdata=as.matrix(batchPCAdata)[,1:argList$nPCs]
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  if(NNmethod=="annoy"){
    idx=Seurat:::AnnoyBuildIndex(data = centroidPCAdata, metric = "euclidean", 
                                 n.trees = n.trees)
    nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = batchPCAdata,k=n.neighbors,include.distance = T,search.k = -1)
  } else {
    nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  }
  
  if(F){
    load(.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
    pd_ps=pd[sl_pseudo$pseudocell,]
    x=matrix("",nrow=nrow(nn.ranked.1$nn.idx),ncol=ncol(nn.ranked.1$nn.idx))
    
    for(i in 1:ncol(nn.ranked.1$nn.idx)){
      x[,i]=pd_ps$anno_orig_cellState[nn.ranked.1$nn.idx[,i]]
    }
    
    tmp=unlist(lapply(1:nrow(x),function(y) sum(x[y,]==pd$anno_orig_cellState[x])/ncol(x)))
    
  }
  
  adj= nn.ranked.1$nn.idx
  j <- as.numeric(t(adj))
  i <- ((1:length(j)) - 1)%/%n.neighbors + 1
  x=as.numeric(t(nn.ranked.1$nn.dists+0.01))
  adj = t(sparseMatrix(i = i, j = j, x = x, dims = c(nrow(batchPCAdata),nrow(centroidPCAdata))))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  
  adj=adj[rowSums(adj)>0.01,]
  
  adj2=adj
  adj2@x=rep(1,length(adj2@x))
  adj=adj[rowSums(adj2)>4,]
  
  tmp_df=data.frame(pseudocells=row.names(centroidPCAdata)[nn.ranked.1$nn.idx[,1]],dist=nn.ranked.1$nn.dists[,1],stringsAsFactors = F)
  tmp_df=aggregate(dist~pseudocells,data=tmp_df,function(x) median(x,na.rm = T))
  nn.centroid.purityIndex=tmp_df$dist
  names(nn.centroid.purityIndex)=as.character(tmp_df$pseudocells)
  
  nn.centroid.purityIndex=nn.centroid.purityIndex[match(row.names(centroidPCAdata),names(nn.centroid.purityIndex))]
  if(sum(is.na(nn.centroid.purityIndex))>0){
    nn.centroid.purityIndex[is.na(nn.centroid.purityIndex)]=0
  }
  tmp_df=quantile(nn.centroid.purityIndex,0.95)
  nn.centroid.purityIndex[nn.centroid.purityIndex>tmp_df]=tmp_df
  rm(tmp_df)
  nn.centroid.purityIndex=nn.centroid.purityIndex[row.names(adj)]
  
  affinities=.extra_matrix_rowNorm(input_mat = adj,rowValues = 1/(nn.centroid.purityIndex+0.000001))#Matrix::Diagonal(x=1/(nn.centroid.purityIndex+0.000001)) %*% adj
  affinities@x=-1*affinities@x^2
  affinities@x=exp(affinities@x)
  
  centroid.affinities <- affinities
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(centroid.affinities) <- row.names(adj)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  prop_mat=centroid.affinities[rowSums(centroid.affinities>0)>3,]
  
  
  if(L2Norm){
    prop_mat2=myL2normFn(inputMat=prop_mat)
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
  } else {
    prop_mat2=prop_mat
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
    c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
    diag(c_c_aff)=1
    c_c_aff@x=pmin(c_c_aff@x,1)
    c_c_aff=sqrt(c_c_aff*t(c_c_aff))
  }
  
  
  
  #c_c_aff[which(c_c_aff<0.1)]=0
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff %*% c_c_aff
  
  
  prop_mat_list=list()
  rowMax_vals=c()
  for(i in seq(1,ncol(prop_mat),50000)){
    tmp=c_c_aff %*% prop_mat[,i:min(i+49999,ncol(prop_mat))]
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
    tmp2=.extra_matrix_rowNorm(input_mat = tmp,rowValues = 1/rowMax_vals)#Matrix::Diagonal(x = 1 / rowMax_vals) %*% tmp
    tmp2=Matrix::drop0(tmp2,tol=0.1)
    tmp2@x=rep(1,length(tmp2@x))
    
    tmp=tmp*tmp2
    tmp=Matrix::drop0(tmp)
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
    rm(tmp)
  }
  
  rm(prop_mat_list)
  gc()
  
  
  matEffectiveSize=.myEffSizePropMat(prop_mat)
  matEffectiveSize=matEffectiveSize$effective_sample_size
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  prop_mat=prop_mat[matEffectiveSize>=argList$min_cluster_size,]
  prop_mat=prop_mat[,colSums(prop_mat)>0]
  
  ##############
  if(mergePseudocells){
    if(L2Norm){
      prop_mat2=myL2normFn(inputMat = prop_mat)
      c_c_aff=t(prop_mat2)
      c_c_aff=prop_mat2 %*% c_c_aff
    } else {
      prop_mat2=prop_mat
      c_c_aff=t(prop_mat2)
      c_c_aff=prop_mat2 %*% c_c_aff
      c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
      diag(c_c_aff)=1
      c_c_aff@x=pmin(c_c_aff@x,1)
      c_c_aff=sqrt(c_c_aff*t(c_c_aff))
    }
    
    
    ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
    ps_merging=cutree(ps_merging,h=merging_strength)
    while(length(unique(ps_merging))<length(ps_merging)){
      ps_merging=data.frame(pseudocell=row.names(prop_mat),id=as.character(ps_merging),stringsAsFactors = F)
      ps_mid=ps_merging[!duplicated(ps_merging$id),]
      colnames(ps_mid)[1]="mid"
      ps_merging=merge(ps_merging,ps_mid,by="id")
      ps_merging=ps_merging[match(row.names(prop_mat),ps_merging$pseudocell),]
      ps_merging=as.matrix(.myOneHotFn(as.character(ps_merging$mid)))
      ps_merging=as(ps_merging,"dgCMatrix")
      ps_merging=t(ps_merging)
      colnames(ps_merging)=row.names(prop_mat)
      prop_mat=ps_merging %*% prop_mat
      rm(ps_merging)
      prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
      
      
      if(L2Norm){
        prop_mat2=myL2normFn(inputMat = prop_mat)
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
      } else {
        prop_mat2=prop_mat
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
        c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
        diag(c_c_aff)=1
        c_c_aff@x=pmin(c_c_aff@x,1)
        c_c_aff=sqrt(c_c_aff*t(c_c_aff))
      }
      
      ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
      ps_merging=cutree(ps_merging,h=merging_strength)
    }
    
    
    prop_mat2=prop_mat
    prop_mat2@x=rep(1,length(prop_mat2@x))
    prop_mat2=myL2normFn(prop_mat2)
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
    ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
    ps_merging=cutree(ps_merging,h=0.1)
    while(length(unique(ps_merging))<length(ps_merging)){
      ps_merging=data.frame(pseudocell=row.names(prop_mat),id=as.character(ps_merging),stringsAsFactors = F)
      ps_mid=ps_merging[!duplicated(ps_merging$id),]
      colnames(ps_mid)[1]="mid"
      ps_merging=merge(ps_merging,ps_mid,by="id")
      ps_merging=ps_merging[match(row.names(prop_mat),ps_merging$pseudocell),]
      ps_merging=as.matrix(.myOneHotFn(as.character(ps_merging$mid)))
      ps_merging=as(ps_merging,"dgCMatrix")
      ps_merging=t(ps_merging)
      colnames(ps_merging)=row.names(prop_mat)
      prop_mat=ps_merging %*% prop_mat
      rm(ps_merging)
      prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
      
      
      if(L2Norm){
        prop_mat2=prop_mat
        prop_mat2@x=rep(1,length(prop_mat2@x))
        prop_mat2=myL2normFn(inputMat = prop_mat2)
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
      } else {
        prop_mat2=prop_mat
        prop_mat2@x=rep(1,length(prop_mat2@x))
        prop_mat2=prop_mat2
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
        c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
        diag(c_c_aff)=1
        c_c_aff@x=pmin(c_c_aff@x,1)
        c_c_aff=sqrt(c_c_aff*t(c_c_aff))
      }
      
      ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
      ps_merging=cutree(ps_merging,h=0.1)
    }
    
    if(F){
      tmp_pd=pd[match(colnames(prop_mat),row.names(pd)),]
      anno=prop_mat %*% as.matrix(.myOneHotFn(tmp_pd$anno_orig_cellState))
      anno_group=apply(anno,1,function(x) which(x==max(x)))
      anno_group=colnames(anno)[anno_group]
      
      matWeights=.myEffSizePropMat(prop_mat)
      matEffectiveSize=matWeights$effective_sample_size
      
      res_check=NULL
      for(i in ps_merging[duplicated(ps_merging)]){
        ii=which(ps_merging==i)
        for(j in 1:(length(ii)-1)){
          for(k in 2:length(ii)){
            Fscore=anno[ii[j],anno_group[ii[j]]]/anno[ii[j],anno_group[ii[k]]]
            Sscore=anno[ii[k],anno_group[ii[k]]]/anno[ii[k],anno_group[ii[j]]]
            res_check=rbind(res_check,data.frame(Fnode=anno_group[ii[j]],Snode=anno_group[ii[k]],effSize1=matEffectiveSize[ii[j]],effsize2=matEffectiveSize[ii[k]],Find=ii[j],Sind=ii[k],Fscore=Fscore,Sscore=Sscore,stringsAsFactors = F))
          }
        }
      }
      table(res_check$Fnode==res_check$Snode)
      length(unique(ps_merging))
      res_check=res_check[res_check$Fnode!=res_check$Snode,]
      table(res_check$Fscore>1.5|res_check$Sscore>1.5)
      res_check[res_check$Fscore>1.5|res_check$Sscore>1.5,]
    }
  }
  
  if(hierarchical_refinement){
    
    x=prop_mat
    x@x=rep(1,length(x@x))
    
    prop_mat2=prop_mat
    prop_mat2@x=rep(1,length(prop_mat2@x))
    prop_mat2=prop_mat2 %*% t(prop_mat2)
    prop_mat2=abs(prop_mat2 - diag(prop_mat2))
    prop_mat2=Matrix::drop0(prop_mat2)
    clust_thr=15
    check_clustering=T
    while(check_clustering&clust_thr>0){
      sl_components=which(prop_mat2<max(clust_thr,1),arr.ind = T)
      sl_components=sl_components[sl_components[,1]!=sl_components[,2],]
      count_singletons=length(setdiff(1:nrow(prop_mat2),c(sl_components[,1],sl_components[,2])))
      sl_components=igraph::graph_from_data_frame(as.data.frame(sl_components), directed = TRUE)
      sl_components=igraph::components(sl_components)
      if(length(unique(sl_components$membership))+count_singletons>200&max(sl_components$csize)/nrow(prop_mat)<0.3){
        check_clustering=F
      } else {
        clust_thr=clust_thr-1
      }
    }
    
    if(length(sl_components$csize)>1){
      included_ps_list=c()
      for(iin in unique(sl_components$membership)){
        tmp_sl=colnames(prop_mat2)[as.numeric(names(sl_components$membership)[which(sl_components$membership==iin)])]
        if(sum(setdiff(row.names(prop_mat),included_ps_list) %in% tmp_sl)>5){
          tmp_sl_cells=prop_mat[row.names(prop_mat) %in% tmp_sl,]
          tmp_sl_cells=colnames(tmp_sl_cells)[colSums(tmp_sl_cells)>0]
          load(.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
          sl_pseudo=sl_pseudo[sl_pseudo$pseudocell %in% tmp_sl_cells,]
          sl_pseudo=intersect(sl_pseudo$cluster,row.names(prop_mat))
          
          tmp_sl_cells=prop_mat[sl_pseudo,]
          tmp_sl_cells=colnames(tmp_sl_cells)[colSums(tmp_sl_cells)>0]
          tmp_sl_cells=setdiff(tmp_sl_cells,included_ps_list)
          if(length(tmp_sl_cells)>2000){
            cluster_count=length(sl_pseudo)
            if(F){
              if(cluster_count<20){
                cluster_count=min(10,cluster_count)
              } else {
                cluster_count=round(cluster_count/2)
              }
            }
            
            #batchPCAdata[tmp_sl_cells,]
            tmp_new_prop=.sconline.subset_purityAnalysis(argList=argList,sl_cells=tmp_sl_cells,inputPCAdata=NULL,batch_variable=batch_variable,cluster_count=min(100,cluster_count))
            tmp_new_prop2=sparseMatrix(i = integer(0), j = integer(0), dims=c(nrow(tmp_new_prop),length(setdiff(colnames(prop_mat),colnames(tmp_new_prop)))))
            colnames(tmp_new_prop2)=setdiff(colnames(prop_mat),colnames(tmp_new_prop))
            row.names(tmp_new_prop2)=row.names(tmp_new_prop)
            tmp_new_prop=cbind(tmp_new_prop,tmp_new_prop2)
            tmp_new_prop=tmp_new_prop[,match(colnames(prop_mat),colnames(tmp_new_prop))]
            included_ps_list=c(included_ps_list,row.names(tmp_new_prop))
            prop_mat=prop_mat[which(!row.names(prop_mat) %in% sl_pseudo),]
            
            prop_mat=rbind(prop_mat,tmp_new_prop)
          }
          
        }
      }
    }
    
    
    if(F){
      tmp_pd=pd[match(colnames(prop_mat),row.names(pd)),]
      anno=prop_mat %*% as.matrix(.myOneHotFn(tmp_pd$anno_orig_cellState))
      anno_group=apply(anno,1,function(x) which(x==max(x)))
      anno_group=colnames(anno)[anno_group]
      names(anno_group)=row.names(prop_mat)
    }
    
  }
  
  if(F){
    rowMax_vals=c()
    for(i in seq(1,nrow(prop_mat),500)){
      tmp_max=as.numeric(qlcMatrix::rowMax(prop_mat[i:min(i+499,nrow(prop_mat)),]))
      rowMax_vals=c(rowMax_vals,as.numeric(tmp_max))
    }
    prop_mat =.extra_matrix_rowNorm(input_mat = prop_mat,rowValues = 1/rowMax_vals)# Matrix::Diagonal(x = 1 / rowMax_vals) %*% prop_mat
    prop_mat=Matrix::drop0(prop_mat,0.1)
    prop_mat = .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / rowSums(prop_mat+0.000001)) %*% prop_mat
  }
  
  cat(paste("Number of retined pseudocells:",nrow(prop_mat),"out of",nrow(centroidPCAdata),"\n"))
  
  ##################
  
  connected_cells=colnames(prop_mat)
  if(argList$singleton.method=="snn"){
    cat("Handling singletons by SNN\n")
    nn.ranked <- Seurat:::NNHelper(query = batchPCAdata[c(connected_cells,setdiff(row.names(batchPCAdata),connected_cells)),], data = batchPCAdata[connected_cells,], k = 20, 
                                   method = "annoy", n.trees = 50, searchtype = "standard", 
                                   eps = 0, metric = "euclidean", cache.index = F, 
                                   index = NULL)
    nn.ranked = Indices(nn.ranked)
    
    snn.matrix = Seurat:::ComputeSNN(nn_ranked = nn.ranked, prune = 1/15)
    rownames(snn.matrix) = c(connected_cells,setdiff(row.names(batchPCAdata),connected_cells))
    colnames(snn.matrix) = c(connected_cells,setdiff(row.names(batchPCAdata),connected_cells))
    
    tst=snn.matrix
    diag(tst)=0
    tst=tst[setdiff(row.names(batchPCAdata),connected_cells),]
    tst=as.numeric(qlcMatrix::rowMax(tst))
    tst=c(rep(1,length(connected_cells)),tst)
    snn.matrix=.extra_matrix_rowNorm(input_mat =snn.matrix,rowValues = 1/(tst+0.000001))#Matrix::Diagonal(x=1/(tst+0.000001))%*% snn.matrix
    diag(snn.matrix)=1
    adj_all=snn.matrix[,colnames(prop_mat)]
  } else if(argList$singleton.method=="fast") {
    cat("Handling singletons by fast method\n")
    if(length(setdiff(row.names(batchPCAdata),connected_cells))>0){
      nn.ranked <- Seurat:::NNHelper(query = batchPCAdata[setdiff(row.names(batchPCAdata),connected_cells),], data = batchPCAdata[connected_cells,], k = 1, 
                                     method = "annoy", n.trees = 50, searchtype = "standard", 
                                     eps = 0, metric = "euclidean", cache.index = F, 
                                     index = NULL)
      nn.ranked = Indices(nn.ranked)
      nn.ranked=connected_cells[nn.ranked[,1]]
      nn.ranked=data.frame(sample=c(connected_cells,setdiff(row.names(batchPCAdata),connected_cells)),map=c(connected_cells,nn.ranked),stringsAsFactors = F)
    } else {
      nn.ranked=data.frame(sample=colnames(prop_mat),map=colnames(prop_mat),stringsAsFactors = F)
    }
    
    nn.ranked=merge(nn.ranked,data.frame(map=colnames(prop_mat),id=1:ncol(prop_mat),stringsAsFactors = F),by="map")
    #nn.ranked=nn.ranked[match(colnames(prop_mat),nn.ranked$sample),]
    adj_all=Matrix::sparseMatrix(i=1:nrow(nn.ranked),j=nn.ranked$id,x=rep(1,nrow(nn.ranked)))
    row.names(adj_all)=nn.ranked$sample
    colnames(adj_all)=colnames(prop_mat)
    rm(nn.ranked)
  } else {
    stop("unrecognized singleton.method")
  }
  
  
  ###############
  #New
  #Excluded the col normalization step
  #prop_mat <- prop_mat %*% Matrix::Diagonal(x = 1 / (colSums(prop_mat)+0.000000000001))
  adj_all=adj_all %*% t(prop_mat)
  #adj_all = Matrix::Diagonal(x = 1 / (rowSums(adj_all)+0.000000000001)) %*% adj_all
  prop_mat=t(adj_all)
  
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  if(colNormalize){
    rowMax_vals=c()
    for(i in seq(1,nrow(prop_mat),500)){
      tmp_max=as.numeric(qlcMatrix::rowMax(prop_mat[i:min(i+499,nrow(prop_mat)),]))
      rowMax_vals=c(rowMax_vals,as.numeric(tmp_max))
    }
    prop_mat =.extra_matrix_rowNorm(input_mat = prop_mat,rowValues = 1/rowMax_vals)# Matrix::Diagonal(x = 1 / rowMax_vals) %*% prop_mat
    rm(rowMax_vals)
    colMax_vals=c()
    for(i in seq(1,ncol(prop_mat),5000)){
      tmp_max=as.numeric(qlcMatrix::colMax(prop_mat[,i:min(i+4999,ncol(prop_mat))]))
      colMax_vals=c(colMax_vals,as.numeric(tmp_max))
    }
    prop_mat = .extra_matrix_colNorm(input_mat = prop_mat,colValues = 1/colMax_vals)#  prop_mat %*% Matrix::Diagonal(x = 1 / colMax_vals)
    prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / rowSums(prop_mat)) %*% prop_mat
  }
  
  
  matWeights=.myEffSizePropMat(prop_mat)
  
  matEffectiveSize=matWeights$effective_sample_size
  matWeights=matWeights$centroid_weights
  matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  
  
  for(i in 1:length(dataArranged)){
    tmp_prop_mat=prop_mat[,match(row.names(dataArranged[[i]]$pcaData),colnames(prop_mat))]
    #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
    tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
    tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
    tmp_weights[tmp_effsize<4]=0
    tmp=list(prop_mat=tmp_prop_mat,data=dataArranged[[i]],pseudocell_sim_mat=c_c_aff,matWeights=tmp_weights,matEffectiveSize=tmp_effsize)
    dataArranged[[i]]=tmp
  }
  
  return(dataArranged)
}

.myConcensusDEFn_step2_detail_newprop3_final_v13=function(dataArranged,centroidPCAdata,argList,exCentroids=NULL,n.neighbors=4,batchPCAdata=NULL,n.trees=50,NNmethod="annoy",L2Norm=T,mergePseudocells=T,batch_variable="anno_batch",merging_strength=0.15,...){
  
  myL2normFn=function(inputMat){
    prop_mat2=rowSums(inputMat^2)
    prop_mat2=sqrt(prop_mat2)
    res=.extra_matrix_rowNorm(input_mat = inputMat,rowValues = 1/(prop_mat2+0.00000001))#Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
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
    idx=Seurat:::AnnoyBuildIndex(data = batchPCAdata, metric = "euclidean", 
                                 n.trees = n.trees)
    nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = batchPCAdata,k=n.neighbors,include.distance = T,search.k = -1)
  } else {
    nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  }
  
  
  affinities=.extra_matrix_rowNorm(input_mat = nn.ranked.1$nn.dists,rowValues = 1/(nn.ranked.1$nn.dists[,2]+0.000001))#Matrix::Diagonal(x=1/(nn.ranked.1$nn.dists[,2]+0.000001)) %*% nn.ranked.1$nn.dists
  affinities@x=-1*affinities@x^2
  affinities@x=exp(affinities@x)
  affinities[,1]=affinities[,2]
  j <- as.numeric(t(nn.ranked.1$nn.idx))
  i <- ((1:length(j)) - 1)%/%n.neighbors + 1
  x=as.numeric(t(affinities))
  adj = sparseMatrix(i = i, j = j, x = x, dims = c(nrow(batchPCAdata),nrow(batchPCAdata)))
  rownames(adj) <- row.names(batchPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  adj= .extra_matrix_rowNorm(adj)#Matrix::Diagonal(x=1/rowSums(adj)) %*% adj
  
  
  
  load(.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
  adj_t=t(adj)
  adj=adj[sl_pseudo$pseudocell,]
  row.names(adj)=sl_pseudo$cluster
  adj=adj %*% adj_t
  
  #adj_t= Matrix::Diagonal(x = 1 / rowSums(adj_t)) %*% adj_t
  
  adj <- .extra_matrix_rowNorm(adj)#Matrix::Diagonal(x = 1 / rowSums(adj)) %*% adj
  
  prop_mat=centroid.affinities=adj
  
  if(F){
    colMax_vals=c()
    for(i in seq(1,ncol(prop_mat),5000)){
      tmp_max=as.numeric(qlcMatrix::colMax(prop_mat[,i:min(i+4999,ncol(prop_mat))]))
      colMax_vals=c(colMax_vals,as.numeric(tmp_max))
    }
    prop_mat = .extra_matrix_colNorm(input_mat = prop_mat,colValues = 1/colMax_vals)#prop_mat %*% Matrix::Diagonal(x = 1 / colMax_vals)
    prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / rowSums(prop_mat)) %*% prop_mat
  }
  
  
  prop_mat2=myL2normFn(inputMat=prop_mat)
  c_c_aff=t(prop_mat2)
  c_c_aff=prop_mat2 %*% c_c_aff
  #c_c_aff=c_c_aff %*%c_c_aff %*% c_c_aff
  
  #c_c_aff=Matrix::drop0(c_c_aff,tol=quantile(c_c_aff@x,0.1))
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff_t=c_c_aff2=c_c_aff
  
  inc_check=T
  inc_count=nrow(prop_mat)
  counter=0
  while((inc_check|counter<10)&counter<15){
    
    counter=counter+1
    
    c_c_aff=c_c_aff %*% c_c_aff_t
    c_c_aff2=(c_c_aff2+c_c_aff)
    c_c_aff=Matrix::drop0(c_c_aff,tol=quantile(c_c_aff@x,0.00001))
    c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
    
    tst=as.numeric(qlcMatrix::rowMax(c_c_aff))
    tst2=as.numeric(diag(c_c_aff))
    tst=as.numeric(qlcMatrix::rowMax(c_c_aff))
    tst=.extra_matrix_rowNorm(input_mat = c_c_aff,rowValues = 1/tst)#Matrix::Diagonal(x=1/tst) %*% c_c_aff
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
    
    length(unique(pd$anno_orig_cellState[pd$sample %in% sl_pseudo$pseudocell]))
    setdiff(pd$anno_orig_cellState,unique(pd$anno_orig_cellState[pd$sample %in% sl_pseudo$pseudocell]))
    
    
    adj2=prop_mat
    
    #adj2=t(adj_all)
    
    adj2=.extra_matrix_rowNorm(adj2)#Matrix::Diagonal(x=1/rowSums(adj2)) %*% adj2
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    pd=pd[colnames(adj2),]
    
    table(is.na(pd$anno_orig_cellState))
    tmp= adj2 %*% as.matrix(.myOneHotFn(pd$anno_orig_cellState))
    
    tmp=apply(tmp,1,function(x){
      y=which(x==max(x))[1]
      tmp2=data.frame(cluster=colnames(tmp)[y],purity=x[y])
      tmp2
    })
    
    tmp=do.call("rbind",tmp)
    tmp2=aggregate(purity~cluster,data=tmp,median)
    tmp$sl="no"
    tmp$sl[unique(sl_ps_list)]="yes"
    aggregate(purity~sl,data=tmp,summary)
  }

  
  c_c_aff=c_c_aff2[sl_ps_list,]
  c_c_aff=.extra_matrix_rowNorm(input_mat = c_c_aff,rowValues = Matrix::Diagonal(x=1/qlcMatrix::rowMax(c_c_aff)))#Matrix::Diagonal(x=1/qlcMatrix::rowMax(c_c_aff)) %*% c_c_aff
  #c_c_aff=Matrix::drop0(c_c_aff,tol=0.05)
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  
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
    tmp2=.extra_matrix_rowNorm(input_mat = tmp,rowValues = 1/rowMax_vals)#Matrix::Diagonal(x = 1 / rowMax_vals) %*% tmp
    tmp2=Matrix::drop0(tmp2,tol=0.1)
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
    rm(tmp)
  }
  
  rm(prop_mat_list)
  gc()
  
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  matEffectiveSize=.myEffSizePropMat(prop_mat)
  matEffectiveSize=matEffectiveSize$effective_sample_size
  summary(matEffectiveSize)
  
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  prop_mat=prop_mat[matEffectiveSize>=argList$min_cluster_size,]
  prop_mat=prop_mat[,colSums(prop_mat)>0]
  
  ##############
  
  if(mergePseudocells){
    
    prop_mat2=myL2normFn(inputMat = prop_mat)
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
    
    
    
    ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
    ps_merging=cutree(ps_merging,h=merging_strength)
    while(length(unique(ps_merging))<length(ps_merging)){
      ps_merging=data.frame(pseudocell=row.names(prop_mat),id=as.character(ps_merging),stringsAsFactors = F)
      ps_mid=ps_merging[!duplicated(ps_merging$id),]
      colnames(ps_mid)[1]="mid"
      ps_merging=merge(ps_merging,ps_mid,by="id")
      ps_merging=ps_merging[match(row.names(prop_mat),ps_merging$pseudocell),]
      ps_merging=as.matrix(.myOneHotFn(as.character(ps_merging$mid)))
      ps_merging=as(ps_merging,"dgCMatrix")
      ps_merging=t(ps_merging)
      colnames(ps_merging)=row.names(prop_mat)
      prop_mat=ps_merging %*% prop_mat
      rm(ps_merging)
      prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
      
      
      if(L2Norm){
        prop_mat2=myL2normFn(inputMat = prop_mat)
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
      } else {
        prop_mat2=prop_mat
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
        c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
        diag(c_c_aff)=1
        c_c_aff@x=pmin(c_c_aff@x,1)
        c_c_aff=sqrt(c_c_aff*t(c_c_aff))
      }
      
      ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
      ps_merging=cutree(ps_merging,h=merging_strength)
    }
    
    
    prop_mat2=prop_mat
    prop_mat2@x=rep(1,length(prop_mat2@x))
    prop_mat2=myL2normFn(prop_mat2)
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
    ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
    ps_merging=cutree(ps_merging,h=0.1)
    while(length(unique(ps_merging))<length(ps_merging)){
      ps_merging=data.frame(pseudocell=row.names(prop_mat),id=as.character(ps_merging),stringsAsFactors = F)
      ps_mid=ps_merging[!duplicated(ps_merging$id),]
      colnames(ps_mid)[1]="mid"
      ps_merging=merge(ps_merging,ps_mid,by="id")
      ps_merging=ps_merging[match(row.names(prop_mat),ps_merging$pseudocell),]
      ps_merging=as.matrix(.myOneHotFn(as.character(ps_merging$mid)))
      ps_merging=as(ps_merging,"dgCMatrix")
      ps_merging=t(ps_merging)
      colnames(ps_merging)=row.names(prop_mat)
      prop_mat=ps_merging %*% prop_mat
      rm(ps_merging)
      prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
      
      
      if(L2Norm){
        prop_mat2=prop_mat
        prop_mat2@x=rep(1,length(prop_mat2@x))
        prop_mat2=myL2normFn(inputMat = prop_mat2)
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
      } else {
        prop_mat2=prop_mat
        prop_mat2@x=rep(1,length(prop_mat2@x))
        prop_mat2=prop_mat2
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
        c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
        diag(c_c_aff)=1
        c_c_aff@x=pmin(c_c_aff@x,1)
        c_c_aff=sqrt(c_c_aff*t(c_c_aff))
      }
      
      ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
      ps_merging=cutree(ps_merging,h=0.1)
    }
    
    if(F){
      tmp_pd=pd[match(colnames(prop_mat),row.names(pd)),]
      anno=prop_mat %*% as.matrix(.myOneHotFn(tmp_pd$anno_orig_cellState))
      anno_group=apply(anno,1,function(x) which(x==max(x)))
      anno_group=colnames(anno)[anno_group]
      
      matWeights=.myEffSizePropMat(prop_mat)
      matEffectiveSize=matWeights$effective_sample_size
      
      res_check=NULL
      for(i in ps_merging[duplicated(ps_merging)]){
        ii=which(ps_merging==i)
        for(j in 1:(length(ii)-1)){
          for(k in 2:length(ii)){
            Fscore=anno[ii[j],anno_group[ii[j]]]/anno[ii[j],anno_group[ii[k]]]
            Sscore=anno[ii[k],anno_group[ii[k]]]/anno[ii[k],anno_group[ii[j]]]
            res_check=rbind(res_check,data.frame(Fnode=anno_group[ii[j]],Snode=anno_group[ii[k]],effSize1=matEffectiveSize[ii[j]],effsize2=matEffectiveSize[ii[k]],Find=ii[j],Sind=ii[k],Fscore=Fscore,Sscore=Sscore,stringsAsFactors = F))
          }
        }
      }
      table(res_check$Fnode==res_check$Snode)
      length(unique(ps_merging))
      res_check=res_check[res_check$Fnode!=res_check$Snode,]
      table(res_check$Fscore>1.5|res_check$Sscore>1.5)
      res_check[res_check$Fscore>1.5|res_check$Sscore>1.5,]
    }
  }
  
  cat(paste("                      Number of retined pseudocells:",nrow(prop_mat),"out of",nrow(centroidPCAdata),"\n"))
  
  if(F){
    prop_mat2=prop_mat %*% adj_t[colnames(prop_mat),]
    prop_mat <- prop_mat2
  }
  
  ##################
  
  connected_cells=colnames(prop_mat)
  
  
  cat("                      Handling singletons by fast method\n")
  if(T){
    if(length(setdiff(row.names(batchPCAdata),connected_cells))>0){
      nn.ranked <- Seurat:::NNHelper(query = batchPCAdata[setdiff(row.names(batchPCAdata),connected_cells),], data = batchPCAdata[connected_cells,], k = 1, 
                                     method = "annoy", n.trees = 50, searchtype = "standard", 
                                     eps = 0, metric = "euclidean", cache.index = F, 
                                     index = NULL)
      nn.ranked = Indices(nn.ranked)
      nn.ranked.arranged=data.frame(sample=connected_cells,map=connected_cells)
      for(icon in 1:ncol(nn.ranked)){
        nn.ranked.arranged=rbind(nn.ranked.arranged,data.frame(sample=setdiff(row.names(batchPCAdata),connected_cells),map=connected_cells[nn.ranked[,icon]],stringsAsFactors = F))
      }
      nn.ranked=nn.ranked.arranged
    } else {
      nn.ranked=data.frame(sample=colnames(prop_mat),map=colnames(prop_mat),stringsAsFactors = F)
    }
    
    nn.ranked=merge(nn.ranked,data.frame(map=colnames(prop_mat),id_map=1:ncol(prop_mat),stringsAsFactors = F),by="map")
    nn.ranked=merge(nn.ranked,data.frame(sample=colnames(adj),id_sample=1:ncol(adj),stringsAsFactors = F),by="sample")
    #nn.ranked=nn.ranked[match(colnames(prop_mat),nn.ranked$sample),]
    adj_all=Matrix::sparseMatrix(i=nn.ranked$id_sample,j=nn.ranked$id_map,x=rep(1,nrow(nn.ranked)))
    row.names(adj_all)=colnames(adj)
    colnames(adj_all)=colnames(prop_mat)
    
    adj_all=.extra_matrix_rowNorm(adj_all)#Matrix::Diagonal(x=1/rowSums(adj_all)) %*% adj_all
    rm(nn.ranked)
    
  } else {
    snn.matrix <-Seurat:::ComputeSNN(nn_ranked = nn.ranked.1$nn.idx, prune = 1/15)
    rownames(x = snn.matrix) <- rownames(batchPCAdata)
    colnames(x = snn.matrix) <- rownames(batchPCAdata)
    tmp.snn=snn.matrix[,connected_cells]
    tmp.snn=.extra_matrix_rowNorm(input_mat = tmp.snn,rowValues = 1/as.numeric(qlcMatrix::rowMax(tmp.snn)))#Matrix::Diagonal(x=1/as.numeric(qlcMatrix::rowMax(tmp.snn))) %*% tmp.snn
    adj_all=Matrix::drop0(tmp.snn,tol=0.9)
    adj_all=.extra_matrix_rowNorm(adj_all)#Matrix::Diagonal(x=1/rowSums(adj_all)) %*% adj_all
  }
  
  
  
  ###############
  #New
  #Excluded the col normalization step
  #prop_mat <- prop_mat %*% Matrix::Diagonal(x = 1 / (colSums(prop_mat)+0.000000000001))
  adj_all=adj_all %*% t(prop_mat)
  #adj_all = Matrix::Diagonal(x = 1 / (rowSums(adj_all)+0.000000000001)) %*% adj_all
  prop_mat=t(adj_all)
  
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  #prop_mat =  Matrix::Diagonal(x = 1 / qlcMatrix::rowMax(prop_mat)) %*% prop_mat
  colMax_vals=c()
  for(i in seq(1,ncol(prop_mat),5000)){
    tmp_max=as.numeric(qlcMatrix::colMax(prop_mat[,i:min(i+4999,ncol(prop_mat))]))
    colMax_vals=c(colMax_vals,as.numeric(tmp_max))
  }
  prop_mat = .extra_matrix_colNorm(input_mat = prop_mat,colValues = 1/colMax_vals)#prop_mat %*% Matrix::Diagonal(x = 1 / colMax_vals)
  #prop_mat=Matrix::drop0(prop_mat,0.1)
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / rowSums(prop_mat)) %*% prop_mat
  
  matWeights=.myEffSizePropMat(prop_mat)
  
  matEffectiveSize=matWeights$effective_sample_size
  matWeights=matWeights$centroid_weights
  matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  
  
  for(i in 1:length(dataArranged)){
    tmp_prop_mat=prop_mat[,match(row.names(dataArranged[[i]]$pcaData),colnames(prop_mat))]
    #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
    tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
    tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
    tmp_weights[tmp_effsize<4]=0
    tmp=list(prop_mat=tmp_prop_mat,data=dataArranged[[i]],pseudocell_sim_mat=c_c_aff,matWeights=tmp_weights,matEffectiveSize=tmp_effsize)
    dataArranged[[i]]=tmp
  }
  
  return(list(dataArranged=dataArranged,prop_mat=prop_mat))
}

.myConcensusDEFn_step2_detail_newprop3_final_v14=function(dataArranged,centroidPCAdata,argList,input_pca_centroids_assignments=NULL,exCentroids=NULL,n.neighbors=4,batchPCAdata=NULL,n.trees=50,k.param = 20,NNmethod="annoy",L2Norm=T,mergePseudocells=T,batch_variable="anno_batch",merging_strength=0.3,hierarchical_mode=F,...){
  require(Matrix)
  
  myL2normFn=function(inputMat){
    prop_mat2=rowSums(inputMat^2)
    prop_mat2=sqrt(prop_mat2)
    res=.extra_matrix_rowNorm(input_mat = inputMat,rowValues = 1/(prop_mat2+0.00000001))#Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
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
    if(!file.exists(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))|hierarchical_mode){
      idx=Seurat:::AnnoyBuildIndex(data = batchPCAdata, metric = "euclidean", 
                                   n.trees = n.trees)
      nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = batchPCAdata,k=k.param,include.distance = T,search.k = -1)
      if(!hierarchical_mode){
        qsave(nn.ranked.1,file=.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
      }
      
    } else {
      nn.ranked.1=qread(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
    }
  } else {
    nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  }
  
  
  affinities=.extra_matrix_rowNorm(input_mat = nn.ranked.1$nn.dists,rowValues = 1/(nn.ranked.1$nn.dists[,2]+0.000001))#Matrix::Diagonal(x=1/(nn.ranked.1$nn.dists[,2]+0.000001)) %*% nn.ranked.1$nn.dists
  affinities@x=-1*affinities@x^2
  affinities@x=exp(affinities@x)
  affinities[,1]=affinities[,2]
  j <- as.numeric(t(nn.ranked.1$nn.idx))
  i <- ((1:length(j)) - 1)%/%n.neighbors + 1
  x=as.numeric(t(affinities))
  adj = sparseMatrix(i = i, j = j, x = x, dims = c(nrow(batchPCAdata),nrow(batchPCAdata)))
  rownames(adj) <- row.names(batchPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  adj=.extra_matrix_rowNorm(input_mat=adj)
  
  
  if(is.null(input_pca_centroids_assignments)){
    load(.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
  } else {
    sl_pseudo=input_pca_centroids_assignments
  }
  
  adj_t=.extra_matrix_rowNorm(t(adj))
  adj=adj[sl_pseudo$pseudocell,]
  row.names(adj)=sl_pseudo$cluster
  adj=adj %*% adj_t
  
  if(ncol(batchPCAdata)<100000){
    adj_sim = sparseMatrix(i = i, j = j, x = 1, dims = c(nrow(batchPCAdata),nrow(batchPCAdata)))
    rownames(adj_sim) <- row.names(batchPCAdata)
    colnames(adj_sim)=c(row.names(batchPCAdata))
    adj_sim= .extra_matrix_rowNorm(adj_sim)#Matrix::Diagonal(x=1/rowSums(adj_sim)) %*% adj_sim
    
    if(is.null(input_pca_centroids_assignments)){
      load(.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
    } else {
      sl_pseudo=input_pca_centroids_assignments
    }
    adj_t_sim=t(adj_sim)
    adj_t_sim=.extra_matrix_rowNorm(adj_t_sim)#Matrix::Diagonal(x=1/rowSums(adj_t_sim)) %*% adj_t_sim
    adj_sim=adj_sim[sl_pseudo$pseudocell,]
    row.names(adj_sim)=sl_pseudo$cluster
    adj_sim=adj_sim %*% adj_t_sim
    adj_sim=.extra_sconline.pseudosim(adj=adj_sim,adj_t=adj_t_sim,argList=argList)
    
  } else {
    adj_sim=.extra_sconline.pseudosim(adj=adj,adj_t=adj_t,argList=argList)
  }
  
 
  adj <- .extra_matrix_rowNorm(adj)#Matrix::Diagonal(x = 1 / rowSums(adj)) %*% adj
  
  prop_mat=centroid.affinities=adj
  
  if(F){
    colMax_vals=c()
    for(i in seq(1,ncol(prop_mat),5000)){
      tmp_max=as.numeric(qlcMatrix::colMax(prop_mat[,i:min(i+4999,ncol(prop_mat))]))
      colMax_vals=c(colMax_vals,as.numeric(tmp_max))
    }
    prop_mat =  .extra_matrix_colNorm(input_mat = prop_mat,colValues = 1/colMax_vals)# prop_mat %*% Matrix::Diagonal(x = 1 / colMax_vals)
    prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / rowSums(prop_mat)) %*% prop_mat
  }
  
  
  prop_mat2=myL2normFn(inputMat=prop_mat)
  c_c_aff=t(prop_mat2)
  c_c_aff=prop_mat2 %*% c_c_aff
  #c_c_aff=c_c_aff %*%c_c_aff %*% c_c_aff
  
  #c_c_aff=Matrix::drop0(c_c_aff,tol=quantile(c_c_aff@x,0.1))
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff_t=c_c_aff2=c_c_aff
  
  inc_check=T
  inc_count=nrow(prop_mat)
  counter=0
  while((inc_check|counter<10)&counter<15){
    
    counter=counter+1
    
    c_c_aff=c_c_aff %*% c_c_aff_t*0.7+c_c_aff_t *0.3
    
    c_c_aff=Matrix::drop0(c_c_aff,tol=quantile(c_c_aff@x,0.00001))
    c_c_aff <- .extra_matrix_rowNorm(c_c_aff)# Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
    
    tst=as.numeric(qlcMatrix::rowMax(c_c_aff))
    tst2=as.numeric(diag(c_c_aff))
    tst=as.numeric(qlcMatrix::rowMax(c_c_aff))
    tst=.extra_matrix_rowNorm(input_mat = c_c_aff,rowValues = 1/tst)#Matrix::Diagonal(x=1/tst) %*% c_c_aff
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
    
    length(unique(pd$anno_orig_cellState[pd$sample %in% sl_pseudo$pseudocell]))
    setdiff(pd$anno_orig_cellState,unique(pd$anno_orig_cellState[pd$sample %in% sl_pseudo$pseudocell]))
    
    
    #adj2=diff2
    
    adj2=prop_mat#=res$prop_mat
    #adj2=t(adj_all)
    
    adj2=.extra_matrix_rowNorm(adj2)#Matrix::Diagonal(x=1/rowSums(adj2)) %*% adj2
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    pd=pd[colnames(adj2),]
    
    table(is.na(pd$anno_orig_cellState))
    tmp= adj2 %*% as.matrix(.myOneHotFn(pd$anno_orig_cellState))
    
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
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)# Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  
  #c_c_aff=Matrix::drop0(c_c_aff,tol=0.05)
  prop_mat_list=list()
  rowMax_vals=c()
  for(i in seq(1,ncol(prop_mat),50000)){
    tmp=c_c_aff %*% prop_mat[colnames(c_c_aff),i:min(i+49999,ncol(prop_mat))]
    prop_mat_list=c(prop_mat_list,list(tmp))
    tmp_max=Matrix::drop0(tmp,tol=0.0001)
    tmp_max=as.numeric(qlcMatrix::rowMax(tmp_max))
    if(length(rowMax_vals)>0){
      rowMax_vals=pmax(rowMax_vals,tmp_max)
    } else {
      rowMax_vals=tmp_max
    }
  }
  
  if(sum(rowMax_vals<=0)>0){
    sl_rows=row.names(prop_mat)[which(rowMax_vals==0)]
    tmp_rowMax_vals=c()
    for(i in seq(1,ncol(prop_mat),50000)){
      tmp=c_c_aff[sl_rows,,drop=F] %*% prop_mat[,i:min(i+49999,ncol(prop_mat))]
      tmp_max=as.numeric(qlcMatrix::rowMax(tmp))
      if(length(tmp_rowMax_vals)>0){
        tmp_rowMax_vals=pmax(tmp_rowMax_vals,tmp_max)
      } else {
        tmp_rowMax_vals=tmp_max
      }
    }
    rowMax_vals[match(sl_rows,colnames(c_c_aff))]=tmp_rowMax_vals
  }
  
  
  for(i in 1:length(prop_mat_list)){
    tmp=prop_mat_list[[i]]
    tmp=Matrix::drop0(tmp,tol=min(rowMax_vals)*0.01)
    #tmp=tmp %*% Matrix::Diagonal(x=1/(qlcMatrix::colMax(tmp)+0.00001))
    tmp2=.extra_matrix_rowNorm(input_mat = tmp,rowValues = 1/rowMax_vals)#Matrix::Diagonal(x = 1 / rowMax_vals) %*% tmp
    tmp2=Matrix::drop0(tmp2,tol=0.05)
    tmp2@x=rep(1,length(tmp2@x))
    row.names(tmp2)=row.names(tmp)
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
    rm(tmp)
  }
  
  rm(prop_mat_list)
  gc()
  
  prop_mat <- .extra_matrix_rowNorm(prop_mat)# Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  matEffectiveSize=.myEffSizePropMat(prop_mat)
  matEffectiveSize=matEffectiveSize$effective_sample_size
  
  
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  prop_mat=prop_mat[matEffectiveSize>=argList$min_cluster_size,]
  prop_mat=prop_mat[,colSums(prop_mat)>0]
  
  ##############
  
  if(mergePseudocells){
    
    prop_mat2=myL2normFn(inputMat = prop_mat)
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
    
    ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
    ps_merging=cutree(ps_merging,h=merging_strength)
    while(length(unique(ps_merging))<length(ps_merging)){
      ps_merging=data.frame(pseudocell=row.names(prop_mat),id=as.character(ps_merging),stringsAsFactors = F)
      ps_mid=ps_merging[!duplicated(ps_merging$id),]
      colnames(ps_mid)[1]="mid"
      ps_merging=merge(ps_merging,ps_mid,by="id")
      ps_merging=ps_merging[match(row.names(prop_mat),ps_merging$pseudocell),]
      ps_merging=as.matrix(.myOneHotFn(as.character(ps_merging$mid)))
      ps_merging=as(ps_merging,"dgCMatrix")
      ps_merging=t(ps_merging)
      colnames(ps_merging)=row.names(prop_mat)
      prop_mat=ps_merging %*% prop_mat
      rm(ps_merging)
      prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
      
      
      if(L2Norm){
        prop_mat2=myL2normFn(inputMat = prop_mat)
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
      } else {
        prop_mat2=prop_mat
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
        c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
        diag(c_c_aff)=1
        c_c_aff@x=pmin(c_c_aff@x,1)
        c_c_aff=sqrt(c_c_aff*t(c_c_aff))
      }
      
      ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
      ps_merging=cutree(ps_merging,h=merging_strength)
    }
    
    
    prop_mat2=prop_mat
    prop_mat2@x=rep(1,length(prop_mat2@x))
    prop_mat2=myL2normFn(prop_mat2)
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
    ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
    ps_merging=cutree(ps_merging,h=0.15)
    while(length(unique(ps_merging))<length(ps_merging)){
      ps_merging=data.frame(pseudocell=row.names(prop_mat),id=as.character(ps_merging),stringsAsFactors = F)
      ps_mid=ps_merging[!duplicated(ps_merging$id),]
      colnames(ps_mid)[1]="mid"
      ps_merging=merge(ps_merging,ps_mid,by="id")
      ps_merging=ps_merging[match(row.names(prop_mat),ps_merging$pseudocell),]
      ps_merging=as.matrix(.myOneHotFn(as.character(ps_merging$mid)))
      ps_merging=as(ps_merging,"dgCMatrix")
      ps_merging=t(ps_merging)
      colnames(ps_merging)=row.names(prop_mat)
      prop_mat=ps_merging %*% prop_mat
      rm(ps_merging)
      prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
      
      
      if(L2Norm){
        prop_mat2=prop_mat
        prop_mat2@x=rep(1,length(prop_mat2@x))
        prop_mat2=myL2normFn(inputMat = prop_mat2)
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
      } else {
        prop_mat2=prop_mat
        prop_mat2@x=rep(1,length(prop_mat2@x))
        prop_mat2=prop_mat2
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
        c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
        diag(c_c_aff)=1
        c_c_aff@x=pmin(c_c_aff@x,1)
        c_c_aff=sqrt(c_c_aff*t(c_c_aff))
      }
      
      ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
      ps_merging=cutree(ps_merging,h=0.15)
    }
    
    if(F){
      tmp_pd=pd[match(colnames(prop_mat),row.names(pd)),]
      anno=prop_mat %*% as.matrix(.myOneHotFn(tmp_pd$anno_orig_cellState))
      anno_group=apply(anno,1,function(x) which(x==max(x)))
      anno_group=colnames(anno)[anno_group]
      
      matWeights=.myEffSizePropMat(prop_mat)
      matEffectiveSize=matWeights$effective_sample_size
      
      res_check=NULL
      for(i in ps_merging[duplicated(ps_merging)]){
        ii=which(ps_merging==i)
        for(j in 1:(length(ii)-1)){
          for(k in 2:length(ii)){
            Fscore=anno[ii[j],anno_group[ii[j]]]/anno[ii[j],anno_group[ii[k]]]
            Sscore=anno[ii[k],anno_group[ii[k]]]/anno[ii[k],anno_group[ii[j]]]
            res_check=rbind(res_check,data.frame(Fnode=anno_group[ii[j]],Snode=anno_group[ii[k]],effSize1=matEffectiveSize[ii[j]],effsize2=matEffectiveSize[ii[k]],Find=ii[j],Sind=ii[k],Fscore=Fscore,Sscore=Sscore,stringsAsFactors = F))
          }
        }
      }
      table(res_check$Fnode==res_check$Snode)
      length(unique(ps_merging))
      res_check=res_check[res_check$Fnode!=res_check$Snode,]
      table(res_check$Fscore>1.5|res_check$Sscore>1.5)
      res_check[res_check$Fscore>1.5|res_check$Sscore>1.5,]
    }
  }
  
  cat(paste("                      Number of retined pseudocells:",nrow(prop_mat),"out of",nrow(centroidPCAdata),"\n"))
  
  if(F){
    prop_mat2=prop_mat %*% adj_t[colnames(prop_mat),]
    prop_mat <- prop_mat2
  }
  
  ##################
  
  if(argList$include.singletons){
    cat("                      Handling singletons by fast method\n")
    connected_cells=colnames(prop_mat)
    counter=0
    snn.matrix <-Seurat:::ComputeSNN(nn_ranked = nn.ranked.1$nn.idx, prune = 1/15)
    rownames(x = snn.matrix) <- rownames(batchPCAdata)
    colnames(x = snn.matrix) <- rownames(batchPCAdata)
    while(length(connected_cells)<nrow(nn.ranked.1$nn.idx)&counter<4){
      counter=counter+1
      if(F){
        if(length(setdiff(row.names(batchPCAdata),connected_cells))>0){
          nn.ranked <- Seurat:::NNHelper(query = batchPCAdata[setdiff(row.names(batchPCAdata),connected_cells),], data = batchPCAdata[connected_cells,], k = 1, 
                                         method = "annoy", n.trees = 50, searchtype = "standard", 
                                         eps = 0, metric = "euclidean", cache.index = F, 
                                         index = NULL)
          nn.ranked = Indices(nn.ranked)
          nn.ranked.arranged=data.frame(sample=connected_cells,map=connected_cells)
          for(icon in 1:ncol(nn.ranked)){
            nn.ranked.arranged=rbind(nn.ranked.arranged,data.frame(sample=setdiff(row.names(batchPCAdata),connected_cells),map=connected_cells[nn.ranked[,icon]],stringsAsFactors = F))
          }
          nn.ranked=nn.ranked.arranged
        } else {
          nn.ranked=data.frame(sample=colnames(prop_mat),map=colnames(prop_mat),stringsAsFactors = F)
        }
        
        nn.ranked=merge(nn.ranked,data.frame(map=colnames(prop_mat),id_map=1:ncol(prop_mat),stringsAsFactors = F),by="map")
        nn.ranked=merge(nn.ranked,data.frame(sample=colnames(adj),id_sample=1:ncol(adj),stringsAsFactors = F),by="sample")
        #nn.ranked=nn.ranked[match(colnames(prop_mat),nn.ranked$sample),]
        adj_all=Matrix::sparseMatrix(i=nn.ranked$id_sample,j=nn.ranked$id_map,x=rep(1,nrow(nn.ranked)))
        row.names(adj_all)=colnames(adj)
        colnames(adj_all)=colnames(prop_mat)
        
        adj_all=.extra_matrix_rowNorm(adj_all)#Matrix::Diagonal(x=1/rowSums(adj_all)) %*% adj_all
        rm(nn.ranked)
        
      } else {
        
        tmp.snn=snn.matrix[,connected_cells]
        tmp.snn=.extra_matrix_rowNorm(input_mat = tmp.snn,rowValues = 1/as.numeric(qlcMatrix::rowMax(tmp.snn)))#Matrix::Diagonal(x=1/as.numeric(qlcMatrix::rowMax(tmp.snn))) %*% tmp.snn
        adj_all=Matrix::drop0(tmp.snn,tol=0.8)
        
        adj_all=.extra_matrix_rowNorm(adj_all)#Matrix::Diagonal(x=1/rowSums(adj_all)) %*% adj_all
      }
      
      adj_all=adj_all %*% t(prop_mat)
      #adj_all = Matrix::Diagonal(x = 1 / (rowSums(adj_all)+0.000000000001)) %*% adj_all
      prop_mat=t(adj_all)
      
      prop_mat <- .extra_matrix_rowNorm(input_mat=prop_mat,rowValues = 1 / (rowSums(prop_mat)+0.000000000001))#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
      
      prop_mat=prop_mat[,colSums(prop_mat)>0]
      
      connected_cells=colnames(prop_mat)
    }
    
    {
      if(length(setdiff(row.names(batchPCAdata),connected_cells))>0){
        
        if(length(setdiff(row.names(batchPCAdata),connected_cells))==1){
          connected_cells=connected_cells[-1]
        }
        
        nn.ranked <- Seurat:::NNHelper(query = batchPCAdata[setdiff(row.names(batchPCAdata),connected_cells),], data = batchPCAdata[connected_cells,], k = 1, 
                                       method = "annoy", n.trees = 50, searchtype = "standard", 
                                       eps = 0, metric = "euclidean", cache.index = F, 
                                       index = NULL)
        nn.ranked = Indices(nn.ranked)
        nn.ranked.arranged=data.frame(sample=connected_cells,map=connected_cells)
        for(icon in 1:ncol(nn.ranked)){
          nn.ranked.arranged=rbind(nn.ranked.arranged,data.frame(sample=setdiff(row.names(batchPCAdata),connected_cells),map=connected_cells[nn.ranked[,icon]],stringsAsFactors = F))
        }
        nn.ranked=nn.ranked.arranged
      } else {
        nn.ranked=data.frame(sample=colnames(prop_mat),map=colnames(prop_mat),stringsAsFactors = F)
      }
      
      nn.ranked=merge(nn.ranked,data.frame(map=colnames(prop_mat),id_map=1:ncol(prop_mat),stringsAsFactors = F),by="map")
      nn.ranked=merge(nn.ranked,data.frame(sample=colnames(adj),id_sample=1:ncol(adj),stringsAsFactors = F),by="sample")
      #nn.ranked=nn.ranked[match(colnames(prop_mat),nn.ranked$sample),]
      adj_all=Matrix::sparseMatrix(i=nn.ranked$id_sample,j=nn.ranked$id_map,x=rep(1,nrow(nn.ranked)))
      row.names(adj_all)=colnames(adj)
      colnames(adj_all)=colnames(prop_mat)
      
      adj_all=.extra_matrix_rowNorm(adj_all)#Matrix::Diagonal(x=1/rowSums(adj_all)) %*% adj_all
      rm(nn.ranked)
      
      adj_all=adj_all %*% t(prop_mat)
      #adj_all = Matrix::Diagonal(x = 1 / (rowSums(adj_all)+0.000000000001)) %*% adj_all
      prop_mat=t(adj_all)
      
      prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
      
    }
  }
  
  
  
  
  
  ###############
  #New
  #Excluded the col normalization step
  #prop_mat <- prop_mat %*% Matrix::Diagonal(x = 1 / (colSums(prop_mat)+0.000000000001))
  
  
  if(argList$colNormalize){
    adj2=adj[row.names(prop_mat),colnames(prop_mat)]
    adj2=.extra_matrix_rowNorm(adj2)#Matrix::Diagonal(x=1/rowSums(adj2)) %*% adj2
    adj2=Matrix::drop0(prop_mat * adj2)
    
    rowMeans_drop0 <- function (dgCMat) {
      RowInd <- dgCMat@i + 1
      sapply(split(dgCMat@x, RowInd), mean)
    }
    
    prop_mat = .extra_matrix_rowNorm(input_mat = prop_mat,rowValues = 1 / (rowMeans_drop0(adj2)))#  Matrix::Diagonal(x = 1 / (rowMeans_drop0(adj2))) %*% prop_mat
    #prop_mat@x=pmin(prop_mat@x,1)
    colMax_vals=c()
    for(i in seq(1,ncol(prop_mat),5000)){
      tmp_max=as.numeric(qlcMatrix::colMax(prop_mat[,i:min(i+4999,ncol(prop_mat))]))
      colMax_vals=c(colMax_vals,as.numeric(tmp_max))
    }
    prop_mat =.extra_matrix_colNorm(input_mat = prop_mat,colValues = 1/colMax_vals)  #prop_mat %*% Matrix::Diagonal(x = 1 / colMax_vals)
    #prop_mat=Matrix::drop0(prop_mat,0.05)
    prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / rowSums(prop_mat)) %*% prop_mat
  }
  
  
  if(ncol(adj)>ncol(prop_mat)){
    tmp_c=adj[row.names(prop_mat),!colnames(adj) %in% colnames(prop_mat)]
    tmp_c@x=rep(0,length(tmp_c@x))
    prop_mat=cbind(prop_mat,tmp_c)
    prop_mat=Matrix::drop0(prop_mat)
  }
  
  matWeights=.myEffSizePropMat(prop_mat)
  
  matEffectiveSize=matWeights$effective_sample_size
  summary(matEffectiveSize)
  
  matWeights=matWeights$centroid_weights
  matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  
  
  for(i in 1:length(dataArranged)){
    tmp_prop_mat=prop_mat[,match(colnames(dataArranged[[i]]$logNorm),colnames(prop_mat))]
    #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
    tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
    tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
    tmp_weights[tmp_effsize<4]=0
    tmp=list(prop_mat=tmp_prop_mat,data=dataArranged[[i]],pseudocell_sim_mat=c_c_aff,matWeights=tmp_weights,matEffectiveSize=tmp_effsize)
    dataArranged[[i]]=tmp
  }
  
  return(list(dataArranged=dataArranged,prop_mat=prop_mat,pseudocell_sim_mat=adj_sim))
}

.myConcensusDEFn_step2_detail_newprop3_final_v15=function(dataArranged,centroidPCAdata,argList,exCentroids=NULL,n.neighbors=4,batchPCAdata=NULL,n.trees=50,k.param = 20,NNmethod="annoy",L2Norm=T,mergePseudocells=T,batch_variable="anno_batch",merging_strength=0.15,...){
  
  myL2normFn=function(inputMat){
    prop_mat2=rowSums(inputMat^2)
    prop_mat2=sqrt(prop_mat2)
    res=.extra_matrix_rowNorm(input_mat = inputMat,rowValues = 1/(prop_mat2+0.00000001))#Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
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
  
  
  affinities=.extra_matrix_rowNorm(input_mat = nn.ranked.1$nn.dists,rowValues =Matrix::Diagonal(x=1/(nn.ranked.1$nn.dists[,2]+0.000001)))#Matrix::Diagonal(x=1/(nn.ranked.1$nn.dists[,2]+0.000001)) %*% nn.ranked.1$nn.dists
  affinities@x=-1*affinities@x^2
  affinities@x=exp(affinities@x)
  affinities[,1]=affinities[,2]
  j <- as.numeric(t(nn.ranked.1$nn.idx))
  i <- ((1:length(j)) - 1)%/%n.neighbors + 1
  x=as.numeric(t(affinities))
  adj = sparseMatrix(i = i, j = j, x = x, dims = c(nrow(batchPCAdata),nrow(batchPCAdata)))
  rownames(adj) <- row.names(batchPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  adj= .extra_matrix_rowNorm(adj)#Matrix::Diagonal(x=1/rowSums(adj)) %*% adj
  
  
  
  load(.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
  adj_t=t(adj)
  adj=adj[,sl_pseudo$pseudocell]
  colnames(adj)=sl_pseudo$cluster
  adj=adj_t %*% adj
  
  #adj_t= Matrix::Diagonal(x = 1 / rowSums(adj_t)) %*% adj_t
  
  adj <- .extra_matrix_rowNorm(adj)#Matrix::Diagonal(x = 1 / rowSums(adj)) %*% adj
  
  prop_mat=t(adj)
  centroid.affinities=prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / rowSums(prop_mat)) %*% prop_mat
  
  if(F){
    colMax_vals=c()
    for(i in seq(1,ncol(prop_mat),5000)){
      tmp_max=as.numeric(qlcMatrix::colMax(prop_mat[,i:min(i+4999,ncol(prop_mat))]))
      colMax_vals=c(colMax_vals,as.numeric(tmp_max))
    }
    prop_mat = .extra_matrix_colNorm(input_mat = prop_mat,colValues = 1/colMax_vals)#prop_mat %*% Matrix::Diagonal(x = 1 / colMax_vals)
    prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / rowSums(prop_mat)) %*% prop_mat
  }
  
  
  prop_mat2=myL2normFn(inputMat=prop_mat)
  c_c_aff=t(prop_mat2)
  c_c_aff=prop_mat2 %*% c_c_aff
  #c_c_aff=c_c_aff %*%c_c_aff %*% c_c_aff
  
  #c_c_aff=Matrix::drop0(c_c_aff,tol=quantile(c_c_aff@x,0.1))
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff_t=c_c_aff2=c_c_aff
  
  inc_check=T
  inc_count=nrow(prop_mat)
  counter=0
  while((inc_check|counter<10)&counter<15){
    
    counter=counter+1
    
    c_c_aff=c_c_aff %*% c_c_aff_t*0.5+c_c_aff_t *0.5
    
    c_c_aff=Matrix::drop0(c_c_aff,tol=quantile(c_c_aff@x,0.00001))
    c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
    
    tst=as.numeric(qlcMatrix::rowMax(c_c_aff))
    tst2=as.numeric(diag(c_c_aff))
    tst=as.numeric(qlcMatrix::rowMax(c_c_aff))
    tst=.extra_matrix_rowNorm(input_mat = c_c_aff,rowValues = 1/tst)#Matrix::Diagonal(x=1/tst) %*% c_c_aff
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
    
    length(unique(pd$anno_orig_cellState[pd$sample %in% sl_pseudo$pseudocell]))
    setdiff(pd$anno_orig_cellState,unique(pd$anno_orig_cellState[pd$sample %in% sl_pseudo$pseudocell]))
    
    
    #adj2=prop_mat
    
    adj2=prop_mat
    #adj2=t(adj_all)
    
    adj2=.extra_matrix_rowNorm(adj2)#Matrix::Diagonal(x=1/rowSums(adj2)) %*% adj2
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    pd=pd[colnames(adj2),]
    
    table(is.na(pd$anno_orig_cellState))
    tmp= adj2 %*% as.matrix(.myOneHotFn(pd$anno_orig_cellState))
    
    tmp=apply(tmp,1,function(x){
      y=which(x==max(x))[1]
      tmp2=data.frame(cluster=colnames(tmp)[y],purity=x[y])
      tmp2
    })
    
    tmp=do.call("rbind",tmp)
    tmp2=aggregate(purity~cluster,data=tmp,median)
    table(tmp2[,2]>0.8)
    table(tmp2[,2]>0.7)
    table(tmp2[,2]>0.6)
    table(tmp2[,2]>0.5)
    dim(tmp2)
    tmp$sl="no"
    tmp$sl[unique(sl_ps_list)]="yes"
    aggregate(purity~sl,data=tmp,summary)
  }
  
  
  c_c_aff=c_c_aff[sl_ps_list,]
  c_c_aff=.extra_matrix_rowNorm(input_mat = c_c_aff,rowValues = 1/qlcMatrix::rowMax(c_c_aff))#Matrix::Diagonal(x=1/qlcMatrix::rowMax(c_c_aff)) %*% c_c_aff
  #c_c_aff=Matrix::drop0(c_c_aff,tol=0.05)
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  
  prop_mat_list=list()
  rowMax_vals=c()
  for(i in seq(1,ncol(prop_mat),50000)){
    tmp=adj[i:min(i+49999,ncol(prop_mat)),colnames(c_c_aff)] %*% c_c_aff
    tmp_max=as.numeric(qlcMatrix::colMax(tmp))
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
    tmp2=.extra_matrix_colNorm(input_mat = tmp,colValues = 1/rowMax_vals) #tmp %*% Matrix::Diagonal(x = 1 / rowMax_vals)
    tmp2=Matrix::drop0(tmp2,tol=0.05)
    tmp2@x=rep(1,length(tmp2@x))
    tmp=tmp*tmp2
    
    #tmp=Matrix::drop0(tmp,tol=0.1)
    prop_mat_list[[i]]=tmp
  }
  
  rm(tmp,tmp2)
  
  prop_mat=prop_mat_list[[1]]
  prop_mat=prop_mat[rowSums(prop_mat)>0,]
  if(length(prop_mat_list)>1){
    for(i in 2:length(prop_mat_list)){
      tmp=prop_mat_list[[i]]
      tmp=tmp[rowSums(tmp)>0,]
      prop_mat=rbind(prop_mat,tmp)
    }
    rm(tmp)
  }
  
  rm(prop_mat_list)
  gc()
  
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  prop_mat_t=t(prop_mat)
  prop_mat_t=.extra_matrix_rowNorm(prop_mat_t)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat_t)+0.000000000001)) %*% prop_mat_t
  matEffectiveSize=.myEffSizePropMat(prop_mat_t)
  matEffectiveSize=matEffectiveSize$effective_sample_size
  summary(matEffectiveSize)
  
  matEffectiveSize=matEffectiveSize[match(colnames(prop_mat),names(matEffectiveSize))]
  prop_mat=prop_mat[,matEffectiveSize>=argList$min_cluster_size]
  prop_mat=prop_mat[rowSums(prop_mat)>0,]
  
  colMax_vals=c()
  for(i in seq(1,nrow(prop_mat),5000)){
    tmp_max=as.numeric(qlcMatrix::rowMax(prop_mat[i:min(i+4999,ncol(prop_mat)),]))
    colMax_vals=c(colMax_vals,as.numeric(tmp_max))
  }
  prop_mat =  .extra_matrix_colNorm(input_mat = prop_mat,colValues = 1/colMax_vals)#Matrix::Diagonal(x = 1 / colMax_vals) %*% prop_mat
  prop_mat=t(prop_mat)
  prop_mat =  .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / rowSums(prop_mat)) %*% prop_mat
  
  ##############
  
  if(mergePseudocells){
    
    prop_mat2=myL2normFn(inputMat = prop_mat)
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
    
    
    
    ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
    ps_merging=cutree(ps_merging,h=merging_strength)
    while(length(unique(ps_merging))<length(ps_merging)){
      ps_merging=data.frame(pseudocell=row.names(prop_mat),id=as.character(ps_merging),stringsAsFactors = F)
      ps_mid=ps_merging[!duplicated(ps_merging$id),]
      colnames(ps_mid)[1]="mid"
      ps_merging=merge(ps_merging,ps_mid,by="id")
      ps_merging=ps_merging[match(row.names(prop_mat),ps_merging$pseudocell),]
      ps_merging=as.matrix(.myOneHotFn(as.character(ps_merging$mid)))
      ps_merging=as(ps_merging,"dgCMatrix")
      ps_merging=t(ps_merging)
      colnames(ps_merging)=row.names(prop_mat)
      prop_mat=ps_merging %*% prop_mat
      rm(ps_merging)
      prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
      
      
      if(L2Norm){
        prop_mat2=myL2normFn(inputMat = prop_mat)
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
      } else {
        prop_mat2=prop_mat
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
        c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
        diag(c_c_aff)=1
        c_c_aff@x=pmin(c_c_aff@x,1)
        c_c_aff=sqrt(c_c_aff*t(c_c_aff))
      }
      
      ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
      ps_merging=cutree(ps_merging,h=merging_strength)
    }
    
    
    prop_mat2=prop_mat
    prop_mat2@x=rep(1,length(prop_mat2@x))
    prop_mat2=myL2normFn(prop_mat2)
    c_c_aff=t(prop_mat2)
    c_c_aff=prop_mat2 %*% c_c_aff
    ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
    ps_merging=cutree(ps_merging,h=0.15)
    while(length(unique(ps_merging))<length(ps_merging)){
      ps_merging=data.frame(pseudocell=row.names(prop_mat),id=as.character(ps_merging),stringsAsFactors = F)
      ps_mid=ps_merging[!duplicated(ps_merging$id),]
      colnames(ps_mid)[1]="mid"
      ps_merging=merge(ps_merging,ps_mid,by="id")
      ps_merging=ps_merging[match(row.names(prop_mat),ps_merging$pseudocell),]
      ps_merging=as.matrix(.myOneHotFn(as.character(ps_merging$mid)))
      ps_merging=as(ps_merging,"dgCMatrix")
      ps_merging=t(ps_merging)
      colnames(ps_merging)=row.names(prop_mat)
      prop_mat=ps_merging %*% prop_mat
      rm(ps_merging)
      prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
      
      
      if(L2Norm){
        prop_mat2=prop_mat
        prop_mat2@x=rep(1,length(prop_mat2@x))
        prop_mat2=myL2normFn(inputMat = prop_mat2)
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
      } else {
        prop_mat2=prop_mat
        prop_mat2@x=rep(1,length(prop_mat2@x))
        prop_mat2=prop_mat2
        c_c_aff=t(prop_mat2)
        c_c_aff=prop_mat2 %*% c_c_aff
        c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
        diag(c_c_aff)=1
        c_c_aff@x=pmin(c_c_aff@x,1)
        c_c_aff=sqrt(c_c_aff*t(c_c_aff))
      }
      
      ps_merging=hclust(as.dist(1-c_c_aff),method = "complete")
      ps_merging=cutree(ps_merging,h=0.15)
    }
    
    if(F){
      tmp_pd=pd[match(colnames(prop_mat),row.names(pd)),]
      anno=prop_mat %*% as.matrix(.myOneHotFn(tmp_pd$anno_orig_cellState))
      anno_group=apply(anno,1,function(x) which(x==max(x)))
      anno_group=colnames(anno)[anno_group]
      
      matWeights=.myEffSizePropMat(prop_mat)
      matEffectiveSize=matWeights$effective_sample_size
      
      res_check=NULL
      for(i in ps_merging[duplicated(ps_merging)]){
        ii=which(ps_merging==i)
        for(j in 1:(length(ii)-1)){
          for(k in 2:length(ii)){
            Fscore=anno[ii[j],anno_group[ii[j]]]/anno[ii[j],anno_group[ii[k]]]
            Sscore=anno[ii[k],anno_group[ii[k]]]/anno[ii[k],anno_group[ii[j]]]
            res_check=rbind(res_check,data.frame(Fnode=anno_group[ii[j]],Snode=anno_group[ii[k]],effSize1=matEffectiveSize[ii[j]],effsize2=matEffectiveSize[ii[k]],Find=ii[j],Sind=ii[k],Fscore=Fscore,Sscore=Sscore,stringsAsFactors = F))
          }
        }
      }
      table(res_check$Fnode==res_check$Snode)
      length(unique(ps_merging))
      res_check=res_check[res_check$Fnode!=res_check$Snode,]
      table(res_check$Fscore>1.5|res_check$Sscore>1.5)
      res_check[res_check$Fscore>1.5|res_check$Sscore>1.5,]
    }
  }
  
  cat(paste("                      Number of retined pseudocells:",nrow(prop_mat),"out of",nrow(centroidPCAdata),"\n"))
  
  if(F){
    prop_mat2=prop_mat %*% adj_t[colnames(prop_mat),]
    prop_mat <- prop_mat2
  }
  
  ##################
  
  connected_cells=colnames(prop_mat)
  
  
  cat("                      Handling singletons by fast method\n")
  if(T){
    if(length(setdiff(row.names(batchPCAdata),connected_cells))>0){
      nn.ranked <- Seurat:::NNHelper(query = batchPCAdata[setdiff(row.names(batchPCAdata),connected_cells),], data = batchPCAdata[connected_cells,], k = 1, 
                                     method = "annoy", n.trees = 50, searchtype = "standard", 
                                     eps = 0, metric = "euclidean", cache.index = F, 
                                     index = NULL)
      nn.ranked = Indices(nn.ranked)
      nn.ranked.arranged=data.frame(sample=connected_cells,map=connected_cells)
      for(icon in 1:ncol(nn.ranked)){
        nn.ranked.arranged=rbind(nn.ranked.arranged,data.frame(sample=setdiff(row.names(batchPCAdata),connected_cells),map=connected_cells[nn.ranked[,icon]],stringsAsFactors = F))
      }
      nn.ranked=nn.ranked.arranged
    } else {
      nn.ranked=data.frame(sample=colnames(prop_mat),map=colnames(prop_mat),stringsAsFactors = F)
    }
    
    nn.ranked=merge(nn.ranked,data.frame(map=colnames(prop_mat),id_map=1:ncol(prop_mat),stringsAsFactors = F),by="map")
    nn.ranked=merge(nn.ranked,data.frame(sample=row.names(adj),id_sample=1:nrow(adj),stringsAsFactors = F),by="sample")
    #nn.ranked=nn.ranked[match(colnames(prop_mat),nn.ranked$sample),]
    adj_all=Matrix::sparseMatrix(i=nn.ranked$id_sample,j=nn.ranked$id_map,x=rep(1,nrow(nn.ranked)))
    row.names(adj_all)=row.names(adj)
    colnames(adj_all)=colnames(prop_mat)
    
    adj_all=.extra_matrix_rowNorm(adj_all)#Matrix::Diagonal(x=1/rowSums(adj_all)) %*% adj_all
    rm(nn.ranked)
    
  } else {
    snn.matrix <-Seurat:::ComputeSNN(nn_ranked = nn.ranked.1$nn.idx, prune = 1/15)
    rownames(x = snn.matrix) <- rownames(batchPCAdata)
    colnames(x = snn.matrix) <- rownames(batchPCAdata)
    tmp.snn=snn.matrix[,connected_cells]
    tmp.snn=.extra_matrix_rowNorm(input_mat = tmp.snn,rowValues = 1/as.numeric(qlcMatrix::rowMax(tmp.snn)))#Matrix::Diagonal(x=1/as.numeric(qlcMatrix::rowMax(tmp.snn))) %*% tmp.snn
    adj_all=Matrix::drop0(tmp.snn,tol=0.9)
    adj_all=.extra_matrix_rowNorm(rowSums)#Matrix::Diagonal(x=1/rowSums(adj_all)) %*% adj_all
  }
  
  
  
  ###############
  #New
  #Excluded the col normalization step
  #prop_mat <- prop_mat %*% Matrix::Diagonal(x = 1 / (colSums(prop_mat)+0.000000000001))
  adj_all=adj_all %*% t(prop_mat)
  #adj_all = Matrix::Diagonal(x = 1 / (rowSums(adj_all)+0.000000000001)) %*% adj_all
  prop_mat=t(adj_all)
  
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  #prop_mat =  Matrix::Diagonal(x = 1 / (qlcMatrix::rowMax(prop_mat)+0.0001)) %*% prop_mat
  colMax_vals=c()
  for(i in seq(1,ncol(prop_mat),5000)){
    tmp_max=as.numeric(qlcMatrix::colMax(prop_mat[,i:min(i+4999,ncol(prop_mat))]))
    colMax_vals=c(colMax_vals,as.numeric(tmp_max))
  }
  prop_mat =  .extra_matrix_colNorm(input_mat = prop_mat,colValues = 1/colMax_vals)#prop_mat %*% Matrix::Diagonal(x = 1 / colMax_vals)
  #prop_mat=Matrix::drop0(prop_mat,0.1)
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / rowSums(prop_mat)) %*% prop_mat
  
  matWeights=.myEffSizePropMat(prop_mat)
  
  matEffectiveSize=matWeights$effective_sample_size
  matWeights=matWeights$centroid_weights
  matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  
  
  for(i in 1:length(dataArranged)){
    tmp_prop_mat=prop_mat[,match(colnames(dataArranged[[i]]$countData),colnames(prop_mat))]
    #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
    tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
    tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
    tmp_weights[tmp_effsize<4]=0
    tmp=list(prop_mat=tmp_prop_mat,data=dataArranged[[i]],pseudocell_sim_mat=c_c_aff,matWeights=tmp_weights,matEffectiveSize=tmp_effsize)
    dataArranged[[i]]=tmp
  }
  
  return(list(dataArranged=dataArranged,prop_mat=prop_mat))
}

.myConcensusDEFn_step2_detail_newprop3_final_v2m=function(dataArranged,centroidPCAdata,argList,exCentroids=NULL,n.neighbors=4,...){
  
  load(.myFilePathMakerFn("harmony-embeddings",argList=.ArgList,pseudoImportant = F))
  batchPCAdata=harmony_embeddings
  myPseudoAffinityMakerFn=function(harmony_embeddings){
    #nn.ranked.1 <- RANN::nn2(harmony_embeddings, k = 10, eps = 0)
    
    graph=Seurat::FindNeighbors(harmony_embeddings,compute.SNN=T)
    graph=as(graph[["snn"]], "dgCMatrix")
    
    #j <- as.numeric(t(nn.ranked.1$nn.idx))
    #i <- ((1:length(j)) - 1)%/%ncol(nn.ranked.1$nn.idx) + 1
    #k=1#as.numeric(t(nn.ranked.1$nn.dists))
    #graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(harmony_embeddings), nrow(harmony_embeddings)))
    
    #graph=graph+t(graph)
    #graph@x=rep(1,length(graph@x))
    
    if(F){
      graph=graph %*% t(graph)
      diag(graph)=0
      
      
      listCols_sparse<-function(X){
        #converts a sparse Matrix into a list of its columns
        #each list item contains only the nonzero elements of the column
        X<-as(X,"CsparseMatrix")
        res<-split(X@x, findInterval(seq_len(nnzero(X)), X@p, left.open=TRUE))
        names(res)<-colnames(X)
        res
      }
      
      colapply_sparse_nonzero<-function(X,FUN,...,mc.cores=1){
        #apply a function FUN to NONZERO elements of each column of sparse Matrix X
        #for an alternative that operates on all values, see colapply_full
        #mc: should parallel processing be used? Only recommended if FUN is slow
        #... additional args passed to mclapply or to FUN
        #this function always returns a list of length ncol(X)
        if(mc.cores>1){
          res=mclapply(listCols_sparse(X),FUN,...,mc.cores=mc.cores)
        } else {
          res=lapply(listCols_sparse(X),FUN,...)
        }
        res=unlist(res)
        X@x=res
        return(X)
      }
      
      affinities=colapply_sparse_nonzero(X=t(graph),FUN=function(x) exp((-3)*((max(x)-x)/(max(x)+1))^2),mc.cores=argList$ncores)
      affinities=sqrt(t(affinities)*affinities)
      diag(affinities)=0
    } else {
      affinities=graph
    }
    
    return(affinities)
  }
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  #batchPCAdata=harmony_embeddings
  batchPCAdata=as.matrix(batchPCAdata)[,1:argList$nPCs]
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  if(T){
    #newPropMethod
    init_affinities=myPseudoAffinityMakerFn(batchPCAdata)
    colnames(init_affinities)=row.names(init_affinities)=row.names(batchPCAdata)
    load(.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
    head(sl_pseudo)
    #init_affinities=init_affinities %*% init_affinities
    diag(init_affinities)=0
    tmp=qlcMatrix::rowMax(init_affinities)
    init_affinities@x[init_affinities@x==0]=tmp
    init_affinities=init_affinities[,sl_pseudo$pseudocell]
    
    colnames(init_affinities)=sl_pseudo$cluster
    init_affinities=init_affinities[,match(row.names(centroidPCAdata),colnames(init_affinities))]
    
    centroid.affinities=t(init_affinities)
    
  } else {
    nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
    
    adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
    rownames(adj) <- row.names(centroidPCAdata)
    colnames(adj)=c(row.names(batchPCAdata))
    for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
      if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0)>0){
        stop("error!")
      }
      adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = nn.ranked.1$nn.dists[,i]
    }
    
    if(F){
      #old way: it doesn't work when pseuodcells are too sparse!
      nn.centroid.purityIndex=apply(adj,1,function(x){
        
        if(sum(x>0,na.rm = T)>0){
          x=x[which(x>0)]
          x=median(x)
        } else {
          x=0
        }
        
        return(x)
      })
    } else {
      tmp_df=data.frame(pseudocells=nn.ranked.1$nn.idx[,1],dist=nn.ranked.1$nn.dists[,1],stringsAsFactors = F)
      tmp_df=aggregate(dist~pseudocells,data=tmp_df,function(x) median(x,na.rm = T))
      nn.centroid.purityIndex=tmp_df$dist
      names(nn.centroid.purityIndex)=as.character(tmp_df$pseudocells)
      
      
      nn.centroid.purityIndex=nn.centroid.purityIndex[match(row.names(centroidPCAdata),names(nn.centroid.purityIndex))]
      if(sum(is.na(nn.centroid.purityIndex))>0){
        nn.centroid.purityIndex[is.na(nn.centroid.purityIndex)]=0
      }
      tmp_df=quantile(nn.centroid.purityIndex,0.95)
      nn.centroid.purityIndex[nn.centroid.purityIndex>tmp_df]=tmp_df
      rm(tmp_df)
    }
    
    
    
    affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
    affinities=do.call("rbind",affinities)
    affinities[which(adj==0)]=0
    centroid.affinities=affinities
    
  }
  
  
  #######################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_adj=adj
    tst_adj[tst_adj>0]=1
    
    tst_adj=tst_adj %*% as.matrix(tst)
    
    .tmp_step1=tst_adj
    
    dataArranged2
  }
  
  #######################
  
  
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  
  
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  prop_mat=centroid.affinities
  
  if(F){
    dataArranged2=dataArranged
    for(i in 1:length(dataArranged2)){
      tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged2[[i]]$pcaData),colnames(prop_mat))],data=dataArranged2[[i]],pseudocell_sim_mat=prop_mat)
      dataArranged2[[i]]=tmp
    }
    p=.sconline.anno2pseudocell_tmp(res_prop=dataArranged,argList=argList,annoCol="anno_orig_cellState",collapse_datasets=T,return_plot=T)
    ggsave(plot=p$plot,file="~/myBucket/torm3.pdf",width=12,height = 12)
    
  }
  
  
  #################################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_affinities=prop_mat %*% as.matrix(tst)
    
    #apply(tst_affinities,2,mean)
    apply(tst_affinities,2,function(x) sum(x>0.5))
    
    
    tst_pd=pd[match(colnames(prop_mat),row.names(pd)),]
    aggregate(prop_mat["8",]~tst_pd$anno_orig_cellState,sum,data=NULL)
    
    eff_size=.myEffSizePropMat(prop_mat)
    
    x=apply(tst_affinities,1,function(x) which(x==max(x))[1])
    x=data.frame(size=eff_size$effective_sample_size,g=x)
    aggregate(size~g,data=x,summary)
    
  }
  
  
  
  ###############################
  
  
  c_c_aff=prop_mat %*% t(prop_mat)
  c_c_aff=c_c_aff/(diag(c_c_aff)+0.0000001)
  diag(c_c_aff)=1
  c_c_aff=pmin(c_c_aff,1)
  c_c_aff=sqrt(c_c_aff*t(c_c_aff))
  #c_c_aff[which(c_c_aff<0.1)]=0
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff
  
  prop_mat=c_c_aff %*% prop_mat
  
  prop_mat2=t(apply(prop_mat,1,function(x) x/max(x)))
  prop_mat[which(prop_mat2<0.1)]=0
  
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  for(i in 1:length(dataArranged)){
    tmp=list(prop_mat=prop_mat[,match(row.names(dataArranged[[i]]$pcaData),colnames(prop_mat))],data=dataArranged[[i]],pseudocell_sim_mat=c_c_aff)
    dataArranged[[i]]=tmp
  }
  
  return(dataArranged)
}

.myConcensusDEFn_step2_detail_newprop3_final_v2_limited=function(data,centroidPCAdata,argList,exCentroids=NULL,runIndx=1,n.neighbors=4,...){
  
  
  if(F){
    pca_centroid[155,]=(pca_centroid[165,]+pca_centroid[145,])/2
  }
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  
  batchPCAdata=as.matrix(data$pcaData)[,1:argList$nPCs]
  dsName=data$dsName
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  
  nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
    if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0)>0){
      stop("error!")
    }
    adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = nn.ranked.1$nn.dists[,i]
  }
  
  
  
  
  #######################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_adj=adj
    tst_adj[tst_adj>0]=1
    
    tst_adj=tst_adj %*% as.matrix(tst)
    
    .tmp_step1=tst_adj
  }
  
  #######################
  
  if(F){
    #old way: it doesn't work when pseuodcells are too sparse!
    nn.centroid.purityIndex=apply(adj,1,function(x){
      
      if(sum(x>0,na.rm = T)>0){
        x=x[which(x>0)]
        x=median(x)
      } else {
        x=0
      }
      
      return(x)
    })
  } else {
    tmp_df=data.frame(pseudocells=nn.ranked.1$nn.idx[,1],dist=nn.ranked.1$nn.dists[,1],stringsAsFactors = F)
    tmp_df=aggregate(dist~pseudocells,data=tmp_df,function(x) median(x,na.rm = T))
    nn.centroid.purityIndex=tmp_df$dist
    names(nn.centroid.purityIndex)=as.character(tmp_df$pseudocells)
    
    
    nn.centroid.purityIndex=nn.centroid.purityIndex[match(row.names(centroidPCAdata),names(nn.centroid.purityIndex))]
    if(sum(is.na(nn.centroid.purityIndex))>0){
      nn.centroid.purityIndex[is.na(nn.centroid.purityIndex)]=0
    }
    tmp_df=quantile(nn.centroid.purityIndex,0.95)
    nn.centroid.purityIndex[nn.centroid.purityIndex>tmp_df]=tmp_df
    rm(tmp_df)
  }
  
  
  
  affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinities[which(adj==0)]=0
  centroid.affinities=affinities
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  prop_mat=centroid.affinities
  
  
  #################################
  #Examination of the propagation performance
  if(F){
    tmp_pd=pd[match(colnames(data$logNormData),row.names(pd)),]
    tst=.myOneHotFn(tmp_pd$anno_orig_cellState)
    tst_affinities=prop_mat %*% as.matrix(tst)
    
    #apply(tst_affinities,2,mean)
    apply(tst_affinities,2,function(x) sum(x>0.5))
    
    
    tst_pd=pd[match(colnames(prop_mat),row.names(pd)),]
    aggregate(prop_mat["8",]~tst_pd$anno_orig_cellState,sum,data=NULL)
    
    eff_size=.myEffSizePropMat(prop_mat)
    
    x=apply(tst_affinities,1,function(x) which(x==max(x))[1])
    x=data.frame(size=eff_size$effective_sample_size,g=x)
    aggregate(size~g,data=x,summary)
    
  }
  
  
  
  ###############################
  
  
  c_c_aff=prop_mat %*% t(prop_mat)
  c_c_aff=c_c_aff/diag(c_c_aff)
  diag(c_c_aff)=1
  c_c_aff=pmin(c_c_aff,1)
  c_c_aff=sqrt(c_c_aff*t(c_c_aff))
  c_c_aff[which(c_c_aff<0.1)]=0
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff
  
  prop_mat=c_c_aff %*% prop_mat
  
  prop_mat2=t(apply(prop_mat,1,function(x) x/max(x)))
  prop_mat[which(prop_mat2<0.1)]=0
  
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  
  return(list(data=data,prop_mat=prop_mat,pseudocell_sim_mat=c_c_aff))
}

.myConcensusDEFn_step2_detail_newprop3v2=function(data,centroidPCAdata,argList,n.adaptiveKernel=20,nPropIter=1,exCentroids=NULL,runIndx=1,n.neighbors=4){
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  batchPCAdata=as.matrix(data$pcaData)[,1:argList$nPCs]
  dsName=data$dsName
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  if(nrow(batchPCAdata)<(3*n.adaptiveKernel+1)){
    n.adaptiveKernel=2
  }
  
  
  nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
    if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0)>0){
      stop("error!")
    }
    adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = nn.ranked.1$nn.dists[,i]
  }
  
  
  neighbor_count=max(nrow(centroidPCAdata)/100,2)
  nn.centroid.purityIndex=RANN::nn2(data = centroidPCAdata, k = neighbor_count, eps = 0)
  nn.centroid.purityIndex=nn.centroid.purityIndex$nn.dists[,neighbor_count]*1.5
  
  
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  nn.dataset <- RANN::nn2(batchPCAdata, k = min(2000,nrow(batchPCAdata)), eps = 0)
  nn.dataset.purityIndx=RANN::nn2(data = centroidPCAdata,query = batchPCAdata, k = 1, eps = 0)
  nn.dataset.purityIndx=as.numeric(nn.dataset.purityIndx$nn.dist[,1])
  
  
  
  affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinities[which(adj==0)]=0
  centroid.affinities=affinities
  
  dists=nn.dataset$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/(nn.dataset.purityIndx[x]+0.01))^2))
  affinities=do.call("rbind",affinities)
  if(ncol(affinities)>(3*n.adaptiveKernel+1)){
    affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
    affinities2=affinities-affinitiesThr
    affinities[affinities2<=0]=0
  }
  nn.dataset$affinities=affinities
  rm(affinities)
  
  j <- as.numeric(t(nn.dataset$nn.idx))
  i <- ((1:length(j)) - 1)%/%ncol(nn.dataset$affinities) + 1
  k=as.numeric(t(nn.dataset$affinities))
  
  graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(batchPCAdata), nrow(batchPCAdata)))
  rm(i,j,k)
  graph <- .extra_matrix_rowNorm(graph)#Matrix::Diagonal(x = 1 / (rowSums(graph))) %*% graph
  
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  if(class(graph)!="dgCMatrix"){
    graph=as(graph,"dgCMatrix")
  }
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(graph) <- row.names(batchPCAdata)
  colnames(graph) <- row.names(batchPCAdata)
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  checkProp=T
  if(nPropIter<1){
    nPropIter=1
  }
  
  #nPropIter=4
  {
    for(i in nPropIter:1){
      if(i==nPropIter){
        prop_mat=centroid.affinities
      } else {
        prop_mat= prop_mat %*% graph
        
      }
    }
    
    
  }
  
  prop_mat=prop_mat*adj
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  c_c_aff=prop_mat %*% t(prop_mat)
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff
  
  prop_mat=c_c_aff %*% prop_mat
  
  centroid_weights=1/rowSums(adj)
  
  centroid_weightsN=1/(1+exp((12-2/(centroid_weights))/2))
  prop_mat2=sweep(prop_mat,1,centroid_weightsN,"*")
  
  final_check=rowSums(prop_mat2)
  if(sum(final_check>1.001)>0){
    stop(paste0("Error in the construction of the propagation matrix for ",dsName," dataset!"))
  }
  
  return(list(data=data,prop_mat=prop_mat2))
}

.myConcensusDEFn_step2_detail_newprop3_cosine=function(data,centroidPCAdata,argList,n.adaptiveKernel=20,nPropIter=1,exCentroids=NULL,runIndx=1,n.neighbors=4){
  
  mymakeCosine=function(inputMatrix){
    inputMatrix_cosine=inputMatrix*inputMatrix
    inputMatrix_cosine=rowSums(inputMatrix_cosine)
    inputMatrix=sweep(inputMatrix,1,inputMatrix_cosine,"/")
    return(inputMatrix)
  }
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  batchPCAdata=as.matrix(data$pcaData)[,1:argList$nPCs]
  batchPCAdata=mymakeCosine(inputMatrix=batchPCAdata)
  
  dsName=data$dsName
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  centroidPCAdata=mymakeCosine(inputMatrix=centroidPCAdata)
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  if(nrow(batchPCAdata)<(3*n.adaptiveKernel+1)){
    n.adaptiveKernel=2
  }
  
  
  nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
    if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0)>0){
      stop("error!")
    }
    adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = nn.ranked.1$nn.dists[,i]
  }
  
  
  nn.centroid.purityIndex=RANN::nn2(data = centroidPCAdata, k = 2, eps = 0)
  nn.centroid.purityIndex=nn.centroid.purityIndex$nn.dists[,2]
  
  
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  nn.dataset <- RANN::nn2(batchPCAdata, k = min(2000,nrow(batchPCAdata)), eps = 0)
  nn.dataset.purityIndx=RANN::nn2(data = centroidPCAdata,query = batchPCAdata, k = 1, eps = 0)
  nn.dataset.purityIndx=as.numeric(nn.dataset.purityIndx$nn.dist[,1])
  
  
  
  affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinities[which(adj==0)]=0
  centroid.affinities=affinities
  
  dists=nn.dataset$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/(nn.dataset.purityIndx[x]+0.01))^2))
  affinities=do.call("rbind",affinities)
  if(ncol(affinities)>(3*n.adaptiveKernel+1)){
    affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
    affinities2=affinities-affinitiesThr
    affinities[affinities2<=0]=0
  }
  nn.dataset$affinities=affinities
  rm(affinities)
  
  j <- as.numeric(t(nn.dataset$nn.idx))
  i <- ((1:length(j)) - 1)%/%ncol(nn.dataset$affinities) + 1
  k=as.numeric(t(nn.dataset$affinities))
  
  graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(batchPCAdata), nrow(batchPCAdata)))
  rm(i,j,k)
  graph <- .extra_matrix_rowNorm(graph)#Matrix::Diagonal(x = 1 / (rowSums(graph))) %*% graph
  
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  if(class(graph)!="dgCMatrix"){
    graph=as(graph,"dgCMatrix")
  }
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(graph) <- row.names(batchPCAdata)
  colnames(graph) <- row.names(batchPCAdata)
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  checkProp=T
  if(nPropIter<1){
    nPropIter=1
  }
  
  #nPropIter=4
  {
    for(i in nPropIter:1){
      if(i==nPropIter){
        prop_mat=centroid.affinities
      } else {
        prop_mat= prop_mat %*% graph
        
      }
    }
    
    
  }
  
  prop_mat=prop_mat*adj
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  c_c_aff=prop_mat %*% t(prop_mat)
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff
  
  prop_mat=c_c_aff %*% prop_mat
  
  centroid_weights=1/rowSums(adj)
  
  centroid_weightsN=1/(1+exp((12-2/(centroid_weights))/2))
  prop_mat2=sweep(prop_mat,1,centroid_weightsN,"*")
  
  final_check=rowSums(prop_mat2)
  if(sum(final_check>1.001)>0){
    stop(paste0("Error in the construction of the propagation matrix for ",dsName," dataset!"))
  }
  
  return(list(data=data,prop_mat=prop_mat2))
}

.myConcensusDEFn_step2_detail_newprop3_scale_v2=function(data,centroidPCAdata,argList,n.adaptiveKernel=20,nPropIter=1,exCentroids=NULL,runIndx=1,n.neighbors=4){
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  batchPCAdata=as.matrix(data$pcaData)[,1:argList$nPCs]
  dsName=data$dsName
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  if(nrow(batchPCAdata)<(3*n.adaptiveKernel+1)){
    n.adaptiveKernel=2
  }
  
  
  nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  
  adj <- matrix(100, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
    if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0& adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]<100)>0){
      stop("Error in constructing the adjancy matrices for propagation!")
    }
    adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = nn.ranked.1$nn.dists[,i]
  }
  
  
  nn.centroid.purityIndex=RANN::nn2(data = centroidPCAdata, k = 2, eps = 0)
  nn.centroid.purityIndex=nn.centroid.purityIndex$nn.dists[,2]
  
  
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  nn.dataset <- RANN::nn2(batchPCAdata, k = min(2000,nrow(batchPCAdata)), eps = 0)
  nn.dataset.purityIndx=RANN::nn2(data = centroidPCAdata,query = batchPCAdata, k = 1, eps = 0)
  nn.dataset.purityIndx=as.numeric(nn.dataset.purityIndx$nn.dist[,1])
  
  
  affinities=lapply(1:nrow(adj),function(x) {
    y=rep(0,ncol(adj))
    y[adj[x,]<100]=exp((-1)*(adj[x,adj[x,]<100]/nn.centroid.purityIndex[x])^2)
    return(y)
  })
  
  #affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  ##zero values in adj don't mean the distance is zero!
  #affinities[which(adj==0)]=0
  centroid.affinities=affinities
  
  dists=nn.dataset$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/(nn.dataset.purityIndx[x]+0.01))^2))
  affinities=do.call("rbind",affinities)
  if(ncol(affinities)>(3*n.adaptiveKernel+1)){
    affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
    affinities2=affinities-affinitiesThr
    affinities[affinities2<=0]=0
  }
  nn.dataset$affinities=affinities
  rm(affinities)
  
  j <- as.numeric(t(nn.dataset$nn.idx))
  i <- ((1:length(j)) - 1)%/%ncol(nn.dataset$affinities) + 1
  k=as.numeric(t(nn.dataset$affinities))
  
  graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(batchPCAdata), nrow(batchPCAdata)))
  rm(i,j,k)
  graph <- .extra_matrix_rowNorm(graph)#Matrix::Diagonal(x = 1 / (rowSums(graph))) %*% graph
  
  #centroid.affinities <-  centroid.affinities %*% Matrix::Diagonal(x = 1 / (colSums(centroid.affinities)+0.000000000001))
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  if(class(graph)!="dgCMatrix"){
    graph=as(graph,"dgCMatrix")
  }
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(graph) <- row.names(batchPCAdata)
  colnames(graph) <- row.names(batchPCAdata)
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  checkProp=T
  if(nPropIter<1){
    nPropIter=1
  }
  
  #nPropIter=4
  {
    for(i in nPropIter:1){
      if(i==nPropIter){
        prop_mat1=prop_mat=centroid.affinities %*% graph
      } else {
        prop_mat= prop_mat %*% graph
        
      }
    }
    
    if(nPropIter>1){
      
      tmp_graph=prop_mat1
      tmp_graph[which(tmp_graph>0)]=1
      prop_mat=prop_mat*tmp_graph
      prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
      
    }
    
  }
  
  #prop_mat=prop_mat*adj
  #prop_mat <- Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  c_c_aff=centroid.affinities %*% t(centroid.affinities)
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff
  
  prop_mat=c_c_aff %*% prop_mat
  prop_mat_scale=t(apply(prop_mat,1,scale))
  prop_mat[which(prop_mat_scale<0.5)]=0
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  #prop_mat2=sweep(prop_mat,1,centroid_weightsN,"*")
  
  #final_check=rowSums(prop_mat2)
  #if(sum(final_check>1.001)>0){
  #  stop(paste0("Error in the construction of the propagation matrix for ",dsName," dataset!"))
  #}
  
  return(list(data=data,prop_mat=prop_mat))
}

.myEffSizePropMat=function(prop_mat){
  require(Matrix)
  
  
  prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / rowSums(prop_mat)+0.0000001) %*% prop_mat
  
  centroid_weights=(rowSums(prop_mat))^2/(rowSums(prop_mat^2)+0.00000001)
  centroid_weightsN=1/(1+exp((6-1/(1/(centroid_weights+0.000000001)))))
  return(list(effective_sample_size=centroid_weights,centroid_weights=centroid_weightsN))
}

.myConcensusDEFn_step2_detail_newprop3_scale=function(data,centroidPCAdata,argList,n.adaptiveKernel=20,nPropIter=1,exCentroids=NULL,runIndx=1,n.neighbors=4){
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=4
  batchPCAdata=as.matrix(data$pcaData)[,1:argList$nPCs]
  dsName=data$dsName
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  if(nrow(batchPCAdata)<(3*n.adaptiveKernel+1)){
    n.adaptiveKernel=2
  }
  
  
  nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  
  adj <- matrix(100, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
    if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0& adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]<100)>0){
      stop("Error in constructing the adjancy matrices for propagation!")
    }
    adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = nn.ranked.1$nn.dists[,i]
  }
  
  
  nn.centroid.purityIndex=RANN::nn2(data = centroidPCAdata, k = 2, eps = 0)
  nn.centroid.purityIndex=nn.centroid.purityIndex$nn.dists[,2]
  
  
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  nn.dataset <- RANN::nn2(batchPCAdata, k = min(2000,nrow(batchPCAdata)), eps = 0)
  nn.dataset.purityIndx=RANN::nn2(data = centroidPCAdata,query = batchPCAdata, k = 1, eps = 0)
  nn.dataset.purityIndx=as.numeric(nn.dataset.purityIndx$nn.dist[,1])
  
  
  affinities=lapply(1:nrow(adj),function(x) {
    y=rep(0,ncol(adj))
    y[adj[x,]<100]=exp((-1)*(adj[x,adj[x,]<100]/nn.centroid.purityIndex[x])^2)
    return(y)
  })
  
  #affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  ##zero values in adj don't mean the distance is zero!
  #affinities[which(adj==0)]=0
  centroid.affinities=affinities
  
  dists=nn.dataset$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/(nn.dataset.purityIndx[x]+0.01))^2))
  affinities=do.call("rbind",affinities)
  if(ncol(affinities)>(3*n.adaptiveKernel+1)){
    affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
    affinities2=affinities-affinitiesThr
    affinities[affinities2<=0]=0
  }
  nn.dataset$affinities=affinities
  rm(affinities)
  
  j <- as.numeric(t(nn.dataset$nn.idx))
  i <- ((1:length(j)) - 1)%/%ncol(nn.dataset$affinities) + 1
  k=as.numeric(t(nn.dataset$affinities))
  
  graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(batchPCAdata), nrow(batchPCAdata)))
  rm(i,j,k)
  graph <- .extra_matrix_rowNorm(graph)#Matrix::Diagonal(x = 1 / (rowSums(graph))) %*% graph
  
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  if(class(graph)!="dgCMatrix"){
    graph=as(graph,"dgCMatrix")
  }
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(graph) <- row.names(batchPCAdata)
  colnames(graph) <- row.names(batchPCAdata)
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  checkProp=T
  if(nPropIter<1){
    nPropIter=1
  }
  
  #nPropIter=4
  {
    for(i in nPropIter:1){
      if(i==nPropIter){
        prop_mat=centroid.affinities
      } else {
        prop_mat= prop_mat %*% graph
        
      }
    }
    
    if(nPropIter>1){
      
      tmp_graph=centroid.affinities
      tmp_graph[which(tmp_graph>0)]=1
      prop_mat=prop_mat*tmp_graph
    }
    
  }
  
  prop_mat=prop_mat*adj
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  c_c_aff=prop_mat %*% t(prop_mat)
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff %*% c_c_aff
  
  prop_mat=c_c_aff %*% prop_mat
  prop_mat_scale=t(apply(prop_mat,1,scale))
  prop_mat[which(prop_mat_scale<0.5)]=0
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  
  centroid_weights= Matrix::Diagonal(x = 1 / (rowMaxs(prop_mat)+0.000000000001)) %*% prop_mat
  centroid_weights=1/rowSums(centroid_weights)
  
  centroid_weights=1/rowSums(adj)
  
  centroid_weightsN=1/(1+exp((12-2/(centroid_weights))/2))
  prop_mat2=sweep(prop_mat,1,centroid_weightsN,"*")
  
  final_check=rowSums(prop_mat2)
  if(sum(final_check>1.001)>0){
    stop(paste0("Error in the construction of the propagation matrix for ",dsName," dataset!"))
  }
  
  return(list(data=data,prop_mat=prop_mat2))
}

.myConcensusDEFn_step2_detail_newprop4=function(data,centroidPCAdata,argList,n.adaptiveKernel=20,nPropIter=3,exCentroids=NULL,runIndx=1,n.neighbors=2){
  
  #data=dataArranged[[1]];centroidPCAdata=pca_centroid;argList=.ArgList;n.adaptiveKernel=20;nPropIter=1;exCentroids=NULL;n.neighbors=2
  batchPCAdata=as.matrix(data$pcaData)[,1:argList$nPCs]
  dsName=data$dsName
  centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),]
  }
  
  if(nrow(batchPCAdata)<(3*n.adaptiveKernel+1)){
    n.adaptiveKernel=2
  }
  
  
  nn.ranked.1 <- RANN::nn2(query = batchPCAdata,data = centroidPCAdata, k = n.neighbors, eps = 0)
  
  adj <- matrix(0, nrow=nrow(centroidPCAdata), ncol=(nrow(batchPCAdata)))
  rownames(adj) <- row.names(centroidPCAdata)
  colnames(adj)=c(row.names(batchPCAdata))
  for(i in seq_len(ncol(nn.ranked.1$nn.idx))) {
    if(sum(adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))]>0)>0){
      stop("error!")
    }
    adj[cbind(nn.ranked.1$nn.idx[,i],1:nrow(nn.ranked.1$nn.idx))] = nn.ranked.1$nn.dists[,i]
  }
  
  
  nn.centroid.purityIndex=RANN::nn2(data = centroidPCAdata, k = 2, eps = 0)
  nn.centroid.purityIndex=nn.centroid.purityIndex$nn.dists[,2]
  
  
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  nn.dataset <- RANN::nn2(batchPCAdata, k = min(2000,nrow(batchPCAdata)), eps = 0)
  nn.dataset.purityIndx=RANN::nn2(data = centroidPCAdata,query = batchPCAdata, k = 1, eps = 0)
  nn.dataset.purityIndx=as.numeric(nn.dataset.purityIndx$nn.dist[,1])
  
  
  
  affinities=lapply(1:nrow(adj),function(x) exp((-1)*(adj[x,]/nn.centroid.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinities[which(adj==0)]=0
  centroid.affinities=affinities
  
  dists=nn.dataset$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/(nn.dataset.purityIndx[x]+0.01))^2))
  affinities=do.call("rbind",affinities)
  if(ncol(affinities)>(3*n.adaptiveKernel+1)){
    affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
    affinities2=affinities-affinitiesThr
    affinities[affinities2<=0]=0
  }
  nn.dataset$affinities=affinities
  rm(affinities)
  
  j <- as.numeric(t(nn.dataset$nn.idx))
  i <- ((1:length(j)) - 1)%/%ncol(nn.dataset$affinities) + 1
  k=as.numeric(t(nn.dataset$affinities))
  
  graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(batchPCAdata), nrow(batchPCAdata)))
  rm(i,j,k)
  graph <- .extra_matrix_rowNorm(graph)#Matrix::Diagonal(x = 1 / (rowSums(graph))) %*% graph
  
  centroid.affinities <- .extra_matrix_rowNorm(centroid.affinities)#Matrix::Diagonal(x = 1 / (rowSums(centroid.affinities)+0.000000000001)) %*% centroid.affinities
  
  if(class(graph)!="dgCMatrix"){
    graph=as(graph,"dgCMatrix")
  }
  
  if(class(centroid.affinities)!="dgCMatrix"){
    centroid.affinities=as(centroid.affinities,"dgCMatrix")
  }
  
  rownames(graph) <- row.names(batchPCAdata)
  colnames(graph) <- row.names(batchPCAdata)
  
  rownames(centroid.affinities) <- row.names(centroidPCAdata)
  colnames(centroid.affinities) <- row.names(batchPCAdata)
  
  checkProp=T
  if(nPropIter<1){
    nPropIter=1
  }
  
  #nPropIter=4
  {
    for(i in nPropIter:1){
      if(i==nPropIter){
        prop_mat=centroid.affinities
      } else {
        prop_mat= prop_mat %*% graph
        
      }
    }
    
    
  }
  
  prop_mat=prop_mat*adj
  prop_mat <- .extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  c_c_aff=prop_mat %*% t(prop_mat)
  c_c_aff <- .extra_matrix_rowNorm(c_c_aff)#Matrix::Diagonal(x = 1 / (rowSums(c_c_aff)+0.000000000001)) %*% c_c_aff
  c_c_aff=c_c_aff %*% c_c_aff
  
  prop_mat=c_c_aff %*% prop_mat
  
  centroid_weights=1/rowSums(adj)
  
  centroid_weightsN=1/(1+exp((12-2/(centroid_weights))/2))
  prop_mat2=sweep(prop_mat,1,centroid_weightsN,"*")
  
  final_check=rowSums(prop_mat2)
  if(sum(final_check>1.001)>0){
    stop(paste0("Error in the construction of the propagation matrix for ",dsName," dataset!"))
  }
  
  return(list(data=data,prop_mat=prop_mat2))
}

.myConcensusDEFn_step2_detail_exp_org=function(inputData,argList,centroidPCAdata,exCentroids=NULL,min_cells=5,resGeneMeanSd=NULL,nperm=15){
  
  #tst=unlist(lapply(res,function(x) ncol(x$data$countData)))
  
  require(data.table,quietly = T)
  normExp=NULL
  
  batchPCAdata=as.matrix(inputData$data$pcaData)[,1:argList$nPCs]
  #expData=inputData$data$expData
  dsName=inputData$data$dsName
  print(dsName)
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),1:argList$nPCs]
  } else {
    centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  }
  
  prop_mat=inputData$prop_mat
  uniformZscore=argList$uniformZscore
  
  #centroid_weights=prop_mat[cbind((ncol(inputData$data$logNormData)+1):nrow(prop_mat),(ncol(inputData$data$logNormData)+1):nrow(prop_mat))]
  #centroid_weights=1-centroid_weights
  #names(centroid_weights)=row.names(prop_mat)[(ncol(inputData$data$logNormData)+1):nrow(prop_mat)]
  
  res2=""
  if(argList$dist_zscore_norm){
    if(uniformZscore){
      load(.myFilePathMakerFn("geneExp_meanSd_weighted",argList=argList))
        resGeneMeanSd=resGeneMeanSd[which(resGeneMeanSd$anno_batch==dsName),]
        resGeneMeanSd=resGeneMeanSd[order(resGeneMeanSd$meanExp_lognorm,decreasing=T),]
        #resGeneMeanSd$w_sd_lognorm[1:150]=pmax(resGeneMeanSd$w_sd_lognorm[1:150],mean(resGeneMeanSd$w_sd_lognorm[101:150]))
        resGeneMeanSd=resGeneMeanSd[match(row.names(inputData$data$logNormData),resGeneMeanSd$gene),]
        if(sum(is.na(resGeneMeanSd$w_sd_lognorm))>0|nrow(resGeneMeanSd)==0){
          stop("Error in matching the centroid gene mean/sd!")
        }
        
        if(argList$regularize){
          reg_sd_data=.myParRegFn(object = resGeneMeanSd,model_pars = c("w_sd_lognorm"),genes_count_mean = "w_mean_lognorm",genes = "gene",doGeneLogTransform=F)
          
          if(!all(resGeneMeanSd$gene==row.names(reg_sd_data))){
            stop("Error in matching the gene names in regularization step")
          }
          resGeneMeanSd$w_sd_lognorm=reg_sd_data[,'w_sd_lognorm']
        }
        
        
        dataN=sweep(inputData$data$logNormData,1,resGeneMeanSd$w_mean_lognorm,"-")
        dataN=sweep(dataN,1,resGeneMeanSd$w_sd_lognorm,"/")
        dataN[which(resGeneMeanSd$w_mean_lognorm==0),]=0
        dataN[which(as.matrix(dataN)>10)]=(10)
        expData=dataN
      
    } else {
      #load(.myFilePathMakerFn("centroids_geneExp_meanSd",argList=argList))
      if(is.null(resGeneMeanSd)){
        load(.myFilePathMakerFn("geneExp_meanSd",argList=argList))
      }
      
      resGeneMeanSd=resGeneMeanSd[which(resGeneMeanSd$anno_batch==dsName),]
      #resGeneMeanSd$sdExp_lognorm[1:150]=pmax(resGeneMeanSd$sdExp_lognorm[1:150],mean(resGeneMeanSd$sdExp_lognorm[101:150]))
      resGeneMeanSd=resGeneMeanSd[match(row.names(inputData$data$logNormData),resGeneMeanSd$gene),]
      if(sum(is.na(resGeneMeanSd$sdExp_lognorm))>0|nrow(resGeneMeanSd)==0){
        stop("Erro in matching the centroid gene mean/sd!")
      }
      
      if(argList$regularize){
        reg_sd_data=.myParRegFn(object = resGeneMeanSd,model_pars = c("sdExp_lognorm"),genes_count_mean = "meanExp_lognorm",genes = "gene",doGeneLogTransform=F)
        
        if(!all(resGeneMeanSd$gene==row.names(reg_sd_data))){
          stop("Error in matching the gene names in regularization step")
        }
        resGeneMeanSd$sdExp_lognorm=pmax(reg_sd_data[,'sdExp_lognorm'],resGeneMeanSd$sdExp_lognorm)
      }
      dataN=sweep(inputData$data$logNormData,1,resGeneMeanSd$meanExp_lognorm,"-")
      dataN=sweep(dataN,1,resGeneMeanSd$sdExp_lognorm,"/")
      dataN[which(resGeneMeanSd$meanExp_lognorm==0),]=0
      tmpInd=which(as.matrix(dataN)>(10))
      if(length(tmpInd)>0){
        dataN[tmpInd]=(10)
      }
      
      expData=dataN
    }
    
    #tmpInd=which(as.matrix(expData)<(-4))
    #if(length(tmpInd)>0){
    #  expData[tmpInd]=(-4)
    #}
    
    tmpInd=apply(expData,1,function(x) sum(abs(x)>2)/length(x))
    tmpInd=which(tmpInd>=0.5)
    if(length(tmpInd)>0){
      expData[tmpInd,]=0
    }
    tmp2=as.matrix(t(expData))
    #tmp2c=matrix(0,ncol=ncol(tmp2),nrow=nrow(centroidPCAdata))
    #tmp2=rbind(tmp2,tmp2c)
    
    #res2=as.matrix(prop_mat)[-c(1:ncol(expData)),] %*% as.matrix(tmp2)
    normExp=tmp2
    res2=eigenMapMatMult(as.matrix(prop_mat), as.matrix(tmp2))
    
    centroid_weights=.myEffSizePropMat(prop_mat)
    res2=sweep(res2,1,centroid_weights$centroid_weights,"*")
    
    regVals=res2
    if(F){
      
      resGeneMeanSd=resGeneMeanSd[match(colnames(tmp2),resGeneMeanSd$gene),]
      geneNets=RANN::nn2(data = resGeneMeanSd[,c("meanExp_count","sdExp_count","geo_mean_count","meanExp_lognorm","sdExp_lognorm","geo_mean_lognorm")], k = 201, eps = 0)
      
      adj <- matrix(0, nrow=nrow(resGeneMeanSd), ncol=nrow(resGeneMeanSd))
      rownames(adj) <- resGeneMeanSd$gene
      colnames(adj)=resGeneMeanSd$gene
      for(i in seq_len(nrow(resGeneMeanSd))) {
        adj[i,geneNets$nn.idx[i,]] = 1
      }
      adj2=as(adj,"dgCMatrix")
      adj2=as.data.frame(summary(adj2))
      
      tmp=t(res2)
      tmp=tmp[adj2[,2],]
      tmp=as.data.frame(tmp)
      
      tmp=as.data.table(tmp)
      tmp$gene=adj2[,1]
      tmp_mean <- tmp[, lapply(.SD, mean), by=gene]
      tmp_mean=as.data.frame(tmp_mean)
      row.names(tmp_mean)=resGeneMeanSd$gene[tmp_mean$gene]
      tmp_mean=tmp_mean[,-1]
      tmp_mean=as.matrix(tmp_mean)
      tmp_sd <- tmp[, lapply(.SD,sd), by=gene]
      tmp_sd=as.data.frame(tmp_sd)
      row.names(tmp_sd)=resGeneMeanSd$gene[tmp_sd$gene]
      tmp_sd=tmp_sd[,-1]
      tmp_sd=as.matrix(tmp_sd)
      
      tmp_mean=t(tmp_mean)
      tmp_mean=tmp_mean[,match(colnames(tmp2),colnames(tmp_mean))]
      tmp_sd=t(tmp_sd)
      tmp_sd=tmp_sd[,match(colnames(tmp2),colnames(tmp_sd))]
      
      
      
      weight_mat=(matrix(centroid_weights,ncol=ncol(regVals),nrow=nrow(regVals)))
      
      regVals=(res2-tmp_mean)/(tmp_sd+0.2)
      regVals[sign(regVals*res2)<0]=0
      regVals2=regVals*weight_mat
      
    }
    
    if(T){
      
      resGeneMeanSd=resGeneMeanSd[match(colnames(tmp2),resGeneMeanSd$gene),]
      for(ic in 1:nrow(res2)){
        df=data.frame(resGeneMeanSd,zscore=res2[ic,],stringsAsFactors = F)
        reg_zscore_data=.myParRegFn(object = df,model_pars = c("zscore"),genes_count_mean = "meanExp_lognorm",genes = "gene",doGeneLogTransform=F)
        summary(reg_zscore_data$zscore)
      }
      adj <- matrix(0, nrow=nrow(resGeneMeanSd), ncol=nrow(resGeneMeanSd))
      rownames(adj) <- resGeneMeanSd$gene
      colnames(adj)=resGeneMeanSd$gene
      for(i in seq_len(nrow(resGeneMeanSd))) {
        adj[i,geneNets$nn.idx[i,]] = 1
      }
      adj2=as(adj,"dgCMatrix")
      adj2=as.data.frame(summary(adj2))
      
      tmp=t(res2)
      tmp=tmp[adj2[,2],]
      tmp=as.data.frame(tmp)
      
      tmp=as.data.table(tmp)
      tmp$gene=adj2[,1]
      tmp_mean <- tmp[, lapply(.SD, mean), by=gene]
      tmp_mean=as.data.frame(tmp_mean)
      row.names(tmp_mean)=resGeneMeanSd$gene[tmp_mean$gene]
      tmp_mean=tmp_mean[,-1]
      tmp_mean=as.matrix(tmp_mean)
      tmp_sd <- tmp[, lapply(.SD,sd), by=gene]
      tmp_sd=as.data.frame(tmp_sd)
      row.names(tmp_sd)=resGeneMeanSd$gene[tmp_sd$gene]
      tmp_sd=tmp_sd[,-1]
      tmp_sd=as.matrix(tmp_sd)
      
      tmp_mean=t(tmp_mean)
      tmp_mean=tmp_mean[,match(colnames(tmp2),colnames(tmp_mean))]
      tmp_sd=t(tmp_sd)
      tmp_sd=tmp_sd[,match(colnames(tmp2),colnames(tmp_sd))]
      
      
      
      weight_mat=(matrix(centroid_weights,ncol=ncol(regVals),nrow=nrow(regVals)))
      
      regVals=(res2-tmp_mean)/(tmp_sd+0.2)
      regVals[sign(regVals*res2)<0]=0
      regVals2=regVals*weight_mat
      
    }
    
    res2=t(res2)
    regVals=t(regVals)
    row.names(res2)=colnames(tmp2)
    colnames(res2)=row.names(prop_mat)[-c(1:(nrow(tmp2)-nrow(tmp2c)))]
    row.names(regVals)=row.names(res2)
    colnames(regVals)=colnames(res2)
    if(T){
      #centroid_weightsN=centroid_weights[match(colnames(regVals),names(centroid_weights))]
      #centroid_weightsN=t(matrix(centroid_weights,nrow=length(centroid_weightsN),ncol=nrow(regVals)))
      #centroid_weightsN=centroid_weightsN*1.96
      #regVals[which(regVals>centroid_weightsN)]=centroid_weightsN[which(regVals>centroid_weightsN)]
      #regVals[which(regVals<(-1.96)*centroid_weightsN)]=((-1.96)*centroid_weightsN)[which(regVals<(-1.96)*centroid_weightsN)]
      
      res2=regVals#*(argList$sensitiveSearch)+(1-argList$sensitiveSearch)*res2
    }
  } else if (argList$dist_zscore_nbinom){
    
    log_umi=apply(inputData$data$countData,2,function(x) log10(sum(x)))
    res=list()
    
    if(!argList$regularize){
      if(uniformZscore){
        load(.myFilePathMakerFn("kmeans_clustering",argList=argList))
        res_clusters=res_clusters[match(colnames(inputData$data$countData),res_clusters$sample),]
        if(sum(is.na(res_clusters$sample))>0){
          stop(paste("Error in kmeans res of",dsName))
        }
        
        for(i in 1:nrow(inputData$data$countData)){
          
          if(sum(inputData$data$countData[i,]>0)>=min_cells){
            df=data.frame(y=inputData$data$countData[i,],log_umi=log_umi)
            fit <- glm(y~log_umi, data = df, family = poisson,weights = res_clusters$cluster_weight)
            theta <- as.numeric(x = MASS::theta.ml(y = df$y, mu = fit$fitted,weights = res_clusters$cluster_weight))
            
            zscores=as.data.frame(matrix(edgeR::zscoreNBinom(df$y,size=theta,mu=fit$fitted.values),nrow = 1))
            colnames(zscores)=colnames(inputData$data$countData)
            row.names(zscores)=row.names(inputData$data$countData)[i]
            res=c(res,list(zscores))
          }
          
        }
      } else {
        for(i in 1:nrow(inputData$data$countData)){
          
          if(sum(inputData$data$countData[i,]>0)>=min_cells){
            #print(i)
            df=data.frame(y=inputData$data$countData[i,],log_umi=log_umi)
            fit <- glm(y~log_umi, data = df, family = poisson)
            theta <- as.numeric(x = MASS::theta.ml(y = df$y, mu = fit$fitted))
            
            zscores=as.data.frame(matrix(edgeR::zscoreNBinom(df$y,size=theta,mu=fit$fitted.values),nrow = 1))
            colnames(zscores)=colnames(inputData$data$countData)
            row.names(zscores)=row.names(inputData$data$countData)[i]
            res=c(res,list(zscores))
          }
          
        }
      }
      
    } else {
      #geoMean=F
      resGeneMeanSd=""
      coefData=NULL
      if(uniformZscore){
        load(.myFilePathMakerFn("kmeans_clustering",argList=argList))
        res_clusters=res_clusters[match(colnames(inputData$data$countData),res_clusters$sample),]
        if(sum(is.na(res_clusters$sample))>0){
          stop(paste("Error in kmeans resesults of",dsName))
        }
        
        
        for(i in 1:nrow(inputData$data$countData)){
          if(sum(inputData$data$countData[i,]>0)>=min_cells){
            df=data.frame(y=inputData$data$countData[i,],log_umi=log_umi)
            fit <- glm(y~log_umi, data = df, family = poisson,weights = res_clusters$cluster_weight)
            theta <- log10(as.numeric(x = MASS::theta.ml(y = df$y, mu = fit$fitted,weights = res_clusters$cluster_weight)))
            
            coefData=rbind(coefData,data.frame(gene=row.names(inputData$data$countData)[i],theta=theta,b0=fit$coefficients[1],b1=fit$coefficients[2]))
          }
        }
        
        load(.myFilePathMakerFn("geneExp_meanSd_weighted",argList=argList))
        colnames(resGeneMeanSd)[colnames(resGeneMeanSd)=="w_mean_count"]="meanExp_count"
        colnames(resGeneMeanSd)[colnames(resGeneMeanSd)=="w_geomean_count"]="geo_mean_count"
      } else {
        
        for(i in 1:nrow(inputData$data$countData)){
          if(sum(inputData$data$countData[i,]>0)>=min_cells){
            df=data.frame(y=inputData$data$countData[i,],log_umi=log_umi)
            fit <- glm(y~log_umi, data = df, family = poisson)
            theta <- log10(as.numeric(x = MASS::theta.ml(y = df$y, mu = fit$fitted)))
            
            coefData=rbind(coefData,data.frame(gene=row.names(inputData$data$countData)[i],theta=theta,b0=fit$coefficients[1],b1=fit$coefficients[2]))
          }
        }
        
        if(is.null(resGeneMeanSd)){
          load(.myFilePathMakerFn("geneExp_meanSd",argList=argList))
        }
      }
      
      
      
      resGeneMeanSd=resGeneMeanSd[which(resGeneMeanSd$anno_batch==dsName),]
      resGeneMeanSd=resGeneMeanSd[,c("gene","meanExp_count","geo_mean_count")]
      coefData=merge(coefData,resGeneMeanSd,by="gene")
      if(argList$geoMean){
        coefData=.myParRegFn(object=coefData,model_pars=c("theta","b0","b1"), genes_count_mean="geo_mean_count",genes="gene",doGeneLogTransform=T, bw_adjust= 3)
      }else {
        coefData=.myParRegFn(object=coefData,model_pars=c("theta","b0","b1"), genes_count_mean="meanExp_count",genes="gene",doGeneLogTransform=T, bw_adjust= 3)
      }
        
      tmpCountData=inputData$data$countData
      tmpCountData=tmpCountData[match(row.names(coefData),row.names(tmpCountData)),]
      
      for(i in 1:nrow(coefData)){
          
        zscores=as.data.frame(matrix(edgeR::zscoreNBinom(tmpCountData[i,],size=10^coefData$theta[i],mu=exp(log_umi*coefData$b1[i]+coefData$b0[i])),nrow = 1))
          
        colnames(zscores)=colnames(tmpCountData)
        row.names(zscores)=row.names(tmpCountData)[i]
        res=c(res,list(zscores))
          
      }
    }
    
    #inputData$data$countData
    res2=data.table::rbindlist(res)
    row.names(res2)=unlist(lapply(res,row.names))
    res2[is.na(res2)]=0
    res3=matrix(0,nrow=length(setdiff(row.names(inputData$data$countData),row.names(res2))),ncol=ncol(inputData$data$countData))
    row.names(res3)=setdiff(row.names(inputData$data$countData),row.names(res2))
    colnames(res3)=colnames(inputData$data$countData)
    res=rbind(res2,res3)
    res=as.data.frame(res)
    row.names(res)=c(row.names(res2),row.names(res3))
    
    res=res[match(row.names(inputData$data$countData),row.names(res)),]
    
    res=as.matrix(res)
    tmpInd=which(res<(-10))
    if(length(tmpInd)>0){
      res[tmpInd]=(-10)
    }
    
    tmpInd=which(res>10)
    if(length(tmpInd)>0){
      res[tmpInd]=10
    }
    
    tmpInd=apply(res,1,function(x) sum(abs(x)>2)/length(x))
    tmpInd=which(tmpInd>=0.5)
    if(length(tmpInd)>0){
      res[tmpInd,]=0
    }
    
    tmp2=t(as.matrix(res))
    tmp2c=matrix(0,ncol=ncol(tmp2),nrow=nrow(centroidPCAdata))
    tmp2=rbind(tmp2,tmp2c)
    normExp=as.matrix(tmp2)
    #res2=as.matrix(prop_mat) %*% as.matrix(tmp2)
    #res2=res2[-c(1:ncol(res)),]
    res2=eigenMapMatMult(as.matrix(prop_mat)[-c(1:ncol(res)),], as.matrix(tmp2))
    
    res2=t(res2)
    row.names(res2)=colnames(tmp2)
    colnames(res2)=row.names(prop_mat)[-c(1:ncol(res))]
  }
  
  if(F){
    
    myWeightedFn=function(inputExpData,weights){
      #res_sd=matrixStats::rowWeightedSds(as.matrix(inputExpData), w = weights)
      #res_m=apply(as.matrix(inputExpData),1,function(x) weighted.mean(x,w=weights))
      inputExpData=as.matrix(inputExpData)
      res_m=apply(inputExpData,1,function(x) sum(x)/sum(weights))
      res_sd=unlist(lapply(1:nrow(inputExpData),function(x) sqrt(sum(weights*(inputExpData[x,]-res_m[x])^2)/sum(weights))))
      #res_sd=apply(as.matrix(inputExpData), 1,sd)
      #res_m=apply(as.matrix(inputExpData),1,mean)
      
      res=data.frame(gene=row.names(inputExpData),w_mean=res_m,w_sd=res_sd)
      
      #reg_sd_data=.myParRegFn(object = res,model_pars = c("w_sd"),genes_count_mean = "w_mean",genes = "gene",doGeneLogTransform=F)
      
      #if(!all(resGeneMeanSd$gene==row.names(reg_sd_data))){
      #  stop("Error in matching the gene names in regularization step")
      #}
      #res$w_sd=pmax(res$w_sd,reg_sd_data[,'w_sd'])
      
      return(res)
    }
    
    centroid_weightsN=centroid_weights[match(colnames(res2),names(centroid_weights))]
    res_w=myWeightedFn(res2,weights=centroid_weightsN)
    res_w=res_w[match(row.names(res2),res_w$gene),]
    
    dataN=sweep(res2,1,res_w$w_mean,"-")
    dataN=sweep(dataN,1,(res_w$w_sd+0.1),"/")
    
    
    
    
  }
  
  
  return(list(res=res2,w_zscore_centroid=centroid_weights,dsName=dsName,normExp=t(normExp[c(1:nrow(batchPCAdata)),]),prop_mat=as.matrix(prop_mat)[-c(1:nrow(batchPCAdata)),c(1:nrow(batchPCAdata))]))
  
}

.myConcensusDEFn_step2_detail_exp_final_archive=function(inputData,argList,sd_offset=0.001,reg_sd=T,zscore_cap=10,minCellCountThr=4){
  
  #inputData=res[[1]];sd_offset=0.005;reg_sd=T
  
  myWeightedTtest=function(inputData,sd_offset,reg_sd,zscore_cap,minCellCountThr){
    require(Matrix,quietly = T)
    
    zscore_cap=abs(zscore_cap)
    if(zscore_cap<2){
      warning("zscore cap was set to below 2. It was changed to 5.")
      zscore_cap=5
    }
    sd_offset=max(sd_offset,0)
    
    geo_mean_count=.extra_gmean(inputData$data$countData)
    resGeneMeanSd=data.frame(gene=names(geo_mean_count),geo_mean_count=geo_mean_count,stringsAsFactors = F)
    
    
    myWeightedVar=function(inputWeight,inputExp){
      require(Matrix,quietly = T)
      inputExp=t(inputExp)
      exp_sq=inputExp*inputExp
      res=(inputWeight%*% exp_sq)
      res=sweep(res,1,rowSums(inputWeight),"*")
      res=res- (inputWeight%*% inputExp)^2
      res=sweep(res,1,rowSums(inputWeight)^2,"/")
      return(res)
    }
    
    prop_mat=inputData$prop_mat
    
    
    prop_mat_c=prop_mat
    if(quantile(apply(prop_mat_c,1,function(x) sum(x>0)),0.25)<(0.85*ncol(prop_mat))){
      prop_mat_c[prop_mat_c>0]=1
    }
    
    prop_mat_c=1-prop_mat_c
    prop_mat_c <- .extra_matrix_rowNorm(prop_mat_c)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat_c))) %*% prop_mat_c
    
    x_exp=prop_mat %*% t(inputData$data$logNormData)
    y_exp=prop_mat_c %*% t(inputData$data$logNormData)
    
    var_x=myWeightedVar(inputWeight = inputData$prop_mat,inputExp = inputData$data$logNormData)
    var_y=myWeightedVar(inputWeight = prop_mat_c,inputExp = inputData$data$logNormData)
    
    n_x=.myEffSizePropMat(prop_mat = prop_mat)
    n_x=n_x$effective_sample_size
    
    n_y=.myEffSizePropMat(prop_mat = prop_mat_c)
    n_y=n_y$effective_sample_size
    
    vxn=sweep(var_x,1,n_x,"/")
    vyn=sweep(var_y,1,n_y,"/")
    
    vxn2=sweep(var_x,1,n_x-1,"*")
    vyn2=sweep(var_y,1,n_y-1,"*")
    
    
    df=(vxn + vyn)^2
    df=df/(sweep(vxn^2,1,(n_x - 1),"/") + 
             sweep(vyn^2,1,n_y - 1,"/"))
    if(sum(vxn==0,na.rm = T)>0){
      df[which(vxn==0)]=matrix(n_y,nrow=nrow(vxn),ncol=ncol(vxn),byrow = F)[which(vxn==0)]
    }
    
    vxn[is.na(as.matrix(vxn))]=1
    sxy=sqrt(vxn+vyn)
    
    d_sxy=sqrt(sweep(vxn2+vyn2,1,n_x+n_y-2,"/"))
    
    sd_reg=sxy
    d_sd_reg=d_sxy
    if(reg_sd){
      sd_reg=lapply(1:nrow(sxy),function(x){
        y=sxy[x,]
        tmp_data=data.frame(exp=resGeneMeanSd$geo_mean_count,sd=sxy[x,])
        slInd=which(tmp_data$exp>0.001&tmp_data$sd>0)
        if(length(slInd)>0){
          tmp_data=tmp_data[slInd,]
          meanExp=tmp_data$exp
          design=splines::ns(log10(meanExp), df = 3, intercept = TRUE)
          fit2=lm.fit(design,log10(tmp_data$sd))
          tmp_data$fit2=pmax(tmp_data$sd,10^fit2$fitted.values)
          
          y[slInd]=tmp_data$fit2
        } else {
          y=pmax(y,1)
        }
        
        return(as.numeric(y))
      })
      sd_reg=do.call("rbind",sd_reg)
      
      d_sd_reg=lapply(1:nrow(d_sxy),function(x){
        y=d_sxy[x,]
        tmp_data=data.frame(exp=resGeneMeanSd$geo_mean_count,sd=d_sxy[x,])
        slInd=which(tmp_data$exp>0.001&tmp_data$sd>0)
        if(length(slInd)>0){
          tmp_data=tmp_data[slInd,]
          meanExp=tmp_data$exp
          design=splines::ns(log10(meanExp), df = 3, intercept = TRUE)
          fit2=lm.fit(design,log10(tmp_data$sd))
          tmp_data$fit2=pmax(tmp_data$sd,10^fit2$fitted.values)
          
          y[slInd]=tmp_data$fit2
        } else {
          y=pmax(y,1)
        }
        
        return(as.numeric(y))
      })
      d_sd_reg=do.call("rbind",d_sd_reg)
    }
    
    t <- (x_exp - y_exp)/(sd_reg+sd_offset)
    cohen_d=(x_exp - y_exp)/d_sd_reg
    hodge_g=cohen_d*(1-3/(4*(n_x+n_y)-9))
    if(sum(is.na(hodge_g))>0){
      hodge_g[is.na(hodge_g)]=0
    }
    
    se.g=0.5*sweep(hodge_g^2,1,n_x+n_y-3.94,"/")
    se.g <- sqrt(sweep(se.g,1, (n_x+n_y)/(n_x*n_y),"+"))
    if(sum(is.na(se.g))>0){
      se.g[is.na(se.g)]=0
    }
    
    #t2 <- (x_exp - y_exp)/(sd_reg)
    
    t_vec=as.numeric(t)
    zscore=qnorm(pt(as.numeric(abs(t_vec)), as.numeric(df),lower.tail = F,log.p = T),lower.tail = F,log.p = T)
    zscore[is.na(zscore)]=0
    zscore=zscore*sign(t_vec)
    zscore=matrix(zscore,nrow=nrow(t),ncol=ncol(t))
    
    colnames(zscore)=row.names(inputData$data$logNormData)
    row.names(zscore)=row.names(inputData$prop_mat)
    
    exp_binary=t(inputData$data$logNormData)
    exp_binary@x=rep(1,length(exp_binary@x))
    
    pct.1=as.matrix(prop_mat %*% exp_binary)
    pct.2=as.matrix(prop_mat_c %*% exp_binary)
    
    exp_norm=t(expm1(inputData$data$logNormData))
    
    exp.1=log2(prop_mat %*% exp_norm+1)
    exp.2=log2(prop_mat_c %*% exp_norm+1)
    
    logFC=as.matrix(exp.1 - exp.2)
    
    n=sqrt(sweep(pct.1,1,n_x,"*")*sweep(pct.2,1,n_y,"*"))
    
    if(sum(zscore>zscore_cap,na.rm = T)>0){
      zscore[which(zscore>zscore_cap)]=zscore_cap
    }
    
    if(sum(zscore<(-1*zscore_cap),na.rm = T)>0){
      zscore[which(zscore<(-1*zscore_cap))]=(-1*zscore_cap)
    }
    
    
    matWeights=.myEffSizePropMat(prop_mat)
    
    matEffectiveSize=matWeights$effective_sample_size
    matWeights=matWeights$centroid_weights
    matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
    matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
    
    
    matWeights=matrix(matWeights,byrow = F,nrow=nrow(zscore),ncol=ncol(zscore))
    matEffectiveSize=matrix(matEffectiveSize,byrow = F,nrow=nrow(zscore),ncol=ncol(zscore))
    row.names(matWeights)=row.names(matEffectiveSize)=row.names(prop_mat)
    colnames(matWeights)=colnames(matEffectiveSize)=colnames(zscore)
    matWeights[which(pct.1*matEffectiveSize<minCellCountThr&pct.2<0.001)]=0
    matWeights[which(zscore==0)]=0
    
    
    return(list(zscore=zscore,pct.1=pct.1,pct.2=pct.2,logFC=logFC,se.g=se.g,hodge_g=hodge_g,n=n,matWeights=matWeights,matEffectiveSize=matEffectiveSize))
  }
  
  tmp_weights=rowSums(as.matrix(inputData$prop_mat))
  
  if(sum(names(argList)=="do.split.prop")==0){
    argList$do.split.prop=T
  }
  
  if(!argList$do.split.prop){
    inputData$prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  }
  
  res=myWeightedTtest(inputData = inputData,sd_offset=sd_offset,reg_sd=reg_sd,zscore_cap = zscore_cap,minCellCountThr=minCellCountThr)
  if(!argList$do.split.prop){
    tmp_weights=matrix(tmp_weights,byrow = F,nrow=nrow(res$zscore),ncol=ncol(res$zscore))
    tmp_weights[res$matWeights==0]=0
    res$matWeights=tmp_weights
  }
  return(c(res,list(dsName=inputData$data$dsName,prop_mat=inputData$prop_mat)))
  
}

.myConcensusDEFn_step2_detail_exp_vLimmaTrend=function(inputData,argList,centroidPCAdata,exCentroids=NULL,min_cells=5,resGeneMeanSd=NULL,nDataset,within_gene_norm=F){
  
  
  #inputData=res[[1]];centroidPCAdata=pca_centroid;exCentroids=exCentroids;argList=argList;nDataset=length(dataArranged);within_gene_norm=F;resGeneMeanSd=NULL
  
  require(data.table,quietly = T)
  normExp=NULL
  
  batchPCAdata=as.matrix(inputData$data$pcaData)[,1:argList$nPCs]
  #expData=inputData$data$expData
  dsName=inputData$data$dsName
  print(dsName)
  if(!is.null(exCentroids)){
    centroidPCAdata=centroidPCAdata[-which(row.names(centroidPCAdata) %in% exCentroids),1:argList$nPCs]
  } else {
    centroidPCAdata=centroidPCAdata[,1:argList$nPCs]
  }
  
  prop_mat=inputData$prop_mat
  uniformZscore=argList$uniformZscore
  
  centroid_weights=rowSums(prop_mat)
  names(centroid_weights)=row.names(prop_mat)
  
  res2=""
  
  expNorm=inputData$data$logNormData
  prop_mat=inputData$prop_mat
  zscoreMat=matrix(0,ncol=nrow(prop_mat),nrow=nrow(expNorm))
  colnames(zscoreMat)=row.names(prop_mat)
  row.names(zscoreMat)=row.names(expNorm)
  
  LFCMat=avgExpMat=zscoreMat
  
  centroid_weights=.myEffSizePropMat(prop_mat)
  
  for(i in which(rowSums(prop_mat)>0)){
    designMat <- model.matrix(~prop_mat[i,])
    # ... after processing with lmFit/eBayes to get fit2 ...
    
    fit <- lmFit(as.matrix(expNorm), design=designMat)
    fit$df.residual=max(pmin(centroid_weights$effective_sample_size[i],10)-2,1)
    fit <- eBayes(fit, trend=TRUE,robust=T)
    resultsFit= topTable(fit, coef=2, n=Inf, sort.by="none")
    resultsFit=resultsFit[match(row.names(zscoreMat),row.names(resultsFit)),]
    zscore=resultsFit$P.Value/2
    zscore=qnorm(zscore,lower.tail = F)*sign(resultsFit$t)
    zscoreMat[,i]=zscore
    LFCMat[,i]=resultsFit$logFC
    avgExpMat[,i]=resultsFit$AveExpr
  }
  
  ###############
  #Local examination of the results
  save(avgExpMat,expNorm,fit,LFCMat,prop_mat,zscoreMat,file="~/myBucket/torm.rda")
  
  rm(list=ls())
  load("~/Desktop/torm.rda")
  load("~/Desktop/torm_prv.rda")
  
  all(row.names(res2)==row.names(zscoreMat))
  
  plot(res2[,25],zscoreMat[,25])
  
  x2=.myEffSizePropMat(prop_mat)
  
  x=sweep(res2,2,x2$centroid_weights,"*")
  
  x[x>10]=10
  plot(x[-2765,141],zscoreMat[-2765,141])
  
  
  
  plot(x[,1],zscoreMat[,1])
  plot(x[,61],zscoreMat[,61])
  
  fit$df.residual=2
  fit <- eBayes(fit, trend=TRUE)
  resultsFit= topTable(fit, coef=2, n=Inf, sort.by="none")
  slzscore=which(zscoreMat>1,arr.ind = T)
  slInd=14195
  slGene=row.names(zscoreMat)[slzscore[slInd,1]]
  slPseudo=colnames(zscoreMat)[slzscore[slInd,2]]
  
  #ENSG00000100412 ENSG00000077522 ENSG00000138107 ENSG00000151694 ENSG00000106948 ENSG00000149925
  slGene="ENSG00000281758"
  slPseudo="25"
  all(colnames(expNorm)==colnames(prop_mat))
  df=data.frame(sample=colnames(expNorm),exp=expNorm[slGene,],prop_score=prop_mat[slPseudo,],stringsAsFactors = F)
  
  plot(df$exp,df$prop_score)
  df=df[order(df$exp,decreasing = T),]
  head(df)
  table(df$exp>0,df$prop_score>0)
  zscoreMat[slGene,slPseudo]
  
  gc()
  return(list(res=zscoreMat,dsName=dsName,normExp=LFCMat,prop_mat=as.matrix(prop_mat)))
  
}

.myConcensusDEFn_step2_FindNeighbors=function (inputCentroids, argList, k.param = 20, compute.SNN = TRUE, jaccard.indx.thr = 1/15, nn.eps = 0, verbose = TRUE) {
  
  nPCs=min(30,ncol(inputCentroids))
  
  inputPCAscores=inputCentroids[,1:nPCs]
  
  n.cells <- nrow(x = inputPCAscores)
  
  if (verbose) {
    message("Computing nearest neighbor graph")
  }
  nn.ranked_org <- RANN::nn2(data = inputPCAscores, k = min(k.param*3,nrow(inputPCAscores)), eps = nn.eps)
  nn.ranked <- nn.ranked_org$nn.idx[,1:min(k.param,nrow(inputPCAscores)/3)]
  
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

.myConcensusDEFn_step2_modified=function(argList,expData=NULL){
  
  reRunCheck=T
  if(file.exists(.myFilePathMakerFn("res_DE_wZscore",argList=argList,uniformImportant=T,propImportant = T))|!file.exists(.myFilePathMakerFn("res_DE_wZscore_pathwayAnalysis",argList=argList,uniformImportant=T,propImportant = T))){
    reRunCheck=tryCatch({load(.myFilePathMakerFn("res_DE_wZscore",argList=argList,uniformImportant=T,propImportant = T));F}, error=function(e) {return(T)})
    if(!reRunCheck){
      reRunCheck=tryCatch({load(.myFilePathMakerFn("res_DE_wZscore_pathwayAnalysis",argList=argList,uniformImportant=T,propImportant = T));F}, error=function(e) {return(T)})
      
    }
  }
  
  if(reRunCheck){
    
    set.seed(1)
    supportingFractionThr=argList$DE_supportingFractionThr
    n.adaptiveKernel=argList$DE_n.adaptiveKernel
    nPropIter=argList$DE_nPropIter
    
    #load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("pca_centroids",argList=argList))
    
    if(!all(row.names(harmony_embeddings)==row.names(pd))){
      stop("Error in matching Names!")
    }
    
    pcaList=split(as.data.frame(harmony_embeddings),pd$anno_batch)
    
    if(is.null(expData)){
      load(.myFilePathMakerFn("exp_merged",argList=argList,expData=T))
      expData=SplitObject(tmp, split.by = "anno_batch")
    }
    
    
    if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename="prevalentGenes"),argList)))){
      cat("           Identifying prevalent genes ...\n")
      prevalenceEstimate=list()
      for(ik in 1:length(expData)){
        tmp=expData[[ik]]@assays$RNA@counts
        tmp=apply(tmp,1,function(x) sum(x>0)/length(x))*100
        prevalenceEstimate=c(prevalenceEstimate,list(data.frame(dsName=expData[[ik]]$anno_batch[1],gene=row.names(expData[[ik]]),prevalence=tmp,stringsAsFactors = F)))
      }
      prevalenceEstimate=do.call("rbind",prevalenceEstimate)
      prevalenceThr=max(1/(2*argList$internal_pseudocell_count),10/median(unlist(lapply(expData,ncol))))
      prevalenceEstimate=aggregate(prevalence~gene,data=prevalenceEstimate,function(x) sum(x>prevalenceThr))
      prevalenceEstimate=prevalenceEstimate[which(prevalenceEstimate$prevalence>=(supportingFractionThr*length(expData))),]
      prevalentGenes=as.character(prevalenceEstimate$gene)
      
      save(prevalentGenes,file=do.call('.myFilePathMakerFn',args=c(list(filename="prevalentGenes"),argList)))
    } else {
      load(do.call('.myFilePathMakerFn',args=c(list(filename="prevalentGenes"),argList)))
    }
    
    
    dataArranged=list()
    for(i in 1:length(expData)){
      dataArranged=c(dataArranged,list(list(logNormData=expData[[i]]@assays$RNA@data,countData=expData[[i]]@assays$RNA@counts,dsName=as.character(expData[[i]]$anno_batch[1]),pcaData=pcaList[[as.character(expData[[i]]$anno_batch[1])]])))
    }
    
    
    if(!file.exists(.myFilePathMakerFn("gene_fractionExp",argList=argList))){
      fractionExpressed=list()
      for(i in 1:length(expData)){
        tmp=apply(expData[[i]]@assays$RNA@counts,1,function(x) sum(x>0)/length(x))
        tmp=data.frame(gene=names(tmp),fraction=tmp,dsName=as.character(expData[[i]]$anno_batch[1]),stringsAsFactors = F)
        fractionExpressed=c(fractionExpressed,list(tmp))
      }
      fractionExpressedList=do.call("rbind",fractionExpressed)
      fractionExpressed=aggregate(fraction~gene,data=fractionExpressedList,mean)
      fractionExpressed$gene=as.character(fractionExpressed$gene)
      save(fractionExpressed,fractionExpressedList,file=.myFilePathMakerFn("gene_fractionExp",argList=argList))
      
    } else {
      load(.myFilePathMakerFn("gene_fractionExp",argList=argList))
    }
    
    #if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename="prop_mat_list_step2",uniformImportant=T),argList)))){
      cat("           Constructing the propagation matrices ...\n")
      #n.adaptiveKernel=5
      #res=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=nPropIter,mc.cores = argList$ncores)
      #res=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=nPropIter,n.neighbors=4,mc.cores = argList$ncores)
      
      #argList$prop.n.neighbors=4
      #centroidPCAdata=pca_centroid;argList=argList;n.adaptiveKernel=n.adaptiveKernel;nPropIter=1;n.neighbors=argList$prop.n.neighbors
      res=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_scale_v2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=argList$prop.n.neighbors,mc.cores = argList$ncores)
      
      .res=res
      myiNPHlbltransferFn=function(inputRes,argList){
        load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
        pd=pd[which(pd$ds_batch=="human_iNPH"),]
        counter=0
        cnames=unique(as.character(pd$anno_orig_cellState))
        res_prop=matrix(0,nrow=argList$internal_pseudocell_count,ncol=length(cnames))
        
        for(ik in 1:length(inputRes)){
          tmp=as.matrix(inputRes[[ik]]$prop_mat)
          if(sum(colnames(tmp) %in% row.names(pd))>0){
            counter=counter+1
            tmp_pd=pd[match(colnames(tmp) , row.names(pd)),]
            tmp_cellState=.myOneHotFn(factor(as.character(tmp_pd$anno_orig_cellState),levels=cnames))
            tmp_prop=tmp %*% as.matrix(tmp_cellState)
            tmp_prop=tmp_prop[,match(cnames,colnames(tmp_prop))]
            if(sum(is.na(tmp_prop))>0){
              tmp_prop[sum(is.na(tmp_prop))>0]=0
            }
            res_prop=res_prop+tmp_prop
            if(counter==1){
              row.names(res_prop)=row.names(tmp_prop)
              colnames(res_prop)=colnames(tmp_prop)
            }
          }
        }
        
        
      }
      
      myBICCNlbltransferFn=function(inputRes,argList){
        load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
        pd=pd[which(pd$ds_batch=="mouse_Yao_34616066_BICCN_sc10Xv2"),]
        pd$cluster_label=gsub(" ","\\.",pd$cluster_label)
        counter=0
        cnames=unique(as.character(pd$cluster_label))
        
        res_prop=matrix(0,nrow=argList$internal_pseudocell_count,ncol=length(cnames))
        
        for(ik in 1:length(inputRes)){
          tmp=as.matrix(inputRes[[ik]]$prop_mat)
          if(sum(colnames(tmp) %in% row.names(pd))>0){
            
            counter=counter+1
            tmp_pd=pd[match(colnames(tmp) , row.names(pd)),]
            tmp_cellState=.myOneHotFn(factor(as.character(tmp_pd$cluster_label),levels=cnames))
            tmp_prop=tmp %*% as.matrix(tmp_cellState)
            tmp_prop=tmp_prop[,match(cnames,colnames(tmp_prop))]
            if(sum(is.na(tmp_prop))>0){
              stop("here")
            }
            res_prop=res_prop+tmp_prop
            if(counter==1){
              row.names(res_prop)=row.names(tmp_prop)
              colnames(res_prop)=colnames(tmp_prop)
            }
          }
        }
        
        save(res_prop,file="~/myBucket/torm_BICCNv2_prop400n3.rda")
      }
      
      mySIMlbltransferFn=function(inputRes,argList){
        load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
        counter=0
        cnames=unique(as.character(pd$anno_orig_cellState))
        
        res_prop=matrix(0,nrow=argList$internal_pseudocell_count,ncol=length(cnames))
        
        for(ik in 1:length(inputRes)){
          tmp=apply(as.matrix(inputRes[[ik]]$prop_mat),2,function(x) {
            if(sum(x)>0){
              return(x/sum(x))
            } else {
              return(x)
            }
            
            })
          if(sum(colnames(tmp) %in% row.names(pd))>0){
            
            counter=counter+1
            tmp_pd=pd[match(colnames(tmp) , row.names(pd)),]
            tmp_cellState=.myOneHotFn(factor(as.character(tmp_pd$anno_orig_cellState),levels=cnames))
            tmp_prop=tmp %*% as.matrix(tmp_cellState)
            tmp_prop=tmp_prop[,match(cnames,colnames(tmp_prop))]
            if(sum(is.na(tmp_prop))>0){
              stop("here")
            }
            res_prop=res_prop+tmp_prop
            if(counter==1){
              row.names(res_prop)=row.names(tmp_prop)
              colnames(res_prop)=colnames(tmp_prop)
            }
          }
        }
        x=res_prop[order(res_prop[,4],decreasing = T),]
        head(x)
        
        save(res_prop,file="~/myBucket/torm_SIM_calcDE3_prop200n4.rda")
      }
      
      
      tmpValCheck=(unlist(lapply(res,length)))
      if(sum(tmpValCheck==1)>0){
        stop(res[[which(tmpValCheck==1)[1]]])
      }
      rm(tmpValCheck)
      
      ###Checking for outlier centroids
      
      
      
      exCentroids=NULL
      if(F){
        for(ipseudo in 1:length(res)){
          #tmp=diag(res[[ipseudo]]$prop_mat)
          #tmp=tmp[(length(tmp)-argList$internal_pseudocell_count+1):length(tmp)]
          #names(tmp)=colnames(res[[ipseudo]]$prop_mat)[(ncol(res[[ipseudo]]$prop_mat)-argList$internal_pseudocell_count+1):ncol(res[[ipseudo]]$prop_mat)]
          
          tmp=1-rowSums(res[[ipseudo]]$prop_mat)
          if(ipseudo==1){
            effSize_list=tmp
          } else {
            effSize_list=effSize_list+tmp
          }
        }
        
        if(length(which(effSize_list<(argList$DE_supportingFractionThr*length(dataArranged))))>0&argList$exclude_non_freq_pseudocells){
          exCentroids=names(effSize_list)[which(effSize_list<(argList$DE_supportingFractionThr*length(dataArranged)))]
          exCentroids=c("172","199")
          cat(paste0("           Excluded ",length(exCentroids)," centroids due to lack of supporting samples\n"))
          res=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=nPropIter,mc.cores = argList$ncores,exCentroids=exCentroids)
          #summary(res[[1]]$prop_mat["15",])
          #table(res[[1]]$prop_mat["15",]>0)
          tmpValCheck=(unlist(lapply(res,length)))
          if(sum(tmpValCheck==1)>0){
            stop(res[[which(tmpValCheck==1)[1]]])
          }
          rm(tmpValCheck)
          
        }
        
      }
      
      
      if(F){
        exCentroids=c("172","199")
        cat(paste0("           Excluded ",length(exCentroids)," centroids due to lack of supporting samples\n"))
        res=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=nPropIter,mc.cores = argList$ncores,exCentroids=exCentroids)
        #summary(res[[1]]$prop_mat["15",])
        #table(res[[1]]$prop_mat["15",]>0)
        tmpValCheck=(unlist(lapply(res,length)))
        if(sum(tmpValCheck==1)>0){
          stop(res[[which(tmpValCheck==1)[1]]])
        }
        rm(tmpValCheck)
      }
      
      cat("           Calculating dataset specific z-scores ...\n")
      #centroidPCAdata=pca_centroid;exCentroids=exCentroids;argList=argList;nDataset=length(dataArranged);within_gene_norm=F
      #108
      require(qs)
      res2=parallel::mclapply(res,.myConcensusDEFn_step2_detail_exp,centroidPCAdata=pca_centroid,exCentroids=exCentroids,argList=argList,nDataset=length(dataArranged),within_gene_norm=F,mc.cores = argList$ncores)#min(8,argList$ncores))
      qsave(res2,.myFilePathMakerFn("prop_res2",argList=argList,qsFormat =T))
      #res2=qread("~/data/results/InN_nPC30/res2.qs")
      #res2=qread(.myFilePathMakerFn("prop_res2",argList=argList,qsFormat =T))
      res=res2
      tmpValCheck=(unlist(lapply(res,length)))
      if(sum(tmpValCheck==1)>0){
        stop(res[[which(tmpValCheck==1)[1]]])
      }
      rm(tmpValCheck)
      
      
      
      #save(res,file=do.call('.myFilePathMakerFn',args=c(list(filename="prop_mat_list_step2",uniformImportant=T),argList)))
      resDistances=data.frame(centroid=names(res[[1]]$w_zscore_centroid),score=res[[1]]$w_zscore_centroid,dataset=res[[1]]$dsName,stringsAsFactors = F)
      if(length(res)>1){
        for(idsScore in 2:length(res)){
          tmp=data.frame(centroid=names(res[[idsScore]]$w_zscore_centroid),score=res[[idsScore]]$w_zscore_centroid,dataset=res[[idsScore]]$dsName,stringsAsFactors = F)
          resDistances=rbind(resDistances,tmp)
        }
      }
      
      resDistances=reshape2::dcast(centroid~dataset,data = resDistances,value.var = "score")
      save(resDistances,file=do.call('.myFilePathMakerFn',args=c(list(filename="centroid_scores",uniformImportant=T,propImportant=T),argList)))
      
      #################
      if(F){
        for(ik in 1:length(res)){
          if(sum(res[[ik]]$res>10)>0){
            res[[ik]]$res[which(res[[ik]]$res>10)]=10
          }
        }
        
        tmp_unique_genes=c()
        for(ik in 1:length(res)){
          tmp_unique_genes=unique(c(tmp_unique_genes,as.character(row.names(res[[ik]]$res))))
        }
        
        for(ik in 1:length(res)){
          res[[ik]]$res=res[[ik]]$res[match(tmp_unique_genes,row.names(res[[ik]]$res)),]
          row.names(res[[ik]]$res)=tmp_unique_genes
          if(sum(is.na(res[[ik]]$res))>0){
            res[[ik]]$res[is.na(res[[ik]]$res)]=0
          }
        }
      }
      
      
      ##########################
      #Counts per dataset
      if(F){
        res_count_binary_ds=list()
        ds_df=data.frame(dsName=pd$ds_batch,batchName=pd$anno_batch,stringsAsFactors = F)
        ds_df=ds_df[!duplicated(ds_df$batchName),]
        ds_df$dsName[which(ds_df$dsName=="mouse_Van_Hove_App_PS1_9mo")]=ds_df$batchName[which(ds_df$dsName=="human_iNPH")]
        for(i in 1:length(res)){
          
          tmp_dsCount=ds_df[which(ds_df$batchName==res[[i]]$dsName),]
          tmp_dsCount=ds_df[which(ds_df$dsName==tmp_dsCount$dsName[1]),]
          
          tmp_w_zscore=res[[i]]$w_zscore_centroid
          
          tmp_fraction=fractionExpressedList[which(fractionExpressedList$dsName ==res[[i]]$dsName),]
          tmp_fraction=tmp_fraction[match(row.names(res[[i]]$res),tmp_fraction$gene),]
          if(sum(is.na(tmp_fraction$fraction))>0){
            tmp_fraction$fraction[is.na(tmp_fraction$fraction)]=0
          }
          tmp_fraction$fraction[which(tmp_fraction$fraction<0.01)]=0
          tmp_fraction$fraction[which(tmp_fraction$fraction>=0.01)]=1
          tmp_fraction=matrix(tmp_fraction$fraction,nrow=nrow(res[[i]]$res),ncol=ncol(res[[i]]$res),byrow = F)
          
          tmpInd=apply(res[[i]]$res,2,function(x) sum(x==0)/length(x))
          
          if(sum(tmpInd>0.95)>0){tmp_w_zscore[which(tmpInd>0.95)]=0
          }
          
          tmp_w_zscore=tmp_fraction
          
          if(sum(is.na(tmp_w_zscore))>0){
            stop(paste("Check",i,"for NA values"))
          }
          
          if(sum(names(res_count_binary_ds)==tmp_dsCount$dsName[1])>0){
            dsInd=which(names(res_count_binary_ds)==tmp_dsCount$dsName[1])
            res_count_binary_ds[[dsInd]]=res_count_binary_ds[[dsInd]]+tmp_w_zscore
          } else {
            res_count_binary_ds=c(res_count_binary_ds,new=list(tmp_w_zscore))
            names(res_count_binary_ds)[which(names(res_count_binary_ds)=="new")]=tmp_dsCount$dsName[1]
          }
          
        }
        
        ds_df=ds_df[grepl("mouse_Van_Hove_App_PS1_9mo",ds_df$dsName),]
        
      }
      
      ##########################
      #merging the data at dataset level
      
      
      if(F){
        load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
        ds_df=data.frame(dsName=pd$ds_batch,batchName=pd$anno_batch,stringsAsFactors = F)
        ds_df=ds_df[!duplicated(ds_df$batchName),]
        nameList=unlist(lapply(res,function(x) x$dsName))
        ds_df=ds_df[match(nameList,ds_df$batchName),]
        ds_df$id=1:nrow(ds_df)
        #mergedResList=list()
        for(ik in unique(ds_df$dsName)){
          tmpInd=which(ds_df$dsName==ik)
          tmp_res=matrix(0,nrow=nrow(res[[1]]$res),ncol=ncol(res[[1]]$res))
          tmp_count=rep(0,length(res[[1]]$w_zscore_centroid))
          for(jk in tmpInd){
            tmp_res=tmp_res+res[[jk]]$res
            tmp_count=tmp_count+res[[1]]$w_zscore_centroid
          }
          tmp_res=tmp_res
          tmp_count=tmp_count
          mergedResList=c(mergedResList,list(list(dsName=ik,res=tmp_res,w_zscore_centroid=tmp_count)))
        }
        
        for(i in 1:length(mergedResList)){
          tmp_count=res_count_w_ds[[mergedResList[[i]]$dsName]]
          mergedResList[[i]]$res=mergedResList[[i]]$res/tmp_count
          mergedResList[[i]]$res[which(tmp_count==0)]=0
        }
        
        qsave(mergedResList,file=do.call('.myFilePathMakerFn',args=c(list(filename="res_datasetlevel",uniformImportant=T,propImportant=T,qsFormat=T),argList)))
        
      }
      
      #.res=res
      #res=mergedResList
      
      if(F){
        #Merging the data at AD/Ctrl level
        ds_df=ds_df[which(ds_df$dsName=="human_iNPH"),]
        anno=read.table("~/myBucket/SC_data/Human/brain/snRNA/iNPH/iNPH_anno_full.txt",sep="\t",header=T,stringsAsFactors=F)
        anno=anno[!duplicated(anno$subjectId),]
        anno$Id=paste0("human_iNPH_",anno$subjectId)
        ds_df=merge(ds_df,anno[,c("Id","Status")],by.x="batchName",by.y="Id")
        ds_df$dx="AD"
        ds_df$dx[ds_df$Status=="Ctrl"]="Ctrl"
        ds_df$dx2=paste0(ds_df$dx,"_2")
        ds_df$dx2[sample(which(ds_df$dx=="Ctrl"),length(which(ds_df$dx=="Ctrl"))/2)]="Ctrl_1"
        ds_df$dx2[sample(which(ds_df$dx=="AD"),length(which(ds_df$dx=="AD"))/2)]="AD_1"
        
        mergedResList=list()
        for(ik in unique(ds_df$dx2)){
          tmpInd=ds_df$id[ds_df$dx2==ik]
          tmp_res=matrix(0,nrow=nrow(res[[1]]$res),ncol=ncol(res[[1]]$res))
          for(jk in tmpInd){
            tmp_res=tmp_res+res[[jk]]$res
          }
          tmp_res=tmp_res/length(tmpInd)
          mergedResList=c(mergedResList,list(list(dsName=ik,res=tmp_res)))
        }
        qsave(mergedResList,file="~/myBucket/InN_nPC30/slRes_iNPH.qs")
        
      }
      
      ##########################
      
      #################
      
   # } else {
    #  load(do.call('.myFilePathMakerFn',args=c(list(filename="prop_mat_list_step2",uniformImportant=T),argList)))
    #}
      
    #Remormalizing the z-score dist of each gene
    
    cat("           Meta-analysis of z-scores ...\n")
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    
    
    res_count=matrix(0,nrow=nrow(res[[1]]$res),ncol(res[[1]]$res))
    res_count_w=matrix(0,nrow=nrow(res[[1]]$res),ncol(res[[1]]$res))
    res_count_sq=matrix(0,nrow=nrow(res[[1]]$res),ncol(res[[1]]$res))
    res_m=matrix(0,nrow=nrow(res[[1]]$res),ncol=ncol(res[[1]]$res))
    res_effectiveSize=matrix(0,nrow=nrow(res[[1]]$res),ncol=ncol(res[[1]]$res))
    res_count_w_ds=list()
    for(i in 1:length(res)){
      
      
      if(length(which(ds_df$batchName==res[[i]]$dsName))>0){
        tmp_dsCountName=ds_df[which(ds_df$batchName==res[[i]]$dsName),]
        #tmp_dsCount=ds_df[which(ds_df$dsName==tmp_dsCount$dsName[1]),]
        
        tmp_dsCount=res_count_binary_ds[[tmp_dsCountName$dsName[1]]]
        tmp_dsCount[which(tmp_dsCount<1)]=1
        tmp_w_zscore=res[[i]]$w_zscore_centroid
        tmp_count_mat=rep(1,ncol(res[[i]]$res))
        
        tmp_fraction=fractionExpressedList[which(fractionExpressedList$dsName ==res[[i]]$dsName),]
        tmp_fraction=tmp_fraction[match(row.names(res[[i]]$res),tmp_fraction$gene),]
        if(sum(is.na(tmp_fraction$fraction))>0){
          tmp_fraction$fraction[is.na(tmp_fraction$fraction)]=0
        }
        tmp_fraction$fraction[which(tmp_fraction$fraction<0.01)]=0
        tmp_fraction$fraction[which(tmp_fraction$fraction>=0.01)]=1
        tmp_fraction=matrix(tmp_fraction$fraction,nrow=nrow(res[[i]]$res),ncol=ncol(res[[i]]$res),byrow = F)
        
        tmp_res_m=res[[i]]$res
        tmp_res_m[which(tmp_fraction==0)]=0
        
        tmpInd=apply(res[[i]]$res,2,function(x) sum(x==0)/length(x))
        
        if(sum(tmpInd>0.95)>0){
          tmp_count_mat[which(tmpInd>0.95)]=0
          tmp_w_zscore[which(tmpInd>0.95)]=0
          tmp_w_zscore_sq[which(tmpInd>0.95)]=0
          tmp_res_m[which(tmpInd>0.95)]=0
        }
        
        tmp_w_zscore=sweep(tmp_fraction,2,tmp_w_zscore,"*")/tmp_dsCount
        
        #tmp_w_zscore=tmp_w_zscore/tmp_dsCount
        
        tmp_w_zscore_sq=tmp_w_zscore^2
        tmp_count_mat=sweep(tmp_fraction,2,tmp_count_mat,"*")
        tmp_count_mat=tmp_count_mat/tmp_dsCount
        tmp_count_mat[which(tmp_count_mat<0.001)]=0
        
        
        if(sum(is.na(tmp_w_zscore))>0){
          stop(paste("Check",i,"for NA values"))
        }
        
        res_m=res_m+tmp_res_m/tmp_dsCount
        
        res_count=res_count+tmp_count_mat
        res_count_w=res_count_w+tmp_w_zscore
        res_count_sq=res_count_sq+tmp_w_zscore_sq
        if(sum(names(res_count_w_ds)==tmp_dsCountName$dsName[1])>0){
          dsInd=which(names(res_count_w_ds)==tmp_dsCountName$dsName[1])
          res_count_w_ds[[dsInd]]=res_count_w_ds[[dsInd]]+tmp_w_zscore
        } else {
          res_count_w_ds=c(res_count_w_ds,new=list(tmp_w_zscore))
          names(res_count_w_ds)[which(names(res_count_w_ds)=="new")]=tmp_dsCountName$dsName[1]
        }
        
      }
      
    }
    gc()
    #res_count=t(matrix(res_count,nrow=ncol(res[[i]]$res),ncol=nrow(res[[i]]$res)))
    #res_count_w=t(matrix(res_count_w,nrow=ncol(res[[i]]$res),ncol=nrow(res[[i]]$res)))
    #res_count_sq=t(matrix(res_count_sq,nrow=ncol(res[[i]]$res),ncol=nrow(res[[i]]$res)))
    
    ds_count_sq=matrix(0,nrow=nrow(res_count_w_ds[[1]]),ncol=ncol(res_count_w_ds[[1]]))
    for(i in 1:length(res_count_w_ds)){
      if(sum(ds_df$dsName==names(res_count_w_ds))>0){
        ds_count_sq=ds_count_sq+(res_count_w_ds[[i]])^2
      }
    }
    
    res_count_sq=matrix(0,nrow=nrow(res_count_w_ds[[1]]),ncol=ncol(res_count_w_ds[[1]]))
    for(i in 1:length(res_count_w_ds)){
      res_count_sq=res_count_sq+res_count_w_ds[[i]]^2
    }
    
    tst=0
    for(i in 1:length(res_count_w_ds)){
      tst=tst+res_count_w_ds[[i]][13507,40]
    }
    
    #res_effectiveSize=res_count_w^2/res_count_sq
    res_effectiveSize=res_count_w^2/res_count_sq
    res_effectiveSize[which(res_count_w==0)]=0
    
    if(length(which(res_effectiveSize<0.001))>0){
      res_effectiveSize2=res_effectiveSize
      res_effectiveSize2[which(res_effectiveSize<0.001)]=0
      res_effectiveSize=res_effectiveSize2
      rm(res_effectiveSize2)
      
    }
    
    #tmp_res_count_sq=sqrt(res_count_sq)
    tmp_res_count_sq=sqrt(ds_count_sq)
    tmp_res_count_sq[which(tmp_res_count_sq<1)]=1
    res_m=res_m/tmp_res_count_sq
    res_m[which(res_count==0)]=0
    
    #Homogeneity analysis
    if(F){
      res_homogeneity_analysis=NULL
      for(ihomogeneity in 1:length(res)){
        tmp_dataN=res[[ihomogeneity]]$normExp
        tmp_prop_mat=res[[ihomogeneity]]$prop_mat
        tmp_zscore=res[[ihomogeneity]]$res
        for(ipseudocell in 1:ncol(res_m)){
          sl_cell_ind=which(tmp_prop_mat[ipseudocell,]>0)
          if(length(sl_cell_ind)>5 & max(tmp_zscore[,ipseudocell])>1){
            pseudo_Ind=res_m[,ipseudocell]
            pseudo_Ind=which(pseudo_Ind>0)
            tmp_cor=cor(res_m[pseudo_Ind,ipseudocell],tmp_dataN[pseudo_Ind,sl_cell_ind])
            tmp_cor=sum(tmp_prop_mat[ipseudocell,sl_cell_ind]*tmp_cor)/sum(tmp_prop_mat[ipseudocell,sl_cell_ind])
            res_homogeneity_analysis=rbind(res_homogeneity_analysis,data.frame(dataset=res[[ihomogeneity]]$dsName,centroid=colnames(res_m)[ipseudocell],cor=tmp_cor,stringsAsFactors = F))
          } else {
            res_homogeneity_analysis=rbind(res_homogeneity_analysis,data.frame(dataset=res[[ihomogeneity]]$dsName,centroid=colnames(res_m)[ipseudocell],cor=NA,stringsAsFactors = F))
          }
          
        }
      }
    }
    
    
    ind=which(res_count!=0,arr.ind = T)
    
    res_arranged=data.frame(gene_index=ind[,1],centroid_index=ind[,2],gene=row.names(res[[1]]$res)[ind[,1]],centroid=colnames(res[[1]]$res)[ind[,2]],zscore=res_m[ind],count=res_count[ind],count_w=res_count_w[ind],effective_size=res_effectiveSize[ind],stringsAsFactors = F)
    fractionExpressed=fractionExpressed[match(res_arranged$gene,fractionExpressed$gene),]
    res_arranged$overall_fractionExpressed=fractionExpressed$fraction
    #res_arranged=res_arranged[which(res_arranged$effective_size>=max(1,supportingFractionThr*length(expData))),]
    #,score_seq=res_indBase[ind]
    
    res_fd=NULL
    for(i in 1:length(expData)){
      fd=as.data.frame(expData[[i]]@assays$RNA@meta.features)
      fd$ensembl_gene_id=gsub("_","-",fd$ensembl_gene_id)
      slCols=c("gene_name","gene_biotype","symbol","gene_short_name","ensembl_gene_id")
      fd=fd[,colnames(fd) %in% slCols]
      if(sum(is.na(fd$ensembl_gene_id))>0){
        fd$ensembl_gene_id[is.na(fd$ensembl_gene_id)]=row.names(fd)[is.na(fd$ensembl_gene_id)]
      }
      
      if(!is.null(res_fd)){
        fd=fd[!fd$ensembl_gene_id %in% res_fd$ensembl_gene_id,]
        if(nrow(fd)>0){
          res_fd=rbind(res_fd,fd)
        }
      } else {
        res_fd=fd
      }
    }
    fd=res_fd
    rm(res_fd)
    
    fd=fd[match(res_arranged$gene,fd$ensembl_gene_id),]
    res_arranged=cbind(res_arranged,fd[,-which(colnames(fd)=="ensembl_gene_id")])
    
    
    
    
    network=.myConcensusDEFn_step2_FindNeighbors(inputCentroids = pca_centroid,argList=argList,verbose = F)
    
    snnNet=as.matrix(network$snn)
    ind=which(snnNet>0,arr.ind = T)
    
    snnNet=data.frame(Fnode=row.names(snnNet)[ind[,1]],Snode=row.names(snnNet)[ind[,2]],weight=snnNet[ind],stringsAsFactors = F)
    snnNet=snnNet[snnNet$Fnode!=snnNet$Snode,]
    snnNet$weight=1-snnNet$weight
    snnNet2=igraph::graph_from_data_frame(snnNet, directed = F, vertices = NULL)
    
    snnNet2=igraph::distances(snnNet2, mode = c("all"), weights = snnNet$weight, algorithm = "dijkstra")
    if(length(which(snnNet2==Inf))>0){
      snnNet2[which(snnNet2==Inf)]=(max(snnNet2[snnNet2!=Inf],na.rm = T)+1)
    }
    
    geneList=unique(as.character(res_arranged$gene))
    
    
    cat("           Calculating cor and MI ...\n")
    resCorMI=NULL
    slGenes=aggregate(zscore~gene,data=res_arranged,function(x) max(x))
    slGenes=slGenes[which(slGenes$zscore>2),]
    
    myCorMIfn=function(inputGenes,res_arranged,snnNet2){
      resCorMI=NULL
      for(i in inputGenes){
        tmp=res_arranged[which(res_arranged$gene==i),]
        if(length(which(tmp$effective_size<1))>0){
          tmp$zscore[which(tmp$effective_size<1)]=0
        }
        
        tmpNet=snnNet2[row.names(snnNet2) %in% tmp$centroid,colnames(snnNet2) %in% tmp$centroid]
        tmp=tmp[match(colnames(tmpNet),tmp$centroid),]
        tmp$zscore[is.na(tmp$zscore)]=0
        if(max(abs(tmp$zscore))>2){
          zOrdered=tmp$zscore[order(tmp$zscore,decreasing = T)]
          
          tmpInd=which(row.names(tmpNet) %in% tmp$centroid[which(tmp$zscore>=max(zOrdered[3],1))])
          if(length(tmpInd)>0){
            if(length(tmpInd)>1){
              distMat=tmpNet[tmpInd,]
              distMat=apply(distMat,2,mean)
            } else {
              distMat=tmpNet[tmpInd,]
            }
            
            
            tmp=data.frame(distance=distMat,heat=tmp$zscore)
            tmpCor=cor(tmp$distance,tmp$heat)
            tmpMI=infotheo::discretize( tmp, disc="equalfreq", nbins=10 )
            tmpMI=infotheo::mutinformation(tmpMI$distance,tmpMI$heat)
            resCorMI=rbind(resCorMI,data.frame(gene=i,cor=tmpCor,MI=tmpMI,stringsAsFactors = F))
            
          } else{
            resCorMI=rbind(resCorMI,data.frame(gene=i,cor=0,MI=0,stringsAsFactors = F))
          }
        }
      }
      return(resCorMI)
    }
    
    if(argList$ncores>1){
      geneList=split(intersect(as.character(slGenes$gene),prevalentGenes), cut(1:length(intersect(as.character(slGenes$gene),prevalentGenes)), argList$ncores, labels = FALSE)) 
    } else {
      geneList=list(intersect(as.character(slGenes$gene),prevalentGenes))
    }
    
    resCorMI=parallel::mclapply(geneList,myCorMIfn,res_arranged=res_arranged,snnNet2=snnNet2,mc.cores = argList$ncores)
    resCorMI=do.call("rbind",resCorMI)
    
    if(length(setdiff(unique(prevalentGenes),resCorMI$gene))>0){
      resCorMI=rbind(resCorMI,data.frame(gene=setdiff(unique(prevalentGenes),resCorMI$gene),cor=0,MI=0,stringsAsFactors = F))
    }
    
    #res_arranged=res_arranged[which(as.character(res_arranged$gene) %in%prevalentGenes),]
    res_arranged=merge(res_arranged,resCorMI,by="gene",all.x=T)
    
    
    res_arranged2=res_arranged
    res_arranged2=res_arranged2[order(res_arranged2$zscore,decreasing = T),]
    save(res_arranged,file=.myFilePathMakerFn("res_DE_wZscore_pathwayAnalysis",argList=argList,uniformImportant=T,propImportant = T))
    
    #zscore_count measures in how many centroids each gene as a nominal significance
    zscore_counts=aggregate(zscore~gene,data=res_arranged,function(x) sum(x>2))
    colnames(zscore_counts)=c("gene","zscore_count")
    res_arranged=merge(res_arranged,zscore_counts,by="gene",all.x=T)
    zscore_counts=aggregate(zscore~gene,data=res_arranged,function(x) sum(abs(x)>2))
    colnames(zscore_counts)=c("gene","zscore_count_abs")
    res_arranged=merge(res_arranged,zscore_counts,by="gene",all.x=T)
    if(sum(res_arranged$zscore_count>0)>0){
      res_arranged=res_arranged[res_arranged$zscore_count>0,]
    }
    
    res_arranged=res_arranged[res_arranged$effective_size/max(res_arranged$effective_size,na.rm = T)>=supportingFractionThr,]
    res_arranged=res_arranged[order(res_arranged$zscore,decreasing = T),]
    
    res_arranged$score_seq="0"
    for(i in 1:length(res)){
      res_arranged$score_seq=paste0(res_arranged$score_seq,",",round(res[[i]]$res[cbind(res_arranged$gene_index,res_arranged$centroid_index)],2))
    }
    res_arranged$score_seq=gsub(",0,",",",res_arranged$score_seq)
    res_arranged$score_seq=gsub("^0,","",res_arranged$score_seq)
    res_arranged$score_seq=gsub(",0$","",res_arranged$score_seq)
    
    save(res_arranged,file=.myFilePathMakerFn("res_DE_wZscore",argList=argList,uniformImportant=T,propImportant = T))
    gc()
    
  }
  
  return("Done!")
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

.myConsensusDEClusteringFn=function(argList,pval_thr,initial_clustering=NULL){
  
  library("gplots")
  library("devtools")
  library(ggnetwork)
  library(network)
  library(sna)
  library(ggplot2)
  
  myCorFn=function(x,zThr=0){
    res=matrix(0,nrow=ncol(x),ncol=ncol(x))
    row.names(res)=colnames(x)
    colnames(res)=colnames(x)
    tmpNames=colnames(res)
    for(i in tmpNames){
      for(j in tmpNames){
        if(sum(x[,i]>zThr|x[,j]>zThr)>5){
          res[i,j]=cor(x[x[,i]>zThr|x[,j]>zThr,i],x[x[,i]>zThr|x[,j]>zThr,j])
        } else {
          if(sum(x[,i]>qnorm(pval_thr,lower.tail = F)|x[,j]>qnorm(pval_thr,lower.tail = F))>0){
            res[i,j]=sum(x[,i]>qnorm(0.01,lower.tail = F)&x[,j]>qnorm(0.01,lower.tail = F))/sum(x[,i]>qnorm(0.01,lower.tail = F)|x[,j]>qnorm(0.01,lower.tail = F))
          } else {
            res[i,j]=0
          }
        }
        
      }
    }
    return(res)
  }
  
  load(.myFilePathMakerFn("res_DE_wZscore",argList=argList,uniformImportant=T))
  res=res_arranged
  rm(res_arranged)
  
  if(file.exists(do.call('.myFilePathMakerFn',args=c(list(filename="centroid_scores",uniformImportant=T),argList)))){
    load(do.call('.myFilePathMakerFn',args=c(list(filename="centroid_scores",uniformImportant=T),argList)))
    load(.myFilePathMakerFn("UMAP_centroid",argList=argList))
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    
    pd_summary=.scatterPlot_summary2d(object=pd,reductionCols=c("UMAP_1","UMAP_2"),n=100)
    resDistances=reshape::melt(resDistances)
    resDistances=merge(resDistances,UMAP_centroid,by="centroid",all.x=T)
    p=ggplot(resDistances,aes(UMAP_1,UMAP_2))+geom_point(color="lightgrey",data=pd_summary,aes(UMAP_1,UMAP_2))+geom_point(shape=21,aes(fill=value,size=value))+scale_fill_gradient2(low = "white",mid = "orange",midpoint = 0.5,high = "red")+facet_wrap(~variable)+theme_bw()+theme_classic()+scale_size(range = c(0.1,3))+theme(panel.grid = element_blank())+scale_x_continuous(limits = c(min(resDistances$UMAP_1)-0.1,max(resDistances$UMAP_1)+0.1))+scale_y_continuous(limits = c(min(resDistances$UMAP_2)-0.1,max(resDistances$UMAP_2)+0.1))
    ggsave(p,file=.myFilePathMakerFn("consensusDE_centroid_contribution",argList=argList,uniformImportant=T,pdf = T,makePlotDir=T),width=24,height=18)
    
  }
  
  
  #constructing the network
  sl_gene=res$gene[which(res$zscore>=2.5)]
  #sl_gene=res$gene[which(res$effective_size>=max(2,(argList$DE_supportingFractionThr*max(res$effective_size)))&res$zscore_count>=max(2,(argList$DE_supportingFractionThr*max(res$effective_size)))&res$cor<(-0.2))]
  sl_gene=unique(sl_gene)
  if(length(sl_gene)==0){
    stop("No marker gene was identified!")
  }
  res=res[res$gene %in% sl_gene,]
  res_arranged=reshape2::dcast(gene~centroid,data = res,value.var = "zscore")
  row.names(res_arranged)=res_arranged[,1]
  res_arranged=res_arranged[,-1]
  res_arranged[is.na(res_arranged)]=0
  
  scale="none"
  colramp = c(colorRampPalette(c("blue","black"))(4),"black",colorRampPalette(c("black","yellow"))(4))
  
  #Define custom dist and hclust functions for use with heatmaps
  mydist=function(c) {as.dist(1-myCorFn(t(c)))}
  myclust=function(d) {hclust(d, method="average")}
  myCentroidNetConstructorFn=function(input_centroid_cor,input_centroid_data,cluster_count=cluster_count){
    require(network)
    require(ggnetwork)
    resNet=NULL
    resNet2=NULL
    for(iclust in 1:(nrow(input_centroid_cor))){
      tmp=as.numeric(input_centroid_cor[iclust,])
      names(tmp)=colnames(input_centroid_cor)
      tmp=tmp[-iclust]
      for(iclust2 in names(tmp)){
        if(sum(input_centroid_data[,iclust2]>2&input_centroid_data[,row.names(input_centroid_cor)[iclust]]>2)>0){
          if(fisher.test(table(input_centroid_data[,iclust2]>2,input_centroid_data[,row.names(input_centroid_cor)[iclust]]>2))$estimate<2){
            tmp[iclust2]=0
          }
        } else {
          tmp[iclust2]=0
        }
        
      }
      tmp=tmp[order(tmp,decreasing = T)]
      tmp=tmp[which(tmp>0.5)]
      
      if(length(tmp)>0){
        tmp=tmp[1:min(length(tmp),10)]
        resNet=rbind(resNet,data.frame(Fnode=row.names(input_centroid_cor)[iclust],Snode=names(tmp),cor=tmp,stringsAsFactors = F))
        tmp=tmp[1:min(length(tmp),5)]
        resNet2=rbind(resNet2,data.frame(Fnode=row.names(input_centroid_cor)[iclust],Snode=names(tmp),cor=tmp,stringsAsFactors = F))
        
      }
    }
    resNet=resNet[!duplicated(.myNetCombination(resNet)),]
    resNet2=resNet2[!duplicated(.myNetCombination(resNet2)),]
    
    resAdj=rbind(resNet,data.frame(Fnode=resNet$Snode,Snode=resNet$Fnode,cor=resNet$cor,stringsAsFactors = F))
    #resAdj$cor=resAdj$cor^3
    resAdj=reshape2::dcast(Fnode~Snode,data=resAdj,value.var = "cor")
    row.names(resAdj)=resAdj[,1]
    resAdj=resAdj[,-1]
    resAdj=as.matrix(resAdj)
    resAdj[is.na(resAdj)]=0
    
    resAdj=resAdj[,match(row.names(resAdj),colnames(resAdj))]
    
    #resClustering=.netFindClusters(inputGraph=resAdj, algorithm = 1,resolution = clustering_res,group.singletons = F,modularity.fxn = 1)
    #resClustering=data.frame(centroid=row.names(resClustering),cluster=as.character(resClustering[,1]),stringsAsFactors = F)
    
    rmlooselyConnected=T
    while(rmlooselyConnected){
      rmlooselyConnected=F
      tormlist=apply(resAdj,1,function(x) sum(x>0)/length(x))
      if(sum(tormlist<0.01)>0){
        rmlooselyConnected=T
        resAdj=resAdj[-which(tormlist<0.01),-which(tormlist<0.01)]
      }
      
    }
    
    #resClustering=kmeans(input_centroid_cor,cluster_count)
    #resClustering=data.frame(centroid=names(resClustering$cluster),cluster=as.character(resClustering$cluster),stringsAsFactors = F)
    resClustering=hclust(as.dist(1-resAdj),method = "average")
    resClustering=cutree(resClustering,cluster_count)
    resClustering=data.frame(centroid=names(resClustering),cluster=as.character(resClustering),stringsAsFactors = F)
    
    resClustering$centroid=as.character(resClustering$centroid)
    
    resNet=resNet2
    
    for(isingletons in unique(resClustering$cluster)){
      if(sum(resClustering$cluster==isingletons)==1){
        if(sum(c(resNet$Fnode,resNet$Snode)==isingletons)>0){
          tmpNet=which(resNet$Fnode==isingletons|resNet$Snode==isingletons)
          tmpNet=resNet[tmpNet,]
          tmpNet=c(tmpNet$Fnode,tmpNet$Snode)
          tmpNet=as.data.frame(table(tmpNet))
          colnames(tmpNet)[1]="Var1"
          tmpNet$Var1=as.character(tmpNet$Var1)
          tmpNet=tmpNet[-which(tmpNet$Var1==isingletons),]
          tmpNet=tmpNet[order(tmpNet$Freq,decreasing = T),]
          tmpNet=tmpNet[tmpNet$Var1 %in% resClustering$centroid,]
          if(nrow(tmpNet)>0){
            resClustering$cluster[resClustering$cluster==isingletons]=resClustering$cluster[resClustering$centroid==tmpNet$Var1[1]]
          } else {
            resClustering=resClustering[-which(resClustering$cluster==isingletons),]
          }
          
        }
      }
    }
    
    clusters_identified=unique(resClustering$cluster)
    
    
    exCentroids=setdiff(row.names(input_centroid_cor),resClustering$centroid)
    if(length(exCentroids)>0){
      resClustering=rbind(resClustering,data.frame(centroid=exCentroids,cluster=max(as.numeric(resClustering$cluster))+1,stringsAsFactors = F))
    }
    
    
    net = network::network(resNet, directed = FALSE,matrix.type="edgelist")
    
    netVerNames=network::network.vertex.names(net)
    net %v% "cluster" = resClustering$cluster[match(netVerNames,as.character(resClustering$centroid))]
    network::set.edge.attribute(net, "weight", resNet$cor^3)
    
    load(.myFilePathMakerFn("UMAP_centroid",argList=argList))
    centroid_layout=UMAP_centroid[match(as.character(netVerNames),as.character(UMAP_centroid$centroid)),-1]
    centroid_layout=as.matrix(centroid_layout)
    colnames(centroid_layout)=c("x","y")
    
    net=ggnetwork:::fortify.network(net,layout = centroid_layout)
    
    resClustering$color=c(hues::iwanthue(length(clusters_identified)),"#808080")[as.numeric(resClustering$cluster)]
    
    clusterCol=resClustering[!duplicated(resClustering$cluster),]
    
    
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    pd_summary=.scatterPlot_summary2d(object=pd,reductionCols=c("UMAP_1","UMAP_2"))
    
    lda.training.data=merge(resClustering,UMAP_centroid,by="centroid")
    #lda.fit=MASS::lda(cluster~UMAP_1+UMAP_2,data=lda.training.data)
    #lda.pred=predict(lda.fit,pd_summary)
    
    knet=RANN::nn2(query=pd_summary[,c("UMAP_1","UMAP_2")],data = lda.training.data[,c("UMAP_1","UMAP_2")],k=12,eps=0)
    affinity_mat=matrix(0,nrow=nrow(knet$nn.idx),ncol=length(unique(lda.training.data$cluster)))
    colnames(affinity_mat)=unique(lda.training.data$cluster)
    
    tmp_affinities=t(apply(knet$nn.dists,1,function(x) exp((-1)*(x/x[3])^2)))
    tmp_affinities=t(apply(tmp_affinities,1,function(x) x/sum(x)))
    for(itr in 1:ncol(knet$nn.idx)){
      tst=.myOneHotFn(inputVector=factor(lda.training.data$cluster[knet$nn.idx[,itr]],levels=unique(lda.training.data$cluster)))
      tst=tst[,match(colnames(affinity_mat),colnames(tst))]
      tst=tmp_affinities[,itr]*tst
      affinity_mat=affinity_mat+tst
      rm(tst)
    }
    
    col_rgb_pallette=resClustering[!duplicated(resClustering$cluster),]
    col_rgb_pallette=col_rgb_pallette[,c("cluster","color")]
    col_rgb_pallette=col_rgb_pallette[match(colnames(affinity_mat),col_rgb_pallette$cluster),]
    col_rgb1=matrix(grDevices::col2rgb(col_rgb_pallette$color)[1,],nrow=nrow(affinity_mat),ncol=nrow(col_rgb_pallette),byrow = T)
    col_rgb2=matrix(grDevices::col2rgb(col_rgb_pallette$color)[2,],nrow=nrow(affinity_mat),ncol=nrow(col_rgb_pallette),byrow = T)
    col_rgb3=matrix(grDevices::col2rgb(col_rgb_pallette$color)[3,],nrow=nrow(affinity_mat),ncol=nrow(col_rgb_pallette),byrow = T)
    
    colnames(col_rgb1)=col_rgb_pallette$cluster
    colnames(col_rgb2)=col_rgb_pallette$cluster
    colnames(col_rgb3)=col_rgb_pallette$cluster
    
    col_rgb1=affinity_mat*col_rgb1
    col_rgb2=affinity_mat*col_rgb2
    col_rgb3=affinity_mat*col_rgb3
    
    col_rgb1=rowSums(col_rgb1)
    col_rgb2=rowSums(col_rgb2)
    col_rgb3=rowSums(col_rgb3)
    
    pd_summary$color=rgb(red=col_rgb1,green=col_rgb2,blue=col_rgb3,maxColorValue = 255)
    
    
    scale_factor=net[!duplicated(net$vertex.names),]
    scale_factor=merge(scale_factor,UMAP_centroid,by.x="vertex.names",by.y="centroid")
    scale_factor1=lm(x~UMAP_1,data=scale_factor)
    pd_summary$UMAP_1=predict(scale_factor1,newdata=pd_summary)
    scale_factor2=lm(y~UMAP_2,data=scale_factor)
    pd_summary$UMAP_2=predict(scale_factor2,newdata=pd_summary)
    
    pd_summary$xend=pd_summary$UMAP_1
    pd_summary$yend=pd_summary$UMAP_2
    
    #predicting the background color
    net=merge(net,clusterCol,by="cluster")
    p=ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) + geom_point(data=pd_summary,aes(UMAP_1,UMAP_2,color=color),alpha=0.5)+
      #geom_rect(aes(xmin=min(c(pd$UMAP_1,net$x)),xmax=max(c(pd$UMAP_1,net$x)),ymin=min(c(pd$UMAP_1,net$x)),ymax=max(c(pd$UMAP_1,net$x)),color="white"),alpha=0.003,fill="white")+
      geom_edges( color = "grey50",aes(size=weight)) +
      geom_point(data=net[is.na(net$cluster),],aes(x,y,fill=color),size=2,color="black",shape=21)+
      geom_nodelabel(aes(color = color, label = as.character(cluster)),fontface = "bold")+
      theme_blank()+scale_size_continuous(range = c(0.04,0.3))+scale_color_identity()+theme(legend.position = "none")
    ggsave(plot=p,.myFilePathMakerFn("consensusDE_clustering_network",argList=argList,uniformImportant=T,pdf = T,makePlotDir=T),width = 9,height = 9)
    return(list(net=resNet,clustering_res=resClustering))
  }
  myCentroidNetConstructorFn_new=function(input_centroid_cor,input_centroid_data){
    require(network)
    require(ggnetwork)
    resNet_10=NULL
    resNet_5=NULL
    for(icentroid in row.names(input_centroid_cor)){
      tmp=as.numeric(input_centroid_cor[icentroid,])
      names(tmp)=colnames(input_centroid_cor)
      tmp=tmp[-which(colnames(input_centroid_cor)==icentroid)]
      
      for(icentroid2 in names(tmp)){
        if(sum(input_centroid_data[,icentroid2]>2&input_centroid_data[,icentroid]>2)>0){
          if(fisher.test(table(input_centroid_data[,icentroid2]>2,input_centroid_data[,icentroid]>2))$estimate<2){
            tmp[icentroid2]=0
          }
        } else {
          tmp[icentroid2]=0
        }
      }
      tmp=tmp[order(tmp,decreasing = T)]
      tmp=tmp[which(tmp>0.5)]
      if(length(tmp)>2){
        tmp_scale=scale(as.numeric(tmp))
        tmp=tmp[which(tmp_scale>(-1.96))]
      }
      
      if(length(tmp)>0){
        tmp=tmp[1:min(length(tmp),10)]
        resNet_10=rbind(resNet_10,data.frame(Fnode=icentroid,Snode=names(tmp),cor=tmp,stringsAsFactors = F))
        tmp=tmp[1:min(length(tmp),5)]
        resNet_5=rbind(resNet_5,data.frame(Fnode=icentroid,Snode=names(tmp),cor=tmp,stringsAsFactors = F))
        
      }
    }
    resNet_5=resNet_5[!duplicated(.myNetCombination(resNet_5)),]
    resNet_10=resNet_10[!duplicated(.myNetCombination(resNet_10)),]
    for(ij in 1:nrow(resNet_10)){
      resNet_10$commonMarkers[ij]=(sum((input_centroid_data[,resNet_10$Fnode[ij]]>2.5&input_centroid_data[,resNet_10$Snode[ij]]>2) | (input_centroid_data[,resNet_10$Fnode[ij]]>2&input_centroid_data[,resNet_10$Snode[ij]]>2.5))+0.5)/(sum(input_centroid_data[,resNet_10$Fnode[ij]]>2.5|input_centroid_data[,resNet_10$Snode[ij]]>2.5)+0.5)
    }
    for(ij in 1:nrow(resNet_5)){
      resNet_5$commonMarkers[ij]=(sum((input_centroid_data[,resNet_5$Fnode[ij]]>2.5&input_centroid_data[,resNet_5$Snode[ij]]>2) | (input_centroid_data[,resNet_5$Fnode[ij]]>2&input_centroid_data[,resNet_5$Snode[ij]]>2.5))+0.5)/(sum(input_centroid_data[,resNet_5$Fnode[ij]]>2.5|input_centroid_data[,resNet_5$Snode[ij]]>2.5)+0.5)
    }
    resNet_10=resNet_10[which(resNet_10$commonMarkers>0.3),]
    resNet_5=resNet_5[which(resNet_5$commonMarkers>0.3),]
    #selecting the connected components
    sg=igraph::graph_from_data_frame(resNet_10[,c("Fnode","Snode")],directed = F)
    sg <- igraph::decompose.graph(sg)
    
    initial_clusters=NULL
    for(icluster in 1:length(sg)){
      tmp_net=igraph::as_data_frame(sg[[icluster]])
      tmp_net$cluster=paste0("C_",icluster)
      initial_clusters=rbind(initial_clusters,tmp_net)
      rm(tmp_net)
    }
    rm(icluster,sg)
    colnames(initial_clusters)[1:2]=c("Fnode","Snode")
    tmp_net=resNet_10[match(.myNetCombination(initial_clusters),.myNetCombination(resNet_10)),]
    initial_clusters$cor=tmp_net$cor
    rm(tmp_net)
    
    resClustering=NULL
    for(iclust in unique(initial_clusters$cluster)){
      tmp_resNet_10=initial_clusters[which(initial_clusters$cluster==iclust),]
      if(nrow(tmp_resNet_10)>4){
        tmp_seeds=tmp_resNet_10
        sg=igraph::graph_from_data_frame(tmp_seeds[,c("Fnode","Snode")],directed = F)
        #tmp=igraph::cluster_louvain(sg, weights = 1/(1+exp(-25*(tmp_resNet_10$cor^0.3-0.8))))
        tmp=igraph::cluster_louvain(sg, weights = tmp_resNet_10$cor^3)
        
        tmp_seeds=NULL
        for(icluster in 1:length(tmp)){
          tmp_node=data.frame(Node=tmp[[icluster]])
          tmp_node$cluster=paste0("S_",icluster)
          tmp_seeds=rbind(tmp_seeds,tmp_node)
          rm(tmp_node)
        }
        tmp_seeds$Node=as.character(tmp_seeds$Node)
        toMergeList=NULL
        if(length(unique(tmp_seeds$cluster))>1){
          for(icluster in 1:(length(unique(tmp_seeds$cluster))-1)){
            for(jcluster in (icluster+1):(length(unique(tmp_seeds$cluster)))){
              slNode1=tmp_seeds$Node[tmp_seeds$cluster==unique(tmp_seeds$cluster)[icluster]]
              slNode2=tmp_seeds$Node[tmp_seeds$cluster==unique(tmp_seeds$cluster)[jcluster]]
              sl_cor_mean=input_centroid_cor[slNode1,slNode2]
              sl_cor_mean=mean(as.numeric(sl_cor_mean))
              if(sl_cor_mean>0.5){
                slNet=tmp_resNet_10[which((tmp_resNet_10$Fnode %in% slNode1&tmp_resNet_10$Snode %in% slNode2)|(tmp_resNet_10$Snode %in% slNode1&tmp_resNet_10$Fnode %in% slNode2)),]
                if(nrow(slNet)>0){
                  slNet1=tmp_resNet_10[which(tmp_resNet_10$Fnode %in% slNode1&tmp_resNet_10$Snode %in% slNode1),]
                  slNet2=tmp_resNet_10[which(tmp_resNet_10$Fnode %in% slNode2 & tmp_resNet_10$Snode %in% slNode2),]
                  if(nrow(slNet1)>1){
                    slSd1=(sd(slNet1$cor)+0.05)
                  } else {
                    slSd1=0.05
                  }
                  if(nrow(slNet2)>1){
                    slSd2=(sd(slNet2$cor)+0.05)
                  } else {
                    slSd2=0.05
                  }
                  if((mean(slNet$cor)-mean(slNet1$cor))/slSd1>(-1.65)&(mean(slNet$cor)-mean(slNet2$cor))/slSd2>(-1.65)){
                    toMergeList=rbind(toMergeList,data.frame(Fc=unique(tmp_seeds$cluster)[icluster],Sc=unique(tmp_seeds$cluster)[jcluster],stringsAsFactors = F))
                  }
                }
                
              }
            }
          }
        }
        
        if(!is.null(toMergeList)){
          toMergeList$mergeId=""
          counter=1
          for(i in unique(c(toMergeList$Fc,toMergeList$Sc))){
            tmp=toMergeList[toMergeList$Fc==i|toMergeList$Sc==i,]
            if(sum(toMergeList$mergeId!="")){
              toMergeList$mergeId=unique(toMergeList$mergeId[toMergeList$mergeId!=""])
            } else {
              toMergeList$mergeId=paste0("M_",counter)
              counter=counter+1
            }
          }
          toMergeList=rbind(data.frame(Fnode=toMergeList$Fc,Snode=toMergeList$mergeId,stringsAsFactors = F),data.frame(Fnode=toMergeList$Sc,Snode=toMergeList$mergeId,stringsAsFactors = F))
          toMergeList=toMergeList[!duplicated(.myNetCombination(toMergeList)),]
          for(i in toMergeList$Fnode){
            tmp_seeds$cluster[which(tmp_seeds$cluster==i)]=toMergeList$Snode[toMergeList$Fnode==i]
          }
        }
        
        
        tmp_seeds$cluster=paste0(iclust,"_",tmp_seeds$cluster)
        resClustering=rbind(resClustering,tmp_seeds)
      } else {
        resClustering=rbind(resClustering,data.frame(Node=unique(c(tmp_resNet_10$Fnode,tmp_resNet_10$Snode)),cluster=iclust,stringsAsFactors = F))
      }
      
    }
    colnames(resClustering)=c("centroid","cluster")
    
    
    initial_clusters=resClustering
    resClustering=NULL
    for(iclust in unique(initial_clusters$cluster)){
      slCentroids=as.character(initial_clusters$centroid)[which(initial_clusters$cluster==iclust)]
      tmp_resNet_10=resNet_10[which(resNet_10$Fnode %in% slCentroids&resNet_10$Snode %in% slCentroids),]
      if(nrow(tmp_resNet_10)>4){
        tmp_seeds=tmp_resNet_10
        sg=igraph::graph_from_data_frame(tmp_seeds[,c("Fnode","Snode")],directed = F)
        #tmp=igraph::cluster_louvain(sg, weights = 1/(1+exp(-25*(tmp_resNet_10$cor^0.3-0.8))))
        tmp=igraph::cluster_louvain(sg, weights = tmp_resNet_10$cor^3)
        
        tmp_seeds=NULL
        for(icluster in 1:length(tmp)){
          tmp_node=data.frame(Node=tmp[[icluster]])
          tmp_node$cluster=paste0("S_",icluster)
          tmp_seeds=rbind(tmp_seeds,tmp_node)
          rm(tmp_node)
        }
        tmp_seeds$Node=as.character(tmp_seeds$Node)
        
        toMergeList=NULL
        if(length(unique(tmp_seeds$cluster))>1){
          for(icluster in 1:(length(unique(tmp_seeds$cluster))-1)){
            for(jcluster in (icluster+1):(length(unique(tmp_seeds$cluster)))){
              slNode1=tmp_seeds$Node[tmp_seeds$cluster==unique(tmp_seeds$cluster)[icluster]]
              slNode2=tmp_seeds$Node[tmp_seeds$cluster==unique(tmp_seeds$cluster)[jcluster]]
              sl_cor_mean=input_centroid_cor[slNode1,slNode2]
              sl_cor_mean=mean(as.numeric(sl_cor_mean))
              if(sl_cor_mean>0.5){
                slNet=tmp_resNet_10[which((tmp_resNet_10$Fnode %in% slNode1&tmp_resNet_10$Snode %in% slNode2)|(tmp_resNet_10$Snode %in% slNode1&tmp_resNet_10$Fnode %in% slNode2)),]
                if(nrow(slNet)>0){
                  slNet1=tmp_resNet_10[which(tmp_resNet_10$Fnode %in% slNode1&tmp_resNet_10$Snode %in% slNode1),]
                  slNet2=tmp_resNet_10[which(tmp_resNet_10$Fnode %in% slNode2 & tmp_resNet_10$Snode %in% slNode2),]
                  if(nrow(slNet1)>1){
                    slSd1=(sd(slNet1$cor)+0.05)
                  } else {
                    slSd1=0.05
                  }
                  if(nrow(slNet2)>1){
                    slSd2=(sd(slNet2$cor)+0.05)
                  } else {
                    slSd2=0.05
                  }
                  if((mean(slNet$cor)-mean(slNet1$cor))/slSd1>(-1.65)&(mean(slNet$cor)-mean(slNet2$cor))/slSd2>(-1.65)){
                    toMergeList=rbind(toMergeList,data.frame(Fc=unique(tmp_seeds$cluster)[icluster],Sc=unique(tmp_seeds$cluster)[jcluster],stringsAsFactors = F))
                  }
                }
                
              }
            }
          }
        }
        if(!is.null(toMergeList)){
          toMergeList$mergeId=""
          counter=1
          for(i in unique(c(toMergeList$Fc,toMergeList$Sc))){
            tmp=toMergeList[toMergeList$Fc==i|toMergeList$Sc==i,]
            if(sum(toMergeList$mergeId!="")){
              toMergeList$mergeId=unique(toMergeList$mergeId[toMergeList$mergeId!=""])
            } else {
              toMergeList$mergeId=paste0("M_",counter)
              counter=counter+1
            }
          }
          toMergeList=rbind(data.frame(Fnode=toMergeList$Fc,Snode=toMergeList$mergeId,stringsAsFactors = F),data.frame(Fnode=toMergeList$Sc,Snode=toMergeList$mergeId,stringsAsFactors = F))
          toMergeList=toMergeList[!duplicated(.myNetCombination(toMergeList)),]
          for(i in toMergeList$Fnode){
            tmp_seeds$cluster[which(tmp_seeds$cluster==i)]=toMergeList$Snode[toMergeList$Fnode==i]
          }
        }
        
        
        tmp_seeds$cluster=paste0(iclust,"_",tmp_seeds$cluster)
        resClustering=rbind(resClustering,tmp_seeds)
      } else {
        resClustering=rbind(resClustering,data.frame(Node=unique(c(tmp_resNet_10$Fnode,tmp_resNet_10$Snode)),cluster=iclust,stringsAsFactors = F))
      }
      
    }
    colnames(resClustering)=c("centroid","cluster")
    
    tmpName=resClustering$cluster
    resClustering$cluster=as.numeric(factor(resClustering$cluster,levels=unique(resClustering$cluster)[order(unique(resClustering$cluster))]))
    resClustering$centroid=as.character(resClustering$centroid)
    for(iclust in unique(resClustering$cluster)){
      tmp=resClustering[resClustering$cluster==iclust,]
      tmp=input_centroid_cor[tmp$centroid,tmp$centroid]
    }
    
    #resClustering=.netFindClusters(inputGraph=resAdj, algorithm = 1,resolution = clustering_res,group.singletons = F,modularity.fxn = 1)
    #resClustering=data.frame(centroid=row.names(resClustering),cluster=as.character(resClustering[,1]),stringsAsFactors = F)
    
    #rmlooselyConnected=T
    #while(rmlooselyConnected){
    #  rmlooselyConnected=F
    #  tormlist=apply(resAdj,1,function(x) sum(x>0)/length(x))
    #  if(sum(tormlist<0.01)>0){
    #    rmlooselyConnected=T
    #    resAdj=resAdj[-which(tormlist<0.01),-which(tormlist<0.01)]
    #  }
    #  
    #}
    
    #resClustering=kmeans(input_centroid_cor,cluster_count)
    #resClustering=data.frame(centroid=names(resClustering$cluster),cluster=as.character(resClustering$cluster),stringsAsFactors = F)
    
    #resClustering=data.frame(centroid=names(resClustering),cluster=as.character(resClustering),stringsAsFactors = F)
    
    
    
    resNet=resNet_5
    exCentroids=setdiff(row.names(input_centroid_cor),unique(c(resNet$Fnode,resNet$Snode)))
    if(length(exCentroids)>0){
      resNet=rbind(resNet,data.frame(Fnode=exCentroids,Snode=exCentroids,cor=1,commonMarkers=1,stringsAsFactors = F))
    }
    
    
    clusters_identified=unique(resClustering$cluster)
    
    exCentroids=setdiff(row.names(input_centroid_cor),resClustering$centroid)
    if(length(exCentroids)>0){
      resClustering=rbind(resClustering,data.frame(centroid=exCentroids,cluster=max(as.numeric(resClustering$cluster))+1,stringsAsFactors = F))
    }
    
    
    net = network::network(resNet, directed = FALSE,matrix.type="edgelist")
    
    netVerNames=network::network.vertex.names(net)
    net %v% "cluster" = resClustering$cluster[match(netVerNames,as.character(resClustering$centroid))]
    network::set.edge.attribute(net, "weight", resNet$cor^3)
    
    load(.myFilePathMakerFn("UMAP_centroid",argList=argList))
    centroid_layout=UMAP_centroid[match(as.character(netVerNames),as.character(UMAP_centroid$centroid)),-1]
    centroid_layout=as.matrix(centroid_layout)
    colnames(centroid_layout)=c("x","y")
    
    net=ggnetwork:::fortify.network(net,layout = centroid_layout)
    
    resClustering$color=c(hues::iwanthue(length(clusters_identified)),"#808080")[as.numeric(resClustering$cluster)]
    
    clusterCol=resClustering[!duplicated(resClustering$cluster),]
    
    
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    pd_summary=.scatterPlot_summary2d(object=pd,reductionCols=c("UMAP_1","UMAP_2"))
    
    lda.training.data=merge(resClustering,UMAP_centroid,by="centroid")
    #lda.fit=MASS::lda(cluster~UMAP_1+UMAP_2,data=lda.training.data)
    #lda.pred=predict(lda.fit,pd_summary)
    
    knet=RANN::nn2(query=pd_summary[,c("UMAP_1","UMAP_2")],data = lda.training.data[,c("UMAP_1","UMAP_2")],k=12,eps=0)
    affinity_mat=matrix(0,nrow=nrow(knet$nn.idx),ncol=length(unique(lda.training.data$cluster)))
    colnames(affinity_mat)=unique(lda.training.data$cluster)
    
    tmp_affinities=t(apply(knet$nn.dists,1,function(x) exp((-1)*(x/x[3])^2)))
    tmp_affinities=t(apply(tmp_affinities,1,function(x) x/sum(x)))
    for(itr in 1:ncol(knet$nn.idx)){
      tst=.myOneHotFn(inputVector=factor(lda.training.data$cluster[knet$nn.idx[,itr]],levels=unique(lda.training.data$cluster)))
      tst=tst[,match(colnames(affinity_mat),colnames(tst))]
      tst=tmp_affinities[,itr]*tst
      affinity_mat=affinity_mat+tst
      rm(tst)
    }
    
    col_rgb_pallette=resClustering[!duplicated(resClustering$cluster),]
    col_rgb_pallette=col_rgb_pallette[,c("cluster","color")]
    col_rgb_pallette=col_rgb_pallette[match(colnames(affinity_mat),col_rgb_pallette$cluster),]
    col_rgb1=matrix(grDevices::col2rgb(col_rgb_pallette$color)[1,],nrow=nrow(affinity_mat),ncol=nrow(col_rgb_pallette),byrow = T)
    col_rgb2=matrix(grDevices::col2rgb(col_rgb_pallette$color)[2,],nrow=nrow(affinity_mat),ncol=nrow(col_rgb_pallette),byrow = T)
    col_rgb3=matrix(grDevices::col2rgb(col_rgb_pallette$color)[3,],nrow=nrow(affinity_mat),ncol=nrow(col_rgb_pallette),byrow = T)
    
    colnames(col_rgb1)=col_rgb_pallette$cluster
    colnames(col_rgb2)=col_rgb_pallette$cluster
    colnames(col_rgb3)=col_rgb_pallette$cluster
    
    col_rgb1=affinity_mat*col_rgb1
    col_rgb2=affinity_mat*col_rgb2
    col_rgb3=affinity_mat*col_rgb3
    
    col_rgb1=rowSums(col_rgb1)
    col_rgb2=rowSums(col_rgb2)
    col_rgb3=rowSums(col_rgb3)
    
    pd_summary$color=rgb(red=col_rgb1,green=col_rgb2,blue=col_rgb3,maxColorValue = 255)
    
    
    scale_factor=net[!duplicated(net$vertex.names),]
    scale_factor=merge(scale_factor,UMAP_centroid,by.x="vertex.names",by.y="centroid")
    scale_factor1=lm(x~UMAP_1,data=scale_factor)
    pd_summary$UMAP_1=predict(scale_factor1,newdata=pd_summary)
    scale_factor2=lm(y~UMAP_2,data=scale_factor)
    pd_summary$UMAP_2=predict(scale_factor2,newdata=pd_summary)
    
    pd_summary$xend=pd_summary$UMAP_1
    pd_summary$yend=pd_summary$UMAP_2
    
    #predicting the background color
    net=merge(net,clusterCol,by="cluster")
    p=ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) + geom_point(data=pd_summary,aes(UMAP_1,UMAP_2,color=color),alpha=0.5)+
      #geom_rect(aes(xmin=min(c(pd$UMAP_1,net$x)),xmax=max(c(pd$UMAP_1,net$x)),ymin=min(c(pd$UMAP_1,net$x)),ymax=max(c(pd$UMAP_1,net$x)),color="white"),alpha=0.003,fill="white")+
      geom_edges( color = "grey50",aes(size=weight)) +
      geom_nodelabel(aes(color = color, label = as.character(cluster)),fontface = "bold")+
      theme_blank()+scale_size_continuous(range = c(0.04,0.3))+scale_color_identity()+theme(legend.position = "none")
    ggsave(plot=p,.myFilePathMakerFn("consensusDE_clustering_network",argList=argList,uniformImportant=T,pdf = T,makePlotDir=T),width = 9,height = 9)
    return(list(net=resNet,clustering_res=resClustering))
  }
  
  myCentroidNetConstructorFn_distanceBased=function(input_centroid_cor,input_centroid_data,argList,initial_clustering){
    require(network)
    require(ggnetwork)
    ORmat=matrix(NA,nrow=ncol(input_centroid_data),ncol=ncol(input_centroid_data))
    
    for(icentroid in 1:(ncol(input_centroid_data)-1)){
      for(jcentroid in (icentroid+1):ncol(input_centroid_data)){
        if(sum(input_centroid_data[,icentroid]>2&input_centroid_data[,jcentroid]>2)>0){
          tmp_estimate=fisher.test(table(input_centroid_data[,icentroid]>2,input_centroid_data[,jcentroid]>2))$estimate
          
        } else {tmp_estimate=0}
        ORmat[icentroid,jcentroid]=tmp_estimate
        ORmat[jcentroid,icentroid]=tmp_estimate
        rm(tmp_estimate)
      }
    }
    
    row.names(ORmat)=colnames(input_centroid_data)
    colnames(ORmat)=colnames(input_centroid_data)
    ORmat=ORmat[match(row.names(input_centroid_cor),row.names(ORmat)),match(colnames(input_centroid_cor),colnames(ORmat))]
    
    if(sum(input_centroid_cor)==0){
      input_centroid_cor[which(ORmat>5)]=1
    }
    
    centroids_w_signal=apply(input_centroid_data,2,function(x) sum(x>3))
    centroids_w_signal=names(centroids_w_signal)[which(centroids_w_signal>3)]
    
    if(length(centroids_w_signal)==0){
      centroids_w_signal=apply(input_centroid_data,2,function(x) sum(x>3))
      centroids_w_signal=names(centroids_w_signal)[which(centroids_w_signal>0)]
    }
    
    OR_counts=apply(ORmat,1,function(x) sum(x>2,na.rm = T))
    OR_counts=OR_counts[names(OR_counts) %in% centroids_w_signal]
    cor_counts=apply(input_centroid_cor,1,function(x) sum(x>0.1))
    cor_counts=cor_counts[names(cor_counts) %in% centroids_w_signal]
    
    centers=NULL
    if(is.null(initial_clustering)){
      if(length(cor_counts)>0){
        centers=c()
        if(length(which(OR_counts==0))>0){
          centers=names(OR_counts)[which(OR_counts==0)]
          OR_counts=OR_counts[-which(names(OR_counts) %in% centers)]
          cor_counts=cor_counts[-which(names(cor_counts) %in% centers)]
        }
        
        cor_counts=cor_counts[order(cor_counts,OR_counts,decreasing = F)]
        
        iterate=T
        while(iterate){
          iterate=F
          slCor=1
          slInd=0
          for(i in 1:length(cor_counts)){
            indItr=cor_counts[i]
            
            if(length(centers)>0){
              if(sum(ORmat[names(indItr),centers]>2&input_centroid_cor[names(indItr),centers]>0.5)==0){
                if(sum(input_centroid_cor[names(indItr),centers])<slCor){
                  slCor=sum(input_centroid_cor[names(indItr),centers])
                  slInd=indItr
                  iterate=T
                }
              }
            } else {
              slInd=indItr
              iterate=T
              break;
            }
          }
          
          if(iterate){
            centers=c(centers,names(slInd))
            cor_counts=cor_counts[-which(names(cor_counts)==names(slInd))]
          }
          if(length(cor_counts)==0){
            iterate=F
          }
          
        }
        
        resClustering=data.frame(centroid=centers,cluster=1:length(centers),stringsAsFactors = F)
        
        if(length(setdiff(colnames(input_centroid_cor),centroids_w_signal))>0){
          resClustering=rbind(resClustering,data.frame(centroid=setdiff(colnames(input_centroid_cor),centroids_w_signal),cluster=0,stringsAsFactors = F))
        }
      } else {
        resClustering=data.frame(centroid=colnames(input_centroid_data),cluster=NA,stringsAsFactors = F)
      }
    } else {
      resClustering=data.frame(centroid=initial_clustering$centroid,cluster=initial_clustering$cluster,stringsAsFactors = F)
      resClustering=resClustering[!is.na(resClustering$cluster),]
    }
    
    
    
    if(nrow(resClustering)>1){
      tmpColors=data.frame(cluster=unique(resClustering$cluster),color=hues::iwanthue(length(unique(unique(resClustering$cluster)))),stringsAsFactors = F)
      resClustering=merge(resClustering,tmpColors,by="cluster",all.x=T)
      
    } else {
      resClustering$color=hues::iwanthue(1)
    }
    
    if(length(setdiff(row.names(input_centroid_cor) , resClustering$centroid))>0){
      
      tmp_sl_clusters=unique(resClustering$cluster)
      tmp_sl_clusters=tmp_sl_clusters[!is.na(tmp_sl_clusters)]
      tmp_sl_clusters=tmp_sl_clusters[tmp_sl_clusters!="0"]
      tmp_cluster=input_centroid_cor[-which(colnames(input_centroid_cor) %in% resClustering$centroid),intersect(colnames(input_centroid_cor),resClustering$centroid[which(resClustering$cluster %in% tmp_sl_clusters)])]
      if(class(tmp_cluster)[1]=="numeric"){
        tmp_cluster2=matrix(tmp_cluster,ncol=length(which(resClustering$cluster %in% tmp_sl_clusters)))
        row.names(tmp_cluster2)=row.names(input_centroid_cor)[-which(colnames(input_centroid_cor) %in% resClustering$centroid)]
        colnames(tmp_cluster2)=resClustering$centroid[which(resClustering$cluster %in% tmp_sl_clusters)]
        tmp_cluster=tmp_cluster2
        rm(tmp_cluster2)
      }
      tmp_cluster=apply(tmp_cluster,1,function(x) colnames(tmp_cluster)[which(x>0.85)])
      if(class(tmp_cluster)=="matrix"){
        tmp_cluster=apply(tmp_cluster,2,list)
        for(ilist in 1:length(tmp_cluster)){
          tmp_cluster[[ilist]]=tmp_cluster[[ilist]][[1]]
        }
      }
      if(length(tmp_cluster)>0){
        for(ic in 1:length(tmp_cluster)){
          if(length(tmp_cluster[[ic]])==1){
            resClustering=rbind(resClustering,data.frame(centroid=names(tmp_cluster)[[ic]],cluster=resClustering$cluster[which(resClustering$centroid== unlist(tmp_cluster[[ic]]))],color=resClustering$color[which(resClustering$centroid== unlist(tmp_cluster[[ic]]))],stringsAsFactors = F))
          }
        }
      }
      
      if(!is.null(centers)){
        if(length(which(row.names(input_centroid_cor) %in% resClustering$centroid))<ncol(input_centroid_cor)){
          otherCentroids=input_centroid_cor[-which(row.names(input_centroid_cor) %in% resClustering$centroid[which(resClustering$cluster!="0")]),centers]
          if(class(otherCentroids)=="numeric"){
            otherCentroids=matrix(otherCentroids,ncol = length(centers))
            colnames(otherCentroids)=centers
            row.names(otherCentroids)=row.names(input_centroid_cor)[-which(row.names(input_centroid_cor) %in% resClustering$centroid[which(resClustering$cluster!="0")])]
          }
          otherOR=ORmat[match(row.names(otherCentroids),row.names(ORmat)),match(colnames(otherCentroids) , colnames(ORmat))]
          otherCentroids[which(otherOR<2)]=0
          otherCentroids[otherCentroids<0.5]=0
          
          affinity_mat=pmax(otherCentroids-0.5,0)
          if(class(affinity_mat)=="numeric"){
            affinity_mat=matrix(affinity_mat,ncol = ncol(otherCentroids))
            colnames(affinity_mat)=colnames(otherCentroids)
            row.names(affinity_mat)=row.names(otherCentroids)
          }
          
          if(length(which(rowSums(affinity_mat)==0))>0){
            affinity_mat=affinity_mat[-which(rowSums(affinity_mat)==0),]
            if(class(affinity_mat)=="numeric"){
              affinity_mat2=matrix(affinity_mat,ncol = ncol(otherCentroids))
              colnames(affinity_mat2)=colnames(otherCentroids)
              row.names(affinity_mat2)=names(affinity_mat)
              affinity_mat=affinity_mat2
              rm(affinity_mat2)
            }
          }
          
          
          if(ncol(affinity_mat)>1){
            tmp_zero=apply(affinity_mat,1,sum)
            affinity_mat=t(apply(affinity_mat,1,function(x) x/sum(x)))
            if(sum(tmp_zero==0)>0){
              affinity_mat[which(tmp_zero==0),]=0
            }
          } else {
            affinity_mat[affinity_mat[,1]>0,1]=1
          }
          
          col_rgb_pallette=resClustering
          col_rgb_pallette=col_rgb_pallette[match(colnames(affinity_mat),col_rgb_pallette$centroid),]
          col_rgb1=matrix(grDevices::col2rgb(col_rgb_pallette$color)[1,],nrow=nrow(affinity_mat),ncol=nrow(col_rgb_pallette),byrow = T)
          col_rgb2=matrix(grDevices::col2rgb(col_rgb_pallette$color)[2,],nrow=nrow(affinity_mat),ncol=nrow(col_rgb_pallette),byrow = T)
          col_rgb3=matrix(grDevices::col2rgb(col_rgb_pallette$color)[3,],nrow=nrow(affinity_mat),ncol=nrow(col_rgb_pallette),byrow = T)
          
          colnames(col_rgb1)=col_rgb_pallette$cluster
          colnames(col_rgb2)=col_rgb_pallette$cluster
          colnames(col_rgb3)=col_rgb_pallette$cluster
          
          col_rgb1=affinity_mat*col_rgb1
          col_rgb2=affinity_mat*col_rgb2
          col_rgb3=affinity_mat*col_rgb3
          
          
          col_rgb1=rowSums(col_rgb1)
          col_rgb2=rowSums(col_rgb2)
          col_rgb3=rowSums(col_rgb3)
          
          for(icolor in 1:length(col_rgb1)){
            if(sum(resClustering$centroid==names(col_rgb1)[icolor])>0){
              resClustering$color[which(resClustering$centroid==names(col_rgb1)[icolor])]=rgb(red=col_rgb1[icolor],green=col_rgb2[icolor],blue=col_rgb3[icolor],maxColorValue = 255)
            } else {
              resClustering=rbind(resClustering,data.frame(centroid=names(col_rgb1)[icolor],cluster=NA,color=rgb(red=col_rgb1[icolor],green=col_rgb2[icolor],blue=col_rgb3[icolor],maxColorValue = 255)))
              
            }
          }
          
          if(length(unique(resClustering$cluster[resClustering$centroid %in% colnames(affinity_mat)]))>1){
            cluster_memberships=.myOneHotFn(inputVector=factor(resClustering$cluster[resClustering$centroid %in% colnames(affinity_mat)]))
            row.names(cluster_memberships)=resClustering$centroid[resClustering$centroid %in% colnames(affinity_mat)]
            cluster_memberships=cluster_memberships[match(colnames(affinity_mat),row.names(cluster_memberships)),]
          } else {
            cluster_memberships=matrix(unique(resClustering$cluster[resClustering$centroid %in% colnames(affinity_mat)]),ncol=1)
            row.names(cluster_memberships)=resClustering$centroid[resClustering$centroid %in% colnames(affinity_mat)]
          }
          
          affinity_mat2=affinity_mat %*% as.matrix(cluster_memberships)
          
          
          tmp_cluster=apply(affinity_mat2,1,function(x) colnames(affinity_mat2)[which(x>0.8)])
          if(length(tmp_cluster)>0){
            for(ic in 1:length(tmp_cluster)){
              if(length(tmp_cluster[[ic]])==1){
                resClustering$cluster[which(resClustering$centroid==names(tmp_cluster)[[ic]])]=unlist(tmp_cluster[[ic]])
              }
            }
          }
          
        }
        
      }
      
      
      tmp_cluster=input_centroid_cor[,resClustering$centroid[!is.na(resClustering$cluster)]]
      tmp_cluster=apply(tmp_cluster,1,function(x) colnames(tmp_cluster)[which(x>(0.9*max(x)))])
      if(length(tmp_cluster)>0){
        for(ic in 1:length(tmp_cluster)){
          if(length(tmp_cluster[[ic]])==1){
            resClustering$cluster[which(resClustering$centroid==names(tmp_cluster)[[ic]])]=resClustering$cluster[which(resClustering$centroid== unlist(tmp_cluster[[ic]]))]
          } else {
            tmp_sl_clusters=unique(resClustering$cluster[which(resClustering$centroid %in% unlist(tmp_cluster[[ic]]))])
            tmp_sl_clusters=tmp_sl_clusters[!is.na(tmp_sl_clusters)]
            tmp_sl_clusters=tmp_sl_clusters[tmp_sl_clusters!="0"]
            if(length(tmp_sl_clusters)>1){
              resClustering$cluster[which(resClustering$centroid==names(tmp_cluster)[[ic]])]=NA
            }
            
          }
        }
      }
    }
    
    
    row.names(resClustering)=resClustering$centroid
    #############################
    
    cor_mat=input_centroid_cor
    cor_mat[which(ORmat<2)]=0
    
    resNet=NULL
    for(icentroid in row.names(input_centroid_cor)){
      tmp=as.numeric(cor_mat[icentroid,])
      names(tmp)=colnames(cor_mat)
      tmp=tmp[-which(colnames(cor_mat)==icentroid)]
      tmp=tmp[order(tmp,decreasing = T)]
      tmp=tmp[which(tmp>0.5)]
      if(length(unique(tmp))>4){
        tmp_scale=scale(as.numeric(tmp))
        tmp=tmp[which(tmp_scale>(-1.96))]
      }
      
      if(length(tmp)>0){
        tmp=tmp[1:min(length(tmp),5)]
        resNet=rbind(resNet,data.frame(Fnode=icentroid,Snode=names(tmp),cor=tmp,stringsAsFactors = F))
        
      }
    }
    resNet=resNet[!duplicated(.myNetCombination(resNet)),]
    
    for(ij in 1:nrow(resNet)){
      resNet$commonMarkers[ij]=(sum((input_centroid_data[,resNet$Fnode[ij]]>2.5&input_centroid_data[,resNet$Snode[ij]]>2) | (input_centroid_data[,resNet$Fnode[ij]]>2&input_centroid_data[,resNet$Snode[ij]]>2.5))+0.5)/(sum(input_centroid_data[,resNet$Fnode[ij]]>2.5|input_centroid_data[,resNet$Snode[ij]]>2.5)+0.5)
    }
    resNet=resNet[which(resNet$commonMarkers>0.3),]
    
    
    #add singletons
    excentroids=setdiff(colnames(input_centroid_cor),c(resNet$Fnode,resNet$Snode))
    if(length(excentroids)>0){
      resNet=rbind(resNet,data.frame(Fnode=excentroids,Snode=excentroids,cor=1,commonMarkers=1,stringsAsFactors = F))
    }
    
    net = network::network(resNet, directed = FALSE,matrix.type="edgelist")
    
    netVerNames=network::network.vertex.names(net)
    resClustering$cluster[which(resClustering$cluster=="0")]=NA
    net %v% "cluster" = resClustering$cluster[match(netVerNames,as.character(resClustering$centroid))]
    network::set.edge.attribute(net, "weight", resNet$cor^3)
    
    if(!is.null(argList)){
      if(file.exists(.myFilePathMakerFn("UMAP_centroid",argList=argList))){
        load(.myFilePathMakerFn("UMAP_centroid",argList=argList))
      }
      
    }
    
    centroid_layout=UMAP_centroid[match(as.character(netVerNames),as.character(UMAP_centroid$centroid)),-1]
    centroid_layout=as.matrix(centroid_layout)
    colnames(centroid_layout)=c("x","y")
    
    net=ggnetwork:::fortify.network(net,layout = centroid_layout)
    
    clusterCol=resClustering[!duplicated(resClustering$cluster),]
    
    if(!is.null(argList)){
      if(file.exists(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))){
        load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
      }
      
    }
    
    pd_summary=.scatterPlot_summary2d(object=pd,reductionCols=c("UMAP_1","UMAP_2"))
    
    lda.training.data=merge(resClustering,UMAP_centroid,by="centroid")
    #lda.fit=MASS::lda(cluster~UMAP_1+UMAP_2,data=lda.training.data)
    #lda.pred=predict(lda.fit,pd_summary)
    
    knet=RANN::nn2(query=pd_summary[,c("UMAP_1","UMAP_2")],data = lda.training.data[,c("UMAP_1","UMAP_2")],k=12,eps=0)
    affinity_mat=matrix(0,nrow=nrow(knet$nn.idx),ncol=length(unique(lda.training.data$centroid)))
    colnames(affinity_mat)=unique(lda.training.data$centroid)
    
    tmp_affinities=t(apply(knet$nn.dists,1,function(x) exp((-1)*(x/x[3])^2)))
    tmp_affinities=t(apply(tmp_affinities,1,function(x) x/sum(x)))
    for(itr in 1:ncol(knet$nn.idx)){
      tst=.myOneHotFn(inputVector=factor(lda.training.data$centroid[knet$nn.idx[,itr]],levels=unique(lda.training.data$centroid)))
      tst=tst[,match(colnames(affinity_mat),colnames(tst))]
      tst=as.matrix(tst)
      tst=tmp_affinities[,itr]*tst
      affinity_mat=affinity_mat+tst
      #print(paste(itr,":",affinity_mat[6151,2]))
      rm(tst)
    }
    
    
    col_rgb_pallette=resClustering
    col_rgb_pallette=col_rgb_pallette[match(colnames(affinity_mat),col_rgb_pallette$centroid),]
    col_rgb1=matrix(grDevices::col2rgb(col_rgb_pallette$color)[1,],nrow=nrow(affinity_mat),ncol=nrow(col_rgb_pallette),byrow = T)
    col_rgb2=matrix(grDevices::col2rgb(col_rgb_pallette$color)[2,],nrow=nrow(affinity_mat),ncol=nrow(col_rgb_pallette),byrow = T)
    col_rgb3=matrix(grDevices::col2rgb(col_rgb_pallette$color)[3,],nrow=nrow(affinity_mat),ncol=nrow(col_rgb_pallette),byrow = T)
    
    colnames(col_rgb1)=col_rgb_pallette$cluster
    colnames(col_rgb2)=col_rgb_pallette$cluster
    colnames(col_rgb3)=col_rgb_pallette$cluster
    
    col_rgb1=affinity_mat*col_rgb1
    col_rgb2=affinity_mat*col_rgb2
    col_rgb3=affinity_mat*col_rgb3
    
    col_rgb1=rowSums(col_rgb1)
    col_rgb2=rowSums(col_rgb2)
    col_rgb3=rowSums(col_rgb3)
    
    pd_summary$color=rgb(red=col_rgb1,green=col_rgb2,blue=col_rgb3,maxColorValue = 255)
    
    
    scale_factor=net[!duplicated(net$vertex.names),]
    scale_factor=merge(scale_factor,UMAP_centroid,by.x="vertex.names",by.y="centroid")
    scale_factor1=lm(x~UMAP_1,data=scale_factor)
    pd_summary$UMAP_1=predict(scale_factor1,newdata=pd_summary)
    scale_factor2=lm(y~UMAP_2,data=scale_factor)
    pd_summary$UMAP_2=predict(scale_factor2,newdata=pd_summary)
    
    pd_summary$xend=pd_summary$UMAP_1
    pd_summary$yend=pd_summary$UMAP_2
    
    #predicting the background color
    net=merge(net,clusterCol,by="cluster")
    
    p=ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) + geom_point(data=pd_summary,aes(UMAP_1,UMAP_2,color=color))+
      #geom_rect(aes(xmin=min(c(pd$UMAP_1,net$x)),xmax=max(c(pd$UMAP_1,net$x)),ymin=min(c(pd$UMAP_1,net$x)),ymax=max(c(pd$UMAP_1,net$x)),color="white"),alpha=0.003,fill="white")+
      geom_edges( color = "grey50",aes(size=weight)) +
      geom_nodelabel(aes(color = color, label = as.character(cluster)),fontface = "bold")+
      theme_blank()+scale_size_continuous(range = c(0.04,0.3))+scale_color_identity()+theme(legend.position = "none")
    #ggsave(plot=p,"~/Desktop/tst.png",width = 9,height = 9)
    if(!is.null(argList)){
      ggsave(plot=p,.myFilePathMakerFn("consensusDE_clustering_network",argList=argList,uniformImportant=T,pdf = T,makePlotDir=T),width = 9,height = 9)
    }
    
    return(list(net=resNet,clustering_res=resClustering))
  }
  
  
  centroid_cor=myCorFn(as.matrix(res_arranged))
  #argList=NULL
  centroid_net=myCentroidNetConstructorFn_distanceBased(input_centroid_cor=centroid_cor,input_centroid_data=as.matrix(res_arranged),initial_clustering=initial_clustering,argList = argList)
  colclusters=centroid_net$clustering_res
  cluster_count=length(unique(centroid_net$clustering_res$cluster))
  
  #selecting the sig genes
  sl_gene=res$gene[which(res$zscore>=2.5),]
  #sl_gene=res$gene[which(res$zscore>qnorm(pval_thr,lower.tail = F)&res$effective_size>=max(2,(argList$DE_supportingFractionThr*max(res$effective_size)))&res$zscore_count>=max(2,(argList$DE_supportingFractionThr*max(res$effective_size)))&res$cor<(-0.2))]
  sl_gene=unique(sl_gene)
  cat(paste0("    Number of selected genes: ",length(sl_gene),"\n"))
  
  res=res[res$gene %in% sl_gene,]
  res_arranged=reshape2::dcast(gene~centroid,data = res,value.var = "zscore")
  row.names(res_arranged)=res_arranged[,1]
  res_arranged=res_arranged[,-1]
  res_arranged[is.na(res_arranged)]=0
  
  myGeneEffectFn=function(icluster,res_arranged,colclusters,input_centroid_cor,pval_thr){
    geneEffectSizes=NULL
    
    if(length(which(row.names(input_centroid_cor) %in% colclusters$centroid[colclusters$cluster==icluster]))>0){
      if(sum(is.na(colclusters$cluster))>0){
        input_centroid_cor=input_centroid_cor[-which(row.names(input_centroid_cor) %in% colclusters$centroid[is.na(colclusters$cluster)]),]
        if(class(input_centroid_cor)=="numeric"){
          input_centroid_cor=matrix(input_centroid_cor,ncol=length(input_centroid_cor))
          row.names(input_centroid_cor)=colclusters$centroid[which(colclusters$cluster==icluster)]
        }
      }
      
      if(nrow(input_centroid_cor)>1){
        input_centroid_cor=input_centroid_cor[-which(row.names(input_centroid_cor) %in% colclusters$centroid[which(colclusters$cluster==icluster)]),colnames(input_centroid_cor) %in% colclusters$centroid[colclusters$cluster==icluster]]
        if(class(input_centroid_cor)[1]!="numeric"){
          input_centroid_cor=rowSums(input_centroid_cor)
        }
        
        input_centroid_cor=input_centroid_cor[order(input_centroid_cor,decreasing = T)]
        input_centroid_cor=input_centroid_cor[1:max(length(which(input_centroid_cor>=input_centroid_cor[sum(colclusters$cluster==icluster,na.rm = T)])),min(max(2*sum(colclusters$cluster==icluster,na.rm = T),10),length(input_centroid_cor)))]
        
        for(ig in 1:nrow(res_arranged)){
          bkg_scores=as.numeric(res_arranged[ig,(colnames(res_arranged) %in% names(input_centroid_cor))])
          bkg_scores=pmax(0,bkg_scores)
          
          tmp_effsize=(mean(as.numeric(res_arranged[ig,colnames(res_arranged) %in% colclusters$centroid[which(colclusters$cluster==icluster)]]))-mean(bkg_scores))/(sd(bkg_scores)+0.1)
          tmp_sig_count=sum(as.numeric(res_arranged[ig,colnames(res_arranged) %in% colclusters$centroid[which(colclusters$cluster==icluster)]])>=qnorm(pval_thr,lower.tail = F))
          geneEffectSizes=rbind(geneEffectSizes,data.frame(gene=row.names(res_arranged)[ig],cluster=icluster,effectSize=tmp_effsize,sigCount=tmp_sig_count,stringsAsFactors = F))
        }
      }
      
    }
    
    return(geneEffectSizes)
  }
  
  slClusters=unique(colclusters$cluster[!is.na(colclusters$cluster)])
  slClusters=slClusters[slClusters!="0"]
  geneEffectSizes=parallel::mclapply(slClusters,myGeneEffectFn,res_arranged=res_arranged,colclusters=colclusters,input_centroid_cor=centroid_cor,pval_thr=pval_thr,mc.cores = 3)
  geneEffectSizes=do.call("rbind",geneEffectSizes)
  if(is.null(geneEffectSizes)){
    stop("No significant cell state was identified!")
  }
  geneEffectSizes=geneEffectSizes[geneEffectSizes$sigCount>0,]
  geneEffectSizes=geneEffectSizes[order(geneEffectSizes$effectSize,decreasing = T),]
  
  rowClusters=NULL
  for(i in 1:ncol(res_arranged)){
    if(sum(res_arranged[,i]>=qnorm(pval_thr,lower.tail = F))>0){
      tmpGenes=row.names(res_arranged)[res_arranged[,i]>=qnorm(pval_thr,lower.tail = F)]
      rowClusters=rbind(rowClusters,data.frame(gene=tmpGenes,centroid=colnames(res_arranged)[i]))
    }
  }
  rowClusters$centroid=as.character(rowClusters$centroid)
  rowClusters=merge(rowClusters,colclusters,by="centroid",all.x=T)
  rowClusters=rowClusters[rowClusters$cluster!="0",]
  tmp_clust_count=aggregate(color~cluster+gene,data=rowClusters,length)
  colnames(tmp_clust_count)[3]="gene_count"
  tmp_clust_count=merge(tmp_clust_count,geneEffectSizes,by=c("gene","cluster"),all.x=T)
  tmp_clust_size=as.data.frame(table(colclusters$cluster))
  tmp_clust_size=tmp_clust_size[tmp_clust_size$Freq>0,]
  tmp_clust_size$Var1=as.character(tmp_clust_size$Var1)
  colnames(tmp_clust_size)=c("cluster","size")
  tmp_clust=merge(tmp_clust_count,tmp_clust_size,by="cluster",all=T)
  tmp_clust$fraction=tmp_clust$gene_count/tmp_clust$size
  tmp_clust=tmp_clust[order(tmp_clust$effectSize,decreasing = T),]
  tmp_clust_white=tmp_clust[tmp_clust$fraction>=(1/3),]
  tmp_clust_white=as.character(tmp_clust$gene)[duplicated(as.character(tmp_clust_white$gene))]
  #tmp_clust=tmp_clust[order(tmp_clust$effectSize,decreasing = T),]
  tmp_clust=tmp_clust[which(tmp_clust$fraction>0),]
  rowClusters=tmp_clust[!duplicated(tmp_clust$gene),]
  
  geneEffectSizes=geneEffectSizes[paste0(geneEffectSizes$gene,"_",geneEffectSizes$cluster) %in% paste0(tmp_clust$gene,"_",tmp_clust$cluster),]
  
  clust_color=aggregate(centroid~cluster+color,data=colclusters,length)
  clust_color=clust_color[order(clust_color$centroid,decreasing = T),]
  clust_color=clust_color[!duplicated(clust_color$cluster),]
  clust_color=clust_color[,c("cluster","color")]
  rowClusters=merge(rowClusters,clust_color,by="cluster",all.x=T)
  
  rlab=rowClusters$color[match(row.names(res_arranged),rowClusters$gene)]
  rlab=matrix(rlab,nrow=1)
  rownames(rlab)=c("Marker")
  if(length(tmp_clust_white)>0){
    rlab[1,match(tmp_clust_white,row.names(res_arranged))]="white"
  }
  rowClusters=data.frame(gene=row.names(res_arranged),cluster=rowClusters$cluster[match(row.names(res_arranged),rowClusters$gene)],isBridge=rlab[1,]=="white",stringsAsFactors = F)
  
  clab=matrix(colclusters$color[match(colnames(res_arranged),colclusters$centroid)],ncol=1)
  colnames(clab)=c("Cluster")
  
  tmpColors=data.frame(cluster=c(paste0("C_",1:cluster_count),"Bridges"),colors=c(hues::iwanthue(cluster_count),"white"),stringsAsFactors = F)
  
  #Create heatmap using custom heatmap.3 source code loaded above
  if(nrow(res_arranged)>1){
    png(file=paste0(gsub("\\.rda","",.myFilePathMakerFn("consensusDE_clustering_heatmap",argList=argList,uniformImportant=T,makePlotDir=T),".pdf")))
    main_title="Cell cluster identification"
    par(cex.main=1)
    res_arranged=as.matrix(res_arranged)
    if(length(which(res_arranged>5))>0){
      res_arranged[which(res_arranged>5)]=5
    }
    if(length(which(res_arranged<(-5)))>0){
      res_arranged[which(res_arranged<(-5))]=(-5)
    }
    
    .myheatmap.3(res_arranged, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(5,10),
                 Rowv=TRUE, Colv=TRUE, ColSideColors=clab, RowSideColors=rlab, symbreaks=FALSE, key=TRUE, symkey=FALSE,
                 density.info="none", trace="none", main=main_title, labCol=FALSE, labRow=FALSE, cexRow=1, col=colramp,
                 ColSideColorsSize=1, RowSideColorsSize=1, KeyValueName="Similarity")
    legend("topright",legend=as.character(tmpColors$cluster),
           fill=as.character(tmpColors$colors), border="black", bty="n", y.intersp = 1.4, cex=0.7)
    dev.off()
    
  }
  
  
  load(.myFilePathMakerFn("exp_merged",argList=argList,expData=T))
  fd=as.data.frame(tmp@assays$RNA@meta.features)
  fd=data.frame(gene=row.names(fd),gene_short_name=fd$gene_short_name,stringsAsFactors = F)
  fd$gene_short_name[is.na(fd$gene_short_name)]=fd$gene[is.na(fd$gene_short_name)]
  rowClusters=merge(rowClusters,fd,by='gene',all.x=T)
  save(rowClusters,colclusters,file=.myFilePathMakerFn("consensusDE_clustering_res",argList=argList,uniformImportant=T))
  
  clustercol=colclusters[!duplicated(colclusters$cluster),c("cluster","color")]
  
  geneEffectSizes=merge(geneEffectSizes,rowClusters[,c("gene","gene_short_name","isBridge")],by="gene",all.x=T)
  geneEffectSizes=merge(geneEffectSizes,clustercol,by="cluster")
  
  
  slClusters=unique(geneEffectSizes$cluster)
  slClusters=slClusters[order(as.numeric(slClusters),decreasing = F)]
  pCloudList=list()
  p=""
  if(length(slClusters)>1){
    for(ic in slClusters){
      tmp=geneEffectSizes[geneEffectSizes$cluster==ic,]
      tmp=tmp[order(tmp$effectSize,decreasing = T),]
      tmp=tmp[which(tmp$effectSize>0.5),]
      if(nrow(tmp)>0){
        tmp=tmp[!is.na(tmp$gene_short_name),]
        tmp=tmp[tmp$gene_short_name!="",]
        tmp=tmp[!duplicated(tmp$gene_short_name),]
        tmp=tmp[1:min(60,nrow(tmp)),]
        #tmp=merge(tmp,rowClusters[,c("gene","gene_short_name","isBridge")],by="gene",all.x=T)
        if(sum(is.na(tmp$gene_short_name))>0){
          tmp$gene_short_name[is.na(tmp$gene_short_name)]=tmp$gene
        }
        
        #tmp=merge(tmp,clustercol,by="cluster")
        tmp$color=sample(c("#C45055","#934FBB","#95B1BB","#83C85F","#B9964B"),size = nrow(tmp),replace = T)
        tmp$color[tmp$isBridge]="black"
        if(nrow(tmp)>0){
          tmp=.myggwordcloud(words = tmp$gene_short_name,freq = tmp$effectSize,colors = tmp$color,scale = c(1.3, 1.5*max(0.3,min(0.8,min(tmp$effectSize)/max(tmp$effectSize)))))+ggtitle(paste("Cluster:",ic))+theme_minimal()+theme(
            plot.title = element_text(face="bold.italic",size = 20,hjust = 0.5))+guides(size = guide_legend(title = "Effect size"))
          pCloudList=c(pCloudList,list(tmp))
        }
      }
      
    }
    if(length(pCloudList)>1){
      p=patchwork::wrap_plots(pCloudList,ncol=max(1,cluster_count/4))
    } else {
      p=pCloudList[[1]]
    }
    
    
  } else {
    ic=slClusters
    tmp=geneEffectSizes[geneEffectSizes$cluster==ic,]
      tmp=tmp[order(tmp$effectSize,decreasing = T),]
      tmp=tmp[which(tmp$effectSize>0.5),]
      if(nrow(tmp)>0){
        tmp=tmp[!is.na(tmp$gene_short_name),]
        tmp=tmp[tmp$gene_short_name!="",]
        tmp=tmp[!duplicated(tmp$gene_short_name),]
        tmp=tmp[1:min(60,nrow(tmp)),]
        #tmp=merge(tmp,rowClusters[,c("gene","gene_short_name","isBridge")],by="gene",all.x=T)
        if(sum(is.na(tmp$gene_short_name))>0){
          tmp$gene_short_name[is.na(tmp$gene_short_name)]=tmp$gene
        }
        
        #tmp=merge(tmp,clustercol,by="cluster")
        tmp$color=sample(c("#C45055","#934FBB","#95B1BB","#83C85F","#B9964B"),size = nrow(tmp),replace = T)
        tmp$color[tmp$isBridge]="black"
        
        tmp=.myggwordcloud(words = tmp$gene_short_name,freq = tmp$effectSize,colors = tmp$color,scale = c(1.3, 1.5*max(0.3,min(0.8,min(tmp$effectSize)/max(tmp$effectSize)))))+ggtitle(paste("Cluster:",ic))+theme_minimal()+theme(
          plot.title = element_text(face="bold.italic",size = 20,hjust = 0.5))+guides(size = guide_legend(title = "Effect size"))
        p=tmp
      }
    
    
  }
  
  if(class(p)[1]!=class("")){
    ggsave(p,file=.myFilePathMakerFn("consensusDE_cluster_markers",argList=argList,uniformImportant=T,pdf2 = T,makePlotDir=T),width=24,height=18)
  }
  
  save(geneEffectSizes,file=.myFilePathMakerFn("consensusDE_markers",argList=argList,uniformImportant=T))
  
  
  return("Done!")
}

.myConsensusDE_markerRepresentationFn=function(argList,inputGenes=NULL,prefix=NULL){
  library(dplyr)
  library(purrr)
  library(cowplot)
  library(patchwork)
  
  load(.myFilePathMakerFn("consensusDE_markers",argList=argList,uniformImportant=T))
  load(.myFilePathMakerFn("UMAP_centroid",argList=argList))
  
  cat("           Generating marker figs (cell based) ...\n")
  load(.myFilePathMakerFn("res_DE_wZscore",argList=argList,uniformImportant=T))
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  load(do.call('.myFilePathMakerFn',args=c(list(filename="exp_merged",expData=T),argList)))
    
  tmp=tmp[,colnames(tmp) %in% row.names(pd)]
  pd2=pd[match(colnames(tmp),row.names(pd)),]
    
  
  plotFn_integrated=function(index,geneNameList,expData,annotationFile,input_cell_bkg_umap,prefix){
    
    genes=unlist(geneNameList[[index]]$gene_short_name)
    
    figWidth=49
    if(length(genes)<21){
      figWidth=49/3*ceiling(length(genes)/7)
    }
    
    figHieght=35
    if(length(genes)<7){
      figHieght=35/7*length(genes)
    }
    
    expData=expData[match(c(geneNameList[[index]]$gene,setdiff(row.names(expData),geneNameList[[index]]$gene)),row.names(expData)),]
    p=list()
    for(i in seq(1,length(genes),7)){
      p1=.myFeaturePlot(inputSeurat = expData,inputDimData = annotationFile,inputGenes = genes[i:min((i+6),length(genes))],combine_figs=F)
      p2=.my2dPlot_counts(inputPCA=annotationFile,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=genes[i:min((i+6),length(genes))],geneNameCol="gene_short_name",expData=expData,combine_figs=F)
      
      tmp_UMAP=list()
      for(ij in i:min((i+6),length(genes))){
        tmp=UMAP_centroid
        tmp$gene_id=geneNameList[[index]]$gene[which(geneNameList[[index]]$gene_short_name==genes[ij])]
        tmp$gene=genes[ij]
        tmp_res=res_arranged[which(res_arranged$gene==unique(tmp$gene_id)),]
        tmp_res=tmp_res[match(as.character(tmp$centroid),as.character(tmp_res$centroid)),]
        tmp$zscore=tmp_res$zscore
        tmp$effectiveSize=tmp_res$effective_size
        tmp$size=tmp_res$count
        if(sum(is.na(tmp$zscore))>0){
          tmp$zscore[is.na(tmp$zscore)]=0
          tmp$effectiveSize[is.na(tmp$effectiveSize)]=0
          tmp$size[is.na(tmp$size)]=0
        }
        tmp$effectiveSize=tmp$effectiveSize+1
        tmp_UMAP=c(tmp_UMAP,list(tmp))
      }
      tmp_UMAP=do.call("rbind",tmp_UMAP)
      
      tmp_UMAP$gene=factor(as.character(tmp_UMAP$gene),levels=genes[i:min((i+6),length(genes))])
      
      tmp_UMAP$effectiveSize=as.numeric(as.character(tmp_UMAP$effectiveSize))
      tmp_UMAP$effectiveSize[is.na(tmp_UMAP$effectiveSize)]=min(min(tmp_UMAP$effectiveSize,na.rm = T),1)
      if(sum(tmp_UMAP$zscore>5)>0){
        tmp_UMAP$zscore[which(tmp_UMAP$zscore>5)]=5
      }
      if(sum(tmp_UMAP$zscore<(-5))>0){
        tmp_UMAP$zscore[which(tmp_UMAP$zscore<(-5))]=(-5)
      }
      p3=tmp_UMAP %>% 
        group_split(gene) %>% 
        map(
          ~ggplot(.,aes(UMAP_1,UMAP_2,color=zscore,size=effectiveSize))+geom_point(data=input_cell_bkg_umap,aes(UMAP_1,UMAP_2),color="lightgrey",size=0.1)+geom_point()+theme_cowplot()+ theme(plot.title = element_text(hjust = 0.5))+ggtitle(toupper(unique(.$gene)))+scale_color_gradientn(colors = c("lightblue","lightblue","black","yellow","orange","red"),breaks=c(-5,-2,0,2,3,5),limits=c(-5,5))+scale_size_continuous(limits=c(0,max(tmp_UMAP$effectiveSize)),breaks=seq(0,max(tmp_UMAP$effectiveSize),5))
        ) #%>% 
        #plot_grid(plotlist = ., align = 'hv',ncol=1)
      
      if(length(genes)<(i+6)&length(genes)>7){
        dfNull=data.frame()
        pNull=ggplot(dfNull)+geom_point()
        for(inull in 1:(i+6-length(genes))){
          p1=c(p1,list(pNull))
          p2=c(p2,list(pNull))
          p3=c(p3,list(pNull))
        }
      }
      
      p=c(p,p1,p3,p2)
      
    }
    p=wrap_plots(p,nrow = min(length(genes),7),byrow = F)
    if(is.null(prefix)){
      ggsave(plot=p,file=.myFilePathMakerFn(filename = paste0("marker_",names(geneNameList)[index]),argList = argList,pdf = T,uniformImportant=T,makePlotDir=T),width = figWidth,height = figHieght)
    } else {
      ggsave(plot=p,file=.myFilePathMakerFn(filename = paste0(prefix,"_marker_",names(geneNameList)[index]),argList = argList,pdf = T,uniformImportant=T,makePlotDir=T),width = figWidth,height = figHieght)
    }
    
    return("Done")
  }
    
  if(!is.null(inputGenes)){
    geneEffectSizes=geneEffectSizes[which(toupper(geneEffectSizes$gene) %in% toupper(inputGenes)),]
  }
  slGenes=geneEffectSizes[order(geneEffectSizes$effectSize,decreasing = T),]
  id=rep(0,nrow(slGenes))
  for(ic in unique(slGenes$cluster)){
    tmpId=which(slGenes$cluster==ic)
    tmpId=split(tmpId,ceiling(seq_along(tmpId)/21))
    counter=0
    for(ij in 1:length(tmpId)){
      counter=counter+1
      slGenes$cluster[tmpId[[ij]]]=paste0("cluster",slGenes$cluster[tmpId[[ij]]],"_",counter)
    }
  }
  slGenes=split(slGenes, slGenes$cluster)
    
  input_cell_bkg_umap=.scatterPlot_summary2d(object=pd2,reductionCols=c("UMAP_1","UMAP_2"))
  
  plotMI=parallel::mclapply(1:length(slGenes),plotFn_integrated,geneNameList=slGenes,expData=tmp[unique(c(1:100,which(row.names(tmp) %in% geneEffectSizes$gene))),],annotationFile=pd2,input_cell_bkg_umap=input_cell_bkg_umap,prefix=prefix,mc.cores = argList$ncores)
  save(plotMI,file="final_errors.rda")
    
  
  return("Done!")
}

.myggwordcloud=function (words, freq, colors = NULL, scale = c(4, 0.5), min.freq = 0, max.words = Inf, rot.per = 0.1, 
                          ordered.colors = FALSE, ...) {
  require(ggwordcloud)
  random.order =F
  words_df <- data.frame(word = words, freq = freq)
  if (min.freq > max(freq)) {
    min.freq <- 0
  }
  ord <- rank(-freq, ties.method = "random")
  words_df <- words_df[ord <= max.words, ]
  
  ord <- order(words_df$freq,words_df$word, decreasing = TRUE)
  
  words_df <- words_df[ord, ]
  words_df <- words_df[words_df$freq >= min.freq, ]
  words_df$normedFreq <- words_df$freq/max(words_df$freq)
  if (!is.null(colors)) {
    words_df$color <- colors
  } else {
    words_df$color <- "black"
    }
  words_df$angle <- 90 * (runif(nrow(words_df)) < rot.per)
  with(words_df, ggplot(data = words_df, aes(label = word, 
                                             size = freq, color = color, angle = angle)) + geom_text_wordcloud(rstep = 0.01, 
                                                                                                               tstep = 0.02, rm_outside = F,show.legend = T,seed = 1, ...) + scale_radius(range = 5 * 
                                                                                                                                                                      c(scale[2], scale[1])) + scale_color_identity())
}

.myDEinvestigationFn=function(slTreatment="AD",argList){
  
  slTreatment="AD"
  supportingFractionThr=argList$DE_supportingFractionThr
  n.adaptiveKernel=argList$DE_n.adaptiveKernel
  nPropIter=argList$DE_nPropIter
  
  #load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  pd$anno_condition[which(pd$anno_condition %in% "Alzheimers_disease")]="AD"
  load(.myFilePathMakerFn("harmony-embeddings",argList=argList))
  load(.myFilePathMakerFn("pca_centroids",argList=argList))
  
  if(!all(row.names(harmony_embeddings)==row.names(pd))){
    stop("Error in matching Names!")
  }
  
  pd_pca=pd
  pd_pca$anno_batch=paste0(pd$anno_batch,"_",pd$anno_condition)
  pcaList=split(as.data.frame(harmony_embeddings),pd_pca$anno_batch)
  rm(pd_pca)
  
  load(.myFilePathMakerFn("exp_merged",argList=argList,expData=T))
  expData=SplitObject(tmp, split.by = "anno_batch")
  
  if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename="prevalentGenes"),argList)))){
    cat("           Identifying prevalent genes ...\n")
    prevalenceEstimate=list()
    for(ik in 1:length(expData)){
      tmp=expData[[ik]]@assays$RNA@counts
      tmp=apply(tmp,1,function(x) sum(x>0)/length(x))*100
      prevalenceEstimate=c(prevalenceEstimate,list(data.frame(dsName=names(expData)[ik],gene=row.names(expData[[ik]]),prevalence=tmp,stringsAsFactors = F)))
    }
    prevalenceEstimate=do.call("rbind",prevalenceEstimate)
    prevalenceEstimate=aggregate(prevalence~gene,data=prevalenceEstimate,function(x) sum(x>1/(2*argList$internal_pseudocell_count)))
    prevalenceEstimate=prevalenceEstimate[which(prevalenceEstimate$prevalence>=(supportingFractionThr*length(expData))),]
    prevalentGenes=as.character(prevalenceEstimate$gene)
    
    save(prevalentGenes,file=do.call('.myFilePathMakerFn',args=c(list(filename="prevalentGenes"),argList)))
  } else {
    load(do.call('.myFilePathMakerFn',args=c(list(filename="prevalentGenes"),argList)))
  }
  
  slDSlist=aggregate(anno_condition~anno_batch,data=pd,function(x) length(which(x == "WT"))>0&length(which(x == slTreatment))>0)
  slDSlist=slDSlist[slDSlist$anno_condition==T,]
  slIds=which(unlist(lapply(expData,function(x) x$anno_batch[1])) %in% as.character(slDSlist$anno_batch))
  expData=expData[slIds]
  
  
  dataArranged=list()
  expData_rearranged=list()
  for(i in 1:length(expData)){
    tmpData=expData[[i]]
    if(length(which(tmpData$anno_condition=="Alzheimers_disease"))>0){
      tmpData$anno_condition[which(tmpData$anno_condition=="Alzheimers_disease")]="AD"
    }
    tmpData=tmpData[,tmpData$anno_condition %in% c("WT","AD")]
    tmpData=SplitObject(tmpData, split.by = "anno_condition")
    names(tmpData)=paste0(unlist(lapply(tmpData,function(x) x$anno_batch[1])),"_",names(tmpData))
    for(ij in 1:length(tmpData)){
      tmpData[[ij]]$anno_batch=names(tmpData)[ij]
    }
    expData_rearranged=c(expData_rearranged,tmpData)
    for(ij in 1:length(tmpData)){
      dataArranged=c(dataArranged,list(list(logNormData=tmpData[[ij]]@assays$RNA@data,countData=tmpData[[ij]]@assays$RNA@counts,dsName=names(tmpData)[ij],pcaData=pcaList[[names(tmpData)[ij]]])))
    }
    rm(tmpData)
    rm(ij)
  }
  
  expData=expData_rearranged
  rm(expData_rearranged)
  resGeneMeanSd=.myGeneMeanSdFn(inputExpList=expData,argList=argList)
  
  res=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=nPropIter,mc.cores = argList$ncores)
  
  cat("           Calculating dataset specific z-scores ...\n")
  res2=parallel::mclapply(res,.myConcensusDEFn_step2_detail_exp,centroidPCAdata=pca_centroid,argList=argList,resGeneMeanSd=resGeneMeanSd,mc.cores = min(8,argList$ncores))
  res=res2
  save(res,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0(slTreatment,"_prop_mat_list_step2"),uniformImportant=T),argList)))
  resDistances=data.frame(centroid=names(res[[1]]$w_zscore_centroid),score=res[[1]]$w_zscore_centroid,dataset=res[[1]]$dsName,stringsAsFactors = F)
  if(length(res)>1){
    for(idsScore in 2:length(res)){
      tmp=data.frame(centroid=names(res[[idsScore]]$w_zscore_centroid),score=res[[idsScore]]$w_zscore_centroid,dataset=res[[idsScore]]$dsName,stringsAsFactors = F)
      resDistances=rbind(resDistances,tmp)
    }
  }
  
  resDistances=reshape2::dcast(centroid~dataset,data = resDistances,value.var = "score")
  save(resDistances,file=do.call('.myFilePathMakerFn',args=c(list(filename=paste0(slTreatment,"_centroid_scores"),uniformImportant=T),argList)))
  
  
  
    
    cat("           Meta-analysis of z-scores ...\n")
    myZscoreArrangerFn=function(res,inputTreatment){
      res_count=matrix(0,nrow=nrow(res[[1]]$res),ncol=ncol(res[[1]]$res))
      res_count_w=matrix(0,nrow=nrow(res[[1]]$res),ncol=ncol(res[[1]]$res))
      res_count_sq=matrix(0,nrow=nrow(res[[1]]$res),ncol=ncol(res[[1]]$res))
      res_m=matrix(0,nrow=nrow(res[[1]]$res),ncol=ncol(res[[1]]$res))
      res_indBase=matrix(0,nrow=nrow(res[[1]]$res),ncol=ncol(res[[1]]$res))
      for(i in which(grepl(paste0("_",inputTreatment,"$"),unlist(lapply(res,function(x)x$dsName))))){
        tmp_w_zscore=t(matrix(res[[i]]$w_zscore_centroid,nrow=ncol(res[[i]]$res),ncol=nrow(res[[i]]$res)))
        tmp_w_zscore_sq=t(matrix(res[[i]]$w_zscore_centroid^2,nrow=ncol(res[[i]]$res),ncol=nrow(res[[i]]$res)))
        tmpInd=which(res[[i]]$res!=0)
        if(length(tmpInd)>0){
          res_m[tmpInd]=res_m[tmpInd]+res[[i]]$res[tmpInd]
          res_count[tmpInd]=res_count[tmpInd]+1
          res_count_w[tmpInd]=res_count_w[tmpInd]+tmp_w_zscore[tmpInd]
          res_count_sq[tmpInd]=res_count_sq[tmpInd]+tmp_w_zscore_sq[tmpInd]
          res_indBase[tmpInd]=paste0(res_indBase[tmpInd],",",round(res[[i]]$res[tmpInd],2))
        } else {
          print(paste("check",i))
        }
      }
      res_indBase=gsub("^0,","",res_indBase)
      
      res_effectiveSize=matrix(0,nrow=nrow(res[[1]]$res),ncol=ncol(res[[1]]$res))
      for(i in which(grepl(paste0("_",inputTreatment,"$"),unlist(lapply(res,function(x)x$dsName))))){
        tmp_w_zscore=t(matrix(res[[i]]$w_zscore_centroid,nrow=ncol(res[[i]]$res),ncol=nrow(res[[i]]$res)))
        res_effectiveSize=res_effectiveSize+(tmp_w_zscore/res_count_w)^2
      }
      res_effectiveSize[which(res_count_w==0)]=0
      res_effectiveSize2=1/res_effectiveSize
      res_effectiveSize2[which(res_effectiveSize==0)]=0
      res_effectiveSize=res_effectiveSize2
      rm(res_effectiveSize2)
      
      tmp_res_count_sq=sqrt(res_count_sq)
      tmp_res_count_sq[which(tmp_res_count_sq<1)]=1
      res_m=res_m/tmp_res_count_sq
      res_m[which(res_count==0)]=0
      
      ind=which(res_count!=0,arr.ind = T)
      
      res_arranged=data.frame(gene=row.names(res[[1]]$res)[ind[,1]],centroid=colnames(res[[1]]$res)[ind[,2]],zscore=res_m[ind],count=res_count[ind],effective_size=res_effectiveSize[ind],score_seq=res_indBase[ind],stringsAsFactors = F)
      colnames(res_arranged)[colnames(res_arranged)=="zscore"]=paste0("zscore_",inputTreatment)
      colnames(res_arranged)[colnames(res_arranged)=="count"]=paste0("count_",inputTreatment)
      colnames(res_arranged)[colnames(res_arranged)=="effective_size"]=paste0("effective_size_",inputTreatment)
      colnames(res_arranged)[colnames(res_arranged)=="score_seq"]=paste0("score_seq_",inputTreatment)
      
      #res_arranged=res_arranged[which(res_arranged$effective_size>=max(1,supportingFractionThr*length(expData))),]
      
      return(list(res_arranged=res_arranged,res_count=res_count,res_count_sq=res_count_sq,res_count_w=res_count_w,res_m=res_m,res_indBase=res_indBase,res_effectiveSize=res_effectiveSize))
    }
    
    res_WT=myZscoreArrangerFn(res=res,inputTreatment = "WT")
    res_treatment=myZscoreArrangerFn(res=res,inputTreatment = slTreatment)
    
    res_arranged=merge(res_WT$res_arranged,res_treatment$res_arranged,by=c("gene","centroid"),all=T)
    
    
    fd=as.data.frame(expData[[1]]@assays$RNA@meta.features)
    slCols=c("gene_name","gene_biotype","symbol","gene_short_name","ensembl_gene_id")
    fd=fd[,colnames(fd) %in% slCols]
    
    res_arranged=merge(res_arranged,fd,by.x="gene",by.y="ensembl_gene_id",all.x=T)
    
    
    tst=res_arranged[which(res_arranged$gene_short_name=="Trem2"),]
    tst=tst[order(tst$zscore_WT,decreasing = T),]
    
  load(.myFilePathMakerFn("consensusDE_clustering_res",argList=argList,uniformImportant=T))
  
  myGeneEffectFn=function(igene,colclusters,res_arranged,slTreatment){
    geneEffectSizes=list()
    for(iclust in unique(colclusters$cluster)){
      tmp=res_arranged[which(res_arranged$centroid %in% colclusters$centroid[which(colclusters$cluster==iclust)]&res_arranged$gene==igene),]
      tmp=effsize::cohen.d(tmp$zscore_WT,tmp[,paste0("zscore_",slTreatment)])$estimate
      geneEffectSizes=c(geneEffectSizes,list(data.frame(gene=igene,cluster=iclust,effectSize=tmp)))
    }
    geneEffectSizes=do.call("rbind",geneEffectSizes)
    return(geneEffectSizes)
  }
  
  geneEffectSizes=parallel::mclapply(unique(res_arranged$gene),myGeneEffectFn,colclusters=colclusters,res_arranged=res_arranged,slTreatment=slTreatment,mc.cores = 20)
  
  geneEffectSizes=do.call("rbind",geneEffectSizes)
  geneEffectSizes$gene=as.character(geneEffectSizes$gene)
  
  ggsave(geneEffectSizes,file="geneEffectSizes_DEanalysis.rda")
  
  
  
}

.extraNewRun=function(argList,newRun=F){
  if(newRun){
    if(file.exists(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))){
      file.remove(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    }
    
    if(file.exists(.myFilePathMakerFn("exp_merged",argList=argList,expData=T))){
      file.remove(.myFilePathMakerFn("exp_merged",argList=argList,expData=T))
    }
    
    if(file.exists(.myFilePathMakerFn("varGenes",argList=argList,varGenes =T))){
      file.remove(.myFilePathMakerFn("varGenes",argList=argList,varGenes =T))
    }
    
    if(file.exists(.myFilePathMakerFn("harmony-embeddings",argList=argList))){
      file.remove(.myFilePathMakerFn("harmony-embeddings",argList=argList))
    }
    
    if(file.exists(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))){
      file.remove(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    }
    
    if(file.exists(.myFilePathMakerFn("UMAP_res",argList=argList))){
      file.remove(.myFilePathMakerFn("UMAP_res",argList=argList))
    }
    
    if(file.exists(.myFilePathMakerFn(paste0("nets10-",argList$nPCs),argList=argList))){
      file.remove(.myFilePathMakerFn(paste0("nets10-",argList$nPCs),argList=argList))
    }
    
    if(file.exists(.myFilePathMakerFn(paste0("res_clustering10_",argList$nPCs),argList = argList))){
      file.remove(.myFilePathMakerFn(paste0("res_clustering10_",argList$nPCs),argList = argList))
    }
    
    if(file.exists(do.call('.myFilePathMakerFn',args=c(list(filename="res_prop_harmony"),argList)))){
      file.remove(do.call('.myFilePathMakerFn',args=c(list(filename="res_prop_harmony"),argList)))
    }
  }
  return("Done")
}

.extraCreateDirFn=function(argList){
  if(!dir.exists(argList$saveDir)){
    dir.create(argList$saveDir,recursive = T)
  }
  
  if(!dir.exists(argList$saveDirGlobal)){
    dir.create(argList$saveDirGlobal,recursive = T)
  }
  
  if(!dir.exists(.myFilePathMakerFn("",argList = argList,returnDir = T,uniformImportant = T))){
    dir.create(.myFilePathMakerFn("",argList = argList,returnDir = T,uniformImportant = T),recursive = T)
  }
  return("Done")
}

.pseudo_sim_metric=function(zscoreMat, weightMat, inputMeta,method="ES",geneset_size=100) {
  #method: ES, AUC
  #zscoreMat=res_cor; weightMat=NULL; inputMeta=res_cor;method="ES"

  
  if(is.null(colnames(zscoreMat))){
    colnames(zscoreMat)=1:ncol(zscoreMat)
  }
  
  if(is.null(colnames(inputMeta))){
    colnames(inputMeta)=1:ncol(inputMeta)
  }
  
  if(!is.null(weightMat)){
    inputDataset=zscoreMat*weightMat
    inputDataset[weightMat<0.5]=NA
  } else {
    inputDataset=zscoreMat
  }
  
  
  tmp_slDs=apply(inputDataset,1,function(x)sum(!is.na(x))/length(x))
  inputDataset=inputDataset[tmp_slDs>0.1,]
  
  res_ES=NULL
  for(idataset in 1:nrow(inputDataset)){
    tmp_ES=NULL
    for(i in 1:nrow(inputMeta)){
      x_test=inputMeta[i,]
      x_test=x_test[order(x_test,decreasing = T)]
      
      #x_test=names(x_test)[1:100]
      if(length(x_test)>10){
        #input=tmp["mouse_Van_Hove_dissociation",!is.na(tmp["mouse_Van_Hove_dissociation",])]
        #input=meta$meta_z[188,!is.na(meta$meta_z[188,])]
        input=inputDataset[idataset,!is.na(inputDataset[idataset,])]
        
        gene.list=order(input,decreasing = T)  
        #gene.set=which(names(input) %in% x_test)
        gene.set=match(names(x_test),names(input))
        gene.set=gene.set[!is.na(gene.set)]
        gene.set=gene.set[1:geneset_size]
        if(method=="AUC"){
          tmp_score=.extra_auc_score(geneExp = as.numeric(input[gene.list]),cluster = as.numeric(gene.list %in% gene.set))
        } else if(method=="ES"){
          weighted.score.type = 1
          correl.vector = NULL
          #gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)
          #gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset)
          
          tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0))  # notice that the sign is 0 (no tag) or 1 (tag)
          no.tag.indicator <- 1 - tag.indicator
          N <- length(gene.list)
          Nh <- length(gene.set)
          Nm <- N - Nh
          correl.vector <- rep(1, N)
          
          alpha <- weighted.score.type
          correl.vector <- abs(correl.vector^alpha)
          sum.correl.tag <- sum(correl.vector[tag.indicator == 1])
          norm.tag <- 1/sum.correl.tag
          norm.no.tag <- 1/Nm
          RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
          tmp_score=max(RES)
          #ES <- RES[max(which(tag.indicator>0))]
          
        } else {
          stop("Unknown error!")
        }
        
        tmp_ES=rbind(tmp_ES,data.frame(score=tmp_score,centroid=row.names(inputMeta)[i],dataset=row.names(inputDataset)[idataset],stringsAsFactors = F))
      }
      
    }
    tmp_ES=tmp_ES[order(tmp_ES$score,decreasing = T),]
    tmp_ES$rank=1:nrow(tmp_ES)
    res_ES=rbind(res_ES,tmp_ES)
  }
  
  return(res_ES)
}

.myMetaQCfn=function(input_res_array, input_meta_scores, pseudo_sim_counts=10, sim_method="ES",fig.panel.size=56){
  
  #input_res_array=res_array;input_meta_scores= meta_z;pseudo_sim_counts=10;sim_method="ES";fig.panel.size=56
  
  #data=qread("~/torm_data_metaQC.qs")
  #input_res_array=data$array;input_meta_scores= data$meta;pseudo_sim_counts=10;sim_method="ES";fig.panel.size=56
  
  pp_scores=.pseudo_sim_metric(zscoreMat=input_meta_scores, weightMat=NULL, inputMeta=input_meta_scores,method=sim_method)
  pp_scores_mat=pp_scores
  pp_scores_mat=reshape2::dcast(centroid~dataset,data=pp_scores_mat,value.var = "score")
  row.names(pp_scores_mat)=pp_scores_mat$centroid
  pp_scores_mat=pp_scores_mat[,-1]
  
  #save(pp_scores_mat,file="~/torm_metaQC_pp_mat.rda")
  #load("~/torm_metaQC_pp_mat.rda")
  
  myPseudoSimRunFn=function(ik,input_res_array,input_meta_scores,method,pseudo_sim_counts){
    res_ds=.pseudo_sim_metric(zscoreMat=as.matrix(input_res_array$zscore[,ik,]), weightMat=as.matrix(input_res_array$matWeights[,ik,]), inputMeta=input_meta_scores,method=sim_method)
    res_ds=res_ds[which(res_ds$rank<=pseudo_sim_counts),]
    tmp_consolidated=NULL
    for(i in unique(res_ds$dataset)){
      tmp=res_ds$centroid[which(res_ds$dataset==i)]
      tmp=mean(as.numeric(pp_scores_mat[dimnames(input_res_array$zscore)[[2]][ik],tmp]))
      tmp_consolidated=rbind(tmp_consolidated,data.frame(dataset=i,score=tmp))
    }
    
    tmp_consolidated$scaled_score=as.numeric(scale(tmp_consolidated$score))
    tmp_consolidated$centroid=dimnames(input_res_array$zscore)[[2]][ik]
    return(list(res_ds=res_ds,res_consolidated=tmp_consolidated))
  }
  
  res=parallel::mclapply(1:dim(input_res_array$zscore)[2],myPseudoSimRunFn,input_res_array=input_res_array[c("matWeights","zscore")],input_meta_scores=input_meta_scores,method=method,pseudo_sim_counts=pseudo_sim_counts,mc.cores = 5)
  
  load(.myFilePathMakerFn("UMAP_centroid",argList=.ArgList))
  save(res,pd,UMAP_centroid,file="~/myBucket/torm_resMetaQC_ExN_d2.rda")
  
  res_consolidated=lapply(res,function(x) x$res_consolidated)
  res_consolidated=do.call("rbind",res_consolidated)
  
  
  save(res_consolidated,file=do.call('.myFilePathMakerFn',args=c(list(filename="metaQC",uniformImportant=T,propImportant=T,qsFormat=T),argList)))
  
  
  load(.myFilePathMakerFn("UMAP_centroid",argList=.ArgList))
  load(.myFilePathMakerFn("UMAP_anno",argList=.ArgList,pseudoImportant = F))
  res_consolidated=merge(res_consolidated,UMAP_centroid,by="centroid")
  
  pd_summary=.scatterPlot_summary2d(object=pd,reductionCols=c("UMAP_1","UMAP_2"),n=300)
  pd_summary$scaled_score=0
  dsNames=unique(res_consolidated$dataset)
  res_consolidated$scaled_score[res_consolidated$scaled_score>(-2)]=0
  
  if(sum(colnames(pd_summary)=="dataset")>0){
    pd_summary=pd_summary[,-which(colnames(pd_summary)=="dataset")]
  }
  
  require(cowplot)
  
  page_counter=0
  for(ipage in seq(1,length(dsNames),fig.panel.size)){
    page_counter=page_counter+1
    tmp_data=res_consolidated[res_consolidated$dataset %in% dsNames[ipage:min(ipage+fig.panel.size-1,length(dsNames))],]
    p=ggplot(data=tmp_data,aes(UMAP_1,UMAP_2,fill=scaled_score))+geom_point(data=pd_summary,color='lightgray')+geom_point(data=tmp_data[which(tmp_data$scaled_score==0),],shape=21,size=1)+
      geom_point(data=tmp_data[which(tmp_data$scaled_score!=0),],shape=21,size=1)+
      facet_wrap(~dataset,ncol = 8)+theme_cowplot()+scale_fill_gradient(low="red",high="white")+theme(legend.position = "none",strip.text = element_text(size = 6))
    
    ggsave(plot=p,file=.myFilePathMakerFn(paste0("metaQC_page",page_counter),argList=argList,uniformImportant=T,pdf = T,makePlotDir=T),width = 10,height = 10/7*ceiling((min(ipage+fig.panel.size,length(dsNames))-ipage)/8))
  }
  
  avg_scores=aggregate(score~centroid+UMAP_1+UMAP_2,data=res_consolidated,mean)
  
  pd_summary$score=0.5
  p=ggplot(data=avg_scores,aes(UMAP_1,UMAP_2,fill=score))+geom_point(data=pd_summary,color='lightgray')+geom_point(data=tmp_data[which(tmp_data$scaled_score==0),],shape=21,size=1)+
    theme_cowplot()+scale_fill_gradient2(low="blue",mid="white",high="red",midpoint = 0.6)+theme(legend.position = "none",strip.text = element_text(size = 6))
  ggsave(plot=p,file=.myFilePathMakerFn(paste0("metaQC_overall"),argList=argList,uniformImportant=T,pdf = T,makePlotDir=T),width = 4,height = 4)
  
  return("Done")
}

.myMainRunFn=function(runIndx,input_highly_var_genes=NULL,ncores=6,newRun=F,includeHuman=F,includeMouse=T,external_DE_path=NULL,external_DE_name=NULL,DE_supportingFractionThr=0.1,DE_n.adaptiveKernel=20,DE_nPropIter=1,uniformZscore=F,dist_zscore_gamma=F,dist_zscore_nbinom=F,regularize=F,geoMean=F,internal_pseudocell_count=200,pval_thr=0.001,exclude_non_freq_pseudocells=F){
  
  if(dist_zscore_gamma){
    dist_zscore_gamma=T
    dist_zscore_norm=F
    dist_zscore_nbinom=F
  } else if(dist_zscore_nbinom){
    dist_zscore_gamma=F
    dist_zscore_norm=F
    dist_zscore_nbinom=T
  } else {
    dist_zscore_gamma=F
    dist_zscore_norm=T
    dist_zscore_nbinom=F
  }
  
  
  .ArgList=.myArgFn(runIndx=runIndx,ncores=ncores,includeHuman=includeHuman,includeMouse=includeMouse,uniformZscore=uniformZscore,DE_supportingFractionThr=DE_supportingFractionThr,DE_n.adaptiveKernel=DE_n.adaptiveKernel,DE_nPropIter=DE_nPropIter,dist_zscore_gamma=dist_zscore_gamma,dist_zscore_norm=dist_zscore_norm,dist_zscore_nbinom=dist_zscore_nbinom,regularize=regularize,geoMean=geoMean,prefix=NULL,newRun = newRun,inputDf=NULL,internal_pseudocell_count=internal_pseudocell_count,external_DE_path=external_DE_path,external_DE_name=external_DE_name)
  .ArgList$exclude_non_freq_pseudocells=exclude_non_freq_pseudocells
  .ArgList$input_highly_var_genes=input_highly_var_genes
  
  
  
  rm(list = setdiff(ls(), lsf.str()))
  
  cat('Loading and preprocessing the datasets\n')
  if(!.ArgList$do.liger){
    reRunCheck=tryCatch({load(.myFilePathMakerFn("pca_anno",argList=.ArgList,pseudoImportant = F));F}, error=function(e) {return(T)})
    
    if(!reRunCheck){
      reRunCheck=tryCatch({load(.myFilePathMakerFn("exp_merged",argList=.ArgList,expData=T));F}, error=function(e) {return(T)})
    }
    
    if(reRunCheck){
      
      data=.myReadDataFn(argList=.ArgList)
      
      
      cat('Selecting the highly variable genes\n')
      if(.ArgList$includeHuman){
        data=.myHighVarGeneSlFn(data,dataorganism="Human",argList = .ArgList)
      } else {
        data=.myHighVarGeneSlFn(data,dataorganism="Mouse",argList = .ArgList)
      }
      
      if(sum(colnames(data$data_m@meta.data)=="QC_Gene_unique_count")>0){
        expData=data$data_m@assays$RNA@data
        expData=expData[row.names(expData) %in% data$varFeatures[,1],]
        expData=as.matrix(expData)
        nUMIdata=data$data_m$QC_Gene_unique_count
        
        cordata=NULL
        for(ik in 1:nrow(expData)){
          cordata=rbind(cordata,data.frame(gene=row.names(expData)[ik],pearson=cor(as.numeric(expData[ik,]),as.numeric(nUMIdata)),stringsAsFactors = F))
        }
        
        print(paste("Removed",sum(abs(cordata$pearson)>.ArgList$UMI_cor_thr,na.rm = T),"out of",nrow(cordata),"Genes due to correlation (abs(cor)>",.ArgList$UMI_cor_thr,") with # unique genes"))
        data$varFeatures=data$varFeatures[which(data$varFeatures[,1] %in% as.character(cordata$gene[which(abs(cordata$pearson)<=.ArgList$UMI_cor_thr)])),]
        
      }
      
      tmp=data[c("varFeatures","allGenes" )]
      if(!file.exists(.myFilePathMakerFn("varGenes",argList=.ArgList,varGenes = T))){
        save(tmp,file=.myFilePathMakerFn("varGenes",argList=.ArgList,varGenes =T))
      }
      
      tmp=data$data_m
      if(!file.exists(.myFilePathMakerFn("exp_merged",argList=.ArgList,expData=T))){
        save(tmp,file=.myFilePathMakerFn("exp_merged",argList=.ArgList,expData=T))
      }
      
      cat('Scaling the data and PCA\n')
      .myPCAfn(data,argList = .ArgList,UMI_cor_thr=.ArgList$UMI_cor_thr)
      
    }
    
    
    cat('Harmony analysis\n')
    #rm(list=ls())
    #source("/data/vahid/data/mySC.R")
    reRunCheck=tryCatch({load(.myFilePathMakerFn("harmony-embeddings",argList=.ArgList));F}, error=function(e) {return(T)})
    
    if(reRunCheck){
      library(future)
      plan("multiprocess", workers = 8)
      #plan()
      options(future.globals.maxSize = 1000 * 1024^4)
      load(.myFilePathMakerFn("pca_anno",argList=.ArgList,pseudoImportant = F))
      
      
      harmony_embeddings <- harmony::HarmonyMatrix(pca_res[,1:.ArgList$nPCs], pd, 'anno_batch', do_pca = FALSE, verbose=FALSE)
      .mySaveFn(harmony_embeddings,file=.myFilePathMakerFn("harmony-embeddings",argList=.ArgList,pseudoImportant = F))
      
      slCols=c("anno_age","anno_batch","anno_condition","anno_orig_cellState","anno_region","anno_sex","anno_genotype","anno_RNAseq_type","anno_RNAseq_type2")
      resSummary=NULL
      for(icol in slCols){
        preEffSizes=c()
        tmpInd=which(!is.na(pd[,icol]))
        if(length(unique(pd[tmpInd,icol]))>1){
          for(i in 1:.ArgList$nPCs){
            fit <- aov(pca_res[tmpInd,i] ~ as.factor(as.character(pd[tmpInd,icol])),data = pd)
            fit=sjstats::eta_sq(fit)
            preEffSizes=c(preEffSizes,fit$etasq)
          }
          postEffSizes=c()
          for(i in 1:.ArgList$nPCs){
            fit <- aov(harmony_embeddings[tmpInd,i] ~ as.factor(as.character(pd[tmpInd,icol])),data = pd)
            fit=sjstats::eta_sq(fit)
            postEffSizes=c(postEffSizes,fit$etasq)
          }
          
          tmp1=data.frame(PC=1:.ArgList$nPCs,eta_sq=preEffSizes)
          tmp1$status="pre_Harmony"
          tmp1$variable=icol
          
          tmp2=data.frame(PC=1:.ArgList$nPCs,eta_sq=postEffSizes)
          tmp2$status="post_Harmony"
          tmp2$variable=icol
          resSummary=rbind(resSummary,tmp1)
          resSummary=rbind(resSummary,tmp2)
        }
        
      }
      save(resSummary,file=.myFilePathMakerFn("harmony_stat_summary",argList=.ArgList))
      
      #kmeans clustering
      set.seed(1)
      res_clust=kmeans(harmony_embeddings[,1:.ArgList$nPCs],.ArgList$internal_pseudocell_count,iter.max = 10000,algorithm = "Lloyd")
      
      #row.names(res_clust$centers)=paste0("C_",1:nrow(res_clust$centers))
      #kcenters=res_clust$centers
      #tormList=c()
      #while(length(tormList)<100){
      #  knet=RANN::nn2(data = kcenters, k = 20, eps = 0)
      #  knet_mutual=rowSums(knet$nn.dist)
      #  tormList=c(tormList,row.names(kcenters)[which(knet_mutual==min(knet_mutual))])
      #  kcenters=kcenters[!row.names(kcenters) %in% tormList,]
      #}
      #kcenters=res_clust$centers[!row.names(res_clust$centers) %in% tormList,]
      #knet=RANN::nn2(data=kcenters,query = harmony_embeddings,k=2,eps=0)
      
      #knet2=RANN::nn2(kcenters,k=10,eps=0)
      #knet2=data.frame(FC=row.names(kcenters),SC=row.names(kcenters)[knet2$nn.idx[,2]],distC=knet2$nn.dists[,2],stringsAsFactors = F)
      #knet3=data.frame(sample=row.names(harmony_embeddings),FC=row.names(kcenters)[knet$nn.idx[,1]],distSample=knet$nn.dists[,1],stringsAsFactors = F)
      #knet2=merge(knet2,knet3,by="FC",all=T)
      
      res_clusters=data.frame(sample=names(res_clust$cluster),cluster_id=res_clust$cluster,stringsAsFactors = F)
      save(res_clusters,file=.myFilePathMakerFn("kmeans_res_clusters",argList=.ArgList))
      
      pca_centroid=res_clust$centers
      save(pca_centroid,file=.myFilePathMakerFn("pca_centroids",argList=.ArgList))
      
    }
    
    if(!file.exists(.myFilePathMakerFn("harmony_stat_summary",argList=.ArgList))){
      cat('Harmony post analysis\n')
      load(.myFilePathMakerFn("pca_anno",argList=.ArgList,pseudoImportant = F))
      load(.myFilePathMakerFn("harmony-embeddings",argList=.ArgList))
      slCols=c("anno_age","anno_batch","anno_condition","anno_orig_cellState","anno_region","anno_sex","anno_genotype","anno_RNAseq_type","anno_RNAseq_type2")
      resSummary=NULL
      for(icol in slCols){
        preEffSizes=c()
        tmpInd=which(!is.na(pd[,icol]))
        if(length(unique(pd[tmpInd,icol]))>1){
          for(i in 1:.ArgList$nPCs){
            fit <- aov(pca_res[tmpInd,i] ~ as.factor(as.character(pd[tmpInd,icol])),data = pd)
            fit=sjstats::eta_sq(fit)
            preEffSizes=c(preEffSizes,fit$etasq)
          }
          postEffSizes=c()
          for(i in 1:.ArgList$nPCs){
            fit <- aov(harmony_embeddings[tmpInd,i] ~ as.factor(as.character(pd[tmpInd,icol])),data = pd)
            fit=sjstats::eta_sq(fit)
            postEffSizes=c(postEffSizes,fit$etasq)
          }
          
          tmp1=data.frame(PC=1:.ArgList$nPCs,eta_sq=preEffSizes)
          tmp1$status="pre_Harmony"
          tmp1$variable=icol
          
          tmp2=data.frame(PC=1:.ArgList$nPCs,eta_sq=postEffSizes)
          tmp2$status="post_Harmony"
          tmp2$variable=icol
          resSummary=rbind(resSummary,tmp1)
          resSummary=rbind(resSummary,tmp2)
        }
        
      }
      save(resSummary,file=.myFilePathMakerFn("harmony_stat_summary",argList=.ArgList))
    }
    
    cat("Running UMAP\n")
    #rm(list = setdiff(ls(), lsf.str()))
    #load(.myFilePathMakerFn("harmony-embeddings",argList=.ArgList))
    reRunCheck=F
    #reRunCheck=tryCatch({load(.myFilePathMakerFn("UMAP_anno",argList=.ArgList,pseudoImportant = F));F}, error=function(e) {return(T)})
    
    if(reRunCheck){
      load(.myFilePathMakerFn("harmony-embeddings",argList=.ArgList))
      load(.myFilePathMakerFn("pca_anno",argList=.ArgList,pseudoImportant = F))
      
      load(.myFilePathMakerFn("pca_centroids",argList=.ArgList))
      
      pca_centroid=pca_centroid[,1:.ArgList$nPCs]
      
      
      
      tst=.reductionUMAPFn(harmony_embeddings,umap.method='umap-learn',testPCAembeddings=pca_centroid)
      resUMAP=tst$embedding
      .mySaveFn(resUMAP,file=.myFilePathMakerFn("UMAP_res",argList=.ArgList))
      
      x=as.data.frame(resUMAP)
      x=cbind(x,pd)
      row.names(x)=row.names(pd)
      pd=x
      .mySaveFn(pd,file=.myFilePathMakerFn("UMAP_anno",argList=.ArgList,pseudoImportant = F))
      
      tst=tst[["test"]]
      UMAP_centroid=data.frame(centroid=row.names(pca_centroid),UMAP_1=tst[,1],UMAP_2=tst[,2],stringsAsFactors = F)
      .mySaveFn(UMAP_centroid,file=.myFilePathMakerFn("UMAP_centroid",argList=.ArgList))
      
    }
  } else {
    
    reRunCheck=tryCatch({load(.myFilePathMakerFn("harmony-embeddings",argList=.ArgList));F}, error=function(e) {return(T)})
    
    if(!reRunCheck){
      reRunCheck=tryCatch({load(.myFilePathMakerFn("exp_merged",argList=.ArgList,expData=T));F}, error=function(e) {return(T)})
    }
    
    if(reRunCheck){
      data=.myReadDataFn(argList=.ArgList)
      
      data_m=.extraExport2SeuratFn(data$data_m)
      
      #################
      #load("~/Mouse_MG_slCells.rda")
      #data_m=data_m[,which(colnames(data_m) %in% row.names(slCells) |grepl("^human_",data_m$anno_batch))]
      ################
      
      data.list=SplitObject(data_m,"anno_batch")
      
      dataList=lapply(data.list,function(x){
        return(x@assays$RNA@counts)
      })
      
      names(dataList)=names(data.list)
      
      data_liger_tmp <- createLiger(raw.data = dataList)
      data_liger_tmp <- liger::normalize(data_liger_tmp)
      data_liger_tmp <- selectGenes(data_liger_tmp,num.genes = 1500,do.plot = F)
      
      data_liger=NULL
      if(!is.null(.ArgList$external_DE_path)){
        data_liger <- createLiger(raw.data = dataList,take.gene.union =T)
        data_liger <- liger::normalize(data_liger)
        data_liger <- selectGenes(data_liger,num.genes = 1500,do.plot = F)
        load(.ArgList$external_DE_path)
        #tst=lapply(data_liger@raw.data,function(x) intersect(row.names(x),external_de_genes))
        data_liger@var.genes=unique(c(data_liger_tmp@var.genes,external_de_genes))
      } else {
        data_liger=data_liger_tmp
      }
      
      if(sum(colnames(data$data_m@meta.data)=="QC_Gene_unique_count")>0){
        expData=data$data_m@assays$RNA@data
        expData=expData[row.names(expData) %in% data$varFeatures[,1],]
        expData=as.matrix(expData)
        nUMIdata=data$data_m$QC_Gene_unique_count
        
        cordata=NULL
        for(ik in 1:nrow(expData)){
          cordata=rbind(cordata,data.frame(gene=row.names(expData)[ik],pearson=cor(as.numeric(expData[ik,]),as.numeric(nUMIdata)),stringsAsFactors = F))
        }
        print(paste("Removed",sum(abs(cordata$pearson)>.ArgList$UMI_cor_thr,na.rm = T),"out of",nrow(cordata),"Genes due to correlation (abs(cor)>",.ArgList$UMI_cor_thr,") with # unique genes"))
        if(sum(abs(cordata$pearson)>.ArgList$UMI_cor_thr)>0){
          data$varFeatures=data$varFeatures[which(data$varFeatures[,1] %in% as.character(cordata$gene[which(abs(cordata$pearson)<=.ArgList$UMI_cor_thr)])),]
          
          data_liger@var.genes=varFeatures
        }
        
      }
      
      
      
      tmp=data_liger@var.genes
      if(!file.exists(.myFilePathMakerFn("varGenes",argList=.ArgList,varGenes = T))){
        save(tmp,file=.myFilePathMakerFn("varGenes",argList=.ArgList,varGenes =T))
      }
      
      tmp=data_m
      if(!file.exists(.myFilePathMakerFn("exp_merged",argList=.ArgList,expData=T))){
        save(tmp,file=.myFilePathMakerFn("exp_merged",argList=.ArgList,expData=T))
      }
      
      cat('Liger analysis\n')
      data_liger <- scaleNotCenter(data_liger,remove.missing=F)
      data_liger <- optimizeALS(data_liger,k = .ArgList$nPCs)
      .data_liger=data_liger
      data_liger <- quantile_norm(data_liger)
      
      harmony_embeddings = data_liger@H.norm
      .mySaveFn(harmony_embeddings,file=.myFilePathMakerFn("harmony-embeddings",argList=.ArgList,pseudoImportant = F))
      
      #kmeans clustering
      set.seed(1)
      res_clust=kmeans(harmony_embeddings[,1:.ArgList$nPCs],.ArgList$internal_pseudocell_count,iter.max = 10000,algorithm = "Lloyd")
      
      res_clusters=data.frame(sample=names(res_clust$cluster),cluster_id=res_clust$cluster,stringsAsFactors = F)
      save(res_clusters,file=.myFilePathMakerFn("kmeans_res_clusters",argList=.ArgList))
      
      pca_centroid=res_clust$centers
      save(pca_centroid,file=.myFilePathMakerFn("pca_centroids",argList=.ArgList))
      
      pd=as.data.frame(data_m@meta.data)
      save(pd,file=.myFilePathMakerFn("pca_anno",argList=.ArgList,pseudoImportant = F))
    }
    
  }
  
  cat("Running UMAP\n")
  #rm(list = setdiff(ls(), lsf.str()))
  #load(.myFilePathMakerFn("harmony-embeddings",argList=.ArgList))
  reRunCheck=tryCatch({load(.myFilePathMakerFn("UMAP_anno",argList=.ArgList,pseudoImportant = F));F}, error=function(e) {return(T)})
  
  if(reRunCheck){
    load(.myFilePathMakerFn("harmony-embeddings",argList=.ArgList))
    load(.myFilePathMakerFn("pca_anno",argList=.ArgList,pseudoImportant = F))
    
    load(.myFilePathMakerFn("pca_centroids",argList=.ArgList))
    
    pca_centroid=pca_centroid[,1:.ArgList$nPCs]
    
    tst=.reductionUMAPFn(harmony_embeddings,umap.method='umap-learn',testPCAembeddings=pca_centroid)
    resUMAP=tst$embedding
    .mySaveFn(resUMAP,file=.myFilePathMakerFn("UMAP_res",argList=.ArgList))
    
    x=as.data.frame(resUMAP)
    x=cbind(x,pd)
    row.names(x)=row.names(pd)
    pd=x
    .mySaveFn(pd,file=.myFilePathMakerFn("UMAP_anno",argList=.ArgList,pseudoImportant = F))
    
    tst=tst[["test"]]
    UMAP_centroid=data.frame(centroid=row.names(pca_centroid),UMAP_1=tst[,1],UMAP_2=tst[,2],stringsAsFactors = F)
    .mySaveFn(UMAP_centroid,file=.myFilePathMakerFn("UMAP_centroid",argList=.ArgList))
    
  }
  
  
  cat("Saving the plots\n")
  
  if(file.exists(.myFilePathMakerFn("UMAP_anno",argList=.ArgList,pseudoImportant = F))&!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename="markers_dsCount_Tushar",pdf=T),.ArgList)))){
    
    if(file.exists(.myFilePathMakerFn("harmony_stat_summary",argList=.ArgList))){
      load(.myFilePathMakerFn("harmony_stat_summary",argList=.ArgList))
      p=ggplot(resSummary,aes(variable,eta_sq,color=status))+geom_boxplot()+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=1,color = "black"))+xlab("")+ylab("variance explained")
      ggsave(plot = p,do.call('.myFilePathMakerFn',args=c(list(filename="plot_Harmony_efficiency",pdf=T),.ArgList)),width = 8,height = 6)
      
    }
    
    
    load(do.call('.myFilePathMakerFn',args=c(list(filename="UMAP_anno",pseudoImportant = F),.ArgList)))
    
    if(sum(grepl("human_Tushar",pd$anno_batch))>0){
      p=.myDimPlotFn(object=pd[which(grepl("human_Tushar",pd$anno_batch)),], dimCols = c("UMAP_1", "UMAP_2"),attCol = "clusters")
      ggsave(plot=p,.myFilePathMakerFn("clustering_Tushar_v1",argList=.ArgList,pdf=T),width = 15,height = 10)
      p=.my2dPlot2(inputPCA = pd,batchCol = "clusters",reductionCols = c('UMAP_1','UMAP_2'))
      ggsave(plot = p,do.call('.myFilePathMakerFn',args=c(list(filename="clustering_Tushar_v2",pdf=T),.ArgList)),width = 12,height = 8)
      
      col_status=rep("other",nrow(pd))
      col_status[which(grepl("human_Tushar",pd$anno_batch))]="Human_Tushar"
      tmp2=pd
      tmp2$col_status=col_status
      p=.myDimPlotFn(object=tmp2, dimCols = c("UMAP_1", "UMAP_2"),attCol = "col_status")
      ggsave(plot=p,.myFilePathMakerFn("clustering_Tushar_v3",argList=.ArgList,pdf=T),width = 15,height = 10)
      p=.my2dPlot2(inputPCA = tmp2,batchCol = "col_status",reductionCols = c('UMAP_1','UMAP_2'))
      ggsave(plot = p,do.call('.myFilePathMakerFn',args=c(list(filename="clustering_Tushar_v4",pdf=T),.ArgList)),width = 12,height = 8)
      
    }
    
    pd$anno_batch[is.na(pd$anno_batch)]=as.character(pd$anno_batch)[is.na(pd$anno_batch)]
    pd$anno_batch[which(grepl("Zywitza_30485812",pd$anno_batch))]="mouse_Zywitza_30485812"
    pd$anno_batch[which(grepl("mouse_Lau_scRNA_32320664",pd$anno_batch))]="mouse_Lau_32320664"
    pd$anno_batch[which(grepl("mouse_Sierksma_31951107",pd$anno_batch))]="mouse_Sierksma_31951107"
    pd$anno_batch[which(grepl("Sousa_30206190",pd$anno_batch))]="mouse_Sousa_30206190"
    pd$anno_batch[which(grepl("mouse_SalaFrigerio",pd$anno_batch))]="mouse_SalaFrigerio_31018141"
    pd$anno_batch[which(grepl("mouse_Masuda",pd$anno_batch))]="mouse_Masuda"
    pd$anno_batch[which(grepl("Chen_28355573",pd$anno_batch))]="mouse_Chen_28355573"
    pd$anno_batch[which(grepl("Chen_30773317",pd$anno_batch))]="mouse_Chen_30773317"
    save(pd,file=do.call('.myFilePathMakerFn',args=c(list(filename="UMAP_anno",pseudoImportant = F),.ArgList)))
                    
    
    p=.my2dPlot2(inputPCA = pd,batchCol = "anno_batch",reductionCols = c('UMAP_1','UMAP_2'))
    ggsave(plot = p,do.call('.myFilePathMakerFn',args=c(list(filename="batchEffect",pdf=T),.ArgList)),width = 36,height = 26)
  
    
    p=ggplot(pd,aes(UMAP_1,UMAP_2,color=anno_batch))+geom_point(alpha=0.5,size=0.5)
    ggsave(plot = p,do.call('.myFilePathMakerFn',args=c(list(filename="plot_batch_v1",pdf=T),.ArgList)),width = 12,height = 8)
    
    
    p=.my2dPlot2(inputPCA = pd,batchCol = "anno_batch",reductionCols = c('UMAP_1','UMAP_2'))
    ggsave(plot = p,do.call('.myFilePathMakerFn',args=c(list(filename="plot_batch_v2",pdf=T),.ArgList)),width = 12,height = 8)
    
    
    tmp=pd[!is.na(pd$anno_orig_cellState),]
    tmp=tmp[tmp$anno_orig_cellState!="",]
    tmp$anno_orig_cellState[grepl("Micro",tmp$anno_orig_cellState)]="Microglia"
    tmp$anno_orig_cellState[tmp$anno_orig_cellState %in% c("MG","MGL","MGL1","MGL2","MGL3","mg","Mic")]="Microglia"
    tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("H1/2M","H1M","H2M","homeostatic mic_s2_c1","homeostatic mic_s2_c2"))]="H1/2M"
    if(F){
      tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("cDC1","cDC2","pDC","migDC"))]="cDC"
      tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("T/NKT cells","T cells","NK cells"))]="T/NKT cells"
      tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("Mnc","Non-cl. monocytes","MCs"))]="Monocyte/Mdc"
      tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("ILC"))]="innate_lymphoid"
      tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("CAM"))]="CNS_assoc_macrophages"
      tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("BAM"))]="BAM"
      
    }
    tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("CPM","CRM"))]="cycling_prolif_microglia"
    tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("IRM"))]="IFN_resp_microglia"
    tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("Proliferation"))]="dissociation_prolif"
    #tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("TRM"))]="transiting_resp_microglia"
    p1=ggplot(tmp,aes(UMAP_1,UMAP_2,color=anno_orig_cellState))+geom_point(alpha=0.5,size=0.5)+scale_color_manual(values=hues::iwanthue(length(unique(as.character(tmp$anno_orig_cellState)))))
    p2=ggplot(tmp,aes(UMAP_1,UMAP_2,color=anno_batch))+geom_point(alpha=0.5,size=0.5)
    p=p1+p2
    ggsave(plot=p,.myFilePathMakerFn("plot_cellStates_v1.1",argList=.ArgList,pdf=T),width = 14,height = 8)
    
    p=.my2dPlot2(inputPCA = tmp,batchCol = "anno_orig_cellState",reductionCols = c('UMAP_1','UMAP_2'))
    ggsave(plot=p,.myFilePathMakerFn("plot_cellStates_v2.1",argList=.ArgList,pdf=T),width = 12,height = 8)
    
    
    tmp=tmp[tmp$anno_orig_cellState!="Microglia",]
    if(nrow(tmp)>0){
      p1=ggplot(tmp,aes(UMAP_1,UMAP_2,color=anno_orig_cellState))+geom_point(alpha=0.5,size=0.5)+scale_color_manual(values=hues::iwanthue(length(unique(as.character(tmp$anno_orig_cellState)))))
      p2=ggplot(tmp,aes(UMAP_1,UMAP_2,color=anno_batch))+geom_point(alpha=0.5,size=0.5)
      p=p1+p2
      ggsave(plot=p,.myFilePathMakerFn("plot_cellStates_v1",argList=.ArgList,pdf=T),width = 14,height = 8)
      
      if(length(which(tmp$anno_orig_cellState %in% c("ARM","IFN_resp_microglia","cycling_prolif_microglia","BAM","DAM","dissociation_prolif","CNS_assoc_macrophages")))>0){
        p=ggplot(tmp[which(tmp$anno_orig_cellState %in% c("ARM","IFN_resp_microglia","cycling_prolif_microglia","BAM","DAM","dissociation_prolif","CNS_assoc_macrophages")),],aes(UMAP_1,UMAP_2,color=anno_orig_cellState))+geom_point(alpha=0.5,size=0.5)+scale_color_manual(values=hues::iwanthue(length(unique(as.character(tmp$anno_orig_cellState)))))
        ggsave(plot=p,.myFilePathMakerFn("plot_cellStates_v3",argList=.ArgList,pdf=T),width = 14,height = 8)
      }
      
      p=.my2dPlot2(inputPCA = tmp,batchCol = "anno_orig_cellState",reductionCols = c('UMAP_1','UMAP_2'))
      ggsave(plot=p,.myFilePathMakerFn("plot_cellStates_v2",argList=.ArgList,pdf=T),width = 12,height = 8)
      
    }
    
    
    
    tmp=pd
    if(length(which(tmp$anno_orig_cellState %in% .ArgList$exCellStates))>0){
      tmp=tmp[-which(tmp$anno_orig_cellState %in% .ArgList$exCellStates),]
    }
    
    if(nrow(tmp)>0){
      tmp=tmp[!grepl("mouse_Sierksma_31951107",tmp$anno_batch),]
      p1=.myDimPlotFn(object=tmp,dimCols = c("UMAP_1","UMAP_2"),attCol = 'anno_RNAseq_type',label = F)
      p2=.my2dPlot2(inputPCA = tmp,batchCol = "anno_RNAseq_type",reductionCols = c('UMAP_1','UMAP_2'))
      p=p1+p2
      ggsave(plot=p,.myFilePathMakerFn("plot_RNAseqType_v1",argList=.ArgList,pdf=T),width = 14,height = 8)
    }
    
    
    tmp=pd[which(pd$anno_batch!="mouse_Wu_29024657"),]
    tmp=tmp[which(!grepl("mouse_Sierksma_31951107",tmp$anno_batch)),]
    if(length(which(tmp$anno_orig_cellState %in% .ArgList$exCellStates))>0){
      tmp=tmp[-which(tmp$anno_orig_cellState %in% .ArgList$exCellStates),]
    }
    
    tmp$anno_RNAseq_type2=as.character(tmp$anno_RNAseq_type2)
    tmp$anno_RNAseq_type2[which(tmp$anno_RNAseq_type2=="Smart-Seq 2")]="Smart-Seq"
    tmp$anno_RNAseq_type2=factor(as.character(tmp$anno_RNAseq_type2),levels=rev(c("Mars-Seq","Smart-Seq","10X Genomics","Drop-Seq")))
    p1=.myDimPlotFn(object=tmp,dimCols = c("UMAP_1","UMAP_2"),attCol = 'anno_RNAseq_type2',label = F)
    p2=.my2dPlot2(inputPCA = tmp,batchCol = "anno_RNAseq_type2",reductionCols = c('UMAP_1','UMAP_2'))
    p=p1+p2
    ggsave(plot=p,.myFilePathMakerFn("plot_RNAseqType_v2",argList=.ArgList,pdf=T),width = 14,height = 8)
    
    
    if(length(which(pd$anno_orig_cellState %in% .ArgList$exCellStates))>0){
      tmp=pd[-which(pd$anno_orig_cellState %in% .ArgList$exCellStates),]
    } else {
      tmp=pd
    }
    
   
    if(length(which(tmp$anno_orig_cellState %in% .ArgList$exCellStates))>0){
      tmp=tmp[-which(tmp$anno_orig_cellState %in% .ArgList$exCellStates),]
    }
    
    p1=.my2dPlot_continuous(inputPCA=tmp,continuous_batch_values="QC_Gene_unique_count",reductionCols=c("UMAP_1","UMAP_2"))+ggtitle("Gene_unique_count")
    p2=.my2dPlot_continuous(inputPCA=tmp,continuous_batch_values="QC_Gene_total_count",reductionCols=c("UMAP_1","UMAP_2"))+ggtitle("Gene_total_count")
    p3=.my2dPlot_continuous(inputPCA=tmp,continuous_batch_values="QC_IEG.pct",reductionCols=c("UMAP_1","UMAP_2"))+ggtitle("QC_IEG.pct")
    p4=.my2dPlot_continuous(inputPCA=tmp,continuous_batch_values="QC_top50_pct",reductionCols=c("UMAP_1","UMAP_2"))+ggtitle("QC_top50_pct")
    p=p1+p2
    ggsave(plot=p,.myFilePathMakerFn("QC_plot_v1",argList=.ArgList,pdf=T),width = 14,height = 8)
    
    p=p3+p4
    ggsave(plot=p,.myFilePathMakerFn("QC_plot_v2",argList=.ArgList,pdf=T),width = 14,height = 8)
    
    tmp=pd[!is.na(pd$anno_condition),]
    tmp=tmp[tmp$anno_condition!="",]
    if(length(which(tmp$anno_condition=="Alzheimers_disease"))>0){
      tmp$anno_condition[which(tmp$anno_condition=="Alzheimers_disease")]="AD"
    }
   
    tmp=tmp[!(tmp$anno_orig_cellState %in% .ArgList$exCellStates),]
    if(nrow(tmp)>0){
      p1=ggplot(tmp,aes(UMAP_1,UMAP_2,color=anno_condition))+geom_point(alpha=0.5,size=0.5)
      p2=ggplot(tmp,aes(UMAP_1,UMAP_2,color=anno_batch))+geom_point(alpha=0.5,size=0.5)
      p=p1+p2
      ggsave(plot=p,.myFilePathMakerFn("plot_treatment_v1",argList=.ArgList,pdf=T),width = 14,height = 8)
      
      p=.my2dPlot2(inputPCA = tmp,batchCol = "anno_condition",reductionCols = c('UMAP_1','UMAP_2'))
      ggsave(plot=p,.myFilePathMakerFn("plot_treatment_v2",argList=.ArgList,pdf=T),width = 12,height = 8)
    }
    
    tmp=pd[!is.na(pd$anno_age),]
    tmp=tmp[which(tmp$anno_condition=="WT"),]
    tmp=tmp[!(tmp$anno_orig_cellState %in% .ArgList$exCellStates),]
    if(nrow(tmp)>0){
      tmp$anno_age[which(tmp$anno_age=="1wk")]="0d-1wk"
      tmp$anno_age=factor(as.character(tmp$anno_age),levels=c("prenatal",'0d-1wk','1wk-1mo','1-3mo','3-6mo','6-9mo','9-12mo','12-24mo'))
      p1=ggplot(tmp,aes(UMAP_1,UMAP_2,color=anno_age))+geom_point(alpha=0.5,size=0.5)+scale_color_manual(values=hues::iwanthue(length(unique(as.character(tmp$anno_age)))))+theme_bw()
      p2=ggplot(tmp,aes(UMAP_1,UMAP_2,color=anno_batch))+geom_point(alpha=0.5,size=0.5)
      p=p1+p2
      ggsave(plot=p,.myFilePathMakerFn("plot_age_v1",argList=.ArgList,pdf=T),width = 14,height = 8)
      
      p1=ggplot(tmp[tmp$anno_age %in% c("prenatal","0d-1wk"),],aes(UMAP_1,UMAP_2,color=anno_age))+geom_point(alpha=0.5,size=0.5)+scale_color_manual(values=hues::iwanthue(length(unique(as.character(tmp$anno_age)))))+theme_bw()
      p2=ggplot(tmp[tmp$anno_age %in% c("prenatal","0d-1wk"),],aes(UMAP_1,UMAP_2,color=anno_batch))+geom_point(alpha=0.5,size=0.5)
      p=p1+p2
      ggsave(plot=p,.myFilePathMakerFn("plot_age_v3",argList=.ArgList,pdf=T),width = 14,height = 8)
      
      if(.ArgList$includeMouse){
        p=.my2dPlot2(inputPCA = tmp,batchCol = "anno_age",reductionCols = c('UMAP_1','UMAP_2'))
        ggsave(plot=p,.myFilePathMakerFn("plot_age_v2",argList=.ArgList,pdf=T),width = 12,height = 8)
        
      }
    }
    
    
    rm(list = setdiff(ls(), lsf.str()))
    #load(.myFilePathMakerFn(paste0("clustering_benchmark_10_",.ArgList$nPCs),argList=.ArgList))
  
  #resAnalysis$cluster_resolution=gsub("res\\.","",resAnalysis$cluster_resolution)
    slResolution=NULL
    if(file.exists(.myFilePathMakerFn(paste0("res_clustering10_",.ArgList$nPCs),argList=.ArgList))){
      load(.myFilePathMakerFn(paste0("res_clustering10_",.ArgList$nPCs),argList=.ArgList))
      
      tmpNames=unlist(lapply(res_clusters,function(x) x$nPCs))
      tmpNames=paste0("nPCs_",tmpNames)
      names(res_clusters)=tmpNames
      
      slResolution=res_clusters[[paste0('nPCs_',.ArgList$nPCs)]]$clusters
      
    }
    
    if(!is.null(slResolution)){
      if(nrow(slResolution)>0){
        
        load(.myFilePathMakerFn("UMAP_anno",argList=.ArgList,pseudoImportant = F))
    
    
        if(!all(row.names(slResolution)==row.names(pd))){
          stop("Error in mapping clustering res and the pd anno file!")
        }
    
        slResolution=cbind(pd,slResolution)
        slResolution$anno_orig_cellState[grepl("Micro",slResolution$anno_orig_cellState)|slResolution$anno_orig_cellState %in% c("MG","MGL","MGL1","MGL2","MGL3")]="Microglia"
        tmp=slResolution[!is.na(slResolution$anno_orig_cellState),]
        
        tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("H1/2M","H1M","H2M","homeostatic mic_s2_c1","homeostatic mic_s2_c2"))]="H1/2M"
        if(F){
          tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("cDC1","cDC2","pDC","migDC"))]="cDC"
          tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("T/NKT cells","T cells","NK cells"))]="T/NKT cells"
          tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("Mnc","Non-cl. monocytes","MCs"))]="Monocyte/Mdc"
          tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("ILC"))]="innate_lymphoid"
          tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("CAM"))]="CNS_assoc_macrophages"
          tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("BAM"))]="BAM"
        }
        
        tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("CPM","CRM"))]="cycl_prolif_microglia"
        tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("IRM"))]="INF_resp_microglia"
        tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("Proliferation"))]="dissociation_prolif"
      
        slClusters=apply(table(tmp$anno_orig_cellState,tmp$res.0.7),2,sum)
        slClusters=names(slClusters)[slClusters>150]
        slClusters=as.numeric(slClusters)
        slClusters=as.character(slClusters[order(slClusters,decreasing = F)])
        
        tmpData=as.data.frame(table(tmp$anno_orig_cellState,tmp$res.0.7))
        tmpData$Var1=as.character(tmpData$Var1)
        tmpData$Var2=as.character(tmpData$Var2)
        tmpData=tmpData[which(tmpData$Var2 %in% slClusters),]
        tmpData$fraction=(-1)
        for(i in 1:nrow(tmpData)){
          tmpData$fraction[i]=tmpData$Freq[i]/sum(tmpData$Freq[which(tmpData$Var1==as.character(tmpData$Var1[i]))])
        }
        
        tmpData$Var2=factor(as.character(tmpData$Var2),levels=slClusters)
        p1=ggplot(tmpData,aes(Var2,Var1,fill=fraction))+geom_tile(color="black")+scale_fill_gradientn(colours=c("white","white","yellow","red"),values=c(0,0.05,0.2,1))+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=1))+ylab("Known cell types")+xlab("Cluster #")
    
        tmp=as.character(slResolution$anno_batch)
        tmp=gsub("Van_Hove","VanHove",tmp)
        tmp=strsplit(tmp,"_")
        tmp=unlist(lapply(tmp,function(x) paste0(x[1],"_",x[2])))
        tmp=data.frame(dataset=tmp,cluster=slResolution$res.0.7,stringsAsFactors = F)
        tmp$cluster=as.character(tmp$cluster)
        tmp=tmp[which(as.character(tmp$cluster) %in% slClusters),]
        tmp=as.data.frame(table(tmp$dataset,tmp$cluster))
        tmp$Var1=as.character(tmp$Var1)
        tmp$Var2=as.character(tmp$Var2)
        tmp$fraction_dataset_norm=(-1)
        for(i in 1:nrow(tmp)){
          tmp$fraction_dataset_norm[i]=tmp$Freq[i]/sum(tmp$Freq[which(tmp$Var1==tmp$Var1[i])])
        }
        tmp$fraction_cluster_norm=(-1)
        for(i in 1:nrow(tmp)){
          tmp$fraction_cluster_norm[i]=tmp$fraction_dataset_norm[i]/sum(tmp$fraction_dataset_norm[which(tmp$Var2==tmp$Var2[i])])
        }
        tmp$Var2=factor(as.character(tmp$Var2),levels=slClusters)
        #tmpCols=c(c("#986800","#2e53d8","#65ca3d","#8e22ab","#c7ff9a","#d6004c","#02a974","#f66626","#02a3eb","#ff906a","#414d98","#ffcb83","#faa5ff"),palette())
        tmpCols=hues::iwanthue(length(unique(tmp$Var1)))
        
        p2=ggplot(tmp,aes(Var2,fraction_cluster_norm,fill=Var1))+geom_bar(stat = "identity",alpha=0.8)+ylab("Cluster composition")+scale_fill_manual(values = tmpCols)+theme_bw()+theme(panel.grid = element_blank())+xlab("Cluster #")
        
        x=as.data.frame(table(slResolution$res.0.7))
        x$Var1=as.character(x$Var1)
        x=x[which(x$Var1 %in% slClusters),]
        x=x[order(x$Freq,decreasing = T),]
        x$Var1=factor(x$Var1,levels=x$Var1)
        p3=ggplot(x,aes(Var1,Freq))+geom_bar(stat='identity')+theme_bw()+scale_y_log10()+ylab("Cluster size")+xlab("Cluster #")
        p=p1/(p2|p3)
    
        ggsave(plot=p,.myFilePathMakerFn("clustering_res_v1",argList=.ArgList,pdf=T),width = 15,height = 10)
  
        slResolution$res.0.7=as.numeric(as.character(slResolution$res.0.7))
        slResolution$res.0.7=factor(slResolution$res.0.7,levels=unique(slResolution$res.0.7[order(slResolution$res.0.7,decreasing = F)]))
        
        p=.myDimPlotFn(object=slResolution[slResolution$res.0.7 %in% slClusters,], dimCols = c("UMAP_1", "UMAP_2"),attCol = "res.0.7")
        ggsave(plot=p,.myFilePathMakerFn("clustering_res_v2",argList=.ArgList,pdf=T),width = 15,height = 10)
        
        
        markerGenes=read.table("~/myBucket/markergenes/macroglia_markers.txt",sep="\t",header = T,stringsAsFactors = F)
        markerGenes=markerGenes[!markerGenes$Gene %in% c("IFIT1","F13A1","MGL2","LYVE1","CENPA","UBE2C","ARG1","FABP5"),]
        markerGenes=markerGenes[which(markerGenes$MainMarker!="Tushar"),]
        load(do.call('.myFilePathMakerFn',args=c(list(filename="exp_merged",expData=T),.ArgList)))
        tmp=tmp[,colnames(tmp) %in% row.names(pd)]
        pd2=pd[match(colnames(tmp),row.names(pd)),]
        p=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd2,inputGenes = markerGenes$Gene)
        ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_microglia_v1",pdf=T),.ArgList)),width = 49,height = 40)
        #########################################################
        
        p=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd2,inputGenes = "MS4A7",inputGroupsCol = 'anno_batch',ncol=8)
        if(any(class(p)!=class(""))){
          ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_MS4A7",pdf=T),.ArgList)),width = 49,height = 49)
        }
        
        p=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd2,inputGenes = "AXL",inputGroupsCol = 'anno_batch',ncol=8)
        if(any(class(p)!=class(""))){
          ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_AXL",pdf=T),.ArgList)),width = 49,height = 49)
        }
        
        
        p=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd2,inputGenes = "APOE",inputGroupsCol = 'anno_batch',ncol=8)
        if(any(class(p)!=class(""))){
          ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_APOE",pdf=T),.ArgList)),width = 49,height = 49)
        }
        
        
        p=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd2,inputGenes = "TMEM119,P2RY12,FCRLS",inputGroupsCol = 'anno_batch',ncol=8)
        if(any(class(p)!=class(""))){
          ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_TMEM119-P2RY12-FCRLS",pdf=T),.ArgList)),width = 49,height = 49)
        }
        
        
        p=.my2dPlot_counts(inputPCA=pd2,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(markerGenes$Gene,function(x) x),geneNameCol="gene_short_name",expData=tmp,ncolumns=6)
        ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_dsCount_v1",pdf=T),.ArgList)),width = 49,height = 40)
        #################################################
        
        markerGenes=read.table("~/myBucket/markergenes/macroglia_markers.txt",sep="\t",header = T,stringsAsFactors = F)
        markerGenes=markerGenes[which(markerGenes$MainMarker=="Yes"),]
        #load(do.call('.myFilePathMakerFn',args=c(list(filename="exp_merged",expData=T),.ArgList)))
        #tmp=tmp[,colnames(tmp) %in% row.names(pd)]
        pd2=pd[match(colnames(tmp),row.names(pd)),]
        p=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd2,inputGenes = markerGenes$Gene)
        ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_microglia_v2",pdf=T),.ArgList)),width = 38,height = 30)
        
        p=.my2dPlot_counts(inputPCA=pd2,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(markerGenes$Gene,function(x) x),geneNameCol="gene_short_name",expData=tmp,ncolumns=6)
        ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_dsCount_v2",pdf=T),.ArgList)),width = 28,height = 20)
        #################################################
        
        markerGenes=read.table("~/myBucket/markergenes/macroglia_markers.txt",sep="\t",header = T,stringsAsFactors = F)
        markerGenes=markerGenes[which(markerGenes$MainMarker=="Tushar"),]
        #load(do.call('.myFilePathMakerFn',args=c(list(filename="exp_merged",expData=T),.ArgList)))
        #tmp=tmp[,colnames(tmp) %in% row.names(pd)]
        pd2=pd[match(colnames(tmp),row.names(pd)),]
        if(length(markerGenes$Gene)>0){
          p=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd2,inputGenes = markerGenes$Gene)
          ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_microglia_Tushar",pdf=T),.ArgList)),width = 38,height = 30)
          
          p=.my2dPlot_counts(inputPCA=pd2,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(markerGenes$Gene,function(x) x),geneNameCol="gene_short_name",expData=tmp,ncolumns=6)
          ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_dsCount_Tushar",pdf=T),.ArgList)),width = 28,height = 20)
          
        }
        #################################################
        
        
        
        missingCells=which(is.na(pd$anno_orig_cellState))
        micCells=which(!pd$anno_orig_cellState %in% .ArgList$exCellStates)
        nonMicCells=which(pd$anno_orig_cellState %in% .ArgList$exCellStates)
        micCells=setdiff(micCells,missingCells)
        nonMicCells=setdiff(nonMicCells,missingCells)
        
        
        #load(do.call('.myFilePathMakerFn',args=c(list(filename="exp_merged",expData=T),.ArgList)))
        
        markerGenes=read.table("~/myBucket/markergenes/macroglia_markers.txt",sep="\t",header = T,stringsAsFactors = F)
        markerGenes=markerGenes[!markerGenes$Gene %in% c("IFIT1","F13A1","MGL2","LYVE1","CENPA","UBE2C","ARG1","FABP5"),]
        #tmp=tmp[,colnames(tmp) %in% row.names(pd)]
        pd2=pd[match(colnames(tmp),row.names(pd)),]
        lbls=rep("microglia",ncol(tmp))
        lbls[is.na(pd2$anno_orig_cellState)]="NA"
        lbls[which(pd2$anno_orig_cellState %in% .ArgList$exCellStates)]="non_microglia"
        p=.myVlnPlot(inputSeurat = tmp,inputGenes = markerGenes$Gene,inputGroups = lbls,pt.size=0,ncol=6)
        ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_microglia_v3",pdf=T),.ArgList)),width = 38,height = 30)
        
        
        
        if(length(nonMicCells)>0){
          pMic=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd,inputGenes = markerGenes$Gene,slIndx=micCells)
          ggsave(plot=pMic,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_micCells",pdf=T),.ArgList)),width=33,height=21)
          
          pNonMic=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd,inputGenes = markerGenes$Gene,slIndx=nonMicCells)
          ggsave(plot=pNonMic,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_nonMicCells",pdf=T),.ArgList)),width=33,height=21)
          
        }
        
      }
    } else {
      load(.myFilePathMakerFn("UMAP_anno",argList=.ArgList,pseudoImportant = F))
      
      markerGenes=read.table("~/myBucket/markergenes/macroglia_markers.txt",sep="\t",header = T,stringsAsFactors = F)
      markerGenes=markerGenes[!markerGenes$Gene %in% c("F13A1","MGL2","CENPA","UBE2C","ARG1","FABP5"),]
      markerGenes=markerGenes[which(markerGenes$MainMarker!="Tushar"),]
      load(do.call('.myFilePathMakerFn',args=c(list(filename="exp_merged",expData=T),.ArgList)))
      tmp=tmp[,colnames(tmp) %in% row.names(pd)]
      pd2=pd[match(colnames(tmp),row.names(pd)),]
      p=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd2,inputGenes = markerGenes$Gene)
      ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_microglia_v1",pdf=T),.ArgList)),width = 49,height = 40)
      #########################################################
      
      
      markerGenes=read.table("~/myBucket/markergenes/macroglia_markers.txt",sep="\t",header = T,stringsAsFactors = F)
      markerGenes=markerGenes[which(markerGenes$MainMarker=="Yes"),]
      #load(do.call('.myFilePathMakerFn',args=c(list(filename="exp_merged",expData=T),.ArgList)))
      #tmp=tmp[,colnames(tmp) %in% row.names(pd)]
      pd2=pd[match(colnames(tmp),row.names(pd)),]
      p=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd2,inputGenes = markerGenes$Gene)
      ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_microglia_v2",pdf=T),.ArgList)),width = 38,height = 30)
      #################################################
      
      markerGenes=read.table("~/myBucket/markergenes/macroglia_markers.txt",sep="\t",header = T,stringsAsFactors = F)
      markerGenes=markerGenes[!markerGenes$Gene %in% c("IFIT1","F13A1","MGL2","LYVE1","CENPA","UBE2C","ARG1","FABP5"),]
      #markerGenes=markerGenes[which(markerGenes$MainMarker=="Tushar"),]
      #load(do.call('.myFilePathMakerFn',args=c(list(filename="exp_merged",expData=T),.ArgList)))
      #tmp=tmp[,colnames(tmp) %in% row.names(pd)]
      pd2=pd[match(colnames(tmp),row.names(pd)),]
      if(length(markerGenes$Gene)>0){
        p=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd2,inputGenes = markerGenes$Gene)
        ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_microglia",pdf=T),.ArgList)),width = 38,height = 30)
        
        p=.my2dPlot_counts(inputPCA=pd2,batch_values="anno_batch",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=lapply(markerGenes$Gene,function(x) x),geneNameCol="gene_short_name",expData=tmp,ncolumns=6)
        ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_dsCount",pdf=T),.ArgList)),width = 28,height = 20)
        
      }
      #################################################
      
      
      
      missingCells=which(is.na(pd$anno_orig_cellState))
      micCells=which(!pd$anno_orig_cellState %in% .ArgList$exCellStates)
      nonMicCells=which(pd$anno_orig_cellState %in% .ArgList$exCellStates)
      micCells=setdiff(micCells,missingCells)
      nonMicCells=setdiff(nonMicCells,missingCells)
      
      
      #load(do.call('.myFilePathMakerFn',args=c(list(filename="exp_merged",expData=T),.ArgList)))
      
      markerGenes=read.table("~/myBucket/markergenes/macroglia_markers.txt",sep="\t",header = T,stringsAsFactors = F)
      markerGenes=markerGenes[!markerGenes$Gene %in% c("F13A1","MGL2","CENPA","UBE2C","ARG1","FABP5"),]
      #tmp=tmp[,colnames(tmp) %in% row.names(pd)]
      pd2=pd[match(colnames(tmp),row.names(pd)),]
      lbls=rep("microglia",ncol(tmp))
      lbls[is.na(pd2$anno_orig_cellState)]="NA"
      lbls[which(pd2$anno_orig_cellState %in% .ArgList$exCellStates)]="non_microglia"
      p=.myVlnPlot(inputSeurat = tmp,inputGenes = markerGenes$Gene,inputGroups = lbls,pt.size=0,ncol=6)
      ggsave(plot=p,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_microglia_v3",pdf=T),.ArgList)),width = 38,height = 30)
      
      
      
      if(length(nonMicCells)>0){
        pMic=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd,inputGenes = markerGenes$Gene,slIndx=micCells)
        ggsave(plot=pMic,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_micCells",pdf=T),.ArgList)),width=33,height=21)
        
        pNonMic=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd,inputGenes = markerGenes$Gene,slIndx=nonMicCells)
        ggsave(plot=pNonMic,file=do.call('.myFilePathMakerFn',args=c(list(filename="markers_nonMicCells",pdf=T),.ArgList)),width=33,height=21)
        
      }
    }
    
  
    
  } else if(!file.exists(.myFilePathMakerFn("UMAP_anno",argList=.ArgList,pseudoImportant = F))){
      print(paste("Error!:",.myFilePathMakerFn("UMAP_anno",argList=.ArgList,pseudoImportant = F)))
  }
  
  cat("Performing consensus marker analysis\n")
  cat("      Consensus marker step 1 ...\n")
  resDE1=.myConcensusDEFn_step1(argList=.ArgList)
  cat("      Consensus marker step 2 ...\n")
  resDE2=.myConcensusDEFn_step2(argList=.ArgList)
  
  cat("Performing consensus marker cluster analysis\n")
  .myConsensusDEClusteringFn(argList=.ArgList,pval_thr=pval_thr,initial_clustering=NULL)
  
  cat("Visualizing the marker data ...\n")
  .myConsensusDE_markerRepresentationFn(argList=.ArgList)
    
    print("Done!")
}

.myMainRunFn_liger=function(inputExpData,liger_embeddings,input_highly_var_genes=NULL,organism,prefix,external_DE_path=NULL,external_DE_name=NULL,runSetting,runIndx,ncores=6,newRun=F,DE_supportingFractionThr=0.1,DE_n.adaptiveKernel=20,DE_nPropIter=1,uniformZscore=F,dist_zscore_gamma=F,dist_zscore_nbinom=F,regularize=F,geoMean=F,internal_pseudocell_count=200,pval_thr=0.001,exclude_non_freq_pseudocells=F){
  
  if(!is.null(liger_embeddings)){
    if(nrow(liger_embeddings)!=ncol(inputExpData)|length(setdiff(colnames(inputExpData),row.names(liger_embeddings)))>0){
      stop("Error! input expData don't match with the provided embeddings!")
    }
  }
  
  if(sum(duplicated(row.names(inputExpData)))>0){
    print(paste(sum(duplicated(row.names(inputExpData))),"Duplicate gene ids were found! duplicates were randomly removed from the data"))
    inputExpData=inputExpData[!duplicated(row.names(inputExpData)),]
  }
  
  if(dist_zscore_gamma){
    dist_zscore_gamma=T
    dist_zscore_norm=F
    dist_zscore_nbinom=F
  } else if(dist_zscore_nbinom){
    dist_zscore_gamma=F
    dist_zscore_norm=F
    dist_zscore_nbinom=T
  } else {
    dist_zscore_gamma=F
    dist_zscore_norm=T
    dist_zscore_nbinom=F
  }
  
  
  .ArgList=.myArgFn(runIndx=runIndx,ncores=ncores,includeHuman=F,includeMouse=F,uniformZscore=uniformZscore,DE_supportingFractionThr=DE_supportingFractionThr,DE_n.adaptiveKernel=DE_n.adaptiveKernel,DE_nPropIter=DE_nPropIter,dist_zscore_gamma=dist_zscore_gamma,dist_zscore_norm=dist_zscore_norm,dist_zscore_nbinom=dist_zscore_nbinom,regularize=regularize,geoMean=geoMean,prefix=prefix,newRun = newRun,inputDf=runSetting,internal_pseudocell_count=internal_pseudocell_count,external_DE_path=external_DE_path,external_DE_name=external_DE_name)
  .ArgList$exclude_non_freq_pseudocells=exclude_non_freq_pseudocells
  .ArgList$input_highly_var_genes=input_highly_var_genes
  
  cat('Loading and preprocessing the datasets\n')
  reRunCheck=tryCatch({load(.myFilePathMakerFn("pca_anno",argList=.ArgList,pseudoImportant = F));F}, error=function(e) {return(T)})
  
  if(!reRunCheck){
    reRunCheck=tryCatch({load(.myFilePathMakerFn("exp_merged",argList=.ArgList,expData=T));F}, error=function(e) {return(T)})
  }
  
  if(reRunCheck){
    
    data=.myReadData_spliterFn(inputData=inputExpData,removeHighExp=.ArgList$excludeHighExp)
    
    
    cat('Selecting the highly variable genes\n')
    data=.myHighVarGeneSlFn(data,dataorganism=organism,argList = .ArgList)
    
    
    if(sum(colnames(data$data_m@meta.data)=="QC_Gene_unique_count")>0){
      expData=data$data_m@assays$RNA@data
      expData=expData[row.names(expData) %in% data$varFeatures[,1],]
      expData=as.matrix(expData)
      nUMIdata=data$data_m$QC_Gene_unique_count
      
      cordata=NULL
      for(ik in 1:nrow(expData)){
        cordata=rbind(cordata,data.frame(gene=row.names(expData)[ik],pearson=cor(as.numeric(expData[ik,]),as.numeric(nUMIdata)),stringsAsFactors = F))
      }
      
      print(paste("Removed",sum(abs(cordata$pearson)>.ArgList$UMI_cor_thr,na.rm = T),"out of",nrow(cordata),"Genes due to correlation (abs(cor)>",.ArgList$UMI_cor_thr,") with # unique genes"))
      data$varFeatures=data$varFeatures[which(data$varFeatures[,1] %in% as.character(cordata$gene[which(abs(cordata$pearson)<=.ArgList$UMI_cor_thr)])),]
      
    }
    
    tmp=data[c("varFeatures","allGenes" )]
    if(!file.exists(.myFilePathMakerFn("varGenes",argList=.ArgList,varGenes = T))){
      save(tmp,file=.myFilePathMakerFn("varGenes",argList=.ArgList,varGenes =T))
    }
    
    tmp=data$data_m
    if(!file.exists(.myFilePathMakerFn("exp_merged",argList=.ArgList,expData=T))){
      save(tmp,file=.myFilePathMakerFn("exp_merged",argList=.ArgList,expData=T))
    }
    
    cat('Scaling the data and PCA\n')
    .myPCAfn(data,argList = .ArgList,UMI_cor_thr=.ArgList$UMI_cor_thr)
    
  }
  
  reRunCheck=tryCatch({load(.myFilePathMakerFn("pca_centroids",argList=.ArgList));F}, error=function(e) {return(T)})
  
  if(reRunCheck){
    if(is.null(liger_embeddings)){
      cat('Harmony analysis\n')
      #rm(list=ls())
      #source("/data/vahid/data/mySC.R")
      
      
      
      library(future)
      plan("multiprocess", workers = 8)
      plan()
      options(future.globals.maxSize = 1000 * 1024^4)
      load(.myFilePathMakerFn("pca_anno",argList=.ArgList,pseudoImportant = F))
      
      
      harmony_embeddings <- harmony::HarmonyMatrix(pca_res[,1:.ArgList$nPCs], pd, 'anno_batch', do_pca = FALSE, verbose=FALSE)
      .mySaveFn(harmony_embeddings,file=.myFilePathMakerFn("harmony-embeddings",argList=.ArgList,pseudoImportant = F))
      
      slCols=c("anno_age","anno_batch","anno_condition","anno_orig_cellState","anno_region","anno_sex","anno_genotype","anno_RNAseq_type","anno_RNAseq_type2")
      resSummary=NULL
      for(icol in slCols){
        preEffSizes=c()
        tmpInd=which(!is.na(pd[,icol]))
        if(length(unique(pd[tmpInd,icol]))>1){
          for(i in 1:.ArgList$nPCs){
            fit <- aov(pca_res[tmpInd,i] ~ as.factor(as.character(pd[tmpInd,icol])),data = pd)
            fit=sjstats::eta_sq(fit)
            preEffSizes=c(preEffSizes,fit$etasq)
          }
          postEffSizes=c()
          for(i in 1:.ArgList$nPCs){
            fit <- aov(harmony_embeddings[tmpInd,i] ~ as.factor(as.character(pd[tmpInd,icol])),data = pd)
            fit=sjstats::eta_sq(fit)
            postEffSizes=c(postEffSizes,fit$etasq)
          }
          
          tmp1=data.frame(PC=1:.ArgList$nPCs,eta_sq=preEffSizes)
          tmp1$status="pre_Harmony"
          tmp1$variable=icol
          
          tmp2=data.frame(PC=1:.ArgList$nPCs,eta_sq=postEffSizes)
          tmp2$status="post_Harmony"
          tmp2$variable=icol
          resSummary=rbind(resSummary,tmp1)
          resSummary=rbind(resSummary,tmp2)
        }
        
      }
      save(resSummary,file=.myFilePathMakerFn("harmony_stat_summary",argList=.ArgList))
      
      #kmeans clustering
      set.seed(1)
      res_clust=kmeans(harmony_embeddings[,1:.ArgList$nPCs],.ArgList$internal_pseudocell_count,iter.max = 10000,algorithm = "Lloyd")
      
      #row.names(res_clust$centers)=paste0("C_",1:nrow(res_clust$centers))
      #kcenters=res_clust$centers
      #tormList=c()
      #while(length(tormList)<100){
      #  knet=RANN::nn2(data = kcenters, k = 20, eps = 0)
      #  knet_mutual=rowSums(knet$nn.dist)
      #  tormList=c(tormList,row.names(kcenters)[which(knet_mutual==min(knet_mutual))])
      #  kcenters=kcenters[!row.names(kcenters) %in% tormList,]
      #}
      #kcenters=res_clust$centers[!row.names(res_clust$centers) %in% tormList,]
      #knet=RANN::nn2(data=kcenters,query = harmony_embeddings,k=2,eps=0)
      
      #knet2=RANN::nn2(kcenters,k=10,eps=0)
      #knet2=data.frame(FC=row.names(kcenters),SC=row.names(kcenters)[knet2$nn.idx[,2]],distC=knet2$nn.dists[,2],stringsAsFactors = F)
      #knet3=data.frame(sample=row.names(harmony_embeddings),FC=row.names(kcenters)[knet$nn.idx[,1]],distSample=knet$nn.dists[,1],stringsAsFactors = F)
      #knet2=merge(knet2,knet3,by="FC",all=T)
      
      res_clusters=data.frame(sample=names(res_clust$cluster),cluster_id=res_clust$cluster,stringsAsFactors = F)
      save(res_clusters,file=.myFilePathMakerFn("kmeans_res_clusters",argList=.ArgList))
      
      pca_centroid=res_clust$centers
      save(pca_centroid,file=.myFilePathMakerFn("pca_centroids",argList=.ArgList))
      
      
    } else {
      
      harmony_embeddings <- liger_embeddings[,1:.ArgList$nPCs]
      .mySaveFn(harmony_embeddings,file=.myFilePathMakerFn("harmony-embeddings",argList=.ArgList,pseudoImportant = F))
      
      #kmeans clustering
      set.seed(1)
      doClustering=T
      itrClustering=0
      while(doClustering&itrClustering<10){
        itrClustering=itrClustering+1
        res_clust=kmeans(harmony_embeddings[,1:.ArgList$nPCs],.ArgList$internal_pseudocell_count,iter.max = 10000,algorithm = "Lloyd")
        if(sum(is.na(res_clust$centers))==0){
          doClustering=F
        }
      }
      
      if(sum(is.na(res_clust$centers))>0){
        stop("Error in identification of the pseudocells")
      }
      
      res_clusters=data.frame(sample=names(res_clust$cluster),cluster_id=res_clust$cluster,stringsAsFactors = F)
      save(res_clusters,file=.myFilePathMakerFn("kmeans_res_clusters",argList=.ArgList))
      
      pca_centroid=res_clust$centers
      save(pca_centroid,file=.myFilePathMakerFn("pca_centroids",argList=.ArgList))
      
    }
  }
  
  
  
  cat("Running UMAP\n")
  #rm(list = setdiff(ls(), lsf.str()))
  #load(.myFilePathMakerFn("harmony-embeddings",argList=.ArgList))
  reRunCheck=tryCatch({load(.myFilePathMakerFn("UMAP_anno",argList=.ArgList,pseudoImportant = F));F}, error=function(e) {return(T)})
  
  if(!reRunCheck){
    reRunCheck=tryCatch({load(.myFilePathMakerFn("UMAP_centroid",argList=.ArgList));F}, error=function(e) {return(T)})
  }
  
  if(reRunCheck){
    load(.myFilePathMakerFn("harmony-embeddings",argList=.ArgList))
    load(.myFilePathMakerFn("pca_anno",argList=.ArgList,pseudoImportant = F))
    
    load(.myFilePathMakerFn("pca_centroids",argList=.ArgList))
    
    pca_centroid=pca_centroid[,1:.ArgList$nPCs]
    
    
    
    tst=.reductionUMAPFn(harmony_embeddings,umap.method='umap-learn',testPCAembeddings=pca_centroid)
    resUMAP=tst$embedding
    .mySaveFn(resUMAP,file=.myFilePathMakerFn("UMAP_res",argList=.ArgList))
    
    x=as.data.frame(resUMAP)
    x=cbind(x,pd)
    row.names(x)=row.names(pd)
    pd=x
    .mySaveFn(pd,file=.myFilePathMakerFn("UMAP_anno",argList=.ArgList,pseudoImportant = F))
    
    tst=tst[["test"]]
    UMAP_centroid=data.frame(centroid=row.names(pca_centroid),UMAP_1=tst[,1],UMAP_2=tst[,2],stringsAsFactors = F)
    .mySaveFn(UMAP_centroid,file=.myFilePathMakerFn("UMAP_centroid",argList=.ArgList))
    
    
  }
  
  cat("Saving the plots\n")
  
  if(file.exists(.myFilePathMakerFn("UMAP_anno",argList=.ArgList,pseudoImportant = F))&!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename="markers_dsCount_Tushar",pdf=T),.ArgList)))){
    
    if(file.exists(.myFilePathMakerFn("harmony_stat_summary",argList=.ArgList))){
      load(.myFilePathMakerFn("harmony_stat_summary",argList=.ArgList))
      p=ggplot(resSummary,aes(variable,eta_sq,color=status))+geom_boxplot()+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=1,color = "black"))+xlab("")+ylab("variance explained")
      ggsave(plot = p,do.call('.myFilePathMakerFn',args=c(list(filename="plot_Harmony_efficiency",pdf=T),.ArgList)),width = 8,height = 6)
      
    }
    
    
    load(do.call('.myFilePathMakerFn',args=c(list(filename="UMAP_anno",pseudoImportant = F),.ArgList)))
    
    if(sum(grepl("human_Tushar",pd$anno_batch))>0){
      p=.myDimPlotFn(object=pd[which(grepl("human_Tushar",pd$anno_batch)),], dimCols = c("UMAP_1", "UMAP_2"),attCol = "clusters")
      ggsave(plot=p,.myFilePathMakerFn("clustering_Tushar_v1",argList=.ArgList,pdf=T),width = 15,height = 10)
      p=.my2dPlot2(inputPCA = pd,batchCol = "clusters",reductionCols = c('UMAP_1','UMAP_2'))
      ggsave(plot = p,do.call('.myFilePathMakerFn',args=c(list(filename="clustering_Tushar_v2",pdf=T),.ArgList)),width = 12,height = 8)
      
    }
    
    pd$anno_batch[is.na(pd$anno_batch)]=as.character(pd$anno_batch)[is.na(pd$anno_batch)]
    pd$anno_batch[which(grepl("Zywitza_30485812",pd$anno_batch))]="mouse_Zywitza_30485812"
    pd$anno_batch[which(grepl("mouse_Lau_scRNA_32320664",pd$anno_batch))]="mouse_Lau_32320664"
    pd$anno_batch[which(grepl("mouse_Sierksma_31951107",pd$anno_batch))]="mouse_Sierksma_31951107"
    pd$anno_batch[which(grepl("Sousa_30206190",pd$anno_batch))]="mouse_Sousa_30206190"
    pd$anno_batch[which(grepl("mouse_SalaFrigerio",pd$anno_batch))]="mouse_SalaFrigerio_31018141"
    pd$anno_batch[which(grepl("mouse_Masuda",pd$anno_batch))]="mouse_Masuda"
    pd$anno_batch[which(grepl("Chen_28355573",pd$anno_batch))]="mouse_Chen_28355573"
    pd$anno_batch[which(grepl("Chen_30773317",pd$anno_batch))]="mouse_Chen_30773317"
    save(pd,file=do.call('.myFilePathMakerFn',args=c(list(filename="UMAP_anno",pseudoImportant = F),.ArgList)))
    
    
    p=.my2dPlot2(inputPCA = pd,batchCol = "anno_batch",reductionCols = c('UMAP_1','UMAP_2'))
    ggsave(plot = p,do.call('.myFilePathMakerFn',args=c(list(filename="batchEffect",pdf=T),.ArgList)),width = 36,height = 26)
    
    
    p=ggplot(pd,aes(UMAP_1,UMAP_2,color=anno_batch))+geom_point(alpha=0.5,size=0.5)
    ggsave(plot = p,do.call('.myFilePathMakerFn',args=c(list(filename="plot_batch_v1",pdf=T),.ArgList)),width = 12,height = 8)
    
    
    p=.my2dPlot2(inputPCA = pd,batchCol = "anno_batch",reductionCols = c('UMAP_1','UMAP_2'))
    ggsave(plot = p,do.call('.myFilePathMakerFn',args=c(list(filename="plot_batch_v2",pdf=T),.ArgList)),width = 12,height = 8)
    
    
    tmp=pd[!is.na(pd$anno_orig_cellState),]
    tmp=tmp[tmp$anno_orig_cellState!="",]
    tmp$anno_orig_cellState[grepl("Micro",tmp$anno_orig_cellState)]="Microglia"
    tmp$anno_orig_cellState[tmp$anno_orig_cellState %in% c("MG","MGL","MGL1","MGL2","MGL3","mg","Mic")]="Microglia"
    tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("H1/2M","H1M","H2M","homeostatic mic_s2_c1","homeostatic mic_s2_c2"))]="H1/2M"
    if(F){
      tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("cDC1","cDC2","pDC","migDC"))]="cDC"
      tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("T/NKT cells","T cells","NK cells"))]="T/NKT cells"
      tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("Mnc","Non-cl. monocytes","MCs"))]="Monocyte/Mdc"
      tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("ILC"))]="innate_lymphoid"
      tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("CAM"))]="CNS_assoc_macrophages"
      tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("BAM"))]="BAM"
      
    }
    tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("CPM","CRM"))]="cycling_prolif_microglia"
    tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("IRM"))]="IFN_resp_microglia"
    tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("Proliferation"))]="dissociation_prolif"
    #tmp$anno_orig_cellState[which(tmp$anno_orig_cellState %in% c("TRM"))]="transiting_resp_microglia"
    p1=ggplot(tmp,aes(UMAP_1,UMAP_2,color=anno_orig_cellState))+geom_point(alpha=0.5,size=0.5)+scale_color_manual(values=hues::iwanthue(length(unique(as.character(tmp$anno_orig_cellState)))))
    p2=ggplot(tmp,aes(UMAP_1,UMAP_2,color=anno_batch))+geom_point(alpha=0.5,size=0.5)
    p=p1+p2
    ggsave(plot=p,.myFilePathMakerFn("plot_cellStates_v1.1",argList=.ArgList,pdf=T),width = 14,height = 8)
    
    p=.my2dPlot2(inputPCA = tmp,batchCol = "anno_orig_cellState",reductionCols = c('UMAP_1','UMAP_2'))
    ggsave(plot=p,.myFilePathMakerFn("plot_cellStates_v2.1",argList=.ArgList,pdf=T),width = 12,height = 8)
    
    tmp=tmp[tmp$anno_orig_cellState!="Microglia",]
    if(nrow(tmp)>0){
      p1=ggplot(tmp,aes(UMAP_1,UMAP_2,color=anno_orig_cellState))+geom_point(alpha=0.5,size=0.5)+scale_color_manual(values=hues::iwanthue(length(unique(as.character(tmp$anno_orig_cellState)))))
      p2=ggplot(tmp,aes(UMAP_1,UMAP_2,color=anno_batch))+geom_point(alpha=0.5,size=0.5)
      p=p1+p2
      ggsave(plot=p,.myFilePathMakerFn("plot_cellStates_v1",argList=.ArgList,pdf=T),width = 14,height = 8)
      
      if(length(which(tmp$anno_orig_cellState %in% c("ARM","IFN_resp_microglia","cycling_prolif_microglia","BAM","DAM","dissociation_prolif","CNS_assoc_macrophages")))>0){
        p=ggplot(tmp[which(tmp$anno_orig_cellState %in% c("ARM","IFN_resp_microglia","cycling_prolif_microglia","BAM","DAM","dissociation_prolif","CNS_assoc_macrophages")),],aes(UMAP_1,UMAP_2,color=anno_orig_cellState))+geom_point(alpha=0.5,size=0.5)+scale_color_manual(values=hues::iwanthue(length(unique(as.character(tmp$anno_orig_cellState)))))
        ggsave(plot=p,.myFilePathMakerFn("plot_cellStates_v3",argList=.ArgList,pdf=T),width = 14,height = 8)
      }
      
      p=.my2dPlot2(inputPCA = tmp,batchCol = "anno_orig_cellState",reductionCols = c('UMAP_1','UMAP_2'))
      ggsave(plot=p,.myFilePathMakerFn("plot_cellStates_v2",argList=.ArgList,pdf=T),width = 12,height = 8)
      
    }
    
    tmp=pd
    if(length(which(tmp$anno_orig_cellState %in% .ArgList$exCellStates))>0){
      tmp=tmp[-which(tmp$anno_orig_cellState %in% .ArgList$exCellStates),]
    }
    
    if(nrow(tmp)>0){
      tmp=tmp[!grepl("mouse_Sierksma_31951107",tmp$anno_batch),]
      p1=.myDimPlotFn(object=tmp,dimCols = c("UMAP_1","UMAP_2"),attCol = 'anno_RNAseq_type',label = F)
      p2=.my2dPlot2(inputPCA = tmp,batchCol = "anno_RNAseq_type",reductionCols = c('UMAP_1','UMAP_2'))
      p=p1+p2
      ggsave(plot=p,.myFilePathMakerFn("plot_RNAseqType_v1",argList=.ArgList,pdf=T),width = 14,height = 8)
    }
    
    
    tmp=pd[which(pd$anno_batch!="mouse_Wu_29024657"),]
    tmp=tmp[which(!grepl("mouse_Sierksma_31951107",tmp$anno_batch)),]
    if(length(which(tmp$anno_orig_cellState %in% .ArgList$exCellStates))>0){
      tmp=tmp[-which(tmp$anno_orig_cellState %in% .ArgList$exCellStates),]
    }
    
    if(length(unique(as.character(tmp$anno_RNAseq_type2)))>0){
      tmp$anno_RNAseq_type2=as.character(tmp$anno_RNAseq_type2)
      tmp$anno_RNAseq_type2[which(tmp$anno_RNAseq_type2=="Smart-Seq 2")]="Smart-Seq"
      tmp$anno_RNAseq_type2=factor(as.character(tmp$anno_RNAseq_type2),levels=rev(c("Mars-Seq","Smart-Seq","10X Genomics","Drop-Seq")))
      p1=.myDimPlotFn(object=tmp,dimCols = c("UMAP_1","UMAP_2"),attCol = 'anno_RNAseq_type2',label = F)
      p2=.my2dPlot2(inputPCA = tmp,batchCol = "anno_RNAseq_type2",reductionCols = c('UMAP_1','UMAP_2'))
      p=p1+p2
      ggsave(plot=p,.myFilePathMakerFn("plot_RNAseqType_v2",argList=.ArgList,pdf=T),width = 14,height = 8)
    }
    
    
    
    if(length(which(pd$anno_orig_cellState %in% .ArgList$exCellStates))>0){
      tmp=pd[-which(pd$anno_orig_cellState %in% .ArgList$exCellStates),]
    } else {
      tmp=pd
    }
    
    #tmp=pd[pd$QC_Gene_unique_count>0,]
    if(length(which(tmp$anno_orig_cellState %in% .ArgList$exCellStates))>0){
      tmp=tmp[-which(tmp$anno_orig_cellState %in% .ArgList$exCellStates),]
    }
    
    p1=.my2dPlot_continuous(inputPCA=tmp,continuous_batch_values="QC_Gene_unique_count",reductionCols=c("UMAP_1","UMAP_2"))+ggtitle("Gene_unique_count")
    p2=.my2dPlot_continuous(inputPCA=tmp,continuous_batch_values="QC_Gene_total_count",reductionCols=c("UMAP_1","UMAP_2"))+ggtitle("Gene_total_count")
    p3=.my2dPlot_continuous(inputPCA=tmp,continuous_batch_values="QC_IEG.pct",reductionCols=c("UMAP_1","UMAP_2"))+ggtitle("QC_IEG.pct")
    p4=.my2dPlot_continuous(inputPCA=tmp,continuous_batch_values="QC_top50_pct",reductionCols=c("UMAP_1","UMAP_2"))+ggtitle("QC_top50_pct")
    p=p1+p2
    ggsave(plot=p,.myFilePathMakerFn("QC_plot_v1",argList=.ArgList,pdf=T),width = 14,height = 8)
    
    p=p3+p4
    ggsave(plot=p,.myFilePathMakerFn("QC_plot_v2",argList=.ArgList,pdf=T),width = 14,height = 8)
    
    tmp=pd[!is.na(pd$anno_condition),]
    tmp=tmp[tmp$anno_condition!="",]
    if(length(which(tmp$anno_condition=="Alzheimers_disease"))>0){
      tmp$anno_condition[which(tmp$anno_condition=="Alzheimers_disease")]="AD"
    }
    
    if(!is.null(.ArgList$exCellStates)){
      tmp=tmp[!(tmp$anno_orig_cellState %in% .ArgList$exCellStates),]
    }
    
    if(nrow(tmp)>0){
      p1=ggplot(tmp,aes(UMAP_1,UMAP_2,color=anno_condition))+geom_point(alpha=0.5,size=0.5)
      p2=ggplot(tmp,aes(UMAP_1,UMAP_2,color=anno_batch))+geom_point(alpha=0.5,size=0.5)
      p=p1+p2
      ggsave(plot=p,.myFilePathMakerFn("plot_treatment_v1",argList=.ArgList,pdf=T),width = 14,height = 8)
      
      p=.my2dPlot2(inputPCA = tmp,batchCol = "anno_condition",reductionCols = c('UMAP_1','UMAP_2'))
      ggsave(plot=p,.myFilePathMakerFn("plot_treatment_v2",argList=.ArgList,pdf=T),width = 12,height = 8)
    }
    
    tmp=pd[!is.na(pd$anno_age),]
    tmp=tmp[which(tmp$anno_condition=="WT"),]
    tmp=tmp[!(tmp$anno_orig_cellState %in% .ArgList$exCellStates),]
    if(nrow(tmp)>0){
      tmp$anno_age[which(tmp$anno_age=="1wk")]="0d-1wk"
      tmp$anno_age=factor(as.character(tmp$anno_age),levels=c("prenatal",'0d-1wk','1wk-1mo','1-3mo','3-6mo','6-9mo','9-12mo','12-24mo'))
      p1=ggplot(tmp,aes(UMAP_1,UMAP_2,color=anno_age))+geom_point(alpha=0.5,size=0.5)+scale_color_manual(values=hues::iwanthue(length(unique(as.character(tmp$anno_age)))))+theme_bw()
      p2=ggplot(tmp,aes(UMAP_1,UMAP_2,color=anno_batch))+geom_point(alpha=0.5,size=0.5)
      p=p1+p2
      ggsave(plot=p,.myFilePathMakerFn("plot_age_v1",argList=.ArgList,pdf=T),width = 14,height = 8)
      
      p1=ggplot(tmp[tmp$anno_age %in% c("prenatal","0d-1wk"),],aes(UMAP_1,UMAP_2,color=anno_age))+geom_point(alpha=0.5,size=0.5)+scale_color_manual(values=hues::iwanthue(length(unique(as.character(tmp$anno_age)))))+theme_bw()
      p2=ggplot(tmp[tmp$anno_age %in% c("prenatal","0d-1wk"),],aes(UMAP_1,UMAP_2,color=anno_batch))+geom_point(alpha=0.5,size=0.5)
      p=p1+p2
      ggsave(plot=p,.myFilePathMakerFn("plot_age_v3",argList=.ArgList,pdf=T),width = 14,height = 8)
      
      if(.ArgList$includeMouse){
        p=.my2dPlot2(inputPCA = tmp,batchCol = "anno_age",reductionCols = c('UMAP_1','UMAP_2'))
        ggsave(plot=p,.myFilePathMakerFn("plot_age_v2",argList=.ArgList,pdf=T),width = 12,height = 8)
        
      }
    }
    
    
  } else if(!file.exists(.myFilePathMakerFn("UMAP_anno",argList=.ArgList,pseudoImportant = F))){
    print(paste("Error!:",.myFilePathMakerFn("UMAP_anno",argList=.ArgList,pseudoImportant = F)))
  }
  
  cat("Performing consensus marker analysis\n")
  cat("      Consensus marker step 1 ...\n")
  resDE1=.myConcensusDEFn_step1(argList=.ArgList)
  cat("      Consensus marker step 2 ...\n")
  resDE2=.myConcensusDEFn_step2(argList=.ArgList)
  
  cat("Performing consensus marker cluster analysis\n")
  .myConsensusDEClusteringFn(argList=.ArgList,pval_thr=pval_thr)
  
  cat("Preparing the shiny object ...\n")
  .myShinyTransferFn(argList=.ArgList)
    
  cat("Visualizing the marker data ...\n")
  .myConsensusDE_markerRepresentationFn(argList=.ArgList)
  
  
  print("Done!")
}

.mySimRunFn=function(inputExpData,liger_embeddings,input_highly_var_genes=NULL,organism,prefix,external_DE_path=NULL,external_DE_name=NULL,runSetting,runIndx,ncores=6,newRun=F,DE_supportingFractionThr=0.1,DE_n.adaptiveKernel=20,DE_nPropIter=1,uniformZscore=F,dist_zscore_gamma=F,dist_zscore_nbinom=F,regularize=F,geoMean=F,internal_pseudocell_count=200,pval_thr=0.001,exclude_non_freq_pseudocells=F){
  
  if(!is.null(liger_embeddings)){
    if(nrow(liger_embeddings)!=ncol(inputExpData)|length(setdiff(colnames(inputExpData),row.names(liger_embeddings)))>0){
      stop("Error! input expData don't match with the provided embeddings!")
    }
  }
  
  if(sum(duplicated(row.names(inputExpData)))>0){
    print(paste(sum(duplicated(row.names(inputExpData))),"Duplicate gene ids were found! duplicates were randomly removed from the data"))
    inputExpData=inputExpData[!duplicated(row.names(inputExpData)),]
  }
  
  if(dist_zscore_gamma){
    dist_zscore_gamma=T
    dist_zscore_norm=F
    dist_zscore_nbinom=F
  } else if(dist_zscore_nbinom){
    dist_zscore_gamma=F
    dist_zscore_norm=F
    dist_zscore_nbinom=T
  } else {
    dist_zscore_gamma=F
    dist_zscore_norm=T
    dist_zscore_nbinom=F
  }
  
  
  .ArgList=.myArgFn(runIndx=runIndx,ncores=ncores,includeHuman=F,includeMouse=F,uniformZscore=uniformZscore,DE_supportingFractionThr=DE_supportingFractionThr,DE_n.adaptiveKernel=DE_n.adaptiveKernel,DE_nPropIter=DE_nPropIter,dist_zscore_gamma=dist_zscore_gamma,dist_zscore_norm=dist_zscore_norm,dist_zscore_nbinom=dist_zscore_nbinom,regularize=regularize,geoMean=geoMean,prefix=prefix,newRun = newRun,inputDf=runSetting,internal_pseudocell_count=internal_pseudocell_count,external_DE_path=external_DE_path,external_DE_name=external_DE_name)
  .ArgList$exclude_non_freq_pseudocells=exclude_non_freq_pseudocells
  .ArgList$input_highly_var_genes=input_highly_var_genes
  
  cat('Loading and preprocessing the datasets\n')
  reRunCheck=tryCatch({load(.myFilePathMakerFn("pca_anno",argList=.ArgList,pseudoImportant = F));F}, error=function(e) {return(T)})
  
  if(!reRunCheck){
    reRunCheck=tryCatch({load(.myFilePathMakerFn("exp_merged",argList=.ArgList,expData=T));F}, error=function(e) {return(T)})
  }
  
  if(reRunCheck){
    
    data=.myReadData_spliterFn(inputData=inputExpData,removeHighExp=.ArgList$excludeHighExp)
    
    
    cat('Selecting the highly variable genes\n')
    data=.myHighVarGeneSlFn(data,dataorganism=organism,argList = .ArgList)
    
    
    if(sum(colnames(data$data_m@meta.data)=="QC_Gene_unique_count")>0){
      expData=data$data_m@assays$RNA@data
      expData=expData[row.names(expData) %in% data$varFeatures[,1],]
      expData=as.matrix(expData)
      nUMIdata=data$data_m$QC_Gene_unique_count
      
      cordata=NULL
      for(ik in 1:nrow(expData)){
        cordata=rbind(cordata,data.frame(gene=row.names(expData)[ik],pearson=cor(as.numeric(expData[ik,]),as.numeric(nUMIdata)),stringsAsFactors = F))
      }
      
      print(paste("Removed",sum(abs(cordata$pearson)>.ArgList$UMI_cor_thr,na.rm = T),"out of",nrow(cordata),"Genes due to correlation (abs(cor)>",.ArgList$UMI_cor_thr,") with # unique genes"))
      data$varFeatures=data$varFeatures[which(data$varFeatures[,1] %in% as.character(cordata$gene[which(abs(cordata$pearson)<=.ArgList$UMI_cor_thr)])),]
      
    }
    
    tmp=data[c("varFeatures","allGenes" )]
    if(!file.exists(.myFilePathMakerFn("varGenes",argList=.ArgList,varGenes = T))){
      save(tmp,file=.myFilePathMakerFn("varGenes",argList=.ArgList,varGenes =T))
    }
    
    tmp=data$data_m
    if(!file.exists(.myFilePathMakerFn("exp_merged",argList=.ArgList,expData=T))){
      save(tmp,file=.myFilePathMakerFn("exp_merged",argList=.ArgList,expData=T))
    }
    
    cat('Scaling the data and PCA\n')
    .myPCAfn(data,argList = .ArgList,UMI_cor_thr=.ArgList$UMI_cor_thr)
    
  }
  
  reRunCheck=tryCatch({load(.myFilePathMakerFn("pca_centroids",argList=.ArgList));F}, error=function(e) {return(T)})
  
  if(reRunCheck){
    if(is.null(liger_embeddings)){
      cat('Harmony analysis\n')
      #rm(list=ls())
      #source("/data/vahid/data/mySC.R")
      
      
      
      library(future)
      plan("multiprocess", workers = 8)
      plan()
      options(future.globals.maxSize = 1000 * 1024^4)
      load(.myFilePathMakerFn("pca_anno",argList=.ArgList,pseudoImportant = F))
      
      
      harmony_embeddings <- harmony::HarmonyMatrix(pca_res[,1:.ArgList$nPCs], pd, 'anno_batch', do_pca = FALSE, verbose=FALSE)
      .mySaveFn(harmony_embeddings,file=.myFilePathMakerFn("harmony-embeddings",argList=.ArgList,pseudoImportant = F))
      
      slCols=c("anno_age","anno_batch","anno_condition","anno_orig_cellState","anno_region","anno_sex","anno_genotype","anno_RNAseq_type","anno_RNAseq_type2")
      resSummary=NULL
      for(icol in slCols){
        preEffSizes=c()
        tmpInd=which(!is.na(pd[,icol]))
        if(length(unique(pd[tmpInd,icol]))>1){
          for(i in 1:.ArgList$nPCs){
            fit <- aov(pca_res[tmpInd,i] ~ as.factor(as.character(pd[tmpInd,icol])),data = pd)
            fit=sjstats::eta_sq(fit)
            preEffSizes=c(preEffSizes,fit$etasq)
          }
          postEffSizes=c()
          for(i in 1:.ArgList$nPCs){
            fit <- aov(harmony_embeddings[tmpInd,i] ~ as.factor(as.character(pd[tmpInd,icol])),data = pd)
            fit=sjstats::eta_sq(fit)
            postEffSizes=c(postEffSizes,fit$etasq)
          }
          
          tmp1=data.frame(PC=1:.ArgList$nPCs,eta_sq=preEffSizes)
          tmp1$status="pre_Harmony"
          tmp1$variable=icol
          
          tmp2=data.frame(PC=1:.ArgList$nPCs,eta_sq=postEffSizes)
          tmp2$status="post_Harmony"
          tmp2$variable=icol
          resSummary=rbind(resSummary,tmp1)
          resSummary=rbind(resSummary,tmp2)
        }
        
      }
      save(resSummary,file=.myFilePathMakerFn("harmony_stat_summary",argList=.ArgList))
      
      #kmeans clustering
      set.seed(1)
      res_clust=kmeans(harmony_embeddings[,1:.ArgList$nPCs],.ArgList$internal_pseudocell_count,iter.max = 10000,algorithm = "Lloyd")
      
      #row.names(res_clust$centers)=paste0("C_",1:nrow(res_clust$centers))
      #kcenters=res_clust$centers
      #tormList=c()
      #while(length(tormList)<100){
      #  knet=RANN::nn2(data = kcenters, k = 20, eps = 0)
      #  knet_mutual=rowSums(knet$nn.dist)
      #  tormList=c(tormList,row.names(kcenters)[which(knet_mutual==min(knet_mutual))])
      #  kcenters=kcenters[!row.names(kcenters) %in% tormList,]
      #}
      #kcenters=res_clust$centers[!row.names(res_clust$centers) %in% tormList,]
      #knet=RANN::nn2(data=kcenters,query = harmony_embeddings,k=2,eps=0)
      
      #knet2=RANN::nn2(kcenters,k=10,eps=0)
      #knet2=data.frame(FC=row.names(kcenters),SC=row.names(kcenters)[knet2$nn.idx[,2]],distC=knet2$nn.dists[,2],stringsAsFactors = F)
      #knet3=data.frame(sample=row.names(harmony_embeddings),FC=row.names(kcenters)[knet$nn.idx[,1]],distSample=knet$nn.dists[,1],stringsAsFactors = F)
      #knet2=merge(knet2,knet3,by="FC",all=T)
      
      res_clusters=data.frame(sample=names(res_clust$cluster),cluster_id=res_clust$cluster,stringsAsFactors = F)
      save(res_clusters,file=.myFilePathMakerFn("kmeans_res_clusters",argList=.ArgList))
      
      pca_centroid=res_clust$centers
      save(pca_centroid,file=.myFilePathMakerFn("pca_centroids",argList=.ArgList))
      
      
    } else {
      
      harmony_embeddings <- liger_embeddings[,1:.ArgList$nPCs]
      .mySaveFn(harmony_embeddings,file=.myFilePathMakerFn("harmony-embeddings",argList=.ArgList,pseudoImportant = F))
      
      #kmeans clustering
      set.seed(1)
      doClustering=T
      itrClustering=0
      while(doClustering&itrClustering<10){
        itrClustering=itrClustering+1
        res_clust=kmeans(harmony_embeddings[,1:.ArgList$nPCs],.ArgList$internal_pseudocell_count,iter.max = 10000,algorithm = "Lloyd")
        if(sum(is.na(res_clust$centers))==0){
          doClustering=F
        }
      }
      
      if(sum(is.na(res_clust$centers))>0){
        stop("Error in identification of the pseudocells")
      }
      
      res_clusters=data.frame(sample=names(res_clust$cluster),cluster_id=res_clust$cluster,stringsAsFactors = F)
      save(res_clusters,file=.myFilePathMakerFn("kmeans_res_clusters",argList=.ArgList))
      
      pca_centroid=res_clust$centers
      save(pca_centroid,file=.myFilePathMakerFn("pca_centroids",argList=.ArgList))
      
    }
  }
  
  
  
  cat("Running UMAP\n")
  #rm(list = setdiff(ls(), lsf.str()))
  #load(.myFilePathMakerFn("harmony-embeddings",argList=.ArgList))
  reRunCheck=tryCatch({load(.myFilePathMakerFn("UMAP_anno",argList=.ArgList,pseudoImportant = F));F}, error=function(e) {return(T)})
  
  if(!reRunCheck){
    reRunCheck=tryCatch({load(.myFilePathMakerFn("UMAP_centroid",argList=.ArgList));F}, error=function(e) {return(T)})
  }
  
  if(reRunCheck){
    load(.myFilePathMakerFn("harmony-embeddings",argList=.ArgList))
    load(.myFilePathMakerFn("pca_anno",argList=.ArgList,pseudoImportant = F))
    
    load(.myFilePathMakerFn("pca_centroids",argList=.ArgList))
    
    pca_centroid=pca_centroid[,1:.ArgList$nPCs]
    
    
    
    tst=.reductionUMAPFn(harmony_embeddings,umap.method='umap-learn',testPCAembeddings=pca_centroid)
    resUMAP=tst$embedding
    .mySaveFn(resUMAP,file=.myFilePathMakerFn("UMAP_res",argList=.ArgList))
    
    x=as.data.frame(resUMAP)
    x=cbind(x,pd)
    row.names(x)=row.names(pd)
    pd=x
    .mySaveFn(pd,file=.myFilePathMakerFn("UMAP_anno",argList=.ArgList,pseudoImportant = F))
    
    tst=tst[["test"]]
    UMAP_centroid=data.frame(centroid=row.names(pca_centroid),UMAP_1=tst[,1],UMAP_2=tst[,2],stringsAsFactors = F)
    .mySaveFn(UMAP_centroid,file=.myFilePathMakerFn("UMAP_centroid",argList=.ArgList))
    
    
  }
  
  
  cat("Performing consensus marker analysis\n")
  cat("      Consensus marker step 1 ...\n")
  resDE1=.myConcensusDEFn_step1(argList=.ArgList)
  cat("      Consensus marker step 2 ...\n")
  resDE2=.myConcensusDEFn_step2(argList=.ArgList)
  #resDE2=.mySimulatedBenchMarkPropagationMethods(argList=.ArgList,saveIntermediate=F)
  
  print(resDE2)
}

.myShinyTransferFn=function(argList){
  save_dir_path=paste0("~/shinyData/",argList$prefix)
  if(dir.exists(save_dir_path)){
    unlink(save_dir_path,recursive = T)
  }
  dir.create(save_dir_path,recursive = T)
  file.copy(.myFilePathMakerFn("UMAP_centroid",argList=argList),save_dir_path)
  file.copy(.myFilePathMakerFn("res_DE_wZscore",argList=argList,uniformImportant=T),save_dir_path)
  file.copy(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F),save_dir_path)
  file.copy(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F),save_dir_path)
  file.copy(.myFilePathMakerFn("pca_centroids",argList=argList),save_dir_path)
  
  load(do.call('.myFilePathMakerFn',args=c(list(filename="exp_merged",expData=T),argList)))
  fd=tmp@assays$RNA@meta.features
  save(fd,file=file.path(save_dir_path,"fd.rda"))
  
  #library(googleCloudStorageR)
  the_files <- file.path(save_dir_path, list.files(path = save_dir_path, all.files = TRUE, recursive = TRUE))
  for(i in 1:length(the_files)){
    gcs_upload(the_files[i], name = paste0("vgazesta/",gsub("~/","",the_files[i])))
  }
  
  
  return("Done")
}

.myBenchMarkPropagationMethods=function(argList){
  
  .myClusterAnalysisFn(argList=argList,clustering_resolution=1.3)
    
  library(qs)
    set.seed(1)
    supportingFractionThr=argList$DE_supportingFractionThr
    n.adaptiveKernel=argList$DE_n.adaptiveKernel
    nPropIter=argList$DE_nPropIter
    
    #load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList))
    load(.myFilePathMakerFn("pca_centroids",argList=argList))
    
    if(!all(row.names(harmony_embeddings)==row.names(pd))){
      stop("Error in matching Names!")
    }
    
    pcaList=split(as.data.frame(harmony_embeddings),pd$anno_batch)
    
    load(.myFilePathMakerFn("exp_merged",argList=argList,expData=T))
    expData=SplitObject(tmp, split.by = "anno_batch")
    
    dataArranged=list()
    for(i in 1:length(expData)){
      dataArranged=c(dataArranged,list(list(logNormData=expData[[i]]@assays$RNA@data,countData=expData[[i]]@assays$RNA@counts,dsName=names(expData)[i],pcaData=pcaList[[names(expData)[i]]])))
    }
    
    #if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename="prop_mat_list_step2",uniformImportant=T),argList)))){
    res=list()
    cat("           Constructing the prop org v1 pro1...\n")
    if(!file.exists(.myFilePathMakerFn("prop_org_v1_p1",argList=argList,qsFormat = T))){
      tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop_org_v1,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,mc.cores = argList$ncores)
      qsave(tmp,file=.myFilePathMakerFn("prop_org_v1_p1",argList=argList,qsFormat = T))
    } else {
      tmp=qread(.myFilePathMakerFn("prop_org_v1_p1",argList=argList,qsFormat = T))
    }
    res=c(res,list(tmp))
    names(res)[length(res)]="prop_org_v1_p1"
    
    cat("           Constructing the prop org v1 pro2...\n")
    if(!file.exists(.myFilePathMakerFn("prop_org_v1_p2",argList=argList,qsFormat = T))){
      tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop_org_v1,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=2,mc.cores = argList$ncores)
      qsave(tmp,file=.myFilePathMakerFn("prop_org_v1_p2",argList=argList,qsFormat = T))
    } else {
      tmp=qread(.myFilePathMakerFn("prop_org_v1_p2",argList=argList,qsFormat = T))
    }
    res=c(res,list(tmp))
    names(res)[length(res)]="prop_org_v1_p2"
    
    cat("           Constructing the prop org v1 pro3...\n")
    if(!file.exists(.myFilePathMakerFn("prop_org_v1_p3",argList=argList,qsFormat = T))){
      tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop_org_v1,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,mc.cores = argList$ncores)
      qsave(tmp,file=.myFilePathMakerFn("prop_org_v1_p3",argList=argList,qsFormat = T))
    } else {
      tmp=qread(.myFilePathMakerFn("prop_org_v1_p3",argList=argList,qsFormat = T))
    }
    res=c(res,list(tmp))
    names(res)[length(res)]="prop_org_v1_p3"
    
    cat("           Constructing the prop org v2 pro3...\n")
    if(!file.exists(.myFilePathMakerFn("prop_org_v2_p3",argList=argList,qsFormat = T))){
      tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop_org_v2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,mc.cores = argList$ncores)
      qsave(tmp,file=.myFilePathMakerFn("prop_org_v2_p3",argList=argList,qsFormat = T))
    } else {
      tmp=qread(.myFilePathMakerFn("prop_org_v2_p3",argList=argList,qsFormat = T))
    }
    res=c(res,list(tmp))
    names(res)[length(res)]="prop_org_v2_p3"
    
    cat("           Constructing the DNN wo AdjPurity pro3...\n")
    if(!file.exists(.myFilePathMakerFn("prop_DNN_noAdjPurity_p3",argList=argList,qsFormat = T))){
      tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop_DNN,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,adjPurityIndx=F,mc.cores = argList$ncores)
      qsave(tmp,file=.myFilePathMakerFn("prop_DNN_noAdjPurity_p3",argList=argList,qsFormat = T))
    } else {
      tmp=qread(.myFilePathMakerFn("prop_DNN_noAdjPurity_p3",argList=argList,qsFormat = T))
    }
    res=c(res,list(tmp))
    names(res)[length(res)]="prop_DNN_noAdjPurity_p3"
    
    cat("           Constructing the DNN w AdjPurity pro3...\n")
    #tst=.myConcensusDEFn_step2_detail_prop_DNN(dataArranged[[1]],centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,adjPurityIndx=T)
    if(!file.exists(.myFilePathMakerFn("prop_DNN_AdjPurity_p3",argList=argList,qsFormat = T))){
      tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop_DNN,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,adjPurityIndx=T,mc.cores = argList$ncores)
      qsave(tmp,file=.myFilePathMakerFn("prop_DNN_AdjPurity_p3",argList=argList,qsFormat = T))
    } else {
      tmp=qread(.myFilePathMakerFn("prop_DNN_AdjPurity_p3",argList=argList,qsFormat = T))
    }
    res=c(res,list(tmp))
    names(res)[length(res)]="prop_DNN_AdjPurity_p3"
    
    cat("           Constructing the runIndx of 1 pro3...\n")
    if(!file.exists(.myFilePathMakerFn("prop_runIndx1_p3",argList=argList,qsFormat = T))){
      tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=1,mc.cores = argList$ncores)
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx1_p3",argList=argList,qsFormat = T))
    } else {
      tmp=qread(.myFilePathMakerFn("prop_runIndx1_p3",argList=argList,qsFormat = T))
    }
    res=c(res,list(tmp))
    names(res)[length(res)]="prop_runIndx1_p3"
    
    cat("           Constructing the runIndx of 2 pro3...\n")
    if(!file.exists(.myFilePathMakerFn("prop_runIndx2_p3",argList=argList,qsFormat = T))){
      tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=2,mc.cores = argList$ncores)
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx2_p3",argList=argList,qsFormat = T))
    } else {
      tmp=qread(.myFilePathMakerFn("prop_runIndx2_p3",argList=argList,qsFormat = T))
    }
    res=c(res,list(tmp))
    names(res)[length(res)]="prop_runIndx2_p3"
    
    cat("           Constructing the runIndx of 3 pro3...\n")
    if(!file.exists(.myFilePathMakerFn("prop_runIndx3_p3",argList=argList,qsFormat = T))){
      tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=3,mc.cores = argList$ncores)
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx3_p3",argList=argList,qsFormat = T))
    } else {
      tmp=qread(.myFilePathMakerFn("prop_runIndx3_p3",argList=argList,qsFormat = T))
    }
    res=c(res,list(tmp))
    names(res)[length(res)]="prop_runIndx3_p3"
    
    cat("           Constructing the runIndx of 4 pro3...\n")
    if(!file.exists(.myFilePathMakerFn("prop_runIndx4_p3",argList=argList,qsFormat = T))){
      tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=4,mc.cores = argList$ncores)
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx4_p3",argList=argList,qsFormat = T))
    } else {
      tmp=qread(.myFilePathMakerFn("prop_runIndx4_p3",argList=argList,qsFormat = T))
    }
    res=c(res,list(tmp))
    names(res)[length(res)]="prop_runIndx4_p3"
    
    
    cat("           Constructing the runIndx of 5 pro3...\n")
    if(!file.exists(.myFilePathMakerFn("prop_runIndx5_p3",argList=argList,qsFormat = T))){
      tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=5,mc.cores = argList$ncores)
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx5_p3",argList=argList,qsFormat = T))
    } else {
      tmp=qread(.myFilePathMakerFn("prop_runIndx5_p3",argList=argList,qsFormat = T))
    }
    res=c(res,list(tmp))
    names(res)[length(res)]="prop_runIndx5_p3"
    
    
    cat("           Constructing the runIndx of 6 pro3...\n")
    if(!file.exists(.myFilePathMakerFn("prop_runIndx6_p3",argList=argList,qsFormat = T))){
      tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=6,mc.cores = argList$ncores)
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx6_p3",argList=argList,qsFormat = T))
    } else {
      tmp=qread(.myFilePathMakerFn("prop_runIndx6_p3",argList=argList,qsFormat = T))
    }
    res=c(res,list(tmp))
    names(res)[length(res)]="prop_runIndx6_p3"
    
    cat("           Constructing the runIndx of 8 pro3...\n")
    if(!file.exists(.myFilePathMakerFn("prop_runIndx8_p3",argList=argList,qsFormat = T))){
      tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=8,mc.cores = argList$ncores)
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx8_p3",argList=argList,qsFormat = T))
    } else {
      tmp=qread(.myFilePathMakerFn("prop_runIndx8_p3",argList=argList,qsFormat = T))
    }
    res=c(res,list(tmp))
    names(res)[length(res)]="prop_runIndx8_p3"
    
    cat("           Constructing the runIndx of 9 pro3...\n")
    if(!file.exists(.myFilePathMakerFn("prop_runIndx9_p3",argList=argList,qsFormat = T))){
      tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=9,mc.cores = argList$ncores)
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx9_p3",argList=argList,qsFormat = T))
    } else {
      tmp=qread(.myFilePathMakerFn("prop_runIndx9_p3",argList=argList,qsFormat = T))
    }
    res=c(res,list(tmp))
    names(res)[length(res)]="prop_runIndx9_p3"
    
    
    res_arranged=list()
    
    for(i in 1:length(res)){
      tmp=res[[i]]
      for(j in 1:length(tmp)){
        tmp_ds=tmp[[j]]
        if(nrow(tmp_ds$prop_mat)>200){
          tmp_prop=tmp_ds$prop_mat[((ncol(tmp_ds$data$logNormData)+1):nrow(tmp_ds$prop_mat)),1:ncol(tmp_ds$data$logNormData)]
        } else {
          tmp_prop=tmp_ds$prop_mat
        }
        
        if(ncol(tmp_prop)>ncol(tmp_ds$data$logNormData)){
          tmp_prop=tmp_prop[,colnames(tmp_prop)%in%colnames(tmp_ds$data$logNormData)]
        }
        
        
        tmp_anno=pd[match(colnames(tmp_prop),row.names(pd)),c("anno_batch",colnames(pd)[grepl("^anno",colnames(pd))])]
        if(sum(is.na(tmp_anno$anno_batch))>0){
          stop("Error in propagation model",i)
        }
        
        #table(tmp_anno$anno_orig_cellState)
        
        if(length(unique(tmp_anno$anno_condition))<2|length(unique(tmp_anno$anno_condition))>5){
          tmp_anno$anno_condition=tmp_anno$anno_condition
        }
        
        .tmp_anno=tmp_anno
        .tmp_prop=tmp_prop
        
        if(sum(!is.na(.tmp_anno$anno_condition))>20 &length(unique(.tmp_anno$anno_condition[!is.na(.tmp_anno$anno_condition)]))>1){
          tmp_anno=.tmp_anno
          tmp_prop=.tmp_prop
          na_ind=which(is.na(tmp_anno$anno_condition))
          if(length(na_ind)>0){
            tmp_anno=tmp_anno[-na_ind,]
            tmp_prop=tmp_prop[,-na_ind]
          }
          
          tmp_cellState=.myOneHotFn(tmp_anno$anno_condition)
          tmp_cellState=tmp_prop %*% (as.matrix(tmp_cellState))
          
          tmp_cellState2=t(apply(tmp_cellState,1,function(x) x/sum(x)))
          tmp_cellStateBkg=matrix(table(tmp_anno$anno_condition)/sum(!is.na(tmp_anno$anno_condition)),nrow=nrow(tmp_cellState2),ncol=ncol(tmp_cellState2),byrow = T)
          tmp_cellState2=tmp_cellState2 - tmp_cellStateBkg
          tmp_cellState2=tmp_cellState2/(1-tmp_cellStateBkg)
          tmp_cellState2_ind_max=rep(NA,nrow(tmp_cellState2))
          for(itr in 1:nrow(tmp_cellState2)){
            if(!is.na(max(tmp_cellState2[itr,]))){
              if(!is.nan(max(tmp_cellState2[itr,]))){
                tmp_cellState2_ind_max[itr]=which(tmp_cellState2[itr,]==max(tmp_cellState2[itr,]))
              }
            }
            
          }
          tmp_cellState2_max_type=colnames(tmp_cellState2)[tmp_cellState2_ind_max]
          tmp_cellState2_max=apply(tmp_cellState2,1,max)
          tmp_res=data.frame(pseudocell=1:nrow(tmp_cellState2),score=tmp_cellState2_max,weigth=rowSums(tmp_prop),dataset=tmp_ds$data$dsName,attribute="anno_condition",method=names(res)[i],cellType=tmp_cellState2_max_type,stringsAsFactors = F)
          res_arranged=c(res_arranged,list(tmp_res))
        }
        
        if(sum(!is.na(.tmp_anno$anno_orig_cellState))>20){
          tmp_anno=.tmp_anno
          tmp_prop=.tmp_prop
          na_ind=which(is.na(tmp_anno$anno_orig_cellState))
          if(length(na_ind)>0){
            tmp_anno=tmp_anno[-na_ind,]
            tmp_prop=tmp_prop[,-na_ind]
          }
          
          na_ind=which(tmp_anno$anno_orig_cellState=="")
          if(length(na_ind)>0){
            tmp_anno=tmp_anno[-na_ind,]
            tmp_prop=tmp_prop[,-na_ind]
          }
          
          if(F){
            if(sum(grepl("Micro",tmp_anno$anno_orig_cellState))>0){
              tmp_anno$anno_orig_cellState[grepl("Micro",tmp_anno$anno_orig_cellState)]="Microglia"
            }
            
            if(sum(tmp_anno$anno_orig_cellState %in% c("MG","MGL","MGL1","MGL2","MGL3","mg","Mic","Micro"))>0){
              tmp_anno$anno_orig_cellState[tmp_anno$anno_orig_cellState %in% c("MG","MGL","MGL1","MGL2","MGL3","mg","Mic","Micro")]="Microglia"
            }
            
            if(length(which(tmp_anno$anno_orig_cellState %in% c("H1/2M","H1M","H2M","homeostatic mic_s2_c1","homeostatic mic_s2_c2")))>0){
              tmp_anno$anno_orig_cellState[which(tmp_anno$anno_orig_cellState %in% c("H1/2M","H1M","H2M","homeostatic mic_s2_c1","homeostatic mic_s2_c2"))]="H1/2M"
            }
            
            if(length(which(tmp_anno$anno_orig_cellState %in% c("cDC1","cDC2","pDC","migDC")))>0){
              tmp_anno$anno_orig_cellState[which(tmp_anno$anno_orig_cellState %in% c("cDC1","cDC2","pDC","migDC"))]="cDC"
            }
            
            if(length(which(tmp_anno$anno_orig_cellState %in% c("T/NKT cells","T cells","NK cells")))>0){
              tmp_anno$anno_orig_cellState[which(tmp_anno$anno_orig_cellState %in% c("T/NKT cells","T cells","NK cells"))]="T/NKT cells"
            }
            
            if(length(which(tmp_anno$anno_orig_cellState %in% c("Mnc","Non-cl. monocytes","MCs")))>0){
              tmp_anno$anno_orig_cellState[which(tmp_anno$anno_orig_cellState %in% c("Mnc","Non-cl. monocytes","MCs"))]="Monocyte/Mdc"
            }
            
            if(length(which(tmp_anno$anno_orig_cellState %in% c("CPM","CRM")))>0){
              tmp_anno$anno_orig_cellState[which(tmp_anno$anno_orig_cellState %in% c("CPM","CRM"))]="cycling_prolif_microglia"
            }
            
            
            na_ind=which(tmp_anno$anno_orig_cellState %in% c("Microglia"))
            if(length(na_ind)>0){
              tmp_anno=tmp_anno[-na_ind,]
              tmp_prop=tmp_prop[,-na_ind]
            }
          }
          

          if(length(unique(tmp_anno$anno_orig_cellState))>1){
            tmp_cellState=.myOneHotFn(tmp_anno$anno_orig_cellState)
            tmp_cellState=tmp_prop %*% (as.matrix(tmp_cellState))
            
            tmp_cellState2=t(apply(tmp_cellState,1,function(x) x/sum(x)))
            tmp_cellStateBkg=matrix(table(tmp_anno$anno_orig_cellState)/sum(!is.na(tmp_anno$anno_orig_cellState)),nrow=nrow(tmp_cellState2),ncol=ncol(tmp_cellState2),byrow = T)
            tmp_cellState2=tmp_cellState2 - tmp_cellStateBkg
            tmp_cellState2=tmp_cellState2/(1-tmp_cellStateBkg)
            tmp_cellState2_ind_max=rep(NA,nrow(tmp_cellState2))
            for(itr in 1:nrow(tmp_cellState2)){
              if(!is.na(max(tmp_cellState2[itr,]))){
                if(!is.nan(max(tmp_cellState2[itr,]))){
                  tmp_cellState2_ind_max[itr]=which(tmp_cellState2[itr,]==max(tmp_cellState2[itr,]))
                }
              }
              
            }
            tmp_cellState2_max_type=colnames(tmp_cellState2)[tmp_cellState2_ind_max]
            tmp_cellState2_max=tmp_cellState2[cbind(1:nrow(tmp_cellState2),tmp_cellState2_ind_max)]
            tmp_res=data.frame(pseudocell=1:nrow(tmp_cellState2),score=tmp_cellState2_max,weigth=rowSums(tmp_prop),dataset=tmp_ds$data$dsName,attribute="anno_orig_cellState",method=names(res)[i],cellType=tmp_cellState2_max_type,stringsAsFactors = F)
            res_arranged=c(res_arranged,list(tmp_res))
          }
          
        }
        
        
        if(sum(!is.na(.tmp_anno$anno_cluster_res))>20){
          tmp_anno=.tmp_anno
          tmp_prop=.tmp_prop
          na_ind=which(is.na(tmp_anno$anno_cluster_res))
          if(length(na_ind)>0){
            tmp_anno=tmp_anno[-na_ind,]
            tmp_prop=tmp_prop[,-na_ind]
          }
          
          na_ind=which(tmp_anno$anno_cluster_res=="")
          if(length(na_ind)>0){
            tmp_anno=tmp_anno[-na_ind,]
            tmp_prop=tmp_prop[,-na_ind]
          }
          
          
          if(length(unique(tmp_anno$anno_cluster_res))>1){
            tmp_cellState=.myOneHotFn(tmp_anno$anno_cluster_res)
            tmp_cellState=tmp_prop %*% (as.matrix(tmp_cellState))
            
            tmp_cellState2=t(apply(tmp_cellState,1,function(x) x/sum(x)))
            tmp_cellStateBkg=matrix(table(tmp_anno$anno_cluster_res)/sum(!is.na(tmp_anno$anno_cluster_res)),nrow=nrow(tmp_cellState2),ncol=ncol(tmp_cellState2),byrow = T)
            tmp_cellState2=tmp_cellState2 - tmp_cellStateBkg
            tmp_cellState2=tmp_cellState2/(1-tmp_cellStateBkg)
            tmp_cellState2_ind_max=rep(NA,nrow(tmp_cellState2))
            for(itr in 1:nrow(tmp_cellState2)){
              if(!is.na(max(tmp_cellState2[itr,]))){
                if(!is.nan(max(tmp_cellState2[itr,]))){
                  tmp_cellState2_ind_max[itr]=which(tmp_cellState2[itr,]==max(tmp_cellState2[itr,]))
                }
              }
              
            }
            tmp_cellState2_max_type=colnames(tmp_cellState2)[tmp_cellState2_ind_max]
            tmp_cellState2_max=tmp_cellState2[cbind(1:nrow(tmp_cellState2),tmp_cellState2_ind_max)]
            tmp_res=data.frame(pseudocell=1:nrow(tmp_cellState2),score=tmp_cellState2_max,weigth=rowSums(tmp_prop),dataset=tmp_ds$data$dsName,attribute="anno_cluster_res",method=names(res)[i],cellType=tmp_cellState2_max_type,stringsAsFactors = F)
            res_arranged=c(res_arranged,list(tmp_res))
          }
          
        }
        
      }
    }
    
    res_arranged=do.call("rbind",res_arranged)
    
    
    return(list(res=res,analysis=res_arranged))
  
  
}

.mySimulatedBenchMarkPropagationMethods=function(argList,saveIntermediate=F,returnDatasetObject=F,path=F){
  
  #argList=.ArgList;saveIntermediate=F;returnDatasetObject=T;path=F
  library(qs)
  set.seed(1)
  supportingFractionThr=argList$DE_supportingFractionThr
  n.adaptiveKernel=argList$DE_n.adaptiveKernel
  nPropIter=argList$DE_nPropIter
  
  #load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  load(.myFilePathMakerFn("harmony-embeddings",argList=argList))
  load(.myFilePathMakerFn("pca_centroids",argList=argList))
  
  if(path){
    pd$anno_orig_cellState="s4"
    pd$anno_orig_cellState[pd$Step<75]="s3"
    pd$anno_orig_cellState[pd$Step<50]="s2"
    pd$anno_orig_cellState[pd$Step<20]="s1"
  }
  
  if(!all(row.names(harmony_embeddings)==row.names(pd))){
    stop("Error in matching Names!")
  }
  
  pcaList=split(as.data.frame(harmony_embeddings),pd$anno_batch)
  
  load(.myFilePathMakerFn("exp_merged",argList=argList,expData=T))
  expData=SplitObject(tmp, split.by = "anno_batch")
  pcaList=pcaList[names(expData)]
  dataArranged=list()
  for(i in 1:length(expData)){
    dataArranged=c(dataArranged,list(list(logNormData=expData[[i]]@assays$RNA@data,countData=expData[[i]]@assays$RNA@counts,dsName=names(expData)[i],pcaData=pcaList[[names(expData)[i]]])))
  }
  
  #if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename="prop_mat_list_step2",uniformImportant=T),argList)))){
  res=list()
  cat("           Constructing the prop org v1 pro1...\n")
  if(!file.exists(.myFilePathMakerFn("prop_org_v1_p1",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop_org_v1,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_org_v1_p1",argList=argList,qsFormat = T))
    }
    
  } else {
    tmp=qread(.myFilePathMakerFn("prop_org_v1_p1",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_org_v1_p1"
  
  cat("           Constructing the prop org v1 pro2...\n")
  if(!file.exists(.myFilePathMakerFn("prop_org_v1_p2",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop_org_v1,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_org_v1_p2",argList=argList,qsFormat = T))
    }
    
  } else {
    tmp=qread(.myFilePathMakerFn("prop_org_v1_p2",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_org_v1_p2"
  
  cat("           Constructing the prop org v1 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_org_v1_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop_org_v1,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_org_v1_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_org_v1_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_org_v1_p3"
  
  cat("           Constructing the prop org v2 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_org_v2_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop_org_v2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_org_v2_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_org_v2_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_org_v2_p3"
  
  cat("           Constructing the DNN wo AdjPurity pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_DNN_noAdjPurity_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop_DNN,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,adjPurityIndx=F,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_DNN_noAdjPurity_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_DNN_noAdjPurity_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_DNN_noAdjPurity_p3"
  
  cat("           Constructing the DNN w AdjPurity pro3...\n")
  #tst=.myConcensusDEFn_step2_detail_prop_DNN(dataArranged[[1]],centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,adjPurityIndx=T)
  if(!file.exists(.myFilePathMakerFn("prop_DNN_AdjPurity_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop_DNN,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,adjPurityIndx=T,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_DNN_AdjPurity_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_DNN_AdjPurity_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_DNN_AdjPurity_p3"
  
  cat("           Constructing the runIndx of 1 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx1_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=1,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx1_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx1_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx1_p3"
  
  cat("           Constructing the runIndx of 1 pro2...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx1_p2",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=2,runIndx=1,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx1_p2",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx1_p2",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx1_p2"
  
  cat("           Constructing the runIndx of 1 pro4...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx1_p4",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=4,runIndx=1,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx1_p4",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx1_p4",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx1_p4"
  
  
  cat("           Constructing the runIndx of 1 pro6...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx1_p6",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=6,runIndx=1,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx1_p6",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx1_p6",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx1_p6"
  
  cat("           Constructing the runIndx of 2 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx2_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx2_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx2_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx2_p3"
  
  cat("           Constructing the runIndx of 3 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx3_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx3_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx3_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx3_p3"
  
  cat("           Constructing the runIndx of 4 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx4_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=4,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx4_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx4_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx4_p3"
  
  
  cat("           Constructing the runIndx of 5 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx5_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=5,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx5_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx5_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx5_p3"
  
  
  cat("           Constructing the runIndx of 6 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx6_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=6,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx6_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx6_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx6_p3"
  
  cat("           Constructing the runIndx of 7 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx7_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=7,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx7_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx7_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx7_p3"
  
  cat("           Constructing the runIndx of 8 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx8_p3",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=8,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx8_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx8_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx8_p3"
  
  cat("           Constructing the runIndx of 9 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx9_p3",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=9,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx9_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx9_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx9_p3"
  
  cat("           Constructing the newProp n1 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newpropn1_p3",argList=argList,qsFormat = T))){
    #centroidPCAdata=pca_centroid;argList=argList;n.adaptiveKernel=n.adaptiveKernel;nPropIter=3;n.neighbors=1
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,n.neighbors=1,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newpropn1_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newpropn1_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newpropn1_p3"
  
  cat("           Constructing the newProp n2 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newpropn2_p3",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,n.neighbors=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newpropn2_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newpropn2_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newpropn2_p3"
  
  cat("           Constructing the newProp n2 pro1...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newpropn2_p1",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newpropn2_p1",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newpropn2_p1",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newpropn2_p1"
  
  cat("           Constructing the newProp n2 pro2...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newpropn2_p2",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=2,n.neighbors=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newpropn2_p2",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newpropn2_p2",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newpropn2_p2"
  
  cat("           Constructing the newProp n2 pro4...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newpropn2_p4",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=4,n.neighbors=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newpropn2_p4",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newpropn2_p4",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newpropn2_p4"
  
  cat("           Constructing the newProp n3 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newpropn3_p3",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,n.neighbors=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newpropn3_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newpropn3_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newpropn3_p3"
  
  cat("           Constructing the newProp2 n3 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop2n3_p3",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,n.neighbors=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop2n3_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop2n3_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop2n3_p3"
  
  cat("           Constructing the newProp2 n3 pro6...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop2n3_p6",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=6,n.neighbors=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop2n3_p6",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop2n3_p6",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop2n3_p6"
  
  cat("           Constructing the newProp2 n2 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop2n2_p3",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,n.neighbors=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop2n2_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop2n2_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop2n2_p3"
  
  cat("           Constructing the newProp3 n2 pro3_3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n2_p3_3",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,n.neighbors=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n2_p3_3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n2_p3_3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n2_p3_3"
  
  cat("           Constructing the newProp3 n2 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n2_p3",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n2_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n2_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n2_p3"
  
  cat("           Constructing the newProp3 n4 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n4_p3",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=4,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n4_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n4_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n4_p3"
  #myConcensusDEFn_step2_detail_newprop3_cosine
  
  cat("           Constructing the newProp3 n3 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n3_p3",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n3_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n3_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n3_p3"
  
  #myConcensusDEFn_step2_detail_newprop3_final
  cat("           Constructing the newProp3-final n4 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n4_p3_final",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_final,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=4,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n4_p3_final",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n4_p3_final",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n4_p3_final"
  
  ##################
  #myConcensusDEFn_step2_detail_newprop3_final
  cat("           Constructing the newProp3-final_v2 n4 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n4_p3_final_v2",argList=argList,qsFormat = T))){
    #centroidPCAdata=pca_centroid;argList=argList;n.adaptiveKernel=n.adaptiveKernel;nPropIter=1;n.neighbors=4
    tmp=.myConcensusDEFn_step2_detail_newprop3_final_v2(dataArranged,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=4)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n4_p3_final_v2",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n4_p3_final_v2",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n4_p3_final_v2"
  #myConcensusDEFn_step2_detail_newprop3_cosine
  
  cat("           Constructing the newProp3-final_v2m n4 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n4_p3_final_v2m",argList=argList,qsFormat = T))){
    #centroidPCAdata=pca_centroid;argList=argList;n.adaptiveKernel=n.adaptiveKernel;nPropIter=1;n.neighbors=4
    tmp=.myConcensusDEFn_step2_detail_newprop3_final_v2m(dataArranged,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=4)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n4_p3_final_v2m",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n4_p3_final_v2m",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n4_p3_final_v2m"
  #myConcensusDEFn_step2_detail_newprop3_cosine
  
  #
  
  
  cat("           Constructing the newProp3-final n3 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n3_p3_final",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_final,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n3_p3_final",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n3_p3_final",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n3_p3_final"
  
  cat("           Constructing the newProp3 n4 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n4_p3",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=4,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n4_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n4_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n4_p3"
  #myConcensusDEFn_step2_detail_newprop3_cosine
  
  cat("           Constructing the newProp3 n3 pro3v2...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n3_p3v2",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3v2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n3_p3v2",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n3_p3v2",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n3_p3v2"
  
  cat("           Constructing the newProp3 n4 pro3v2...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n4_p3v2",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3v2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=4,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n4_p3v2",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n4_p3v2",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n4_p3v2"
  
  cat("           Constructing the newProp3 n2 pro3_cosine...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n2_p3_cosine",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_cosine,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n2_p3_cosine",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n2_p3_cosine",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n2_p3_cosine"
  
  cat("           Constructing the newProp3 n4 pro3_cosine...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n4_p3_cosine",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_cosine,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=4,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n4_p3_cosine",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n4_p3_cosine",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n4_p3_cosine"
  
  cat("           Constructing the newProp3 n6 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n6_p3",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=6,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n6_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n6_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n6_p3"
  
  
  cat("           Constructing the newProp3 n8 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n8_p3",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=8,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n8_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n8_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n8_p3"
  
  cat("           Constructing the newProp4 n4 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop4n4_p3",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop4,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=4,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop4n4_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop4n4_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop4n4_p3"
  
  cat("           Constructing the newProp4 n3 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop4n3_p3",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop4,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop4n3_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop4n3_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop4n3_p3"
  
  cat("           Constructing the newProp4 n6 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop4n6_p3",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop4,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=6,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop4n6_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop4n6_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop4n6_p3"
  
  
  #myConcensusDEFn_step2_detail_newprop3_scale_v2
  cat("           Constructing the newProp3scalev2 n4 pro1...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3scalev2n4_p1",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_scale_v2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=4,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3scalev2n4_p1",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3scalev2n4_p1",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3scalev2n4_p1"
  
  cat("           Constructing the newProp3scalev2 n3 pro1...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3scalev2n3_p1",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_scale_v2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3scalev2n3_p1",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3scalev2n3_p1",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3scalev2n3_p1"
  
  cat("           Constructing the newProp3scalev2 n3 p2...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3scalev2n3_p2",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_scale_v2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=2,n.neighbors=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3scalev2n3_p2",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3scalev2n3_p2",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3scalev2n3_p2"
  
  cat("           Constructing the newProp3scale n4 pro1...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3scalen4_p1",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_scale,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=4,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3scalen4_p1",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3scalen4_p1",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3scalen4_p1"
  
  cat("           Constructing the newProp3scale n3 pro1...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3scalen3_p1",argList=argList,qsFormat = T))){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_scale,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=3,mc.cores = argList$ncores)
    qsave(tmp,file=.myFilePathMakerFn("prop_newprop3scalen3_p1",argList=argList,qsFormat = T))
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3scalen3_p1",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3scalen3_p1"
  
  res_arranged=list()
  
  for(icounter in 1:length(res)){
    tmp=res[[icounter]]
    
    for(j in 1:length(tmp)){
      tmp_ds=tmp[[j]]
      if(nrow(tmp_ds$prop_mat)>200){
        tmp_prop=tmp_ds$prop_mat[((ncol(tmp_ds$data$logNormData)+1):nrow(tmp_ds$prop_mat)),1:ncol(tmp_ds$data$logNormData)]
      } else {
        tmp_prop=tmp_ds$prop_mat
      }
      
      if(ncol(tmp_prop)>ncol(tmp_ds$data$logNormData)){
        tmp_prop=tmp_prop[,colnames(tmp_prop)%in%colnames(tmp_ds$data$logNormData)]
      }
      
      
      tmp_anno=pd[match(colnames(tmp_prop),row.names(pd)),c("anno_batch",colnames(pd)[grepl("^anno",colnames(pd))])]
      if(sum(is.na(tmp_anno$anno_batch))>0){
        stop("Error in propagation model",names(res)[icounter])
      }
      
      #table(tmp_anno$anno_orig_cellState)
      
      if(length(unique(tmp_anno$anno_condition))<2|length(unique(tmp_anno$anno_condition))>5){
        tmp_anno$anno_condition=tmp_anno$anno_condition
      }
      
      .tmp_anno=tmp_anno
      .tmp_prop=tmp_prop
      
      if(sum(!is.na(.tmp_anno$anno_orig_cellState))>20){
        tmp_anno=.tmp_anno
        tmp_prop=.tmp_prop
        na_ind=which(is.na(tmp_anno$anno_orig_cellState))
        if(length(na_ind)>0){
          tmp_anno=tmp_anno[-na_ind,]
          tmp_prop=tmp_prop[,-na_ind]
        }
        
        na_ind=which(tmp_anno$anno_orig_cellState=="")
        if(length(na_ind)>0){
          tmp_anno=tmp_anno[-na_ind,]
          tmp_prop=tmp_prop[,-na_ind]
        }
        
        if(length(unique(tmp_anno$anno_orig_cellState))>1){
          tmp_cellState=.myOneHotFn(tmp_anno$anno_orig_cellState)
          tmp_cellState=tmp_prop %*% (as.matrix(tmp_cellState))
          
          effective_sample_size=.myEffSizePropMat(tmp_prop)
          effective_sample_size=effective_sample_size$effective_sample_size
        
          tmp_cellState2=t(apply(tmp_cellState,1,function(x) x/sum(x)))
          tmp_cellStateBkg=matrix(table(tmp_anno$anno_orig_cellState)/sum(!is.na(tmp_anno$anno_orig_cellState)),nrow=nrow(tmp_cellState2),ncol=ncol(tmp_cellState2),byrow = T)
          tmp_cellState2=tmp_cellState2 - tmp_cellStateBkg
          tmp_cellState2=tmp_cellState2/(1-tmp_cellStateBkg)
          tmp_cellState2_ind_max=rep(NA,nrow(tmp_cellState2))
          for(itr in 1:nrow(tmp_cellState2)){
            if(!is.na(max(tmp_cellState2[itr,]))){
              if(!is.nan(max(tmp_cellState2[itr,]))){
                tmp_cellState2_ind_max[itr]=which(tmp_cellState2[itr,]==max(tmp_cellState2[itr,]))
              }
            }
            
          }
          tmp_cellState2_max_type=colnames(tmp_cellState2)[tmp_cellState2_ind_max]
          tmp_cellState2_max=tmp_cellState2[cbind(1:nrow(tmp_cellState2),tmp_cellState2_ind_max)]
          tmp_res=data.frame(pseudocell=1:nrow(tmp_cellState2),score=tmp_cellState2_max,weight=rowSums(tmp_prop),dataset=tmp_ds$data$dsName,attribute="anno_orig_cellState",method=names(res)[icounter],cellType=tmp_cellState2_max_type,effective_sample_size=effective_sample_size,stringsAsFactors = F)
          res_arranged=c(res_arranged,list(tmp_res))
        }
        
      }
     
    }
  }
  
  res_arranged=do.call("rbind",res_arranged)
  
  library(qs)
  if(saveIntermediate){
    qsave(res_arranged,file=.myFilePathMakerFn("sim_prop_valuation",argList=argList,qsFormat = T))
  }
  
  
  res_summary=NULL
  for(igroup in unique(res_arranged$cellType)){
    tmp_main=res_arranged[which(res_arranged$cellType==igroup),]
    for(imethod in unique(tmp_main$method)){
      tmp_test=tmp_main[which(tmp_main$method==imethod),]
      tmp_bkg=tmp_main[which(tmp_main$method!=imethod),]
      
      all=effsize::cohen.d(tmp_test$score,tmp_bkg$score)$estimate
      top0.5=effsize::cohen.d(tmp_test$score[which(tmp_test$weight>0.5)],tmp_bkg$score[which(tmp_test$weight>0.5)])$estimate
      top0.8=effsize::cohen.d(tmp_test$score[which(tmp_test$weight>0.8)],tmp_bkg$score[which(tmp_test$weight>0.8)])$estimate
      
      sample_all=effsize::cohen.d(tmp_test$effective_sample_size,tmp_bkg$effective_sample_size)$estimate
      sample_top0.5=effsize::cohen.d(tmp_test$effective_sample_size[which(tmp_test$weight>0.5)],tmp_bkg$effective_sample_size[which(tmp_test$weight>0.5)])$estimate
      sample_top0.8=effsize::cohen.d(tmp_test$effective_sample_size[which(tmp_test$weight>0.8)],tmp_bkg$effective_sample_size[which(tmp_test$weight>0.8)])$estimate
      
      res_summary=rbind(res_summary,data.frame(method=imethod,group=igroup,median_score=median(tmp_test$score),median_score.5=median(tmp_test$score[which(tmp_test$weight>0.5)]),median_score.8=median(tmp_test$score[which(tmp_test$weight>0.8)]),effSize=all,effSize.5=top0.5,effSize.8=top0.8,median_sample_size=sample_all,median_sample_size.5=sample_top0.5,median_sample_size.8=sample_top0.8,stringsAsFactors=F))
    }
  }
  
  if(T){
    qsave(res_summary,file=.myFilePathMakerFn("sim_prop_summary_valuation",argList=argList,qsFormat = T))
  }
  
  
  if(returnDatasetObject){
    return(list(res_summary=res_summary,res=res_arranged))
  } else {
    return(res_summary)
  }
}

.mySimulatedBenchMarkPropagationMethods_excludeGroup=function(argList,returnDatasetObject=F,path=F,num_excluded_groups=10,...){
  
  saveIntermediate=F
  #argList=.ArgList;saveIntermediate=T
  library(qs)
  set.seed(1)
  supportingFractionThr=argList$DE_supportingFractionThr
  n.adaptiveKernel=argList$DE_n.adaptiveKernel
  nPropIter=argList$DE_nPropIter
  
  #load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  load(.myFilePathMakerFn("harmony-embeddings",argList=argList))
  load(.myFilePathMakerFn("pca_centroids",argList=argList))
  
  if(path){
    pd$anno_orig_cellState="s4"
    pd$anno_orig_cellState[pd$Step<75]="s3"
    pd$anno_orig_cellState[pd$Step<50]="s2"
    pd$anno_orig_cellState[pd$Step<20]="s1"
  }
  
  if(!all(row.names(harmony_embeddings)==row.names(pd))){
    stop("Error in matching Names!")
  }
  
  pd_org=pd
  for(icounter in 1:num_excluded_groups){
    r_sample=sample(unique(pd$anno_batch),1)
    r_group=as.character(sample(unique(pd$anno_orig_cellState[pd$anno_batch==r_sample]),1))
    pd=pd[-which(pd$anno_orig_cellState==r_group&pd$anno_batch==r_sample),]
  }
  harmony_embeddings=harmony_embeddings[match(row.names(pd),row.names(harmony_embeddings)),]
  load(.myFilePathMakerFn("exp_merged",argList=argList,expData=T))
  tmp=tmp[,match(row.names(pd),colnames(tmp))]
  
  
  pcaList=split(as.data.frame(harmony_embeddings),pd$anno_batch)
  
  
  {
    fd=tmp@assays$RNA@meta.features
    fd$dePattern="none"
    dwn=apply(fd[,grepl("DEFacGroup",colnames(fd))],1,function(x) sum(x<1))
    fd$dePattern[which(dwn>0)]="down_others"
    fd$dePattern[which(fd$DEFacGroup1>1)]="Group1"
    fd$dePattern[which(fd$DEFacGroup2>1)]="Group2"
    fd$dePattern[which(fd$DEFacGroup3>1)]="Group3"
    fd$dePattern[which(fd$DEFacGroup4>1)]="Group4"
    fd$dePattern[which(fd$DEFacGroup5>1)]="Group5"
    upmix=apply(fd[,grepl("DEFacGroup",colnames(fd))],1,function(x) sum(x>1))
    fd$dePattern[which(upmix>1)]="upmix"
    fdPattern=fd
  }
  
  expData=SplitObject(tmp, split.by = "anno_batch")
  pcaList=pcaList[names(expData)]
  dataArranged=list()
  for(i in 1:length(expData)){
    dataArranged=c(dataArranged,list(list(logNormData=expData[[i]]@assays$RNA@data,countData=expData[[i]]@assays$RNA@counts,dsName=names(expData)[i],pcaData=pcaList[[names(expData)[i]]])))
  }
  
  #if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename="prop_mat_list_step2",uniformImportant=T),argList)))){
  res=list()
  cat("           Constructing the prop org v1 pro1...\n")
  if(!file.exists(.myFilePathMakerFn("prop_org_v1_p1",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop_org_v1,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_org_v1_p1",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_org_v1_p1",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_org_v1_p1"
  
  cat("           Constructing the prop org v1 pro2...\n")
  if(!file.exists(.myFilePathMakerFn("prop_org_v1_p2",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop_org_v1,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_org_v1_p2",argList=argList,qsFormat = T))
    }
    
  } else {
    tmp=qread(.myFilePathMakerFn("prop_org_v1_p2",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_org_v1_p2"
  
  cat("           Constructing the prop org v1 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_org_v1_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop_org_v1,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_org_v1_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_org_v1_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_org_v1_p3"
  
  cat("           Constructing the prop org v2 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_org_v2_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop_org_v2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_org_v2_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_org_v2_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_org_v2_p3"
  
  cat("           Constructing the DNN wo AdjPurity pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_DNN_noAdjPurity_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop_DNN,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,adjPurityIndx=F,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_DNN_noAdjPurity_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_DNN_noAdjPurity_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_DNN_noAdjPurity_p3"
  
  cat("           Constructing the DNN w AdjPurity pro3...\n")
  #tst=.myConcensusDEFn_step2_detail_prop_DNN(dataArranged[[1]],centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,adjPurityIndx=T)
  if(!file.exists(.myFilePathMakerFn("prop_DNN_AdjPurity_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop_DNN,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,adjPurityIndx=T,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_DNN_AdjPurity_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_DNN_AdjPurity_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_DNN_AdjPurity_p3"
  
  cat("           Constructing the runIndx of 1 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx1_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=1,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx1_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx1_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx1_p3"
  
  cat("           Constructing the runIndx of 1 pro2...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx1_p2",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=2,runIndx=1,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx1_p2",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx1_p2",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx1_p2"
  
  cat("           Constructing the runIndx of 1 pro4...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx1_p4",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=4,runIndx=1,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx1_p4",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx1_p4",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx1_p4"
  
  
  cat("           Constructing the runIndx of 1 pro6...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx1_p6",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=6,runIndx=1,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx1_p6",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx1_p6",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx1_p6"
  
  cat("           Constructing the runIndx of 2 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx2_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx2_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx2_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx2_p3"
  
  cat("           Constructing the runIndx of 3 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx3_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx3_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx3_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx3_p3"
  
  cat("           Constructing the runIndx of 4 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx4_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=4,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx4_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx4_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx4_p3"
  
  
  cat("           Constructing the runIndx of 5 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx5_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=5,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx5_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx5_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx5_p3"
  
  
  cat("           Constructing the runIndx of 6 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx6_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=6,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx6_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx6_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx6_p3"
  
  cat("           Constructing the runIndx of 7 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx7_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=7,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx7_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx7_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx7_p3"
  
  cat("           Constructing the runIndx of 8 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx8_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=8,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx8_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx8_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx8_p3"
  
  cat("           Constructing the runIndx of 9 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_runIndx9_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,runIndx=9,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_runIndx9_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_runIndx9_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_runIndx9_p3"
  
  ###################
  cat("           Constructing the newProp n1 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newpropn1_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,n.neighbors=1,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newpropn1_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newpropn1_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newpropn1_p3"
  
  cat("           Constructing the newProp n2 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newpropn2_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,n.neighbors=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newpropn2_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newpropn2_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newpropn2_p3"
  
  cat("           Constructing the newProp n2 pro1...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newpropn2_p1",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newpropn2_p1",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newpropn2_p1",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newpropn2_p1"
  
  cat("           Constructing the newProp n2 pro2...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newpropn2_p2",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=2,n.neighbors=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newpropn2_p2",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newpropn2_p2",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newpropn2_p2"
  
  cat("           Constructing the newProp n2 pro4...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newpropn2_p4",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=4,n.neighbors=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newpropn2_p4",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newpropn2_p4",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newpropn2_p4"
  
  cat("           Constructing the newProp n3 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newpropn3_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,n.neighbors=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newpropn3_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newpropn3_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newpropn3_p3"
  
  cat("           Constructing the newProp2 n3 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop2n3_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,n.neighbors=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop2n3_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop2n3_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop2n3_p3"
  
  cat("           Constructing the newProp2 n3 pro6...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop2n3_p6",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=6,n.neighbors=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop2n3_p6",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop2n3_p6",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop2n3_p6"
  
  cat("           Constructing the newProp2 n2 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop2n2_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,n.neighbors=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop2n2_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop2n2_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop2n2_p3"
  
  cat("           Constructing the newProp3 n2 pro3_3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n2_p3_3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=3,n.neighbors=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n2_p3_3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n2_p3_3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n2_p3_3"
  
  cat("           Constructing the newProp3 n2 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n2_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n2_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n2_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n2_p3"
  
  cat("           Constructing the newProp3 n4 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n4_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=4,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n4_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n4_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n4_p3"
  #myConcensusDEFn_step2_detail_newprop3_cosine
  
  cat("           Constructing the newProp3 n3 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n3_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n3_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n3_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n3_p3"
  
  #myConcensusDEFn_step2_detail_newprop3_final
  cat("           Constructing the newProp3-final n4 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n4_p3_final",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_final,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=4,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n4_p3_final",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n4_p3_final",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n4_p3_final"
  #myConcensusDEFn_step2_detail_newprop3_cosine
  
  cat("           Constructing the newProp3-final n3 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n3_p3_final",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_final,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n3_p3_final",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n3_p3_final",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n3_p3_final"
  
  cat("           Constructing the newProp3 n4 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n4_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=4,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n4_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n4_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n4_p3"
  #myConcensusDEFn_step2_detail_newprop3_cosine
  
  cat("           Constructing the newProp3 n3 pro3v2...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n3_p3v2",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3v2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n3_p3v2",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n3_p3v2",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n3_p3v2"
  
  cat("           Constructing the newProp3 n4 pro3v2...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n4_p3v2",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3v2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=4,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n4_p3v2",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n4_p3v2",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n4_p3v2"
  
  cat("           Constructing the newProp3 n2 pro3_cosine...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n2_p3_cosine",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_cosine,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=2,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n2_p3_cosine",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n2_p3_cosine",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n2_p3_cosine"
  
  cat("           Constructing the newProp3 n4 pro3_cosine...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n4_p3_cosine",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_cosine,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=4,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n4_p3_cosine",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n4_p3_cosine",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n4_p3_cosine"
  
  cat("           Constructing the newProp3 n6 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n6_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=6,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n6_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n6_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n6_p3"
  
  
  cat("           Constructing the newProp3 n8 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3n8_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=8,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3n8_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3n8_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3n8_p3"
  
  cat("           Constructing the newProp4 n4 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop4n4_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop4,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=4,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop4n4_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop4n4_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop4n4_p3"
  
  cat("           Constructing the newProp4 n3 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop4n3_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop4,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop4n3_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop4n3_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop4n3_p3"
  
  cat("           Constructing the newProp4 n6 pro3...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop4n6_p3",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop4,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=6,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop4n6_p3",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop4n6_p3",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop4n6_p3"
  
  
  #myConcensusDEFn_step2_detail_newprop3_scale_v2
  cat("           Constructing the newProp3scalev2 n4 pro1...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3scalev2n4_p1",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_scale_v2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=4,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3scalev2n4_p1",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3scalev2n4_p1",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3scalev2n4_p1"
  
  cat("           Constructing the newProp3scalev2 n3 pro1...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3scalev2n3_p1",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_scale_v2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3scalev2n3_p1",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3scalev2n3_p1",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3scalev2n3_p1"
  
  cat("           Constructing the newProp3scalev2 n3 p2...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3scalev2n3_p2",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_scale_v2,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=2,n.neighbors=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3scalev2n3_p2",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3scalev2n3_p2",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3scalev2n3_p2"
  
  cat("           Constructing the newProp3scale n4 pro1...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3scalen4_p1",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_scale,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=4,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3scalen4_p1",argList=argList,qsFormat = T))
    }
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3scalen4_p1",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3scalen4_p1"
  
  cat("           Constructing the newProp3scale n3 pro1...\n")
  if(!file.exists(.myFilePathMakerFn("prop_newprop3scalen3_p1",argList=argList,qsFormat = T))|!saveIntermediate){
    tmp=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_scale,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=3,mc.cores = argList$ncores)
    if(saveIntermediate){
      qsave(tmp,file=.myFilePathMakerFn("prop_newprop3scalen3_p1",argList=argList,qsFormat = T))
    }
    
  } else {
    tmp=qread(.myFilePathMakerFn("prop_newprop3scalen3_p1",argList=argList,qsFormat = T))
  }
  res=c(res,list(tmp))
  names(res)[length(res)]="prop_newprop3scalen3_p1"
  
  res_arranged=list()
  
  for(icounter in 1:length(res)){
    tmp=res[[icounter]]
    
    for(j in 1:length(tmp)){
      tmp_ds=tmp[[j]]
      if(nrow(tmp_ds$prop_mat)>ncol(tmp_ds$data$logNormData)){
        tmp_prop=tmp_ds$prop_mat[((ncol(tmp_ds$data$logNormData)+1):nrow(tmp_ds$prop_mat)),1:ncol(tmp_ds$data$logNormData)]
      } else {
        tmp_prop=tmp_ds$prop_mat
      }
      
      if(ncol(tmp_prop)>ncol(tmp_ds$data$logNormData)){
        tmp_prop=tmp_prop[,colnames(tmp_prop)%in%colnames(tmp_ds$data$logNormData)]
      }
      
      tmp[[j]]$prop_mat=tmp_prop
      
      tmp_anno=pd[match(colnames(tmp_prop),row.names(pd)),c("anno_batch",colnames(pd)[grepl("^anno",colnames(pd))])]
      if(sum(is.na(tmp_anno$anno_batch))>0){
        stop("Error in propagation model ",names(res)[icounter])
      }
      
      
      #table(tmp_anno$anno_orig_cellState)
      
      if(length(unique(tmp_anno$anno_condition))<2|length(unique(tmp_anno$anno_condition))>5){
        tmp_anno$anno_condition=tmp_anno$anno_condition
      }
      
      .tmp_anno=tmp_anno
      .tmp_prop=tmp_prop
      
      if(sum(!is.na(.tmp_anno$anno_orig_cellState))>20){
        tmp_anno=.tmp_anno
        tmp_prop=.tmp_prop
        na_ind=which(is.na(tmp_anno$anno_orig_cellState))
        if(length(na_ind)>0){
          tmp_anno=tmp_anno[-na_ind,]
          tmp_prop=tmp_prop[,-na_ind]
        }
        
        na_ind=which(tmp_anno$anno_orig_cellState=="")
        if(length(na_ind)>0){
          tmp_anno=tmp_anno[-na_ind,]
          tmp_prop=tmp_prop[,-na_ind]
        }
        
        if(length(unique(tmp_anno$anno_orig_cellState))>1){
          tmp_cellState=.myOneHotFn(tmp_anno$anno_orig_cellState)
          tmp_cellState=tmp_prop %*% (as.matrix(tmp_cellState))
          
          effective_sample_size=.myEffSizePropMat(tmp_prop)
          effective_sample_size=effective_sample_size$effective_sample_size
          
          tmp_cellState2=t(apply(tmp_cellState,1,function(x) x/sum(x)))
          tmp_cellStateBkg=matrix(table(tmp_anno$anno_orig_cellState)/sum(!is.na(tmp_anno$anno_orig_cellState)),nrow=nrow(tmp_cellState2),ncol=ncol(tmp_cellState2),byrow = T)
          tmp_cellState2=tmp_cellState2 - tmp_cellStateBkg
          tmp_cellState2=tmp_cellState2/(1-tmp_cellStateBkg)
          tmp_cellState2_ind_max=rep(NA,nrow(tmp_cellState2))
          for(itr in 1:nrow(tmp_cellState2)){
            if(!is.na(max(tmp_cellState2[itr,]))){
              if(!is.nan(max(tmp_cellState2[itr,]))){
                tmp_cellState2_ind_max[itr]=which(tmp_cellState2[itr,]==max(tmp_cellState2[itr,]))
              }
            }
            
          }
          tmp_cellState2_max_type=colnames(tmp_cellState2)[tmp_cellState2_ind_max]
          tmp_cellState2_max=tmp_cellState2[cbind(1:nrow(tmp_cellState2),tmp_cellState2_ind_max)]
          tmp_res=data.frame(pseudocell=1:nrow(tmp_cellState2),score=tmp_cellState2_max,weight=rowSums(tmp_prop),dataset=tmp_ds$data$dsName,attribute="anno_orig_cellState",method=names(res)[icounter],cellType=tmp_cellState2_max_type,effective_sample_size=effective_sample_size,stringsAsFactors = F)
          res_arranged=c(res_arranged,list(tmp_res))
        }
        
      }
      
    }
    res[[icounter]]=tmp
  }
  
  res_arranged=do.call("rbind",res_arranged)
  
  library(qs)
  if(saveIntermediate){
    qsave(res_arranged,file=.myFilePathMakerFn("sim_prop_valuation",argList=argList,qsFormat = T))
  }
  
  myPropIndx=function(icounter,res,fdPattern,pd,prop_method){
    res_DE_summary=list()
    tmp=res[[icounter]]$zscore
    if(sum(is.na(tmp))>0){
      tmp[is.na(tmp)]=0
    }
    tmp_fd=fdPattern[match(colnames(tmp),row.names(fdPattern)),]
    
    tmp_pd=pd[match(colnames(res[[icounter]]$prop_mat),row.names(pd)),]
    anno_pro=res[[icounter]]$prop_mat %*% as.matrix(.myOneHotFn(tmp_pd$anno_orig_cellState))
    anno_pro_df=apply(anno_pro,1,function(x) which(x==max(x))[1])
    anno_pro_df=data.frame(pseudocell=1:nrow(anno_pro),Group=colnames(anno_pro)[anno_pro_df],stringsAsFactors = F)
    
    for(z_score_thr in c(1,1.5,2,2.5,3)){
      
      if(sum(tmp>z_score_thr)>0){
        comp_res=lapply(1:nrow(tmp),function(x){
          y=as.data.frame(table(tmp[x,]>z_score_thr,tmp_fd$dePattern))
          y2=as.data.frame(table(tmp_fd$dePattern))
          colnames(y2)[2]="total_sum"
          y=merge(y,y2,by.x="Var2",by.y="Var1")
          y$fraction=y$Freq/y$total_sum
          if(length(which(y$Var2=="none"&y$Var1==T))>0){
            y$fraction[which(y$Var2=="none"&y$Var1==T)]=y$Freq[which(y$Var2=="none"&y$Var1==T)]/sum(y$Freq[which(y$Var1==T)])
          }
          y=y[y$Var1==T,]
          y=y[,c("Var2","Freq","total_sum","fraction")]
          colnames(y)[1]="Group"
          y$Group=as.character(y$Group)
          if(nrow(y)>0){
            y$pseudocell=x
          }
          return(y)
        })
        comp_res=do.call("rbind",comp_res)
        
        comp_res2=merge(comp_res,anno_pro_df,by=c("Group","pseudocell"))
        comp_res=rbind(comp_res2,comp_res[comp_res$Group %in% c("down_others","none","upmix"),])
        aggregate(fraction~Group,data=comp_res,summary)
        comp_res$prop_method=prop_method
        comp_res$dataset=icounter
        comp_res$zscore_thr=z_score_thr
        res_DE_summary=c(res_DE_summary,list(comp_res))
      }
      
    }
    res_DE_summary=do.call("rbind",res_DE_summary)
    return(res_DE_summary)
  }
  
  res_summary=NULL
  for(igroup in unique(res_arranged$cellType)){
    tmp_main=res_arranged[which(res_arranged$cellType==igroup),]
    for(imethod in unique(tmp_main$method)){
      tmp_test=tmp_main[which(tmp_main$method==imethod),]
      tmp_bkg=tmp_main[which(tmp_main$method!=imethod),]
      
      all=effsize::cohen.d(tmp_test$score,tmp_bkg$score)$estimate
      top0.5=effsize::cohen.d(tmp_test$score[which(tmp_test$weight>0.5)],tmp_bkg$score[which(tmp_test$weight>0.5)])$estimate
      top0.8=effsize::cohen.d(tmp_test$score[which(tmp_test$weight>0.8)],tmp_bkg$score[which(tmp_test$weight>0.8)])$estimate
      
      sample_all=effsize::cohen.d(tmp_test$effective_sample_size,tmp_bkg$effective_sample_size)$estimate
      sample_top0.5=effsize::cohen.d(tmp_test$effective_sample_size[which(tmp_test$weight>0.5)],tmp_bkg$effective_sample_size[which(tmp_test$weight>0.5)])$estimate
      sample_top0.8=effsize::cohen.d(tmp_test$effective_sample_size[which(tmp_test$weight>0.8)],tmp_bkg$effective_sample_size[which(tmp_test$weight>0.8)])$estimate
      
      res_summary=rbind(res_summary,data.frame(method=imethod,group=igroup,median_score=median(tmp_test$score),median_score.5=median(tmp_test$score[which(tmp_test$weight>0.5)]),median_score.8=median(tmp_test$score[which(tmp_test$weight>0.8)]),effSize=all,effSize.5=top0.5,effSize.8=top0.8,median_sample_size=sample_all,median_sample_size.5=sample_top0.5,median_sample_size.8=sample_top0.8,stringsAsFactors=F))
    }
  }
  
  res_DE_summary=list()
  for(iprop in 1:length(res)){
    tmp_res=res[[iprop]]
    gc()
    
    
    exCentroids=NULL
    
    tmp_res=parallel::mclapply(tmp_res,.myConcensusDEFn_step2_detail_exp_final,mc.cores = 10)
    gc()
    tmpValCheck=(unlist(lapply(tmp_res,length)))
    if(sum(tmpValCheck==1)>0){
      stop(tmp_res[[which(tmpValCheck==1)[1]]])
    }
    rm(tmpValCheck)
    
    
    
    tmp_res=parallel::mclapply(1:length(tmp_res),myPropIndx,res=tmp_res,fdPattern=fdPattern,pd=pd,prop_method=names(res)[iprop],mc.cores=10)
   
    res_DE_summary=c(res_DE_summary,tmp_res)
  }
  res_DE_summary=do.call("rbind",res_DE_summary)
  
  tmp_pd=as.data.frame(table(pd$anno_batch,pd$anno_orig_cellState))
  tmp_pd$isPresent="yes"
  tmp_pd$isPresent[tmp_pd$Freq==0]="no"
  tmp_pd=tmp_pd[,c("Var1","Var2","isPresent")]
  colnames(tmp_pd)[1:2]=c("dataset","Group")
  res_DE_summary=merge(res_DE_summary,tmp_pd,by=c("dataset","Group"),all.x=T)
  res_DE_summary$isPresent[res_DE_summary$Group %in% c("none","down_others","upmix")]="yes"
  
  if(F){
    qsave(res_summary,file=.myFilePathMakerFn("sim_prop_summary_valuation",argList=argList,qsFormat = T))
  }
  
  
  if(returnDatasetObject){
    return(list(res_summary=res_summary,res=res_arranged,res_DE_summary=res_DE_summary))
  } else {
    return(list(res_summary=res_summary,res_DE_summary=res_DE_summary))
  }
}

.myClusterAnalysisFn=function(argList,clustering_resolution=1.3){
  
  options(future.globals.maxSize = 100000 * 1024^4)
  library(future)
  plan("multiprocess", workers = 10)
  plan()
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  if(sum(colnames(pd)=="anno_cluster_res")==0){
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList))
    
    harmony_net=.netEuclidean(inputPCAscores=harmony_embeddings, input.is.distance.matrix = FALSE, k.param = 20, compute.SNN = TRUE, 
                              jaccard.indx.thr = 1/15, nn.eps = 0, verbose = TRUE)
    clusters=.netFindClusters(inputGraph=harmony_net[["snn"]], algorithm = 1,resolution = clustering_resolution,group.singletons = TRUE,modularity.fxn = 1)
    
    
    clusters$anno=1
    clusters=clusters[match(row.names(pd),row.names(clusters)),]
    pd$anno_cluster_res=paste0("C",as.character(clusters[,1]))
    
    save(pd,file=.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  }
  
  
}

.myZscoreMergerFn=function(argList,inputPropRes,expData){
  
  
  
    
    set.seed(1)
    supportingFractionThr=argList$DE_supportingFractionThr
    n.adaptiveKernel=argList$DE_n.adaptiveKernel
    nPropIter=argList$DE_nPropIter
    
    #load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList))
    load(.myFilePathMakerFn("pca_centroids",argList=argList))
    
    if(!all(row.names(harmony_embeddings)==row.names(pd))){
      stop("Error in matching Names!")
    }
    
    pcaList=split(as.data.frame(harmony_embeddings),pd$anno_batch)
    
    if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename="prevalentGenes"),argList)))){
      cat("           Identifying prevalent genes ...\n")
      prevalenceEstimate=list()
      for(ik in 1:length(expData)){
        tmp=expData[[ik]]@assays$RNA@counts
        tmp=apply(tmp,1,function(x) sum(x>0)/length(x))*100
        prevalenceEstimate=c(prevalenceEstimate,list(data.frame(dsName=expData[[ik]]$anno_batch[1],gene=row.names(expData[[ik]]),prevalence=tmp,stringsAsFactors = F)))
      }
      prevalenceEstimate=do.call("rbind",prevalenceEstimate)
      prevalenceThr=max(1/(2*argList$internal_pseudocell_count),10/median(unlist(lapply(expData,ncol))))
      prevalenceEstimate=aggregate(prevalence~gene,data=prevalenceEstimate,function(x) sum(x>prevalenceThr))
      prevalenceEstimate=prevalenceEstimate[which(prevalenceEstimate$prevalence>=(supportingFractionThr*length(expData))),]
      prevalentGenes=as.character(prevalenceEstimate$gene)
      
      save(prevalentGenes,file=do.call('.myFilePathMakerFn',args=c(list(filename="prevalentGenes"),argList)))
    } else {
      load(do.call('.myFilePathMakerFn',args=c(list(filename="prevalentGenes"),argList)))
    }
    
    
    dataArranged=list()
    for(i in 1:length(expData)){
      dataArranged=c(dataArranged,list(list(logNormData=expData[[i]]@assays$RNA@data,countData=expData[[i]]@assays$RNA@counts,dsName=as.character(expData[[i]]$anno_batch[1]),pcaData=pcaList[[as.character(expData[[i]]$anno_batch[1])]])))
    }
    
    
    if(!file.exists(.myFilePathMakerFn("gene_fractionExp",argList=argList))){
      fractionExpressed=list()
      for(i in 1:length(expData)){
        tmp=apply(expData[[i]]@assays$RNA@counts,1,function(x) sum(x>0)/length(x))
        tmp=data.frame(gene=names(tmp),fraction=tmp,dsName=as.character(expData[[i]]$anno_batch[1]),stringsAsFactors = F)
        fractionExpressed=c(fractionExpressed,list(tmp))
      }
      fractionExpressedList=do.call("rbind",fractionExpressed)
      fractionExpressed=aggregate(fraction~gene,data=fractionExpressedList,mean)
      fractionExpressed$gene=as.character(fractionExpressed$gene)
      save(fractionExpressed,fractionExpressedList,file=.myFilePathMakerFn("gene_fractionExp",argList=argList))
      
    } else {
      load(.myFilePathMakerFn("gene_fractionExp",argList=argList))
    }
    
    
    require(qs)
    
    res=inputPropRes
    
    ds_df=data.frame(dsName=pd$anno_batch,batchName=pd$anno_batch,stringsAsFactors = F)
    
    cat("           Meta-analysis of z-scores ...\n")
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    
    res_count=matrix(0,nrow=nrow(res[[1]]$res),ncol(res[[1]]$res))
    res_count_w=matrix(0,nrow=nrow(res[[1]]$res),ncol(res[[1]]$res))
    res_count_sq=matrix(0,nrow=nrow(res[[1]]$res),ncol(res[[1]]$res))
    res_m=matrix(0,nrow=nrow(res[[1]]$res),ncol=ncol(res[[1]]$res))
    res_effectiveSize=matrix(0,nrow=nrow(res[[1]]$res),ncol=ncol(res[[1]]$res))
    res_count_w_ds=list()
    for(i in 1:length(res)){
      
      
      if(length(which(ds_df$batchName==res[[i]]$dsName))>0){
        tmp_dsCountName=ds_df[which(ds_df$batchName==res[[i]]$dsName),]
        #tmp_dsCount=ds_df[which(ds_df$dsName==tmp_dsCount$dsName[1]),]
        
        tmp_w_zscore=res[[i]]$w_zscore_centroid
        tmp_count_mat=rep(1,ncol(res[[i]]$res))
        
        tmp_fraction=fractionExpressedList[which(fractionExpressedList$dsName ==res[[i]]$dsName),]
        tmp_fraction=tmp_fraction[match(row.names(res[[i]]$res),tmp_fraction$gene),]
        if(sum(is.na(tmp_fraction$fraction))>0){
          tmp_fraction$fraction[is.na(tmp_fraction$fraction)]=0
        }
        tmp_fraction$fraction[which(tmp_fraction$fraction<0.001)]=0
        tmp_fraction$fraction[which(tmp_fraction$fraction>=0.001)]=1
        tmp_fraction=matrix(tmp_fraction$fraction,nrow=nrow(res[[i]]$res),ncol=ncol(res[[i]]$res),byrow = F)
        
        tmp_res_m=res[[i]]$res
        tmp_res_m[which(tmp_fraction==0)]=0
        
        #removing pseudocells that are not assigned to the dataset
        tmpInd=apply(res[[i]]$res,2,function(x) sum(x==0)/length(x))
        
        if(sum(tmpInd>0.95)>0){
          tmp_count_mat[which(tmpInd>0.95)]=0
          tmp_w_zscore[which(tmpInd>0.95)]=0
          tmp_res_m[which(tmpInd>0.95)]=0
        }
        
        #tmp_fraction is binary indicating if the gene is expressed in the dataset
        tmp_w_zscore=sweep(tmp_fraction,2,tmp_w_zscore,"*")
        
        #tmp_w_zscore=tmp_w_zscore/tmp_dsCount
        
        tmp_w_zscore_sq=tmp_w_zscore^2
        tmp_count_mat=sweep(tmp_fraction,2,tmp_count_mat,"*")
        tmp_count_mat[which(tmp_count_mat<0.001)]=0
        
        
        if(sum(is.na(tmp_w_zscore))>0){
          stop(paste("Check",i,"for NA values"))
        }
        
        res_m=res_m+tmp_res_m
        
        res_count=res_count+tmp_count_mat
        res_count_w=res_count_w+tmp_w_zscore
        res_count_sq=res_count_sq+tmp_w_zscore_sq
        if(sum(names(res_count_w_ds)==tmp_dsCountName$dsName[1])>0){
          dsInd=which(names(res_count_w_ds)==tmp_dsCountName$dsName[1])
          res_count_w_ds[[dsInd]]=res_count_w_ds[[dsInd]]+tmp_w_zscore
        } else {
          res_count_w_ds=c(res_count_w_ds,new=list(tmp_w_zscore))
          names(res_count_w_ds)[which(names(res_count_w_ds)=="new")]=tmp_dsCountName$dsName[1]
        }
        
      }
      
    }
    gc()
    #res_count=t(matrix(res_count,nrow=ncol(res[[i]]$res),ncol=nrow(res[[i]]$res)))
    #res_count_w=t(matrix(res_count_w,nrow=ncol(res[[i]]$res),ncol=nrow(res[[i]]$res)))
    #res_count_sq=t(matrix(res_count_sq,nrow=ncol(res[[i]]$res),ncol=nrow(res[[i]]$res)))
    
    ds_count_sq=matrix(0,nrow=nrow(res_count_w_ds[[1]]),ncol=ncol(res_count_w_ds[[1]]))
    for(i in 1:length(res_count_w_ds)){
      if(sum(ds_df$dsName==names(res_count_w_ds))>0){
        ds_count_sq=ds_count_sq+(res_count_w_ds[[i]])^2
      }
    }
    
    res_count_sq=matrix(0,nrow=nrow(res_count_w_ds[[1]]),ncol=ncol(res_count_w_ds[[1]]))
    for(i in 1:length(res_count_w_ds)){
      res_count_sq=res_count_sq+res_count_w_ds[[i]]^2
    }
    
    res_effectiveSize=res_count_w^2/res_count_sq
    res_effectiveSize[which(res_count_w==0)]=0
    
    if(length(which(res_effectiveSize<0.001))>0){
      res_effectiveSize2=res_effectiveSize
      res_effectiveSize2[which(res_effectiveSize<0.001)]=0
      res_effectiveSize=res_effectiveSize2
      rm(res_effectiveSize2)
      
    }
    
    #tmp_res_count_sq=sqrt(res_count_sq)
    tmp_res_count_sq=sqrt(ds_count_sq)
    tmp_res_count_sq[which(tmp_res_count_sq<1)]=1
    res_m=res_m/tmp_res_count_sq
    res_m[which(res_count==0)]=0
    
    
    ind=which(res_count!=0,arr.ind = T)
    
    res_arranged=data.frame(gene_index=ind[,1],centroid_index=ind[,2],gene=row.names(res[[1]]$res)[ind[,1]],centroid=colnames(res[[1]]$res)[ind[,2]],zscore=res_m[ind],count=res_count[ind],count_w=res_count_w[ind],effective_size=res_effectiveSize[ind],stringsAsFactors = F)
    fractionExpressed=fractionExpressed[match(res_arranged$gene,fractionExpressed$gene),]
    res_arranged$overall_fractionExpressed=fractionExpressed$fraction
    #res_arranged=res_arranged[which(res_arranged$effective_size>=max(1,supportingFractionThr*length(expData))),]
    #,score_seq=res_indBase[ind]
    
    res_fd=NULL
    for(i in 1:length(expData)){
      fd=as.data.frame(expData[[i]]@assays$RNA@meta.features)
      fd$ensembl_gene_id=gsub("_","-",fd$ensembl_gene_id)
      slCols=c("gene_name","gene_biotype","symbol","gene_short_name","ensembl_gene_id")
      fd=fd[,colnames(fd) %in% slCols]
      if(sum(is.na(fd$ensembl_gene_id))>0){
        fd$ensembl_gene_id[is.na(fd$ensembl_gene_id)]=row.names(fd)[is.na(fd$ensembl_gene_id)]
      }
      
      if(!is.null(res_fd)){
        fd=fd[!fd$ensembl_gene_id %in% res_fd$ensembl_gene_id,]
        if(nrow(fd)>0){
          res_fd=rbind(res_fd,fd)
        }
      } else {
        res_fd=fd
      }
    }
    fd=res_fd
    rm(res_fd)
    
    fd=fd[match(res_arranged$gene,fd$ensembl_gene_id),]
    res_arranged=cbind(res_arranged,fd[,-which(colnames(fd)=="ensembl_gene_id")])
    
    network=.myConcensusDEFn_step2_FindNeighbors(inputCentroids = pca_centroid,argList=argList,verbose = F)
    
    snnNet=as.matrix(network$snn)
    ind=which(snnNet>0,arr.ind = T)
    
    snnNet=data.frame(Fnode=row.names(snnNet)[ind[,1]],Snode=row.names(snnNet)[ind[,2]],weight=snnNet[ind],stringsAsFactors = F)
    snnNet=snnNet[snnNet$Fnode!=snnNet$Snode,]
    snnNet$weight=1-snnNet$weight
    snnNet2=igraph::graph_from_data_frame(snnNet, directed = F, vertices = NULL)
    
    snnNet2=igraph::distances(snnNet2, mode = c("all"), weights = snnNet$weight, algorithm = "dijkstra")
    if(length(which(snnNet2==Inf))>0){
      snnNet2[which(snnNet2==Inf)]=(max(snnNet2[snnNet2!=Inf],na.rm = T)+1)
    }
    
    geneList=unique(as.character(res_arranged$gene))
    
    
    cat("           Calculating cor and MI ...\n")
    resCorMI=NULL
    slGenes=aggregate(zscore~gene,data=res_arranged,function(x) max(x))
    slGenes=slGenes[which(slGenes$zscore>2),]
    
    myCorMIfn=function(inputGenes,res_arranged,snnNet2){
      resCorMI=NULL
      for(i in inputGenes){
        tmp=res_arranged[which(res_arranged$gene==i),]
        if(length(which(tmp$effective_size<1))>0){
          tmp$zscore[which(tmp$effective_size<1)]=0
        }
        
        tmpNet=snnNet2[row.names(snnNet2) %in% tmp$centroid,colnames(snnNet2) %in% tmp$centroid]
        tmp=tmp[match(colnames(tmpNet),tmp$centroid),]
        tmp$zscore[is.na(tmp$zscore)]=0
        if(max(abs(tmp$zscore))>2){
          zOrdered=tmp$zscore[order(tmp$zscore,decreasing = T)]
          
          tmpInd=which(row.names(tmpNet) %in% tmp$centroid[which(tmp$zscore>=max(zOrdered[3],1))])
          if(length(tmpInd)>0){
            if(length(tmpInd)>1){
              distMat=tmpNet[tmpInd,]
              distMat=apply(distMat,2,mean)
            } else {
              distMat=tmpNet[tmpInd,]
            }
            
            
            tmp=data.frame(distance=distMat,heat=tmp$zscore)
            tmpCor=cor(tmp$distance,tmp$heat)
            tmpMI=infotheo::discretize( tmp, disc="equalfreq", nbins=10 )
            tmpMI=infotheo::mutinformation(tmpMI$distance,tmpMI$heat)
            resCorMI=rbind(resCorMI,data.frame(gene=i,cor=tmpCor,MI=tmpMI,stringsAsFactors = F))
            
          } else{
            resCorMI=rbind(resCorMI,data.frame(gene=i,cor=0,MI=0,stringsAsFactors = F))
          }
        }
      }
      return(resCorMI)
    }
    
    if(argList$ncores>1){
      geneList=split(intersect(as.character(slGenes$gene),prevalentGenes), cut(1:length(intersect(as.character(slGenes$gene),prevalentGenes)), argList$ncores, labels = FALSE)) 
    } else {
      geneList=list(intersect(as.character(slGenes$gene),prevalentGenes))
    }
    
    resCorMI=parallel::mclapply(geneList,myCorMIfn,res_arranged=res_arranged,snnNet2=snnNet2,mc.cores = argList$ncores)
    resCorMI=do.call("rbind",resCorMI)
    
    if(length(setdiff(unique(prevalentGenes),resCorMI$gene))>0){
      resCorMI=rbind(resCorMI,data.frame(gene=setdiff(unique(prevalentGenes),resCorMI$gene),cor=0,MI=0,stringsAsFactors = F))
    }
    
    #res_arranged=res_arranged[which(as.character(res_arranged$gene) %in%prevalentGenes),]
    res_arranged=merge(res_arranged,resCorMI,by="gene",all.x=T)
    
    
    res_arranged2=res_arranged
    res_arranged2=res_arranged2[order(res_arranged2$zscore,decreasing = T),]
    res_w_path=res_arranged
    #zscore_count measures in how many centroids each gene as a nominal significance
    zscore_counts=aggregate(zscore~gene,data=res_arranged,function(x) sum(x>2))
    colnames(zscore_counts)=c("gene","zscore_count")
    res_arranged=merge(res_arranged,zscore_counts,by="gene",all.x=T)
    zscore_counts=aggregate(zscore~gene,data=res_arranged,function(x) sum(abs(x)>2))
    colnames(zscore_counts)=c("gene","zscore_count_abs")
    res_arranged=merge(res_arranged,zscore_counts,by="gene",all.x=T)
    if(sum(res_arranged$zscore_count>0)>0){
      res_arranged=res_arranged[res_arranged$zscore_count>0,]
    }
    
    res_arranged=res_arranged[res_arranged$effective_size/max(res_arranged$effective_size,na.rm = T)>=supportingFractionThr,]
    res_arranged=res_arranged[order(res_arranged$zscore,decreasing = T),]
    
    res_arranged$score_seq="0"
    for(i in 1:length(res)){
      res_arranged$score_seq=paste0(res_arranged$score_seq,",",round(res[[i]]$res[cbind(res_arranged$gene_index,res_arranged$centroid_index)],2))
    }
    res_arranged$score_seq=gsub(",0,",",",res_arranged$score_seq)
    res_arranged$score_seq=gsub("^0,","",res_arranged$score_seq)
    res_arranged$score_seq=gsub(",0$","",res_arranged$score_seq)
    
    res=res_arranged
  
  
  return(list(res=res_arranged,res_w_path=res_w_path))
}

.myDE2stepFn_BN_temp=function(argList){
  
  if(!file.exists(.myFilePathMakerFn("DE_analysis_res",argList=argList))){
    set.seed(1)
    supportingFractionThr=argList$DE_supportingFractionThr
    n.adaptiveKernel=argList$DE_n.adaptiveKernel
    nPropIter=argList$DE_nPropIter
    
    #load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("pca_centroids",argList=argList))
    
    if(!all(row.names(harmony_embeddings)==row.names(pd))){
      stop("Error in matching Names!")
    }
    
    pcaList=split(as.data.frame(harmony_embeddings),pd$anno_batch)
    
    if(T){
      load(.myFilePathMakerFn("exp_merged",argList=argList,expData=T))
      fd=tmp@assays$RNA@meta.features
      fd$dePattern="none"
      dwn=apply(fd[,grepl("DEFacGroup",colnames(fd))],1,function(x) sum(x<1))
      fd$dePattern[which(dwn>0)]="down_others"
      fd$dePattern[which(fd$DEFacGroup1>1)]="up1"
      fd$dePattern[which(fd$DEFacGroup2>1)]="up2"
      fd$dePattern[which(fd$DEFacGroup3>1)]="up3"
      fd$dePattern[which(fd$DEFacGroup4>1)]="up4"
      fd$dePattern[which(fd$DEFacGroup5>1)]="up5"
      upmix=apply(fd[,grepl("DEFacGroup",colnames(fd))],1,function(x) sum(x>1))
      fd$dePattern[which(upmix>1)]="upmix"
      fdPattern=fd
      
      expData=SplitObject(tmp, split.by = "anno_batch")
    }
    
    
    if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename="prevalentGenes"),argList)))){
      cat("           Identifying prevalent genes ...\n")
      prevalenceEstimate=list()
      for(ik in 1:length(expData)){
        tmp=expData[[ik]]@assays$RNA@counts
        tmp=apply(tmp,1,function(x) sum(x>0)/length(x))*100
        prevalenceEstimate=c(prevalenceEstimate,list(data.frame(dsName=expData[[ik]]$anno_batch[1],gene=row.names(expData[[ik]]),prevalence=tmp,stringsAsFactors = F)))
      }
      prevalenceEstimate=do.call("rbind",prevalenceEstimate)
      prevalenceThr=max(1/(2*argList$internal_pseudocell_count),10/median(unlist(lapply(expData,ncol))))
      prevalenceEstimate=aggregate(prevalence~gene,data=prevalenceEstimate,function(x) sum(x>prevalenceThr))
      prevalenceEstimate=prevalenceEstimate[which(prevalenceEstimate$prevalence>=(supportingFractionThr*length(expData))),]
      prevalentGenes=as.character(prevalenceEstimate$gene)
      
      save(prevalentGenes,file=do.call('.myFilePathMakerFn',args=c(list(filename="prevalentGenes"),argList)))
    } else {
      load(do.call('.myFilePathMakerFn',args=c(list(filename="prevalentGenes"),argList)))
    }
    
    
    dataArranged=list()
    for(i in 1:length(expData)){
      dataArranged=c(dataArranged,list(list(logNormData=expData[[i]]@assays$RNA@data,countData=expData[[i]]@assays$RNA@counts,dsName=as.character(expData[[i]]$anno_batch[1]),pcaData=pcaList[[as.character(expData[[i]]$anno_batch[1])]])))
    }
    
    
    if(!file.exists(.myFilePathMakerFn("gene_fractionExp",argList=argList))){
      fractionExpressed=list()
      for(i in 1:length(expData)){
        tmp=apply(expData[[i]]@assays$RNA@counts,1,function(x) sum(x>0)/length(x))
        tmp=data.frame(gene=names(tmp),fraction=tmp,dsName=as.character(expData[[i]]$anno_batch[1]),stringsAsFactors = F)
        fractionExpressed=c(fractionExpressed,list(tmp))
      }
      fractionExpressedList=do.call("rbind",fractionExpressed)
      fractionExpressed=aggregate(fraction~gene,data=fractionExpressedList,mean)
      fractionExpressed$gene=as.character(fractionExpressed$gene)
      save(fractionExpressed,fractionExpressedList,file=.myFilePathMakerFn("gene_fractionExp",argList=argList))
      
    } else {
      load(.myFilePathMakerFn("gene_fractionExp",argList=argList))
    }
    
    res_summary=NULL
    DE_res=NULL
    prop_mat_fileList=c("prop_org_v1_p1","prop_org_v1_p2","prop_org_v1_p3","prop_org_v2_p3","prop_DNN_noAdjPurity_p3","prop_DNN_AdjPurity_p3","prop_runIndx1_p3","prop_runIndx1_p2","prop_runIndx1_p4","prop_runIndx1_p6","prop_runIndx2_p3","prop_runIndx3_p3","prop_runIndx4_p3","prop_runIndx5_p3","prop_runIndx6_p3","prop_runIndx7_p3","prop_runIndx8_p3","prop_runIndx9_p3","prop_newpropn1_p3","prop_newpropn2_p3","prop_newpropn2_p1","prop_newpropn2_p2","prop_newpropn2_p4","prop_newpropn3_p3","prop_newprop2n3_p3","prop_newprop2n3_p6","prop_newprop2n2_p3","prop_newprop3n2_p3_3","prop_newprop3n2_p3","prop_newprop3n4_p3","prop_newprop3n3_p3","prop_newprop3n6_p3","prop_newprop3n8_p3","prop_newprop4n4_p3","prop_newprop4n3_p3","prop_newprop4n6_p3","prop_newprop3scalev2n4_p1","prop_newprop3scalev2n3_p1","prop_newprop3scalev2n3_p2","prop_newprop3scalen4_p1","prop_newprop3scalen3_p1")
    for(iprop in prop_mat_fileList){
      res=qread(.myFilePathMakerFn(iprop,argList=argList,qsFormat = T))
      
      tmp=res
      
      for(j in 1:length(tmp)){
        tmp_ds=tmp[[j]]
        if(nrow(tmp_ds$prop_mat)>200){
          tmp_prop=tmp_ds$prop_mat[((ncol(tmp_ds$data$logNormData)+1):nrow(tmp_ds$prop_mat)),1:ncol(tmp_ds$data$logNormData)]
        } else {
          tmp_prop=tmp_ds$prop_mat
        }
        
        if(ncol(tmp_prop)>ncol(tmp_ds$data$logNormData)){
          tmp_prop=tmp_prop[,colnames(tmp_prop)%in%colnames(tmp_ds$data$logNormData)]
        }
        
        res[[j]]$prop_mat=tmp_prop
        
      }
      
      exCentroids=NULL
      
      #inputData=res[[1]];centroidPCAdata=pca_centroid;exCentroids=exCentroids;argList=argList;nDataset=length(dataArranged);within_gene_norm=F
      res2=parallel::mclapply(res,.myConcensusDEFn_step2_detail_exp,centroidPCAdata=pca_centroid,exCentroids=exCentroids,argList=argList,nDataset=length(dataArranged),within_gene_norm=F,mc.cores = argList$ncores)
      #qsave(res2,.myFilePathMakerFn("prop_res2",argList=argList,qsFormat =T))
      #res2=qread("~/data/results/InN_nPC30/res2.qs")
      #res2=qread(.myFilePathMakerFn("prop_res2",argList=argList,qsFormat =T))
      res=res2
      tmpValCheck=(unlist(lapply(res,length)))
      if(sum(tmpValCheck==1)>0){
        stop(res[[which(tmpValCheck==1)[1]]])
      }
      rm(tmpValCheck)
      
      res_m=.myZscoreMergerFn(argList=argList,inputPropRes=res,expData=expData)
      
      tmp_res=res_m$res[order(res_m$res$zscore,decreasing = T),]
      tmp_res=tmp_res[!duplicated(tmp_res$gene),]
      tmp_res=merge(tmp_res,fdPattern,by.x="gene",by.y="Gene")
      DE_res=rbind(DE_res,tmp_res)
      for(z_score_thr in c(2,2.5,3,4,5)){
        if(sum(res_m$res$zscore>z_score_thr)>0){
          sl_genes=res_m$res
          sl_genes=sl_genes[which(sl_genes$zscore>z_score_thr),]
          sl_genes=unique(sl_genes$gene)
          
          tmp_fd=fdPattern
          tmp_fd$sl="no"
          tmp_fd$sl[tmp_fd$Gene %in% sl_genes]="Yes"
          tmp_fd=as.data.frame(table(tmp_fd$dePattern,tmp_fd$sl))
          tmp_fd$zThr=z_score_thr
          tmp_fd$fraction=0
          for(i in unique(tmp_fd$Var1)){
            tmp_fd$fraction[which(tmp_fd$Var1==i)]=tmp_fd$Freq[which(tmp_fd$Var1==i)]/sum(tmp_fd$Freq[which(tmp_fd$Var1==i)])
          }
          if(length(which(tmp_fd$Var1=="none"&tmp_fd$Var2=="Yes"))>0){
            tmp_fd$fraction[which(tmp_fd$Var1=="none"&tmp_fd$Var2=="Yes")]=tmp_fd$Freq[which(tmp_fd$Var1=="none"&tmp_fd$Var2=="Yes")]/sum(tmp_fd$Freq[which(tmp_fd$Var2=="Yes")])
          }
          tmp_fd=tmp_fd[tmp_fd$Var2=="Yes",]
          tmp_fd$prop_method=iprop
          res_summary=rbind(res_summary,tmp_fd)
        }
        
      }
      gc()
      
    }
    
    save(res_summary,file=.myFilePathMakerFn("DE_analysis_summary",argList=argList))
    save(DE_res,file=.myFilePathMakerFn("DE_analysis_res",argList=argList))
    
    
  }
  
  
  
  
  return("Done")
}

.myDE2stepFn_BN=function(argList){
  
  if(!file.exists(.myFilePathMakerFn("DE_analysis_summary",argList=argList))){
    set.seed(1)
    supportingFractionThr=argList$DE_supportingFractionThr
    n.adaptiveKernel=argList$DE_n.adaptiveKernel
    nPropIter=argList$DE_nPropIter
    
    #load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("pca_centroids",argList=argList))
    
    if(!all(row.names(harmony_embeddings)==row.names(pd))){
      stop("Error in matching Names!")
    }
    
    pcaList=split(as.data.frame(harmony_embeddings),pd$anno_batch)
    
    if(T){
      load(.myFilePathMakerFn("exp_merged",argList=argList,expData=T))
      fd=tmp@assays$RNA@meta.features
      fd$dePattern="none"
      dwn=apply(fd[,grepl("DEFacGroup",colnames(fd))],1,function(x) sum(x<1))
      fd$dePattern[which(dwn>0)]="down_others"
      fd$dePattern[which(fd$DEFacGroup1>1)]="Group1"
      fd$dePattern[which(fd$DEFacGroup2>1)]="Group2"
      fd$dePattern[which(fd$DEFacGroup3>1)]="Group3"
      fd$dePattern[which(fd$DEFacGroup4>1)]="Group4"
      fd$dePattern[which(fd$DEFacGroup5>1)]="Group5"
      upmix=apply(fd[,grepl("DEFacGroup",colnames(fd))],1,function(x) sum(x>1))
      fd$dePattern[which(upmix>1)]="upmix"
      fdPattern=fd
      
      expData=SplitObject(tmp, split.by = "anno_batch")
    }
    
    dataArranged=list()
    for(i in 1:length(expData)){
      dataArranged=c(dataArranged,list(list(logNormData=expData[[i]]@assays$RNA@data,countData=expData[[i]]@assays$RNA@counts,dsName=as.character(expData[[i]]$anno_batch[1]),pcaData=pcaList[[as.character(expData[[i]]$anno_batch[1])]])))
    }
    
    res_summary=list()
    
    prop_mat_fileList=c("prop_newprop3n3_p3_final","prop_newprop3n4_p3_final","prop_org_v1_p1","prop_org_v1_p2","prop_org_v1_p3","prop_org_v2_p3","prop_DNN_noAdjPurity_p3","prop_DNN_AdjPurity_p3","prop_runIndx1_p3","prop_runIndx1_p2","prop_runIndx1_p4","prop_runIndx1_p6","prop_runIndx2_p3","prop_runIndx3_p3","prop_runIndx4_p3","prop_runIndx5_p3","prop_runIndx6_p3","prop_runIndx7_p3","prop_runIndx8_p3","prop_runIndx9_p3","prop_newpropn1_p3","prop_newpropn2_p3","prop_newpropn2_p1","prop_newpropn2_p2","prop_newpropn2_p4","prop_newpropn3_p3","prop_newprop2n3_p3","prop_newprop2n3_p6","prop_newprop2n2_p3","prop_newprop3n4_p3","prop_newprop3n3_p3","prop_newprop3n6_p3","prop_newprop3n8_p3","prop_newprop4n4_p3","prop_newprop4n3_p3","prop_newprop4n6_p3","prop_newprop3scalev2n4_p1","prop_newprop3scalev2n3_p1","prop_newprop3scalev2n3_p2","prop_newprop3scalen4_p1","prop_newprop3scalen3_p1")
    for(iprop in prop_mat_fileList){
      res=qread(.myFilePathMakerFn(iprop,argList=argList,qsFormat = T))
      gc()
      
      
      for(j in 1:length(res)){
        tmp_ds=res[[j]]
        if(nrow(tmp_ds$prop_mat)>200){
          tmp_prop=tmp_ds$prop_mat[((ncol(tmp_ds$data$logNormData)+1):nrow(tmp_ds$prop_mat)),1:ncol(tmp_ds$data$logNormData)]
        } else {
          tmp_prop=tmp_ds$prop_mat
        }
        
        if(ncol(tmp_prop)>ncol(tmp_ds$data$logNormData)){
          tmp_prop=tmp_prop[,colnames(tmp_prop)%in%colnames(tmp_ds$data$logNormData)]
        }
        
        res[[j]]$prop_mat=tmp_prop
        
      }
      
      exCentroids=NULL
      
      res=parallel::mclapply(res,.myConcensusDEFn_step2_detail_exp_final,mc.cores = 10)
      gc()
      tmpValCheck=(unlist(lapply(res,length)))
      if(sum(tmpValCheck==1)>0){
        stop(res[[which(tmpValCheck==1)[1]]])
      }
      rm(tmpValCheck)
      
      myPropIndx=function(icounter,res,fdPattern,pd){
        res_summary=list()
        tmp=res[[icounter]]$zscore
        if(sum(is.na(tmp))>0){
          tmp[is.na(tmp)]=0
        }
        tmp_fd=fdPattern[match(colnames(tmp),row.names(fdPattern)),]
        
        tmp_pd=pd[match(colnames(res[[icounter]]$prop_mat),row.names(pd)),]
        anno_pro=res[[icounter]]$prop_mat %*% as.matrix(.myOneHotFn(tmp_pd$anno_orig_cellState))
        anno_pro_df=apply(anno_pro,1,function(x) which(x==max(x))[1])
        anno_pro_df=data.frame(pseudocell=1:nrow(anno_pro),Group=colnames(anno_pro)[anno_pro_df],stringsAsFactors = F)
        
        for(z_score_thr in c(1,1.5,2,2.5,3)){
          
          if(sum(tmp>z_score_thr)>0){
            comp_res=lapply(1:nrow(tmp),function(x){
              y=as.data.frame(table(tmp[x,]>z_score_thr,tmp_fd$dePattern))
              y2=as.data.frame(table(tmp_fd$dePattern))
              colnames(y2)[2]="total_sum"
              y=merge(y,y2,by.x="Var2",by.y="Var1")
              y$fraction=y$Freq/y$total_sum
              if(length(which(y$Var2=="none"&y$Var1==T))>0){
                y$fraction[which(y$Var2=="none"&y$Var1==T)]=y$Freq[which(y$Var2=="none"&y$Var1==T)]/sum(y$Freq[which(y$Var1==T)])
              }
              y=y[y$Var1==T,]
              y=y[,c("Var2","Freq","total_sum","fraction")]
              colnames(y)[1]="Group"
              y$Group=as.character(y$Group)
              if(nrow(y)>0){
                y$pseudocell=x
              }
              return(y)
            })
            comp_res=do.call("rbind",comp_res)
            
            comp_res2=merge(comp_res,anno_pro_df,by=c("Group","pseudocell"))
            comp_res=rbind(comp_res2,comp_res[comp_res$Group %in% c("down_others","none","upmix"),])
            aggregate(fraction~Group,data=comp_res,summary)
            comp_res$prop_method=iprop
            comp_res$dataset=icounter
            comp_res$zscore_thr=z_score_thr
            res_summary=c(res_summary,list(comp_res))
          }
          
        }
        res_summary=do.call("rbind",res_summary)
        return(res_summary)
      }
      
      tmp_res=parallel::mclapply(1:length(res),myPropIndx,res=res,fdPattern=fdPattern,pd=pd,mc.cores=10)
      res_summary=c(res_summary,tmp_res)
    }
    res_summary=do.call("rbind",res_summary)
    save(res_summary,file=.myFilePathMakerFn("DE_analysis_summary",argList=argList))
    
    
  }
  
  
  
  
  return("Done")
}

