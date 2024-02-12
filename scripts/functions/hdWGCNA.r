############################
# datExpr
###########################


#' SetDatExpr
#'
#' This function specifies the gene expression matrix for co-expression network analysis.
#'
#' @param seurat_obj A Seurat object
#' @param group_name A string containing a group present in the provided group.by column or in the Seurat Idents. A character vector can be provided to select multiple groups at a time.
#' @param use_metacells A logical determining if we use the metacells (TRUE) or the full expression matrix (FALSE)
#' @param group.by A string containing the name of a column in the Seurat object with cell groups (clusters, cell types, etc). If NULL (default), hdWGCNA uses the Seurat Idents as the group.
#' @param multi.group.by A string containing the name of a column in the Seurat object with groups for consensus WGCNA (dataset, sample, condition, etc)
#' @param multi_group_name A string containing the name of a group present in the multi.group.by column.
#' @param wgcna_name A string containing the name of the WGCNA slot in seurat_obj@misc. Default = NULL which retrieves the currently active WGCNA data
#' @keywords scRNA-seq
#' @export
SetDatExpr_Alt <- function(
    seurat_obj,
    group_name,
    use_metacells=TRUE,
    group.by=NULL,
    multi.group.by = NULL,
    multi_group_name = NULL,
    return_seurat = TRUE,
    wgcna_name=NULL,
    assay=NULL,
    slot = 'data',
    ...
){
  
  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  
  # get parameters from seurat object
  params <- GetWGCNAParams(seurat_obj, wgcna_name)
  genes_use <- GetWGCNAGenes(seurat_obj, wgcna_name)
  
  # get metacell object
  m_obj <- GetMetacellObject(seurat_obj, wgcna_name)
  
  # use metacells or whole seurat object?
  if(use_metacells & !is.null(m_obj)){
    s_obj <- m_obj
  } else{
    if(is.null(m_obj)){warning("Metacell Seurat object not found. Using full Seurat object instead.")}
    s_obj <- seurat_obj
  }
  
  # get the metadata from the seurat object:
  seurat_meta <- s_obj@meta.data
  
  if(is.null(assay)){
    assay <- DefaultAssay(s_obj)
    warning(paste0('assay not specified, trying to use assay ', assay))
  }
  
  # check the assay:
  if(!(assay %in% names(s_obj@assays))){
    stop("Assay not found. Check names(seurat_obj@assays) or names(GetMetacellObject(seurat_obj)@assays)")
  }
  
  if(!is.null(group.by)){
    
    # check that group.by is in the Seurat object & in the metacell object:
    if(!(group.by %in% colnames(s_obj@meta.data))){
      m_cell_message <- ""
      if(use_metacells){m_cell_message <- "metacell"}
      stop(paste0(group.by, ' not found in the meta data of the ', m_cell_message, ' Seurat object'))
    }
    
    # check that the selected groups are in the Seurat object:
    if(!all(group_name %in% s_obj@meta.data[[group.by]])){
      groups_not_found <- group_name[!(group_name %in% s_obj@meta.data[[group.by]])]
      stop(
        paste0("Some groups in group_name are not found in the seurat_obj: ", paste(groups_not_found, collapse=', '))
      )
    }
    
  }
  
  # columns to group by for cluster/celltype
  if(!is.null(group.by)){
    seurat_meta <- seurat_meta %>% subset(get(group.by) %in% group_name)
  }
  
  # check that the group names are actually in the group.by column:
  
  # subset further if multiExpr:
  if(!is.null(multi.group.by)){
    seurat_meta <- seurat_meta %>% subset(get(multi.group.by) %in% multi_group_name)
  }
  
  # get list of cells to use
  cells <- rownames(seurat_meta)
  
  # get expression data from seurat obj
  datExpr <- as.data.frame(
    Seurat::GetAssayData(
      s_obj,
      assay=assay,
      slot=slot
    )
  )[genes_use,cells]
  
  # transpose data
  datExpr <- as.data.frame(t(datExpr))
  
  # only get good genes:
  # if(is.null(multi.group.by)){
  #   gene_list = genes_use[WGCNA::goodGenes(datExpr, ...)]
  #   datExpr <- datExpr[,gene_list]
  # }
  
  if(return_seurat){
    
    gene_list <- genes_use[WGCNA::goodGenes(datExpr, ...)]
    datExpr <- datExpr[,gene_list]
    
    # update the WGCNA gene list:
    seurat_obj <- SetWGCNAGenes(seurat_obj, gene_list, wgcna_name)
    
    # set the datExpr in the Seurat object
    seurat_obj@misc[[wgcna_name]]$datExpr <- datExpr
    out <- seurat_obj
  } else{
    out <- datExpr
  }
  out
}