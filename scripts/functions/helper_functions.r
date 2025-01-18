##############################################################################################################################
filter_lowly_exp_pathways = function(expressed, all_paths, pathway_gene_threshold=0.33){
    gsets_per_celltype = list()

    for(i in names(expressed)){
        index = unlist(lapply(names(all_paths), function(x) (sum(all_paths[[x]]%in%expressed[[i]])/length(all_paths[[x]]))>pathway_gene_threshold))
        p = all_paths[(index)]
        x = lapply(names(p), function(y) intersect(expressed[[i]], p[[y]]))
        names(x) = names(p)
        x = x[!duplicated(x)]
        gsets_per_celltype[[i]] = x
    }
    return(gsets_per_celltype)
}

##############################################################################################################################
filter_genesets <- function(expressed, all_paths, pathway_gene_threshold = 0.33, min_genes = 5, max_genes = inf) {
  # filter by absolute gene count thresholds
  gene_count_filter <- sapply(seq_along(all_paths), function(x) {
    path_length <- length(all_paths[[x]])
    path_length >= min_genes && path_length <= max_genes
  })
  
  size_filtered_paths <- all_paths[gene_count_filter]
  
  # calculate percentage of pathway genes that are expressed
  index <- sapply(seq_along(size_filtered_paths), function(x) {
    (sum(size_filtered_paths[[x]] %in% expressed) / length(size_filtered_paths[[x]])) > pathway_gene_threshold
  })
  
  # flter pathways meeting threshold
  filtered_paths <- size_filtered_paths[index]
  
  # get intersections and remove duplicates
  x <- lapply(seq_along(filtered_paths), function(y) {
    intersect(expressed, filtered_paths[[y]])
  })
  names(x) <- names(filtered_paths)
  x <- x[!duplicated(x)]
  
  return(x)
}

##############################################################################################################################
read.geneset = function(path_to_gset){
    invisible(capture.output({
        bp = GSA.read.gmt(path_to_gset)
    }))
    out = bp$genesets
    out = lapply(1:length(out), function(x) out[[x]][out[[x]]!=''])
    names(out) = bp$geneset.names
    return(out)
}

####################################################################################################################
get_gset_names_by_category = function(cat, gsets){
  gset = unlist(lapply(gsets, function(x) unlist(sum(sapply(cat, grepl, x))>0)))
  gset = (gsets[gset])
  return(gset)
}
####################################################################################################################

##############################################################################################################################
get_pathway_fits = function(order, av_expression, pathways, top_20, summary){
    # run GSVA on GO pathways
    all_data = list()
    print('running GSVA...')
    out_bp_terms = lapply(order, function(x) t(gsva(as.matrix(av_expression[[x]]), pathways[[x]], mx.diff=TRUE, verbose=F, kcdf=c("Gaussian"), min.sz=5, max.sz = 150, parallel.sz=16)))
    names(out_bp_terms) = order
    all_data[['gsva_out']] = out_bp_terms

    # get linear model fits
    print('getting linear model fits...')
    fits = get_fits(out_bp_terms, summary)
    all_data[['fits_all']] = fits

    # get matrix of scores for heatmap
    print('get matrix of scores')
    scores = get_scores(fits)
    all_data[['scores_all']] = scores

    print('filter by score 1.3')
    names = unique(unname(unlist(lapply(names(scores$all), function(x) rownames(scores$all[[x]])))))
    mat = get_matrix(scores$all, names)
    mat = mat[unname(rowSums(abs(mat)>1.3)>0),]

    if(top_20==TRUE){
        index = unique(unname(unlist(lapply(colnames(mat), function(x) order(abs(mat[[x]]),decreasing = T)[1:20]))))
        mat = mat[index,]
    }
    all_data[['scores_filtered']] = mat
    return(all_data)
}
####################################################################################################################
get_fits = function(gsva_out, meta){
    fits = list()
    for(i in names(gsva_out)){
        predict = meta[as.character(rownames(gsva_out[[i]])),]
        mod = model.matrix(~APOE4 + amyloid + nft + age_death + msex + pmi, data=predict)
        fits[[i]] = fit.gsva(mod, i, gsva_out, 'APOE4')
    }
    return(fits)
}

####################################################################################################################
fit.gsva = function(mod1, i, gsva.per.celltype, coef){
fit <- lmFit(t(gsva.per.celltype[[i]]), design=mod1)
fit <- eBayes(fit)
allgenesets <- topTable(fit, coef=coef, number=Inf, confint = T) %>% .[order(.$P.Value, decreasing = F),]
allgenesets$celltype = i
allgenesets$names = rownames(allgenesets)
return(allgenesets)

}

####################################################################################################################
get_scores = function(fits){
    outs = list()
    all = list()
    for(i in names(fits)){
        df = fits[[i]]
        df$celltype = i
        df = df[order(abs(df$logFC),decreasing = T),]
        all[[i]] = df[,c('celltype', 'logFC','P.Value', 'names')]
    }
    return(list('all' = all))
}

####################################################################################################################
get_matrix = function(scores, top_paths){
    df = do.call('rbind',scores)
    df$score = sign(df$logFC) * -log10(df$P.Value)
    df = as.data.frame(pivot_wider(df[,c('celltype', 'names', 'score')], values_from = 'score', names_from = 'celltype'))
    df[is.na(df)] = 0
    rownames(df) = df$names
    df$names = NULL
    return(df[top_paths,])
}
####################################################################################################################
get_stratified_fits = function(stratified_summary, gsva_out_full, model, coefficient){
    out = lapply(names(gsva_out_full), function(x) gsva_out_full[[x]][rownames(gsva_out_full[[x]])%in%as.character(rownames(stratified_summary)),])
    names(out) = names(gsva_out_full)
    fit <- lmFit(t(out[['Oli']]), design=model)
    fit <- eBayes(fit)
    allgenesets <- topTable(fit, coef=coefficient, number=Inf, confint = T) %>% .[order(.$P.Value, decreasing = F),]
    return(allgenesets)
}

####################################################################################################################
# Define function to get linear model fits for each gene set
get_fits = function(gsva_out, meta, covariates, random_effect){


    # Function to get linear model fits for each gene set

    # Parameters:
    
    # gsva_out: A list of gene set variation analysis (GSVA) scores for each gene set
    # meta: A data frame containing metadata for the samples corresponding to the gene set scores
    # covariates: A string specifying the covariates to include in the linear model

    # Returns:
    # A list of linear model fits for each gene set

    # Usage:
    # fits <- get_fits(gsva_out, meta, covariates)

    fits = list()
    
    for(i in names(gsva_out)){

        # Get metadata for the samples corresponding to the gene set scores
        predict = meta[[i]]
        groups_in_meta = unique(predict$pathology.group)

        if (length(unique(predict$pathology.group))>=2){

            # Define the formula based on the selected covariates
            formula="~ 0"
            for(icov in covariates){
                #icov <- gsub(" ", ".", icov)
                if(length(unique(predict[, icov]))>1){
                    if(class(predict[, icov]) == class(factor)){
                        predict[, icov] = as.character(predict[, icov])
                    }
                    formula = paste0(formula, " + ", icov)
                } else {
                    warning(paste0("Excluding ", icov," covariate as it has only one level!"))
                }
            }
            # Generate the design matrix based on the formula and data
            mod <- model.matrix(as.formula(formula), data = predict)

            if ('early' %in% groups_in_meta && 'late' %in% groups_in_meta && 'no' %in% groups_in_meta) {
            contrasts <- makeContrasts(
                pathology.groupearly - pathology.groupno,
                pathology.grouplate - pathology.groupearly,
                pathology.grouplate - pathology.groupno,
                (pathology.groupearly + pathology.grouplate)/2 - pathology.groupno,
                levels=colnames(mod)
            )
            } else if ('early' %in% groups_in_meta && 'no' %in% groups_in_meta && !('late' %in% groups_in_meta)) {
            contrasts <- makeContrasts(pathology.groupearly - pathology.groupno, levels=colnames(mod))
            } else if ('late' %in% groups_in_meta && 'early' %in% groups_in_meta && !('no' %in% groups_in_meta)) {
            contrasts <- makeContrasts(pathology.grouplate - pathology.groupearly, levels=colnames(mod))
            } else if ('late' %in% groups_in_meta && 'no' %in% groups_in_meta && !('early' %in% groups_in_meta)) {
            contrasts <- makeContrasts(pathology.grouplate - pathology.groupno, levels=colnames(mod))
            }

            # Fit the linear model using the gene set scores as the outcome
            fits[[i]] = fit.gsva(mod, i, gsva_out, 'pathology.group', contrasts, meta, random_effect, groups_in_meta)

      } else {
            print('WARNING****************')
            print(paste0("Excluding cell type ",i," as it has only one covariate level!"))
      }
      
    }
    return(fits)
}
   
####################################################################################################################
# Define function to fit linear models and get gene set scores for a single gene set

fit.gsva = function(mod1, i, gsva.per.celltype, coef, contrasts, meta, random_effect, groups_in_meta) {
    # This function fits linear models and obtains gene set scores for a single gene set.

    # Args:
    # mod1: A model matrix specifying the linear model design, created using the metadata for the samples corresponding to the gene set scores.

    # i: A character string indicating the name of the gene set being analyzed.
    # gsva.per.celltype: A list containing gene set scores for all cell types.
    # coef: A character string specifying the linear model coefficient to be used for ranking gene sets based on P-values.
    # contrasts: A contrast matrix to be used in the linear model.

    # Returns:
    # A list containing a table of gene set results for each contrast and each cell type.

    # Fit the linear model with the gene set scores as the outcome
    expression <- gsva.per.celltype[[i]]
    meta <- meta[[i]]
    
    if (random_effect != 'None') {
      cor <- duplicateCorrelation(expression, mod1, block = as.character(meta[, random_effect]))
      fit <- if (!is.nan(cor$consensus.correlation) && abs(cor$consensus.correlation) > 0.01) 
                  lmFit(expression, design = mod1, block = as.character(meta[, random_effect]), correlation = cor$consensus.correlation)
              else
                  lmFit(expression, design = mod1)
    } else {
        fit <- lmFit(expression, design = mod1)
    }
    
    fit <- contrasts.fit(fit, contrasts)
    fit <- eBayes(fit)

    # Extract all genesets, ranked by their P-values
    allgenesets <- list()
    
    if ('early' %in% groups_in_meta && 'late' %in% groups_in_meta && 'no' %in% groups_in_meta) {

        allgenesets[['early_vs_no']] <- topTable(fit, 
                                                adjust.method = "BH",
                                                coef = "pathology.groupearly - pathology.groupno", 
                                                number = Inf,
                                                confint = TRUE) %>% arrange(P.Value)

        allgenesets[['late_vs_early']] <- topTable(fit, 
                                                    adjust.method = "BH", 
                                                    coef = "pathology.grouplate - pathology.groupearly", 
                                                    number = Inf, 
                                                    confint = TRUE) %>% arrange(P.Value)

        allgenesets[['late_vs_no']] <- topTable(fit, 
                                                adjust.method = "BH", 
                                                coef = "pathology.grouplate - pathology.groupno", 
                                                number = Inf, 
                                                confint = TRUE) %>% arrange(P.Value)

        allgenesets[['ad_vs_no']] <- topTable(fit, 
                                              adjust.method = "BH", 
                                              coef = "(pathology.groupearly + pathology.grouplate)/2 - pathology.groupno", 
                                              number = Inf, 
                                              confint = TRUE) %>% arrange(P.Value)

    } else if ('early' %in% groups_in_meta && 'no' %in% groups_in_meta && !('late' %in% groups_in_meta)) {

        allgenesets[['early_vs_no']] <- topTable(fit, 
                                                adjust.method = "BH", 
                                                coef = "pathology.groupearly - pathology.groupno", 
                                                number = Inf, 
                                                confint = TRUE) %>% arrange(P.Value)

    } else if ('late' %in% groups_in_meta && 'early' %in% groups_in_meta && !('no' %in% groups_in_meta)) {

        allgenesets[['late_vs_early']] <- topTable(fit, 
                                                    adjust.method = "BH", 
                                                    coef = "pathology.grouplate - pathology.groupearly", 
                                                    number = Inf, 
                                                    confint = TRUE) %>% arrange(P.Value)

    } else if ('late' %in% groups_in_meta && 'no' %in% groups_in_meta && !('early' %in% groups_in_meta)) {

        allgenesets[['late_vs_no']] <- topTable(fit,
                                                adjust.method = "BH", 
                                                coef = "pathology.grouplate - pathology.groupno", 
                                                number = Inf, 
                                                confint = TRUE) %>% arrange(P.Value)
    }
    
    # Add the celltype and gene set names to the output table
    for (j in names(allgenesets)) {
        allgenesets[[j]]$celltype = i
        allgenesets[[j]]$names = rownames(allgenesets[[j]])
    }
    
    return(allgenesets)
}



####################################################################################################################
#' Fit Linear Models for Pathway Activity Scores
#' 
#' @description
#' Analyzes pathway activity scores across cell types using linear modeling, with 
#' optional contrast analysis between early and late stages.
#'
#' @param pa_scores Dataframe containing pathway activity scores
#' @param meta Dataframe with sample metadata
#' @param covariates Vector of covariate names to include in model
#' @param random_effect Name of random effect variable (or 'None')
#' @param cell_type Current cell type being analyzed
#' @param use_contrasts Boolean indicating whether to perform contrast analysis
#'
#' @return List of linear model fits and results for the cell type
get_fits_cps = function(pa_scores, meta, covariates, random_effect, cell_type, use_contrasts=FALSE) {

    i = cell_type
    predict = meta

    # Build formula dynamically from covariates
    formula = "~ "
    for(icov in covariates) {
        # Skip covariates with only one level
        if(length(unique(predict[, icov])) > 1) {
            if(is.factor(class(predict[, icov]))) {
                predict[, icov] = as.character(predict[, icov])
            }
            formula = paste0(formula, " + ", icov)
        } else {
            warning(paste0("Excluding ", icov," covariate as it has only one level!"))
        }
    }


    # Create design matrix and fit model
    mod <- model.matrix(as.formula(formula), data = predict)

    if (use_contrasts) {
        contrasts <- makeContrasts(CPSlate - CPSearly, levels=colnames(mod))
        fits = fit.cps(mod, i, pa_scores, meta, contrasts, random_effect)
    } else {
        fits = fit.cps(mod, i, pa_scores, meta, contrasts=NULL, random_effect)
    }
    return(fits)
}

####################################################################################################################
#' Fit Linear Model for Single Cell Type
#' 
#' @description
#' Performs detailed linear modeling for a single cell type, including empirical
#' Bayes moderation and optional random effects/contrasts analysis.
#'
#' @param mod1 Design matrix for linear model
#' @param i Current cell type name
#' @param pa_scores Matrix of pathway activity scores
#' @param meta Sample metadata
#' @param contrasts Contrast matrix for differential analysis
#' @param random_effect Random effect variable name
#'
#' @return List containing:
#'   - CPS: Main analysis results
#'   - late_vs_early: Contrast analysis (if requested)
fit.cps = function(mod1, i, pa_scores, meta, contrasts, random_effect) {
    expression <- pa_scores

    # handle random effects if specified
    if (!is.null(random_effect)) {
        # calculate correlation for repeated measurements
        cat("       Fetching duplicate correlation\n")
        cor <- duplicateCorrelation(expression, mod1, 
                                  block = as.character(meta[, random_effect]))
        
        # Only use correlation if significant
        fit <- if (!is.nan(cor$consensus.correlation) && 
                abs(cor$consensus.correlation) > 0.01) {
            cat("       Fitting model with random effect correction with correlation:", round(cor$consensus.correlation, 3), "\n")
            lmFit(expression, design = mod1, 
                block = as.character(meta[, random_effect]), 
                correlation = cor$consensus.correlation)
        } else {
            cat("       Correlation not significant, using standard fit\n")
            lmFit(expression, design = mod1)
        }
    } else {
        cat("       Running Standard fit\n")
        fit <- lmFit(expression, design = mod1)
    }

    # standard errors
    coef_summary <- data.frame(StdErr = sqrt(fit$stdev.unscaled))

    # empirical Bayes moderation
    fit_ebayes <- eBayes(fit)

    # get results
    results <- list()
    results[['CPS']] <- topTable(fit_ebayes, adjust.method = "BH", 
                                number = Inf, confint = TRUE) %>% 
                        arrange(P.Value)

    # add metadata
    results[['CPS']]$celltype = i
    results[['CPS']]$names = rownames(results[['CPS']])
    results[['CPS']] <- merge(results[['CPS']], 
                             coef_summary,
                             by = "row.names", 
                             all = FALSE,
                             all.x = TRUE) %>% 
                        arrange(P.Value)

    # contrast analysis if requested
    if (!is.null(contrasts)) {
        fit_contrasts <- contrasts.fit(fit, contrasts)
        fit_contrasts <- eBayes(fit_contrasts)
        results[['late_vs_early']] <- topTable(fit_contrasts, 
                                             adjust.method = "BH", 
                                             coef = "CPSlate - CPSearly",  
                                             number = Inf, 
                                             confint = TRUE) %>% 
                                    arrange(P.Value)
        results[['late_vs_early']]$celltype = i
        results[['late_vs_early']]$names = rownames(results[['late_vs_early']])
    }

    return(results)
}


####################################################################################################################
# Define function to get a data frame of gene set scores for all gene sets
get_scores = function(fits){

    # This function takes a list of fitted linear models for each gene set and returns a data frame of gene set scores for all gene sets.

    # Args:
    # fits: A list of fitted linear models for each gene set generated by the function get_fits().

    # Returns:

    # A list containing the table of gene set scores for all gene sets. The output is a list with one element, 
    # 'all', which is a list of lists. Each nested list corresponds to a gene set and contains a data frame of gene set scores for all cell types.
    
    # The data frame has the following columns:
    #     - 'celltype': A character string representing the cell type.
    #     - 'logFC': A numeric value representing the log-fold change of the gene set score for the given cell type and gene set.
    #     - 'P.Value': A numeric value representing the P-value of the moderated t-test for the given cell type and gene set.
    #     - 'names': A character string representing the gene set name.


    outs = list()
    all = list()
    for(i in names(fits)){
        all[[i]] = list()
        for(j in names(fits[[i]])){
            # Extract gene set scores for a single gene set
            df = fits[[i]][[j]]
            # Add the celltype to the output table
            df$celltype = i
            # Sort the table by absolute log fold-change, in descending order
            df = df[order(abs(df$logFC),decreasing = T),]
            all[[i]][[j]] = df[,c('celltype', 'logFC', 'adj.P.Val', 'P.Value', 'names')]
        }
    }
    # Return a list containing the table of gene set scores for all gene sets
    return(list('all' = all))
}

####################################################################################################################
# Define function to get a matrix of gene set scores for the top gene sets, in a format suitable for heatmap plotting
get_matrix = function(scores, top_paths){


    # Define function to get a matrix of gene set scores for the top gene sets, in a format suitable for heatmap plotting.

    # Parameters:

    # - scores: list
    #     List containing the table of gene set scores for all gene sets.
    #     Each entry in the list corresponds to one cell type, and contains a list of data frames, with one data frame for each gene set.

    # - top_paths: integer
    #     Number of top gene sets to include in the output matrix, sorted by their absolute log fold-change values in descending order.

    # Returns:
    # - matrix: matrix
    #     Matrix of gene set scores for the top gene sets, with one row per gene set and one column per cell type. 
    #     The rows are sorted by their absolute log fold-change values in descending order, 
    #     and the columns are sorted alphabetically by cell type name.


    # Combine the gene set scores for all gene sets into a single data frame
    df = do.call('rbind', scores)
    # Transform the scores into a matrix, with one row per gene set and one column per cell type
    df$score = sign(df$logFC) * -log10(df$P.Value)
    
    # Check if there's only one unique celltype
    if(length(unique(df$celltype)) == 1) {
    unique_celltype <- unique(df$celltype)[[1]]
    df[[unique_celltype]] <- df$score
    df$score <- NULL
    df$celltype <- NULL
    } else {
    # If there's more than one celltype, use pivot_wider as before
    df = as.data.frame(pivot_wider(df[,c('celltype', 'names', 'score')], values_from = 'score', names_from = 'celltype'))
    }
    # Replace missing values with 0

    df[is.na(df)] = 0
    # Use the gene set names as row names
    rownames(df) = df$names
    df$names = NULL
    
    return(df[top_paths,])
}


# Function to generate a heatmap with row annotations based on gene expression data

# Arguments:
# - matrix: Numeric matrix containing gene expression data
# - names: Character vector of gene names
# - path_names_go: Data frame mapping gene names to Gene Ontology (GO) terms
# - path_slot: Column index specifying the GO term to be used for annotation
# - cell_types: Column indices specifying the cell types to be included in the heatmap

get_go_hmap = function(matrix, names, path_names_go, path_slot, cell_type_order, reorder=FALSE) {
  
    # Create a temporary data frame based on the absolute value of matrix elements, 
    # with values set to 1 or -1 depending on whether they exceed a threshold of 1.3
    temp = as.data.frame((abs(matrix) > 1.3) * sign(matrix))
    
    # Reorder the rows of the matrix based on the counts of different cell types

    ordering = 'order('
    for (ct in cell_type_order){
        ordering <- paste0(ordering, '(temp$', ct, '),')
    }

    ordering <- paste0(ordering, 'decreasing = T)')
    ordering_expr <- parse(text = ordering)

    temp_order = rownames(temp[eval(ordering_expr), ])

    # temp_order = rownames(temp[order((temp$Excitatory),
    #                             (temp$Inhibitory), 
    #                             (temp$Astrocyte), 
    #                             (temp$Microglia), 
    #                             (temp$Oligodendrocyte), 
    #                             (temp$OPC), decreasing = T),])


    matrix = matrix[temp_order,]
    
    if (reorder){
        matrix = matrix[sort(rownames(matrix)), ]
    }
    
    # Find the index of each pathway name to be annotated in the reordered matrix
    index = unname(unlist(lapply(names, function(x) which(rownames(matrix) == x))))
    
    # Create row annotation for the heatmap using gene names and corresponding GO terms
    ha = rowAnnotation(foo = anno_mark(at = index, labels = path_names_go[rownames(matrix), path_slot][index]), 
                       gp = gpar(fontsize = 30))


    # Generate the heatmap using the specified matrix and cell types, with the row annotation
    h = Heatmap(as.matrix(matrix[, cell_type_order]),
                right_annotation = ha, 
                col = colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red')), 
                show_row_names = F, 
                border = T, 
                cluster_rows = F, 
                cluster_columns = F, 
                row_km = 0,
                rect_gp = gpar(col = 'black', lwd = .5)
               )
    
    # Return the generated heatmap
    return(h)
}

#' Convert Seurat object to SingleCellExperiment with merged metadata
#' @param seurat_obj Seurat object to convert
#' @param adata_annot AnnData object with annotations
#' @param cell_type_column Name of the cell type column
#' @param subject_id Name of the subject ID column
#' @return SingleCellExperiment object with merged metadata
process_metacell_data <- function(seurat_obj, adata_annot, cell_type_column, subject_id) {

  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required. Please install it.")
  }
  
  if (missing(seurat_obj)) stop("seurat_obj is required")
  if (missing(adata_annot)) stop("adata_annot is required")
  if (missing(cell_type_column)) stop("cell_type_column is required")
  if (missing(subject_id)) stop("subject_id is required")

  tryCatch({
    sce <- as.SingleCellExperiment(GetMetacellObject(seurat_obj))
    
    metadata_cols <- setdiff(colnames(colData(adata_annot)), 
                           cell_type_column)
    
    metadata_raw <- colData(adata_annot)[, metadata_cols] %>%
      as.data.frame() %>%
      subset(!duplicated(.[, subject_id]), )
    
    metacell_coldata <- as.data.frame(colData(sce))
    
    merged_rowdata <- merge(
      x = metacell_coldata,
      y = metadata_raw,
      by.x = subject_id,
      by.y = subject_id,
      all.x = TRUE
    )
    
    if (nrow(merged_rowdata) != nrow(colData(sce))) {
      warning("Merge resulted in different number of rows. Check your subject IDs.")
    }
    
    rownames(merged_rowdata) <- rownames(colData(sce))
    colData(sce) <- DataFrame(merged_rowdata)
    
    return(sce)
    
  }, error = function(e) {
    stop(paste("Error processing metacell data:", e$message))
  })
}


find_parent_key <- function(cell_type, subclass) {
  for (key in names(subclass)) {
    if (cell_type %in% subclass[[key]]) {
      return(key)
    }
  }
  return(NA)  # Return NA if no match found
}

clean_strings <- function(covariates, separator = "_", preserve_case = FALSE, replace_hyphen = TRUE) {
  clean_string <- function(text) {
    # Choose pattern based on replace_hyphen flag
    pattern <- if(replace_hyphen) {
      "[^[:alnum:][:space:]]"  # Replace all special chars including hyphens
    } else {
      "[^[:alnum:][:space:]-]"  # Replace all special chars except hyphens
    }
    
    # Replace special chars and spaces with separator
    cleaned <- gsub(pattern, separator, text)
    cleaned <- gsub("\\s+", separator, cleaned)
    cleaned <- gsub(paste0(separator, "+"), separator, cleaned)
    
    # Remove leading/trailing separators
    cleaned <- trimws(cleaned, whitespace = separator)
    
    if (!preserve_case) {
      cleaned <- tolower(cleaned)
    }
    return(cleaned)
  }
  
  if (length(covariates) == 1 && is.character(covariates)) {
    return(clean_string(covariates))
  }
  return(sapply(covariates, clean_string))
}


#' Run differential analysis for a single gene set, method and cell type combination
#' @description Performs differential analysis using pathway activity scores, accounting for random effects
#' and covariates while handling potential errors gracefully.
#'
#' @param pas_object SingleCellExperiment object containing pathway activity scores
#' @param gene_set_db Named list of gene set collections (e.g., 'gabitto', KEGG)
#' @param gene_set String indicating which gene set to analyze (e.g., "gabitto")
#' @param method String specifying the pathway activity scoring method used (e.g 'AUCell')
#' @param cell_type String specifying the cell type to analyze
#' @param test_categories Named list of test category filters
#' @param test_cat String indicating which test category filter to apply
#' @param cell_type_column String specifying the column name containing cell type information
#' @param random_effect String specifying the column name for random effect (e.g., "subject_id")
#' @param factor String specifying the column name for the continuous factor of interest
#' @param covariates Character vector of column names to include as covariates
#' @param numeric_covariates Character vector of covariates that should be treated as numeric
#' @param data_dir String specifying path to data directory containing expressed genes
#' @param pathway_gene_threshold Numeric threshold for filtering pathways (default: 0.33)
#'
#' @return List containing fitted model results for each pathway
#' @export


run_differential_analysis <- function(
    pas_object,
    gene_set_db,
    gene_set,
    method,
    cell_type,
    test_categories,
    test_cat,
    cell_type_column,
    random_effect,
    factor = "continuous_pseudoprogression_score",
    covariates,
    numeric_covariates,
    data_dir,
    pathway_gene_threshold = 0.33
) {
    # Input validation
    if (!is.character(gene_set) || length(gene_set) != 1) {
        stop("gene_set must be a single character string")
    }
    if (!gene_set %in% names(gene_set_db)) {
        stop("gene_set not found in gene_set_db")
    }
    
    # tryCatch({
        # Read gene set terms
        pw_terms <- read.geneset(gene_set_db[[gene_set]])
        
        # Extract summary data
        col_data <- colData(pas_object[[gene_set]][[method]])
        pa_scores <- as.data.frame(Matrix::as.matrix(assay(pas_object[[gene_set]][[method]], 'X')))
        
        # Get expressed genes
        expressed_genes_path <- file.path(
            data_dir,
            "anndata",
            find_parent_key(cell_type, subclass),
            clean_strings(cell_type, preserve_case=TRUE),
            paste0(clean_strings(cell_type, preserve_case=TRUE), "_pseudobulked_expressed_genes.rds")
        )
        
        if (!file.exists(expressed_genes_path)) {
            stop(paste("Expressed genes file not found:", expressed_genes_path))
        }
        
        expressed_genes <- readRDS(expressed_genes_path)
        
        # Filter data based on test category
        if (!test_cat %in% names(test_categories)) {
            stop(paste("Test category", test_cat, "not found in test_categories"))
        }
        
        full_expression <- paste0("col_data$", test_categories[[test_cat]])
        summary <- col_data[
            col_data[[cell_type_column]] == cell_type &
            eval(parse(text = full_expression)),
        ]
        
        if (nrow(summary) == 0) {
            stop("No samples remain after filtering")
        }
        
        # Process random effects and factors
        if(!is.null(random_effect)){
            summary[[random_effect]] <- as.factor(summary[[random_effect]])
        }
        summary[['CPS']] <- summary[[factor]]
        
        # Convert numeric covariates
        if (!identical(numeric_covariates, "None")) {
            for (cov in numeric_covariates) {
                if (!cov %in% colnames(summary)) {
                    stop(paste("Covariate not found in data:", cov))
                }
                summary[[cov]] <- as.numeric(summary[[cov]])
            }
        }
        
        # Filter pathways and get scores
        filtered_pathways <- names(filter_genesets(expressed_genes, pw_terms, pathway_gene_threshold))
        
        if (length(filtered_pathways) == 0) {
            stop("No pathways remain after filtering")
        }
        
        scores <- pa_scores[filtered_pathways, rownames(summary)]
        
        # Fit models
        fits <- get_fits_cps(scores, summary, covariates, random_effect, cell_type, use_contrasts=FALSE)

        return(fits)
        
    # }, error = function(e) {
    #     if (grepl("       F-statistics not found in fit", e$message)) {
    #         stop(paste("       F-statistics error for", cell_type, ":", e$message))
    #     } else if (grepl("       undefined columns selected", e$message)) {
    #         stop(paste("       Column matching error for", cell_type, ":", e$message))
    #     } else {
    #         stop(paste("       Error in analysis for", cell_type, ":", e$message))
    #     }
    # })
}


get.fits <- function(pas_object, curated_gene_set_db, pas_methods, cell_types, 
                    test_categories, cell_type_column, factor, covariates, 
                    numeric_covariates, data_dir, pathway_gene_threshold = 0.33) {
  fit_result <- list()
  
  safe_run_differential <- purrr::safely(run_differential_analysis)
  
  for (gs in names(curated_gene_set_db)) {
    fit_result[[gs]] <- purrr::map(pas_methods, function(method) {
      cat(sprintf("\nProcessing %s with method %s\n", gs, method))
      
      purrr::map(cell_types, function(cell_type) {
        cat(sprintf("  Cell type: %s\n", cell_type))
        
        purrr::map(names(test_categories), function(test_cat) {
          cat(sprintf("    Category: %s\n", test_cat))
          
          result <- safe_run_differential(
            pas_object = pas_object,
            gene_set_db = curated_gene_set_db,
            gene_set = gs,
            method = method,
            cell_type = cell_type,
            test_categories = test_categories,
            test_cat = test_cat,
            cell_type_column = cell_type_column,
            random_effect = NULL,
            factor = factor,
            covariates = covariates,
            numeric_covariates = numeric_covariates,
            data_dir = data_dir,
            pathway_gene_threshold = pathway_gene_threshold
          )
          
          if (!is.null(result$error)) {
            cat(sprintf("    Error in %s-%s-%s: %s\n", 
                       method, cell_type, test_cat, result$error))
            return(NULL)
          }
          
          result$result
        }) %>% setNames(names(test_categories))
      }) %>% setNames(cell_types)
    }) %>% setNames(pas_methods)
  }
  
  fit_result
}


convert_pathway_list_to_df <- function(pathway_list) {
  all_genes <- vector()
  all_modules <- vector()
  
  for(module_name in names(pathway_list)) {
    genes <- pathway_list[[module_name]]
    # sdd each gene-pathway pair
    all_genes <- c(all_genes, genes)
    all_modules <- c(all_modules, rep(module_name, length(genes)))
  }
  
  result_df <- data.frame(
    gene_name = all_genes,
    module = factor(all_modules),
    stringsAsFactors = FALSE
  )
}
# 
