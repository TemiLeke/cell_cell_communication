
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
      fit <- if (!is.nan(cor$consensus.correlation) && abs(cor$consensus.correlation) < 0.01) 
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
