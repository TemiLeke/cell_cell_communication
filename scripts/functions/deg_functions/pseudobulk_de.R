#' Run pseudobulk differential expression methods
#' 
#' Run pseudobulk differential expression methods on single-cell data
#' 
#' @param input a single-cell matrix to be converted, with features (genes) in rows
#'   and cells in columns. Alternatively, a \code{Seurat}, \code{monocle3}, or 
#'   or \code{SingleCellExperiment} object can be directly input.
#' @param meta the accompanying meta data whereby the rownames match the column
#'   names of \code{input}.
#' @param replicate_col the vector in \code{meta} containing the replicate 
#'   information. Defaults to \code{replicate}.
#' @param cell_type_col the vector in \code{meta} containing the cell type 
#'   information. Defaults to \code{cell_type}.
#' @param label_col the vector in \code{meta} containing the experimental
#'   label. Defaults to \code{label}. 
#' @param min_cells the minimum number of cells in a cell type to retain it.
#'   Defaults to \code{3}.
#' @param min_reps the minimum number of replicates in a cell type to retain it.
#'   Defaults to \code{2}.
#' @param min_features the minimum number of expressing cells (or replicates) 
#'   for a gene to retain it. Defaults to \code{0}.
#' @param de_method the specific differential expression testing method to use.
#'   Defaults to edgeR.
#' @param de_type the specific parameter of the differential expression testing
#'   method. Defaults to LRT for edgeR, LRT for DESeq2, and trend for limma.
#' @return a data frame containing differential expression results.
#'  
#' @importFrom magrittr %<>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr %>% mutate n_distinct
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmQLFit glmQLFTest
#'   glmFit glmLRT topTags cpm
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom purrr map 
#' @importFrom stats model.matrix
#' @importFrom methods new
#' 
#' 



pseudobulk_de = function(input, 
                         meta = NULL, 
                         replicate_col = 'replicate',
                         cell_type_col = 'cell_type',
                         label_col = 'label',
                         test_name = 'ad_vs_no',
                         ref_level = 'no',
                         latent_vars = NULL,
                         min_cells = 3,
                         min_reps = 2,
                         min_features = 0,
                         de_family = 'pseudobulk',
                         de_method = 'edgeR',
                         de_type = 'LRT') {
  # check args
  if (de_method == 'limma') {
    if (de_type != 'voom') {
      # change default type to use
      de_type = 'trend'  
    }
  }

  
  # get the pseudobulks list
  pseudobulks = to_pseudobulk(
    input = input,
    meta = meta,
    replicate_col = replicate_col,
    cell_type_col = cell_type_col,
    label_col = label_col,
    latent_vars = latent_vars,
    min_cells = min_cells,
    min_reps = min_reps,
    min_features = min_features,
    external = F
  )


  results = map(pseudobulks, function(x) {
    # create targets matrix

    targets = data.frame(group_sample = colnames(x)) %>%
      mutate(group = gsub(".*\\:", "", group_sample)) %>%
      mutate(replicate = as.character(gsub(":.*", "", group_sample)))

    
    # # Ensure the 'replicate' column in 'meta' is also character
    # meta$replicate <- as.character(meta$replicate)

    # # Merge the targets with meta based on 'replicate'
    # targets <- merge(targets, meta, by = "replicate", all.x = TRUE)

    # Map values from meta into targets for each variable in latent_vars
    for (var in latent_vars) {
      mapping_dict <- setNames(meta[[var]], meta[[replicate_col]])
      targets[var] <- map_chr(targets$replicate, ~ mapping_dict[.x])
    }
    # Convert latent_vars columns to factors
    targets <- targets %>% 
      mutate(across(all_of(latent_vars), as.factor))

    # temp_meta <- distinct(meta, replicate, .keep_all = TRUE)
    # temp_meta$replicate = temp_meta[, replicate_col]
    # # temp_meta <- mutate_all(temp_meta, as.factor)
    # targets <- merge(targets, temp_meta, by = "replicate", all.x = FALSE)

    # optionally, carry over factor levels from entire dataset
    if (is.factor(meta$label)){
      targets$group %<>% factor(levels = levels(meta$label))
    }
    # if (n_distinct(targets$group) > 2)
    #   return(NULL)
    
    targets$group <- relevel(targets$group, ref=ref_level)
  
    formula <- as.formula(paste0("~ ", paste(c('group', latent_vars), collapse = " + ")))

    # create design
    design <- model.matrix(formula, data = targets)

    # design = model.matrix(~ group, data = targets)
    
    DE = switch(de_method,
                edgeR = {
                  tryCatch({
                    y = DGEList(counts = x, group = targets$group) %>%
                      calcNormFactors(method = 'TMM') %>%
                      estimateDisp(design)

                    if (test_name == 'ad_vs_no'){
                      contrasts <- makeContrasts(0.5*(grouplate + groupearly), levels=design)
                    } else if (test_name == 'late_vs_no'){
                      contrasts <- makeContrasts(grouplate, levels=design)
                    } 
                     else if (test_name == 'late_vs_early'){
                      contrasts <- makeContrasts(grouplate - groupearly, levels=design)
                    }
                    else if (test_name == 'early_vs_no'){
                      contrasts <- makeContrasts(groupearly, levels=design)
                    }
                    else if (test_name == 'ad_vs_no_binary'){
                      contrasts <- makeContrasts(groupad, levels=design)
                    }
                    test = switch(de_type,
                                  QLF = {
                                    fit = glmQLFit(y, design)
                                    test = glmQLFTest(fit, contrast = contrasts) #coef = -1)
                                  },
                                  LRT = {
                                    fit = glmFit(y, design = design)
                                    test = glmLRT(fit, contrast = contrasts)
                                  })
                    res = topTags(test, n = Inf) %>%
                      as.data.frame() %>%
                      rownames_to_column('gene') %>%
                      # flag metrics in results
                      mutate(de_family = 'pseudobulk',
                             de_method = de_method,
                             de_type = de_type)
                  }, error = function(e) {
                    message(e)
                    data.frame()
                  })
                },
                DESeq2 = {
                  tryCatch({
                    dds = DESeqDataSetFromMatrix(countData = x,
                                                 colData = targets,
                                                 design = ~ group)
                    dds = switch(de_type,
                                 Wald = {
                                   dds = try(DESeq(dds,
                                                   test = 'Wald',
                                                   fitType = 'parametric',
                                                   sfType = 'poscounts',
                                                   betaPrior = F))
                                 },
                                 LRT = {
                                   dds = try(DESeq(dds,
                                                   test = 'LRT',
                                                   reduced = ~ 1,
                                                   fitType = 'parametric',
                                                   sfType = 'poscounts',
                                                   betaPrior = F))
                                 }
                    )

                    mod_mat <- model.matrix(design(dds), colData(dds))
               
                    if (test_name == 'ad_vs_no'){
                            ref_level <- colMeans(mod_mat[dds$group == 'no', ])
                            test_level <- colMeans(mod_mat[dds$group %in% c("early", "late"),])
                    } else if (test_name == 'ad_vs_no_binary'){
                            ref_level <- colMeans(mod_mat[dds$group == 'no', ])
                            test_level <- colMeans(mod_mat[dds$group == "ad", ])
                    } else if (test_name %in% c('late_vs_no', 'late_vs_early', 'early_vs_no')){
                            ref_level <- colMeans(mod_mat[dds$group == strsplit(test_name, "_vs_")[[1]][2], ])
                            test_level <- colMeans(mod_mat[dds$group == strsplit(test_name, "_vs_")[[1]][1], ])
                    }
                    # calculate coefficient vectors for each group
                    res = results(dds, contrast = test_level - ref_level)
                    # write
                    res = as.data.frame(res) %>%
                      mutate(gene = rownames(x)) %>%
                      # flag metrics in results
                      mutate(de_family = 'pseudobulk',
                             de_method = de_method,
                             de_type = de_type)
                  }, error = function(e) {
                    message(e)
                    data.frame()
                  })
                },
                limma = {
                  tryCatch({
                    x = switch(de_type,
                               trend = {
                                 trend_bool = T
                                 dge = DGEList(as.matrix(x), group = targets$group)
                                 dge = calcNormFactors(dge)
                                 x = new("EList")
                                 x$E = cpm(dge, log = TRUE, prior.count = 3)
                                 x
                               },
                               voom = {
                                 counts = all(as.matrix(x) %% 1 == 0)
                                 if (counts) {
                                   trend_bool = F
                                   x = voom(as.matrix(x), design)
                                   x
                                 }
                               })


                    if (test_name == 'ad_vs_no'){
                      contrasts <- makeContrasts((grouplate + groupearly)/2, levels=design)
                      coef = "(grouplate + groupearly)/2"
                    } else if (test_name == 'ad_vs_no_binary'){
                      contrasts <- makeContrasts(groupad, levels=design)
                      coef = "groupad"
                    } else if (test_name == 'late_vs_no'){
                      contrasts <- makeContrasts(grouplate, levels=design)
                      coef = "grouplate"
                    } 
                     else if (test_name == 'late_vs_early'){
                      contrasts <- makeContrasts(grouplate - groupearly, levels=design)
                      coef = "grouplate - groupearly"
                    }
                    else if (test_name == 'early_vs_no'){
                      contrasts <- makeContrasts(groupearly, levels=design)
                      coef <- "groupearly"
                    }

                    # get fit
                    fit = lmFit(x, design)
                    fit = contrasts.fit(fit, contrasts)
                    fit = eBayes(fit, trend = trend_bool, robust = trend_bool)

                    # format the results
                    res = fit %>%
                      # extract all coefs except intercept
                      topTable(number = Inf, coef = coef) %>%
                      rownames_to_column('gene') %>%
                      # flag metrics in results
                      mutate(
                        de_family = 'pseudobulk',
                        de_method = de_method,
                        de_type = de_type)
                  }, error = function(e) {
                    message(e)
                    data.frame()
                  })
                }
    )
  })
  results %<>% bind_rows(.id = 'cell_type')
}
