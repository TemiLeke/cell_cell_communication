import scanpy as sc
import pandas as pd
import glob
import re
import anndata as ad
import warnings
import os
import numpy as np
from joblib import Parallel, delayed
import warnings
from datetime import datetime
from igraph import *
import rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import argparse
from functions.helper_functions import clean_strings

sc.settings.n_jobs = 40

# parser = argparse.ArgumentParser()
# parser.add_argument("--threads", help="threads")
# parser.add_argument("--adata", help="adata")
# parser.add_argument("--target", help="target")
# parser.add_argument("--region", help="region")
# parser.add_argument("--split_key", help="split_key")
# parser.add_argument("--outgroup_de", help="outgroup_de")
# parser.add_argument("--covariates", help="covariates")
# parser.add_argument("--random_effect", help="random_effect")
# parser.add_argument("--covariate_formula", help="covariate_formula")
# parser.add_argument("--tests", help="tests")
# parser.add_argument("--layer", help="layer")
# parser.add_argument("--offset_variable", help="offset_variable")
# parser.add_argument("--group_cell", help="group_cell")
# parser.add_argument("--save_dir", help="save_dir")

# args = parser.parse_args()
# pwd = os.getcwd()
warnings.filterwarnings("ignore")
warnings.filterwarnings('ignore', category=FutureWarning, module='anndata')
warnings.filterwarnings('ignore', message='.*reticulate.*')


def run_de_tests(sub, region, target, split_key, outgroup_de, outgroup, covariates, 
                 random_effect, covariate_formula, tests, layer, offset_variable, group_cell, save_dir):
    
    """
    Performs differential expression analysis using the NEBULA method, which is designed for single-cell RNA-seq data 
    with complex experimental designs.

    Parameters
    ----------
    sub : AnnData
        AnnData object containing the single-cell RNA-seq data
    region : str
        Name of the region/area to analyze, used for output file organization
    target : str
        Target identifier/name for the analysis
    split_key : str
        Column name in sub.obs that contains the grouping variable for comparison
    outgroup_de : bool
        Whether to compare differential expression between groups of cells/supertypes
    outgroup : str
        The specific group to compare against others
    covariates : list
        List of covariate names to include in the model
    random_effect : str
        Column name in sub.obs specifying the random effect (e.g., subject ID)
    covariate_formula : str
        Formula specifying how covariates should be included in the model
    tests : list
        List of additional tests to perform
    layer : str
        Layer in AnnData object to use for analysis. Use "X" for main layer
    offset_variable : str
        Column name in sub.obs containing the offset variable (e.g., size factors)
    group_cell: bool
        Whether to group cells by subject id using the nebula function
    save_dir: str
        Parent folder/directory for saving results.

    Returns
    -------
    None
        Results are saved to disk in CSV format

    File Outputs
    -----------
    - {outgroup}_vs_all_DE.csv : Results of comparing outgroup vs all other groups
    - {outgroup}_across_{test}_DE.csv : Results of additional specified tests
    - Associated RDS files containing full NEBULA model results

    Notes
    -----
    - Automatically handles data preprocessing including:
      - Removal of low-count genes
      - Min-max normalization of numeric covariates
      - Addition of pseudocounts where necessary
    - Implements quality control by filtering results based on:
      - Model convergence
      - Overdispersion parameters
    - Uses R's NEBULA package through rpy2 interface
    - Supports parallel processing through parent script

    Example Usage
    ------------
    run_de_tests(
        sub=adata,
        region="cortex",
        target="Astrocyte",
        split_key="cell_type",
        outgroup_de=False,
        outgroup="excitatory",
        covariates=["age", "sex"],
        random_effect="subject",
        covariate_formula="age + sex",
        tests=["condition"],
        layer="counts",
        offset_variable="size_factors",
        group_cell=True,
        save_dir = '/documents/'

    )
    """
    
    r_Matrix = importr("Matrix")
    ro.r("library(stats)")
    ro.r("library(nebula)")
    ro.r("library(here)")
    
    if outgroup_de:
        if os.path.exists(os.path.join(save_dir, outgroup, "versus_all", clean_strings(outgroup, "_",  True) + "_vs_all_DE.csv")) is False:
        
            print(str(datetime.now()) + " -- Starting " + clean_strings(outgroup, "_",  True) + " versus all", flush=True)
            sub.obs["comparison"] = "0"
            sub.obs.loc[sub.obs[split_key] == outgroup, "comparison"] = "1"

            if layer != "X":
                tmp = sub.layers[layer].T.tocoo()
            else:
                tmp = sub.X.T.tocoo()

            
            counts = r_Matrix.sparseMatrix(
                i=ro.IntVector(tmp.row + 1),
                j=ro.IntVector(tmp.col + 1),
                x=ro.FloatVector(tmp.data),
                dims=ro.IntVector(tmp.shape)
            )
            
            with localconverter(ro.default_converter + pandas2ri.converter):
                obs = ro.conversion.py2rpy(sub.obs)
            
            ro.r.assign("obs", obs)
            ro.r.assign("counts", counts)
            ro.r("covariates <- c('comparison', " + str(covariates).replace("[", "").replace("]", "") + ")")
            ro.r("df <- model.matrix(~" + covariate_formula + "comparison, data=obs[,covariates])")

            with localconverter(ro.default_converter + pandas2ri.converter):
                df = ro.r("as.data.frame(df)")
                    
            ro.r("data_g <- group_cell(count=counts, id=obs$" + random_effect + ", offset=obs$" + offset_variable + ", pred=df)")
            ordered = ro.r("data_g")
            if (group_cell is True) and (isinstance(ordered, rpy2.rinterface_lib.sexp.NULLType) is True):
                ro.r("re <- nebula(counts, obs$" + random_effect + ", offset=obs$" + offset_variable + ", pred=df, covariance=TRUE)")
                ro.r("saveRDS(re, '" + os.path.join(save_dir, outgroup, "versus_all", clean_strings(outgroup, "_",  True) + "_vs_all_DE.rds") + "')")
            else:
                ro.r("re <- nebula(data_g$count, data_g$id, offset=data_g$offset, pred=data_g$pred, covariance=TRUE)")
                ro.r("saveRDS(re, '" + os.path.join(save_dir, outgroup, "versus_all", clean_strings(outgroup, "_",  True) + "_vs_all_DE.rds") + "')")

            with localconverter(ro.default_converter + pandas2ri.converter):
                results = ro.r("re$summary")
                covariance = ro.r("re$covariance")
                overdispersion = ro.r("re$overdispersion")
                convergence = ro.r("re$convergence")
                random_effect = ro.r("re$random_effect")

            results = results.loc[convergence == 1, :]
            results = results.loc[overdispersion["Subject"] < 1, :]
            gene_index = [j - 1 for j in results["gene_id"].to_list()]
            results.index = adata.var_names[gene_index]

            results.to_csv(os.path.join(save_dir, outgroup, "versus_all", clean_strings(outgroup, "_",  True) + "_versus_all_DE.csv"))
            print(str(datetime.now()) + " -- " + clean_strings(outgroup, "_",  True) + " versus all was written to disk", flush=True)

        else:
            print(str(datetime.now()) + " -- Skipping " + clean_strings(outgroup, "_",  True) + " versus all (already exists)", flush=True)

    else: 
        print(str(datetime.now()) + " -- No outgroup specified, running analysis for all groups and tests", flush=True)

    sub = sub[sub.obs[split_key] == outgroup]
    starting_features = sub.shape[1]
    
    if layer != "X":
        counts_per_cell = sub.layers[layer].sum(axis=0) / sub.shape[0]
    else:
        counts_per_cell = sub.X.sum(axis=0) / sub.shape[0]
    
    sub = sub[:, counts_per_cell > 0.005].copy()
    ending_features = sub.shape[1]
    if starting_features != ending_features:
        print(str(datetime.now()) + " -- Removing " + str(starting_features - ending_features) + " features from " + clean_strings(outgroup, "_",  True) + " for low numbers of counts per cell", flush=True)
        
    if layer != "X":
        global_nonzero_values = (sub.layers[layer] > 0).sum(axis=0)
    else:
        global_nonzero_values = (sub.X > 0).sum(axis=0)
    to_remove = []
    for z in covariates:
        try:
            sub.obs[z] = (sub.obs[z] - sub.obs[z].min()) / sub.obs[z].max()
            print(str(datetime.now()) + " -- Detected " + z + " as an integer or float, applying a min-max normalization", flush=True)
        except:
            pass
            
        if sub.obs[z].unique().shape[0] < 2:
            to_remove.append(z)
            covariate_formula = covariate_formula.replace(z + " + ", "")
            print(str(datetime.now()) + " -- Removing " + z + " from the covariate formula", flush=True)
            continue
        
        if sub.obs[z].unique().shape[0] > 10:
            continue
        
        for y in sub.obs[z].unique():
            if layer != "X":
                spec_nonzero_values = (sub[sub.obs[z] == y, :].layers[layer] > 0).sum(axis=0)
            else:
                spec_nonzero_values = (sub[sub.obs[z] == y, :].X > 0).sum(axis=0)
            if len(sub[:, (spec_nonzero_values == 0) & (global_nonzero_values != 0)].var_names) > 0:
                    print(str(datetime.now()) + " -- Adding 3 pseudocounts to " + z + "=" + str(y) + " for " + str(len(sub[:, (spec_nonzero_values == 0) & (global_nonzero_values != 0)].var_names)) + " genes in " + clean_strings(outgroup, "_",  True) + " because all values are 0 for that covariate", flush=True)
                    for x in sub[:, (spec_nonzero_values == 0) & (global_nonzero_values != 0)].var_names:
                        random_three = np.random.choice(sub[sub.obs[z] == y].obs_names, 3)
                        if layer != "X":
                            sub.layers[layer][[i in random_three for i in sub.obs_names], np.where(sub.var_names == x)[0][0]] = 1
                        else:
                            sub.X[[i in random_three for i in sub.obs_names], np.where(sub.var_names == x)[0][0]] = 1
    
    covariates = np.setdiff1d(covariates, to_remove).tolist()

    if layer != "X":
        tmp = sub.layers[layer].T.tocoo()
    else:
        tmp = sub.X.T.tocoo()
        
    counts = r_Matrix.sparseMatrix(
        i=ro.IntVector(tmp.row + 1),
        j=ro.IntVector(tmp.col + 1),
        x=ro.FloatVector(tmp.data),
        dims=ro.IntVector(tmp.shape)
    )
    with localconverter(ro.default_converter + pandas2ri.converter):
        obs = ro.conversion.py2rpy(sub.obs)
    
    ro.r.assign("obs", obs)
    ro.r.assign("counts", counts)
    
    if tests == [""]:
        return
    
    for i in tests:
        if os.path.exists(os.path.join(save_dir, outgroup, i)) == False:
            os.makedirs(os.path.join(save_dir, outgroup, i))
            
        if os.path.exists(os.path.join(save_dir, outgroup, i, clean_strings(target, "_",  True) + "_" + clean_strings(outgroup, "_",  True) + "_across_" + i + "_DE.csv")) == False:
            print(str(datetime.now()) + " -- Starting " + clean_strings(outgroup, "_",  True) + " across " + i)
            ro.r("covariates <- c('" + i + "', " + str(covariates).replace("[", "").replace("]", "") + ")")
            ro.r("df <- model.matrix(~" + covariate_formula + i + ", data=obs[,covariates])")
            ro.r("data_g <- group_cell(count=counts, id=obs$" + random_effect + ", offset=obs$" + offset_variable + ", pred=df)")

            with localconverter(ro.default_converter + pandas2ri.converter):
                df = ro.r("as.data.frame(df)")
            
            ordered = ro.r("data_g")
            if (group_cell is False) and (isinstance(ordered, rpy2.rinterface_lib.sexp.NULLType) is True):
                ro.r("re <- nebula(counts, obs$" + random_effect + ", pred=df, offset=obs$" + offset_variable + ", covariance=TRUE)")
                ro.r("saveRDS(re, '" + os.path.join(save_dir, outgroup, i, clean_strings(target, "_",  True) + "_" + clean_strings(outgroup, "_",  True) + "_across_" + i + "_DE.rds") + "')")
            else:
                ro.r("re <- nebula(data_g$count, data_g$id, pred=data_g$pred, offset=data_g$offset, covariance=TRUE)")
                ro.r("saveRDS(re, '" + os.path.join(save_dir, outgroup, i, clean_strings(target, "_",  True) + "_" + clean_strings(outgroup, "_",  True) + "_across_" + i + "_DE.rds") + "')")

            with localconverter(ro.default_converter + pandas2ri.converter):
                results = ro.r("re$summary")
                overdispersion = ro.r("re$overdispersion")
                convergence = ro.r("re$convergence")
                random_effect = ro.r("re$random_effect")

            results = results.loc[convergence == 1, :]
            gene_index = [j - 1 for j in results["gene_id"].to_list()]
            results.index = sub.var_names[gene_index]
            results.to_csv(os.path.join(save_dir, outgroup, i, clean_strings(target, "_",  True) + "_" + clean_strings(outgroup, "_",  True) + "_across_" + i + "_DE.csv"))
            print(str(datetime.now()) + " -- " + clean_strings(outgroup, "_",  True) + " along " + i + " was written to disk", flush=True)
        
        else:
            print(str(datetime.now()) + " -- Skipping " + clean_strings(outgroup, "_",  True) + " along " + i + " (already exists)", flush=True)
    return


def build_effect_size_anndata(
    results_dir,
    glob_pattern,
    file_pattern,
    test,
    adata,
    subclass,
    celltype,
    blacklisted_genes=[],
    filter_genes=True,
    normalize=True,
    effect_size_cutoff=10,
    ):

    cell_types = [d for d in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir, d))]
    results_files = sorted([glob.glob(os.path.join(results_dir, cell_type, test, glob_pattern))[0] for cell_type in cell_types])
    subclasses = cell_types
    pattern = r'^(' + '|'.join(map(re.escape, sorted(cell_types, key=len, reverse=True))) + r')_'
    supertypes = [re.sub(pattern, '', os.path.basename(i).replace(file_pattern, '')) 
                for i in results_files]
    # subclasses = [re.sub("([^_]+)_(.*)(_[0-9]{1,2})?(-SEAAD)?$", "\\1", os.path.basename(i).replace(file_pattern, "").replace("Lamp5_Lhx6", "Lamp5 Lhx6")) for i in results_files]
    # supertypes = [re.sub("([^_]+)_(.*)$", "\\2", os.path.basename(i).replace(file_pattern, "").replace("Lamp5_Lhx6", "Lamp5 Lhx6")) for i in results_files]
    # supertypes = [i.replace("Lamp5 Lhx6", "Lamp5_Lhx6").replace("L2 3 IT", "L2/3 IT").replace("L5 6 NP", "L5/6 NP") for i in supertypes]
        
    effect_sizes = pd.DataFrame(np.zeros((len(adata.obs[celltype].cat.categories), adata.shape[1])), 
                                columns=adata.var_names, 
                                index=clean_strings(adata.obs[celltype].cat.categories, preserve_case=True))
    pvalues = pd.DataFrame(np.ones((len(adata.obs[celltype].cat.categories), adata.shape[1])), 
                            columns=adata.var_names, 
                            index=clean_strings(adata.obs[celltype].cat.categories, preserve_case=True))
    std_errors = pd.DataFrame(np.ones((len(adata.obs[celltype].cat.categories), adata.shape[1])), 
                                columns=adata.var_names, 
                                index=clean_strings(adata.obs[celltype].cat.categories, preserve_case=True))
    
    for i,j in enumerate(results_files):
        print(os.path.basename(j))
        results = pd.read_csv(j, index_col=0)
        effect_sizes.loc[supertypes[i], results.index] = results.loc[:, "logFC_" + test]
        pvalues.loc[supertypes[i], results.index] = results.loc[:, "p_" + test]
        std_errors.loc[supertypes[i], results.index] = results.loc[:, "se_" + test]
    
    effect_sizes = ad.AnnData(X=effect_sizes)

    if subclass!=celltype:
        subclasses = adata.obs.loc[:, [subclass, celltype]].drop_duplicates().sort_values(by=celltype).loc[:, subclass].to_list()
    else:
        subclasses = adata.obs[subclass].drop_duplicates().sort_values().to_list()
    effect_sizes.obs[subclass] = subclasses
    effect_sizes.obs[subclass] = effect_sizes.obs[subclass].astype("category")
    
    pvalues = ad.AnnData(X=pvalues)
    pvalues.obs[subclass] = subclasses
    std_errors = ad.AnnData(X=std_errors)
    std_errors.obs[subclass] = subclasses

    effect_sizes = effect_sizes.T
    pvalues = pvalues.T
    std_errors = std_errors.T

    effect_sizes.var["Class"] = "Non-neuronal and non-neural"
    effect_sizes.var.loc[[i in ["Sst", "Sst Chodl", "Pvalb", "Chandelier", "Vip", "Sncg", "Pax6", "Lamp5", "Lamp5 Lhx6"] for i in effect_sizes.var["Subclass"]], "Class"] = "Neuronal: GABAergic"
    effect_sizes.var.loc[[i in ["L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", "L5/6 NP", "L6 CT", "L6b"] for i in effect_sizes.var["Subclass"]], "Class"] = "Neuronal: Glutamatergic"

    pvalues.var["Class"] = "Non-neuronal and non-neural"
    pvalues.var.loc[[i in ["Sst", "Sst Chodl", "Pvalb", "Chandelier", "Vip", "Sncg", "Pax6", "Lamp5", "Lamp5 Lhx6"] for i in pvalues.var["Subclass"]], "Class"] = "Neuronal: GABAergic"
    pvalues.var.loc[[i in ["L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", "L5/6 NP", "L6 CT", "L6b"] for i in pvalues.var["Subclass"]], "Class"] = "Neuronal: Glutamatergic"

    std_errors.var["Class"] = "Non-neuronal and non-neural"
    std_errors.var.loc[[i in ["Sst", "Sst Chodl", "Pvalb", "Chandelier", "Vip", "Sncg", "Pax6", "Lamp5", "Lamp5 Lhx6"] for i in std_errors.var["Subclass"]], "Class"] = "Neuronal: GABAergic"
    std_errors.var.loc[[i in ["L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", "L5/6 NP", "L6 CT", "L6b"] for i in std_errors.var["Subclass"]], "Class"] = "Neuronal: Glutamatergic"


    effect_sizes.var["Subclass"] = effect_sizes.var["Subclass"].cat.reorder_categories(
        adata.obs["Subclass"].cat.categories,
    )
    
    if "Supertype" in adata.obs.columns:
        effect_sizes = effect_sizes[:, adata.obs["Supertype"].cat.categories].copy()
        effect_sizes.var["Supertype"] = effect_sizes.var.index.copy()

    if normalize == True:
        pvalues.layers["pvalues"] = pvalues.X.copy()
        pvalues.X = -np.log10(pvalues.X) * np.sign(effect_sizes.X)
        
        effect_sizes.X = np.nan_to_num(effect_sizes.X)
        effect_sizes.X[effect_sizes.X > effect_size_cutoff] = 0
        effect_sizes.X[effect_sizes.X < -1 * effect_size_cutoff] = 0
        effect_sizes.layers["effect_sizes"] = effect_sizes.X.copy()
        effect_sizes.X = effect_sizes.X / std_errors.X
        effect_sizes.X = np.nan_to_num(effect_sizes.X)

    if len(blacklisted_genes) > 0:
        effect_sizes = effect_sizes[[i not in blacklisted_genes for i in effect_sizes.obs_names], :].copy()
        pvalues = pvalues[[i not in blacklisted_genes for i in pvalues.obs_names], :].copy()
        std_errors = std_errors[[i not in blacklisted_genes for i in std_errors.obs_names], :].copy()
        
    return effect_sizes, pvalues, std_errors


def get_standardized_effects(results_path, factor, filter_genes=False, thr=0.05):
    results = pd.read_csv(results_path)
   
    cols = {
        'effect': f"logFC_{factor}",
        'pval': f"p_{factor}",
        'se': f"se_{factor}",
        'z_scores': f'z_scores_{factor}'
    }
    
    if filter_genes:
        results = results[results[cols['pval']] < thr]
    
    results[cols['sign_log_pval']] = -np.log10(results[cols['pval']]) * np.sign(results[cols['effect']])
    results[cols['z_scores']] = np.nan_to_num(results[cols['effect']]) / np.nan_to_num(results[cols['se']])
    results['z_scores*pval'] = results[cols['z_scores']] * np.abs(results[cols['sign_log_pval']])
    
    results.set_index('Unnamed: 0', inplace=True)
    results.index.name = 'gene_name'
    
    results = results[[cols['effect'], cols['pval'], cols['se'], cols['z_scores'], 'sign_log_pval', 'z_scores*pval']]
    
    results.rename(columns={
        cols['effect']: 'logFC',
        cols['pval']: 'pval',
        cols['se']: 'se',
        cols['z_scores']: 'z_scores'},
        inplace=True)

    return results

# threads = args.threads
# threads = int(threads)
# target = args.target
# region = args.region
# split_key = args.split_key
# outgroup_de = args.outgroup_de
# covariates = args.covariates
# covariates = [a for a in covariates.split(",")]
# covariate_formula = args.covariate_formula
# tests = args.tests
# tests = [a for a in tests.split(",")]
# layer = args.layer
# adata = args.adata
# adata.obs[split_key] = adata.obs[split_key].astype("category")
# offset_variable = args.offset_variable
# group_cell = args.group_cell
# save_dir = args.save_dir

# Parallel(n_jobs=threads)(
#     delayed(run_de_tests)(
#         adata, region, target, split_key, outgroup_de, b, covariates, covariate_formula, tests,
#          layer, offset_variable, group_cell, save_dir) for b in adata.obs[split_key].cat.categories
# )