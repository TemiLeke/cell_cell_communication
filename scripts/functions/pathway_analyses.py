
import scipy
import warnings
import anndata2ri
import pandas as pd
import scanpy as sc
import numpy as np
import seaborn as sb
import decoupler as dc
from scipy import sparse
from anndata import AnnData
from tabnanny import verbose
import matplotlib.pyplot as plt
from gsva_prep import prep_gsva
from typing import Optional, Union
from matplotlib.pyplot import rcParams
from statsmodels.stats.multitest import multipletests
from sklearn.model_selection import train_test_split
from rpy2.robjects.conversion import localconverter




def rescale_matrix(S, log_scale=False):
    """
    Sums cell-level counts by factors in label vector

    Parameters
    ----------
    S : np.ndarray, scipy.sparse.csr_matrix or pandas.DataFrame
        Matrix with read counts (gene x cell)
    log_scale : bool, optional (default: False)
        Whether to log-transform the rescaled matrix

    Returns
    -------
    B : np.ndarray or scipy.sparse.csr_matrix
        Scaled and log-transformed matrix
    """
    if isinstance(S, pd.DataFrame):
        S = S.values
    elif isinstance(S, np.ndarray):
        pass
    elif isinstance(S, scipy.sparse.csr_matrix):
        S = S.toarray()
    else:
        raise ValueError('Input S must be a pandas.DataFrame, numpy.ndarray or scipy.sparse.csr_matrix')
        
    cs = np.sum(S, axis=0)
    cs[cs == 0] = 1
    B = np.median(cs) * (S / cs)
    if log_scale:
        B = np.log1p(B)
    return B

def normalize_default(adata, log_scale=True):
    """
    Normalizes gene expression matrix by total count and scales by median

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    log_scale : bool, optional (default: True)
        Whether to log-transform the rescaled matrix.

    Returns
    -------
    adata : AnnData
        Annotated data matrix with normalized and scaled expression values.
    """
    if 'counts' in adata.layers.keys():
        print('normalizaing data using count data in .layers["counts] ')
        S = adata.layers['counts']
    else:
        print('normaling data using count data in .X')
        S = adata.X 
    B = rescale_matrix(S, log_scale=log_scale)
    adata.X = B
    return adata


def normalize_matrix(
    X: Union[np.ndarray, sparse.spmatrix],
    top_features_frac: float = 1.0,
    scale_factor: Union[str, float, int, np.ndarray, None] = "median",
    transformation: Union[str, None] = "log",
    anchor_features: Union[np.ndarray, None] = None,
) -> Union[np.ndarray, sparse.spmatrix]:

    X = X.astype(dtype=np.float64)

    # Which features (i.e. genes) should we use to compute library sizes?
    if anchor_features is not None:
        lib_sizes = np.array(np.mean(X[:, anchor_features], axis=1))
    else:
        if top_features_frac < 1.0:
            universality = np.array(np.mean(X > 0, axis=0))
            selected_features = np.flatnonzero(universality > (1 - top_features_frac))
            lib_sizes = np.array(np.mean(X[:, selected_features], axis=1))
        else:
            lib_sizes = np.array(np.mean(X, axis=1))

    # Note: mean as opposed to sum

    # Normalize library sizes
    if isinstance(X, sparse.spmatrix):
        X_scaled = X.multiply(1 / lib_sizes)
    else:
        try:
            X_scaled = X / lib_sizes
        except ValueError:
            lib_sizes = np.reshape(lib_sizes, (-1, 1))
            X_scaled = X / lib_sizes

    # scale normalized columns
    if scale_factor == "median":
        kappa = np.median(np.array(np.sum(X, axis=1) / np.sum(X_scaled, axis=1)))
        X_scaled_norm = X_scaled * kappa
    elif isinstance(scale_factor, (int, float)):
        X_scaled_norm = X_scaled * scale_factor
    elif isinstance(scale_factor, np.ndarray):
        if sparse.issparse(X_scaled):
            X_scaled_norm = X_scaled.multiply(scale_factor)
        else:
            X_scaled_norm = X_scaled / scale_factor

    # For compatibility with C
    if sparse.issparse(X_scaled_norm):
        X_scaled_norm = sparse.csc_matrix(X_scaled_norm)

    # Post-transformation
    if transformation == "log":
        X_scaled_norm_trans = np.log1p(X_scaled_norm)
    elif transformation == "tukey":
        if sparse.issparse(X_scaled_norm):
            nnz_idx = X_scaled_norm.nonzero()
            ii = nnz_idx[0]
            jj = nnz_idx[1]
            vv = X_scaled_norm[ii, jj]
            vv_transformed = np.sqrt(vv) + np.sqrt(1 + vv)
            X_scaled_norm[ii, jj] = vv_transformed
        else:
            X_scaled_norm[X_scaled_norm < 0] = 0
            vv = X_scaled_norm[X_scaled_norm != 0]
            vv_transformed = np.sqrt(vv) + np.sqrt(1 + vv)
            X_scaled_norm[X_scaled_norm != 0] = vv_transformed
        
    # elif transformation == "lsi":
    #     if sparse.issparse(X_scaled_norm):
    #         X_scaled_norm_trans = _an.LSI(X_scaled_norm)
    #     else:
    #         X_scaled_norm_sp = sparse.csc_matrix(X_scaled_norm)
    #         X_scaled_norm_trans = _an.LSI(X_scaled_norm_sp).toarray()
    else:
        X_scaled_norm_trans = X_scaled_norm

    return X_scaled_norm_trans


def normalize_actionet(
    adata: AnnData,
    layer_key: Optional[str] = None,
    layer_key_out: Optional[str] = None,
    top_features_frac: float = 1.0,
    scale_factor: Union[str, float, int, np.ndarray, None] = "median",
    transformation: Union[str, None] = "log",
    anchor_features: Union[np.ndarray, None] = None,
    copy: Optional[bool] = False,
) -> Optional[AnnData]:
    adata = adata.copy() if copy else adata

    if "metadta" in adata.uns.keys():
        if "norm_method" in adata.uns["metadata"].keys():  # Already normalized? leave it alone!
            # return adata if copy else None
            warnings.warn("AnnData object is prenormalized. Please make sure to use the right assay.")

    if layer_key is None and "input_assay" in adata.uns["metadata"].keys():
        layer_key = adata.uns["metadata"]["input_assay"]

    if layer_key is not None:
        if layer_key not in adata.layers.keys():
            raise ValueError("Did not find adata.layers['" + layer_key + "']. ")
        S = adata.layers[layer_key]
    else:
        S = adata.X

    if sparse.issparse(S):
        UE = set(S.data)
    else:
        UE = set(S.flatten())

    nonint_count = len(UE.difference(set(np.arange(0, max(UE) + 1))))
    if 0 < nonint_count:
        warnings.warn("Input [count] assay has non-integer values, which looks like a normalized matrix. Please make sure to use the right assay.")

    S = normalize_matrix(
        S,
        anchor_features=anchor_features,
        top_features_frac=top_features_frac,
        scale_factor=scale_factor,
        transformation=transformation,
    )

    adata.uns["metadata"] = {}
    adata.uns["metadata"]["norm_method"] = "default_top%.2f_%s" % (
        top_features_frac,
        transformation,
    )

    if layer_key_out is not None:
        adata.uns["metadata"]["default_assay"] = layer_key_out
        adata.layers[layer_key_out] = S
    else:
        adata.uns["metadata"]["default_assay"] = None
        adata.X = S

    return adata if copy else None

def read_pathways(filename):
    with open(filename, 'r') as temp_f:
        col_count = [ len(l.split("\t")) for l in temp_f.readlines() ]
    column_names = [i for i in range(0, max(col_count))]
    ### Read csv
    return pd.read_csv(filename, header=None, delimiter="\t", names=column_names)



def filter_expressed_genes_by_celltype(adata: AnnData, 
                                      threshold: float=0.05,
                                      filter_genes_from: str='singlecell', 
                                      subject_id: str='Subject', 
                                      groups_col = 'cell_type'):
    """

        Function to filter expressed genes by cell type based on a threshold

    Parameters:
    -----------
    adata : AnnData object
        Annotated Data matrix with rows representing genes and columns representing cells.
    threshold : float, optional (default=0.05)
        The threshold to use for filtering expressed genes based on the minimum number of cells they are detected in.
    filter_genes_from: str, optional (default=`singlecell`)
        Whether to filter genes that meet threshold in pseudobulk data or singlecell data.
    subject_id (str): a string indicating the column containing individual identifiers.


    Returns:
    --------
    expressed_genes_per_celltype : pandas DataFrame
        A dataframe where the rows are the gene names and columns are the cell types, 
        containing only the genes that are expressed in at least the specified percentage of cells for each cell type.


    """

    # Initialize empty dictionaries to store the expressed genes and gene sets per cell type
    expressed_genes_per_celltype = {}
    gene_set_per_celltype = {}

    if filter_genes_from=='pseudobulk':
        # Get pseudo-bulk profile
        adata = dc.get_pseudobulk(adata,
                                sample_col=subject_id,
                                groups_col=groups_col,
                                layer='counts',
                                mode='sum',
                                min_cells=0,
                                min_counts=0
                                )
        # Loop through each unique cell type in the input AnnData object

        for cell_type in adata.obs[groups_col].unique():

            expressed_genes_per_celltype[cell_type] = dc.filter_by_prop(adata[adata.obs[groups_col]==cell_type],
                                                                      min_prop=threshold)
    
    elif filter_genes_from=='singlecell':
        # Loop through each unique cell type in the input AnnData object

        for cell_type in adata.obs[groups_col].unique():

            # Calculate the number of cells based on the specified threshold
            percent = threshold
            num_cells = round(percent*len(adata[adata.obs[groups_col]==cell_type]))

            # Filter genes based on minimum number of cells and store the resulting gene names
            expressed_genes_per_celltype[cell_type], _ = sc.pp.filter_genes(adata[adata.obs[groups_col]==cell_type].layers['counts'],
                                                                            min_cells=num_cells, inplace=False)
            expressed_genes_per_celltype[cell_type] = list(adata.var_names[expressed_genes_per_celltype[cell_type]])

    # Convert the dictionary of expressed genes per cell type to a Pandas DataFrame
    expressed_genes_per_celltype = pd.DataFrame.from_dict(expressed_genes_per_celltype, orient='index').transpose()

    return expressed_genes_per_celltype


def filter_lowly_exp_genes(expressed: pd.DataFrame,
                           all_paths: pd.DataFrame,
                           threshold: float = 0.33):

    """
    Filters lowly expressed gene sets based on a threshold and pathway membership.

    Parameters:
    -----------
    expressed: pandas.DataFrame
        A DataFrame of expressed genes with cell types as columns and gene IDs as rows.
    all_paths: pandas.DataFrame
        A DataFrame of gene sets with pathways as columns and gene IDs as rows.
    threshold: float, optional (default=0.33)
        A proportion threshold used to filter gene sets based on their expression in each cell type.

    Returns:
    --------
    gene_set_per_celltype: dict of pandas.DataFrame
        A dictionary of gene sets per cell type, with cell type names as keys and gene set dataframes as values.
        Each gene set dataframe has three columns: 'description', 'member', and 'name'.
    """

    # Initialize empty dictionaries to store the gene sets and gene sets per cell type
    gene_set = {}
    gene_set_per_celltype = {}

    # Loop through each cell type in the input Pandas DataFrame of expressed genes
    for cell_type in expressed.columns:
        # Determine which pathways have a proportion of genes above the specified threshold
        index = [sum(all_paths[x].isin(expressed[cell_type]))/\
                 (sum(~all_paths[x].isin([None]))) > threshold for x in all_paths.columns]
        # Filter pathways based on threshold and store the resulting gene sets
        p = all_paths.loc[:, index]
        x = {y: pd.Series(list(set(expressed[cell_type]).intersection(set(p[y])))) for y in p.columns}
        x = {k: v for k, v in x.items() if not v.empty}
        gene_set[cell_type] = x

        # Convert the gene sets to Pandas DataFrames and store them in a dictionary by cell type
        gene_set_per_celltype[cell_type] = pd.DataFrame(columns=['description', 'member', 'name'])    
        for pathway, gene_list in gene_set[cell_type].items():

            df = pd.DataFrame(columns=['description', 'member', 'name'])  
            df['member'] = gene_list
            df['name'] = pathway
            df['description'] = pathway.split(" ")[-1]                                                          
            gene_set_per_celltype[cell_type] = pd.concat([gene_set_per_celltype[cell_type], df], join='outer', ignore_index=True)

        # Sort the resulting gene sets by description and member
        gene_set_per_celltype[cell_type].sort_index(axis=1, inplace=True)
        gene_set_per_celltype[cell_type].sort_index(axis=0, inplace=True)


    return gene_set_per_celltype


def get_ind_level_ave(adata: AnnData, subject_id: str = 'Subject', method: str = "agg_x_num",
                     expressed_genes_per_celltype: dict = {}, filter_genes_at_threshold: bool = True, groups_col='cell_type'):
    """
    Get averaged expression data for each cell type and individual in an AnnData object.
    

    Args:

        adata (AnnData): An AnnData object with read counts (gene x cell).
        subject_id (str): a string indicating the column containing individual identifiers.
        method (str): a string indicating the method to be used. The default is "agg_x_num".
        filter_genes_at_threshold (bool): A boolean indicating whether to filter genes based on threshold. The default is True.
        expressed_genes_per_celltype (float): A dictionary of the genes to be filtered for each celltype.
        
    Returns:

        Dictionary: A dictionary of data frames with averaged expression data for each cell type and individual.

    """

    if method == "agg_x_norm":

        avs_logcounts_cellxind = {}
        # loop over each unique cell type in the annotation metadata
        for cell_type in adata.obs[groups_col].unique():

            # filter genes based on threshold
            if filter_genes_at_threshold:
                adata_temp = adata[adata.obs[groups_col]==cell_type].copy()
                # sc.pp.filter_genes(adata_temp, min_cells=gene_celltype_threshold*adata_temp.n_obs)
                adata_temp = adata_temp[:, adata_temp.var_names.isin(expressed_genes_per_celltype[cell_type].tolist())]
            else:
                adata_temp = adata.copy()

            # Get pseudo-bulk profile
            pdata = dc.get_pseudobulk(adata_temp, sample_col=subject_id, groups_col=groups_col, layer='counts', mode='sum',
                                    min_cells=0, min_counts=0)
            
            # genes = dc.filter_by_prop(pdata, min_prop=0.05, min_smpls=1)
            # pdata = pdata[:, genes].copy()

            # Normalize and log transform

            # sc.pp.normalize_total(pdata, 1e06)
            # sc.pp.log1p(pdata)

            pdata.layers['counts'] = pdata.X
            pdata = normalize_actionet(pdata, layer_key = 'counts', layer_key_out = None, 
                                top_features_frac = 1.0, scale_factor = "median", 
                                transformation = "log", anchor_features = None, copy = True)

            # Store the log-normalized, averaged expression data for each individual and cell type
            avs_logcounts_cellxind[cell_type] = pd.DataFrame(pdata.X.T, columns=pdata.obs[subject_id], index=pdata.var_names)
            
            del adata_temp, pdata
            
    elif method == 'norm_x_agg':

        def sum_counts(counts, label, cell_labels, gene_labels):

            """
            Sums cell-level counts by factors in label vector.
            
            Args:
                counts (AnnData): An AnnData object with read counts (gene x cell).
                label (pd.DataFrame): Variable of interest by which to sum counts.
                cell_labels (pd.Index): Vector of cell labels.
                gene_labels (pd.Index): Vector of gene labels.
                
            Returns:
                Dictionary: A dictionary with the following keys:
                    - 'summed_counts': A data frame with summed counts.
                    - 'ncells': A data frame with the number of cells used per summation.
            """
            # Create a data frame with the label vector and add a column of 1s for counting.
            label_df = pd.DataFrame(label)
            label_df.columns = ['ID']
            label_df['index'] = 1

            # Add a column for cell type and pivot the data frame to create a matrix of counts.
            label_df['celltype'] = cell_labels
            label_df = label_df.pivot_table(index='celltype', columns='ID', values='index', aggfunc=np.sum, fill_value=0)
            label_df = label_df.astype(float)

            # Multiply the counts matrix by the gene expression matrix to get summed counts.
            summed_counts = pd.DataFrame(counts.X.T @ label_df.values, index = gene_labels, columns= label_df.columns)

            # Sum the number of cells used for each summation.
            ncells = label_df.sum()

            # Return the summed counts and number of cells as a dictionary.
            return {'summed_counts': summed_counts, 'ncells': ncells}

        
        # Get metadata from the AnnData object.
        meta = adata.obs # Get metadata


        # Create a data frame of labels by combining cell type and individual metadata fields.
        # Sum counts by individual
        labels = pd.DataFrame(meta[groups_col].astype(str) + '_' + meta[subject_id].astype(str), columns=['individual'])

        # Sum counts by individual and store the results in a dictionary.
        summed_logcounts_cellxind = sum_counts(adata, labels, adata.obs_names, adata.var_names)

        # Calculate averages for each cell type and individual and store the results in a dictionary.
        # Get averages corresponding to both count matrices
        avs_logcounts = np.array(summed_logcounts_cellxind['summed_counts'].values) / np.array(summed_logcounts_cellxind['ncells'].values)
        # avs_logcounts = np.array(summed_logcounts_cellxind['summed_counts'].values)
        avs_logcounts = pd.DataFrame(avs_logcounts, index = summed_logcounts_cellxind['summed_counts'].index,
                                                columns=summed_logcounts_cellxind['summed_counts'].columns)
        

        # Split the averages by cell type and individual and store the results in a dictionary.
        # Split column names into two parts: cell type and individual
        x = [col.split('_') for col in avs_logcounts.columns]
        celltype = [col[0] for col in x]
        individual = [col[1] for col in x]

        # Get unique cell types in the dataset
        celltype_unique = np.unique(celltype)

        # Create an empty dictionary to store the average counts for each cell type and individual    
        avs_by_ind_out = {}

        # Loop over the unique cell types and subset the average counts for each cell type and individual
        for i in celltype_unique:
            index = np.array(celltype)==i
            df = avs_logcounts.loc[:, index]
            df.columns = np.array(individual)[index]
            avs_by_ind_out[i] = df
            
            if filter_genes_at_threshold:
                # num_cells = round(gene_celltype_threshold*len(adata[adata.obs[groups_col]==cell_type]))
                # # Filter genes based on minimum number of cells and store the resulting gene names
                # gene_mask, _ = sc.pp.filter_genes(adata[adata.obs[groups_col]==cell_type].layers['counts'],
                #                                                                 min_cells=num_cells, 
                #                                                                 inplace=False)
                # genes = list(adata.var_names[gene_mask])
                avs_by_ind_out[i] = avs_by_ind_out[i].loc[expressed_genes_per_celltype[i], :]
            else:
                adata = adata.copy()
        # Store the dictionary of average counts for each cell type and individual    
        avs_logcounts_cellxind = avs_by_ind_out

        # Return the dictionary of average counts for each cell type and individual

    return  avs_logcounts_cellxind


def plot_and_select_top_deps(all_pathways: pd.DataFrame(), 
                            list_of_paths_to_annotate: list = [],
                            linewidths: float = 0.07,
                            cell_types: list = ["Excitatory", "Inhibitory", "Astrocyte",
                                                     "Microglia", "Oligodendrocyte", "OPC", "Endothelial"],
                            save_path='cell_type_specific',
                            save_prefix: str = 'mathys_pfc', 
                            filter: bool=False,
                            cell_type_specific: bool = True,
                            test_name: str = ''): 

    if cell_type_specific:
        # Plot certain cell_type specific pathways 
        collated_df = pd.DataFrame(all_pathways.groupby(all_pathways.index).agg({'score_adj': list, 'celltype': list, 
                                    'logFC': list, 'P.Value': list, 'shortened': list, 'highlight': list}))
        # filter pathways only expressed in one cell type
        mask = collated_df["celltype"].apply(len) == 1
        df = collated_df[mask]

        # create pathway by cell type pivot table
        scores_table = pd.pivot_table(all_pathways, values='score_adj', index='pathway', columns='celltype')
        scores_table = scores_table.loc[df.index]
        scores_table['shortened'] = df.shortened.apply(lambda x: x[0])
        scores_table['highlight'] = df.highlight.apply(lambda x: x[0])
        scores_table.sort_values(by=[cell_type for cell_type in all_pathways.celltype.unique()], inplace=True)

        # drop pathways with same shortened names ??
        scores_table = scores_table.drop_duplicates(subset='shortened', keep='first')

        ###### Plot Cell type specific data

        if filter:
            xticks = cell_types

            # select only pathways that should be visualized
            shortened_names = scores_table[scores_table.shortened.isin(list_of_paths_to_annotate)]['shortened']
            scores_table = scores_table[scores_table.shortened.isin(list_of_paths_to_annotate)]

            n_rows = len(scores_table)

            fig, ax1 = plt.subplots(1, 1, figsize=(0.5, n_rows*0.095), sharex=False, layout='constrained')
            fig.tight_layout()
    
            # order table by cell type name
            # scores_table = scores_table.reindex(columns=['Excitatory', 'Inhibitory', 'Astrocyte', 'Oligodendrocyte',
            #                                                     'OPC', 'Microglia'])
            scores_table = scores_table[xticks].apply(pd.to_numeric, errors='coerce')

            g1 = sb.heatmap(scores_table, cmap='bwr', center=0, vmin=-2.5, vmax=2.5, robust=False, annot=None, fmt='.1g', 
                            linewidths=linewidths, linecolor='black', annot_kws=None, cbar_kws={'shrink': 0.2},
                            cbar_ax=None, square=False,ax=ax1, xticklabels=xticks, yticklabels=shortened_names, mask=None,) 


            cax = g1.figure.axes[-1]

            g1.set_title(f'Select Cell-type-specific Pathways in {test_name.split("_")[0]}- vs {test_name.split("_")[-1]}-pathology',
                          fontsize=3)           
            g1.set_ylabel('')
            g1.set_xlabel('')

            ax1.tick_params(axis='both', which='major', labelsize=4, length=1.5, width=0.5)
            cax.tick_params(labelsize=4, length=1.5, width=0.5, which="major")

            plt.tight_layout()
            plt.savefig(save_path, bbox_inches='tight')
            plt.show(block=False)

        else:
            xticks = cell_types


            scores_table = scores_table[scores_table.shortened!='None']     
            yticklabels = scores_table['shortened']
            # order table by cell type name
            
            scores_table = scores_table[xticks].apply(pd.to_numeric, errors='coerce')

            n_rows = len(scores_table)

            fig, ax1 = plt.subplots(1, 1, figsize=(0.5, n_rows*0.095), sharex=False, layout='constrained')
            fig.tight_layout()

            g1 = sb.heatmap(scores_table, cmap='bwr', center=0, vmin=-2.5, vmax=2.5, robust=False, annot=None, fmt='.1g', 
                            linewidths=linewidths, linecolor='black', annot_kws=None, cbar_kws={'shrink': 0.1},
                            cbar_ax=None, square=False, ax=ax1, xticklabels=xticks, yticklabels=yticklabels, mask=None,) 


            cax = g1.figure.axes[-1]

            g1.set_title(f'All Cell-type-specific Pathways in {test_name.split("_")[0]}- vs {test_name.split("_")[-1]}-pathology', 
                         fontsize=3)           
            g1.set_ylabel('')
            g1.set_xlabel('')

            ax1.tick_params(axis='both', which='major', labelsize=2, length=1.5, width=0.25)
            cax.tick_params(labelsize=4, length=1.5, width=0.25, which="major")


            plt.tight_layout()
            plt.savefig(save_path, bbox_inches='tight')
            plt.show(block=False)


    else:
        # Plot certain cell_type specific pathways 
        collated_df = pd.DataFrame(all_pathways.groupby(all_pathways.index).agg({'score_adj': list, 'celltype': list, 
                                    'logFC': list, 'P.Value': list, 'shortened': list, 'highlight': list}))
        # filte pathways only expressed in one cell type
        mask = collated_df["celltype"].apply(len) > 1
        df = collated_df[mask]

        # create pathway by cell type pivot table
        scores_table = pd.pivot_table(all_pathways, values='score_adj', index='pathway', columns='celltype')
        scores_table = scores_table.loc[df.index]
        scores_table['shortened'] = df.shortened.apply(lambda x: x[0])
        scores_table['highlight'] = df.highlight.apply(lambda x: x[0])
        scores_table.sort_values(by=[cell_type for cell_type in all_pathways.celltype.unique()], inplace=True)

        # drop pathways with same shortened names ??
        scores_table = scores_table.drop_duplicates(subset='shortened', keep='first')

        ###### Plot Cell type specific data

        if filter:
            xticks = cell_types

            # select only pathways that should be visualized
            shortened_names = scores_table[scores_table.shortened.isin(list_of_paths_to_annotate)]['shortened']
            scores_table = scores_table[scores_table.shortened.isin(list_of_paths_to_annotate)]
        
            # order table by cell type name
            scores_table = scores_table[xticks].apply(pd.to_numeric, errors='coerce')

            n_rows = len(scores_table)

            fig, ax1 = plt.subplots(1, 1, figsize=(0.5, n_rows*0.095), sharex=False, layout='constrained')
            fig.tight_layout()

            g1 = sb.heatmap(scores_table, cmap='bwr', center=0, vmin=-2.5, vmax=2.5, robust=False, annot=None, fmt='.1g', 
                            linewidths=0.15, linecolor='black', annot_kws=None, cbar_kws={'shrink': 0.2},
                            cbar_ax=None, square=False,ax=ax1, xticklabels=xticks, yticklabels=shortened_names, mask=None,) 

            cax = g1.figure.axes[-1]

            g1.set_title(f'Select Shared Pathways in {test_name.split("_")[0]}- vs {test_name.split("_")[-1]}-pathology', fontsize=3)           
            g1.set_ylabel('')
            g1.set_xlabel('')

            ax1.tick_params(axis='both', which='major', labelsize=4, length=1.5, width=0.5)
            cax.tick_params(labelsize=4, length=1.5, width=0.5, which="major")

            plt.tight_layout()
            plt.savefig(save_path, bbox_inches='tight')
            plt.show(block=False)

        else:
            xticks = cell_types

            scores_table = scores_table[scores_table.shortened!='None']     
            yticklabels = scores_table['shortened']
            # order table by cell type name
            
            scores_table = scores_table[xticks].apply(pd.to_numeric, errors='coerce')

            n_rows = len(scores_table)

            fig, ax1 = plt.subplots(1, 1, figsize=(0.5, n_rows*0.095), sharex=False, layout='constrained')
            fig.tight_layout()

            g1 = sb.heatmap(scores_table, cmap='bwr', center=0, vmin=-2.5, vmax=2.5, robust=False, annot=None, fmt='.1g', 
                            linewidths=linewidths, linecolor='black', annot_kws=None, cbar_kws={'shrink': 0.1},
                            cbar_ax=None, square=False, ax=ax1, xticklabels=xticks, yticklabels=yticklabels, mask=None,) 

            cax = g1.figure.axes[-1]

            g1.set_title(f'All Broad Pathways in {test_name.split("_")[0]}- vs {test_name.split("_")[-1]}-pathology', fontsize=3)           
            g1.set_ylabel('')
            g1.set_xlabel('')

            ax1.tick_params(axis='both', which='major', labelsize=2, length=1.5, width=0.25)
            cax.tick_params(labelsize=4, length=1.5, width=0.25, which="major")

            plt.tight_layout()
            plt.savefig(save_path, bbox_inches='tight')
            plt.show(block=False)

    return 


def multi_study_pathway_overlap(pathway_scores: dict = {},
                                p_thresh: float = 0.05,
                                filtered_pathways: list = [],
                                cell_types: list = ["Excitatory", "Inhibitory", "Astrocyte",
                                                     "Microglia", "Oligodendrocyte", "OPC", 
                                                     'Endothelial',],
                                test_name: str = 'ad_vs_no',
                                linewidths: float = 0.015,
                                top_n: int = 10,
                                pathways: list = [],
                                filter: bool = False,
                                save_path: str = 'ad_vs_no/',
                                method: str = 'cell_type_overlap'):

    """
    This function generates a heatmap of the overlapping pathways across multiple studies. The heatmap displays the adjusted
    pathway scores across different cell types for each pathway in each study. The function also returns a dictionary of
    filtered scores that contain only the overlapping pathways across the studies.

    Parameters:
    -----------
    pathway_scores : dict
        A dictionary of pathway scores for different studies.
    filtered_pathways : list, optional
        A list of pathways to be used as a filter.
    cell_types : list, optional
        A list of cell types to be included in the heatmap. Default is ["Excitatory", "Inhibitory", "Astrocyte",
        "Microglia", "Oligodendrocyte", "OPC", "Endothelial"].
    test_name : str, optional
        The name of the test being compared. Default is 'ad_vs_no'.
    top_n : int, optional
        The number of top pathways to be included in the heatmap. Default is 10.
    pathways : list, optional
        A list of pathways to be included in the heatmap. If not empty, only these pathways will be included in the
        heatmap. Default is [].
    filter : bool, optional
        If True, the function will filter out pathways that are not present in the filtered_pathways list. Default is
        False.
    save_suffix : str, optional
        A suffix to be added to the output file name. Default is 'ad_vs_no'.
    method : str, optional
        The method used to generate the overlap. 'cell_type_overlap' will generate the overlap based on cell type.
        'global_overlap' will generate the overlap based on all pathways in the studies. Default is 'cell_type_overlap'.

    Returns:
    --------
    filtered_scores : dict
        A dictionary of pathway scores for the overlapping pathways across the studies.

    Examples:
    ---------
    >>> multi_study_pathway_overlap(pathway_scores, filtered_pathways=['pathway1', 'pathway2'],
                                        cell_types=['Excitatory', 'Astrocyte'], test_name='ad_vs_no', filter=True)
    """

    scores = {}
    for i, study in enumerate(pathway_scores.keys()):
        scores[study] = pathway_scores[study][test_name][pathway_scores[study][test_name]['celltype'].isin(cell_types)]
        scores[study] = scores[study][scores[study]['P.Value'] < p_thresh]

    if method == "cell_type_overlap":
        overlap = []
        for cell_type in cell_types:
            eval_string = []
            for i, study in enumerate(pathway_scores.keys()):
                eval_string.append(f'set(scores["{study}"][scores["{study}"].celltype=="{cell_type}"].pathway)')

            eval_string = '&'.join(eval_string)
            overlap.extend(list(eval(eval_string)))

    elif method == "global_overlap":
        overlap = []
        eval_string = []
        for i, study in enumerate(pathway_scores.keys()):
            eval_string.append(f'set(scores["{study}"].pathway)')

        eval_string = '&'.join(eval_string)
        overlap.extend(list(eval(eval_string)))
            

    if filter:
        n_rows = len(set(filtered_pathways) & set(overlap))
    else:
        n_rows = len(overlap)

    if n_rows > 0:
            
        fig, axs = plt.subplots(1, 3, figsize=(3.5, n_rows*0.095), gridspec_kw={'width_ratios':[0.85, 0.85, 1]}, sharex=False,
                                    sharey=True, layout='constrained')
        fig.tight_layout()

        filtered_scores = {}
        shortened_names = {}

        for i, study in enumerate(pathway_scores.keys()):

            filtered_scores[study] = scores[study][scores[study].pathway.isin(overlap)]

            sorted_cell_types = sorted(filtered_scores[study].celltype.unique(), 
                                       key=lambda x: cell_types.index(x))
            
            filtered_scores[study] = pd.pivot_table(filtered_scores[study], values='score_adj', index='pathway', columns='celltype')
            filtered_scores[study] = filtered_scores[study][sorted_cell_types]

            if filter:
                filtered_scores[study] = filtered_scores[study].loc[filtered_scores[study].index.isin(filtered_pathways)]

            shortened_names[study] = [' '.join(name.split(" ")[:-1]) for name in filtered_scores[study].index]
            # shortened_names[study] = filtered_scores[study].index

            cbar=True if study==list(pathway_scores.keys())[-1] else False
            g1 = sb.heatmap(filtered_scores[study], cmap='bwr', center=0, vmin=-2.5, vmax=2.5, robust=False, annot=None, fmt='.1g', 
                            linewidths=linewidths, linecolor='black', annot_kws=None, cbar_kws={'shrink': 0.2}, cbar=cbar,
                            cbar_ax=None, square=False, ax=axs[i], xticklabels=sorted_cell_types, yticklabels=shortened_names[study], mask=None,) 

            axs[i].tick_params(axis='both', which='major', labelsize=2.5, length=1.5, width=0.5)

            g1.set_title(study.split('_')[-1].upper(), fontsize=3)           
            g1.set_ylabel('', fontsize=4)
            g1.set_xlabel('')

        cax = g1.figure.axes[-1]
        cax.tick_params(labelsize=4, length=1.5, width=0.5, which="major")

        # plt.tight_layout()
        # if filter:  
        #     plt.savefig(f'../results/pathway_meta_analysis/filtered_overlap_pathway_diff_exp_patterns_{save_suffix}.pdf', bbox_inches='tight')
        # else:

        plt.suptitle(f"{test_name.split('_')[0].capitalize()}- vs {test_name.split('_')[-1]}-pathology", fontsize=4)

        plt.savefig(save_path, bbox_inches='tight')       
        plt.show(block=False)  

    else:
        filtered_scores = {}
        for i, study in enumerate(pathway_scores.keys()):
            filtered_scores[study] = None
        overlap = []

        print(f'There are no overlapping pathways for {test_name} at current P.Value threshold of {p_thresh}...')
        print('Moving to next test...')


    return overlap


def save_plot(fig, ax, save):
    if save is not None:
        if ax is not None:
            if fig is not None:
                fig.savefig(save, bbox_inches='tight')
            else:
                raise ValueError("fig is None, cannot save figure.")
        else:
            raise ValueError("ax is None, cannot save figure.")
        

def check_if_matplotlib(return_mpl=False):
    if not return_mpl:
        try:
            import matplotlib.pyplot as plt
        except Exception:
            raise ImportError('matplotlib is not installed. Please install it with: pip install matplotlib')
        return plt
    else:
        try:
            import matplotlib as mpl
        except Exception:
            raise ImportError('matplotlib is not installed. Please install it with: pip install matplotlib')
        return mpl


def check_if_seaborn():
    try:
        import seaborn as sns
    except Exception:
        raise ImportError('seaborn is not installed. Please install it with: pip install seaborn')
    return sns


def check_if_adjustText():
    try:
        import adjustText as at
    except Exception:
        raise ImportError('adjustText is not installed. Please install it with: pip install adjustText')
    return at


def filter_limits(df, sign_limit=None, lFCs_limit=None):
        
    """
    Filters a DataFrame by limits of the absolute value of the columns pvals and logFCs.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame to be filtered.
    sign_limit : float, None
        The absolute value limit for the p-values. If None, defaults to infinity.
    lFCs_limit : float, None
        The absolute value limit for the logFCs. If None, defaults to infinity.

    Returns
    -------
    pd.DataFrame
        The filtered DataFrame.
    """

    # Define limits if not defined
    if sign_limit is None:
        sign_limit = np.inf
    if lFCs_limit is None:
        lFCs_limit = np.inf

    # Filter by absolute value limits
    msk_sign = df['pvals'] < np.abs(sign_limit)
    msk_lFCs = np.abs(df['logFCs']) < np.abs(lFCs_limit)
    df = df.loc[msk_sign & msk_lFCs]

    return df


def plot_volcano(data, x, y, x_label, y_label='-log10(pvals)', annotate=True, 
                    annot_by='top', names=[], 
                    top=5, sign_thr=0.05, lFCs_thr=0.5, sign_limit=None, lFCs_limit=None,
                    figsize=(7, 5), dpi=100, ax=None, return_fig=False, save=None, 
                    fontsizes={"on_plot": 4}):
    """
    Plot logFC and p-values from a long formated data-frame.

    Parameters
    ----------
    data : pd.DataFrame
        Results of DEA in long format.
    x : str
        Column name of data storing the logFCs.
    y : str
        Columns name of data storing the p-values.
    x_label: str
        Aternate name for LogFC to be included in plot. If None, defaults to x
    y_label: str
        Aternate name for p-values to be included in plot. If None, defaults to y
    annotate: bool
        Whether to annotate labels.
    annot_by: str
        Determines how to annotate the plot for top features. It can be either 'top' or 'name'. 
        If set to 'top', the top top differentially expressed features will be annotated. If set to 'name',
        only the features specified in names will be annotated.
    names: list[]:
        A list of feature names to be annotated in the plot. Only used if annot_by is set to 'name'.
    top : int
        Number of top differentially expressed features to show.
    sign_thr : float
        Significance threshold for p-values.
    lFCs_thr : float
        Significance threshold for logFCs.
    sign_limit : float
        Limit of p-values to plot in -log10.
    lFCs_limit : float
        Limit of logFCs to plot in absolute value.
    figsize : tuple
        Figure size.
    dpi : int
        DPI resolution of figure.
    ax : Axes, None
        A matplotlib axes object. If None returns new figure.
    return_fig : bool
        Whether to return a Figure object or not.
    save : str, None
        Path to where to save the plot. Infer the filetype if ending on {`.pdf`, `.png`, `.svg`}.

    Returns
    -------
    fig : Figure, None
        If return_fig, returns Figure object.
    """


    if x_label is None:
        x_label = x

    if y_label is None:
        y_label = y

    # Load plotting packages
    plt = check_if_matplotlib()
    at = check_if_adjustText()

    # Transform sign_thr
    sign_thr = -np.log10(sign_thr)

    # Extract df
    df = data.copy()
    df['logFCs'] = df[x]
    df['pvals'] = -np.log10(df[y])

    # Filter by limits
    df = filter_limits(df, sign_limit=sign_limit, lFCs_limit=lFCs_limit)

    # Define color by up or down regulation and significance
    df['weight'] = 'gray'
    up_msk = (df['logFCs'] >= lFCs_thr) & (df['pvals'] >= sign_thr)
    dw_msk = (df['logFCs'] <= -lFCs_thr) & (df['pvals'] >= sign_thr)
    df.loc[up_msk, 'weight'] = '#D62728'
    df.loc[dw_msk, 'weight'] = '#1F77B4'

    # Plot
    fig = None
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)

    n = df.shape[0]
    size = 120000 / (100*n)

    df.plot.scatter(x='logFCs', y='pvals', c='weight', sharex=False, ax=ax, s=size)

    # Draw sign lines
    ax.axhline(y=sign_thr, linestyle='--', color="black")
    ax.axvline(x=lFCs_thr, linestyle='--', color="black")
    ax.axvline(x=-lFCs_thr, linestyle='--', color="black")

    # Plot top sign features
    signs = df[up_msk | dw_msk].sort_values('pvals', ascending=False)

    # Add labels
    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)
    
    if annotate:
        if annot_by == 'top':
            signs = signs.iloc[:top]
        elif annot_by == 'name':
            signs = signs.loc[signs.index.isin(names)]
        
        texts = []
        for x, y, s in zip(signs['logFCs'], signs['pvals'], signs.index):
            texts.append(ax.text(x, y, s, fontsize=fontsizes['on_plot']))
        if len(texts) > 0:
            at.adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black'), ax=ax)

    save_plot(fig, ax, save)

    if return_fig:
        return fig