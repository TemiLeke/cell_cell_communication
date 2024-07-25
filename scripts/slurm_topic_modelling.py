# %%
import os
import sys
import pickle
import requests
import warnings
import logging
import itertools
import numpy as np
import pandas as pd
import pandas as pd
import scanpy as sc
import pyranges as pr
import pybiomart as pbm
from pycisTopic.qc import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import rcParams
from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import *
from pycisTopic.clust_vis import *
from pycisTopic.topic_binarization import *
from pycisTopic.diff_features import *
from pycisTopic.iterative_peak_calling import *
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling


#supress warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

_stderr = sys.stderr                                                         
null = open(os.devnull,'wb')

# %% [markdown]
# # **scATAC-seq preprocessing using pycisTopic**
# 
# Now we will preprocess the scATAC-seq side of the multiome data.
# 
# Most importantly we will:
# 
# 1. generate pseudobulk ATAC-seq profiles per cell type/subclass and call peaks
# 2. merge these peaks into a consensus peak-set
# 3. do quality control on the scATAC-seq barcodes
# 4. run topic modeling to find sets of co-accessible regions and to impute chromatin accessibility resolving the issue of drop outs
# 
# For this we will use the python package [pycisTopic](https://pycistopic.readthedocs.io/en/latest/).

# %% [markdown]
# 
# ### **Data Prep Parameters**
# 
# Specify aditional information and importantly, the location of the ATAC fragments file, whihc is the main input into `pycisTopic`

# %%
save_prefix = 'seaad_mtg' # this takes the format '{StudyName}_{ThreeLetterAccronymForBrainRegion}'
load_cell_data = True           # whether to load cell type annotations or to use dataframe containing metadata along with filtered cell type barcode and sample information
load_pseudobulk_profiles = True # Whether to split and (or) merge fragments into cell type specific bed or bigwig files or just load ones already split and merged
stage = 'split_n_merge'
load_peaks = True              # whether to load dictionary already called peaks
run_consensus_peaks = False      # whether to generate a consensus peak set
run_qc = False
run_cistopic_sample_processing = False
load_cistopic_object = True
use_mallet = True

exclude_cells = ['L4 IT', 'L5 IT', 'Vip', 'Pvalb', 'Oligodendrocyte', 'Sst', 'L6 IT', 'Lamp5', 'Astrocyte', 'L6 IT Car3', 'Sncg', 'Microglia-PVM', 'OPC', 'Lamp5 Lhx6', 'L5/6 NP', 'L6 CT', 'L6b', 'Pax6', 'Chandelier', 'VLMC', 'L5 ET', 'Sst Chodl', 'Endothelial']
include_cells = None

cell_type_column = 'Subclass' # 'Supertype (non-expanded)', 'Subclass'

subclass = {'seaad_mtg': {'excitatory': ['L5 IT', 'L2/3 IT', 'L4 IT', 'L6 IT', 'L6 IT Car3', 'L5/6 NP', 'L6b', 'L6 CT', 'L5 ET'],
                          'inhibitory': ['Pvalb', 'Sst', 'Lamp5 Lhx6', 'Vip', 'Lamp5', 'Sncg', 'Chandelier', 'Sst Chodl', 'Pax6'],
                          'astrocyte': ['Astrocyte'],
                          'microglia': ['Microglia-PVM'],
                          'opc': ['OPC'],
                          'oligodendrocyte': ['Oligodendrocyte'],
                          'endothelial': ['Endothelial'],
                           'vlmc': ['VLMC'],
                         },

            'seaad_pfc': {'excitatory': ['L5 IT', 'L2/3 IT', 'L4 IT', 'L6 IT', 'L6 IT Car3', 'L5/6 NP', 'L6b', 'L6 CT', 'L5 ET'],
                          'inhibitory': ['Pvalb', 'Sst', 'Lamp5 Lhx6', 'Vip', 'Lamp5', 'Sncg', 'Chandelier', 'Sst Chodl', 'Pax6'],
                          'astrocyte': ['Astrocyte'],
                          'microglia': ['Microglia-PVM'],
                          'opc': ['OPC'],
                          'oligodendrocyte': ['Oligodendrocyte'],
                          'endothelial': ['Endothelial'],
                           'vlmc': ['VLMC'],
                         },
            
            'gazestani_pfc': {'excitatory': ['FEZF2', 'RORB', 'THEMIS', 'LINC00507', 'CTGF'],
                              'inhibitory': ['PVALB', 'SST', 'LHX6', 'LAMP5', 'VIP', 'NDNF'],
                              'astrocyte': ['WIF1', 'CHI3L1', 'PTCSC3', 'GRIA1'],
                              'microglia': ['Myeloid', 'CX3CR1', 'GPNMB', 'Prolif', 'IRM', 'Macrophage', 'CRM'],
                              'opc': ['PBX3', 'ANKRD55', 'BRCA2', 'OPC', 'PLP1'],
                              'oligodendrocyte': ['GRIK2', 'BACE2', 'PLXDC2', 'SLC38A1'],
                              'endothelial': ['ABCC9', 'TMEM45B', 'GPR126', 'C7', 'HRK', 'IGFBP6', 'DOCK8', 'G0S2', 'APLN', 'COL8A1'],
                        }
            }

cell_group = 'astrocyte'
cell_supertype = subclass[save_prefix][cell_group][0]

map_meta = True

get_cell_types = False   # whether to reformat cell_type annotation

subject_ids_for_study = {'leng_sfg': 'PatientID',
                        'leng_etc': 'PatientID',
                        'seaad_mtg': 'Donor ID', 
                        'seaad_pfc': 'Donor ID', 
                        'gazestani_pfc': 'individualID'}

subject_id = subject_ids_for_study[save_prefix]     # for leng this is `PatientID` for mathys is 'Subject', and allen is 'individualID'
metadata = f'../data/raw/{save_prefix}/{save_prefix}_metadata.csv' # Metatdata location
meta = pd.read_csv(metadata, encoding_errors='ignore')

region_name = save_prefix.split('_')[-1].upper()

save_dir = f'../data/SEA-AD/{region_name}/ATACseq/results'
tmp_dir = f'../data/SEA-AD/{region_name}/ATACseq/temp_files'

if not os.path.exists(save_dir):
    os.makedirs(save_dir)


fragments_dir = f'../data/SEA-AD/{region_name}/ATACseq/fragments/'
fragments_dict = {fragment.split('-')[0]: fragments_dir + fragment for fragment in os.listdir(fragments_dir) if fragment.endswith('.tsv.gz')}
fragments_metadata = pd.read_csv(fragments_dir+f'{save_prefix}_fragments_metadata.csv')
fragments_metadata['fileID'] = fragments_metadata['fileID'].astype(str)


# %%

## **Creating a cisTopic object**

# Now that we have good quality barcodes we will generate a binary count matrix of ATAC-seq fragments over consensus peaks. 
# This matrix, along with metadata, will be stored in a cisTopic object and be used for topic modeling.
# For SCENIC+ we will only keep cells which passed quality metrics in the scATA-seq assay. Also, **[Blacklist regions](https://www.nature.com/articles/s41598-019-45839-z)** will be removed from this count matrix.

path_to_regions = os.path.join(save_dir, "consensus_peak_calling/consensus_regions.bed")
path_to_blacklist= '../scripts/functions/pycisTopic/blacklist/hg38-blacklist.v2.bed'
pycistopic_qc_output_dir = os.path.join(save_dir, "qc")

os.makedirs(f"{save_dir}/cistopic_objects", exist_ok=True)

from pycisTopic.cistopic_class import parallel_cistopic_processing

if run_cistopic_sample_processing:

    parallel_cistopic_processing(
                            fragments_dict = fragments_dict,
                            pycistopic_qc_output_dir = pycistopic_qc_output_dir,
                            sample_id_to_barcodes_passing_filters = sample_id_to_barcodes_passing_filters,
                            cell_data = cell_data,
                            path_to_regions = path_to_regions,
                            path_to_blacklist = path_to_blacklist,
                            save_dir = save_dir,
                            n_jobs = 8,
                            use_partition = False
                            )

if load_cistopic_object:

    print('Loading cistopic object!')
    
    atac = sc.read_h5ad(os.path.join(save_dir, 'cistopic_objects/all_sample_cistopic_anndata.h5ad'))
    
    cistopic_obj = create_cistopic_object(fragment_matrix=atac.layers['counts'].T.astype(int),
                                    cell_names=atac.obs_names.to_list(),
                                    region_names=atac.var_names.to_list(),
                                    tag_cells=False
                                    )
    cistopic_obj.add_cell_data(atac.obs)   

    del atac

else:

    cistopic_obj_list = []
    for sample_id in fragments_dict:
        cistopic_obj_list.append(pickle.load(open(os.path.join(save_dir, f'cistopic_objects/{sample_id}_cistopic.pkl'), 'rb')))

    del cistopic_obj_list
        
    cistopic_obj = merge(cistopic_obj_list, split_pattern='.')

    cistopic_metadata = cistopic_obj.cell_data.copy()
    cistopic_metadata = cistopic_metadata.merge(cell_data, how='left', left_index=True, right_index=True)

    cistopic_metadata = cistopic_metadata.drop(columns=[col for col in cistopic_metadata.columns if col.endswith("_y")])
    cistopic_metadata = cistopic_metadata.rename(columns={col: col[:-2] for col in cistopic_metadata.columns if col.endswith("_x")})

    cistopic_metadata = cistopic_metadata.loc[:, ~cistopic_metadata.columns.duplicated(keep='first')]

    cistopic_obj.add_cell_data(cistopic_metadata)

    pickle.dump(cistopic_obj, open(os.path.join(save_dir, 'cistopic_objects/all_sample_cistopic.pkl'), 'wb')) 

    del cistopic_metadata
    
    adata = sc.AnnData(cistopic_obj.binary_matrix.T)
    adata.layers['counts'] = cistopic_obj.fragment_matrix.T
    adata.obs = cistopic_obj.cell_data
    adata.var = cistopic_obj.region_data

    atac = sc.read_h5ad('/media/tadeoye/Volume1/SEA-AD/MTG/ATACseq/anndata/SEAAD_MTG_ATACseq_final-nuclei.2024-02-13.h5ad')

    atac.obs['cell_barcode'] = atac.obs['sample_id'].apply(lambda x: "-".join([x.split('-')[0], '1']))
    atac.obs['sample_id'] = atac.obs['sample_id'].apply(lambda x: x.split('-')[-1]).astype(str)
    atac.obs['cell_barcode'] = atac.obs['cell_barcode'] + '.' + atac.obs['sample_id']
    atac.obs.set_index('cell_barcode', inplace=True)
    atac = atac[adata.obs.index]

    adata.uns = atac.uns.copy()
    adata.obsm = atac.obsm.copy()
    adata.obsp = atac.obsp.copy()

    del atac

    for col in list(adata.obs.columns):
        if adata.obs[col].dtype=='O':
            if pd.isna(adata.obs[col][0]):
                adata.obs[col] = adata.obs[col].astype('category')
            elif adata.obs[col][0].dtype in ['int64', 'int32', 'float32', 'float64']:
                adata.obs[col] = adata.obs[col].astype(float)
            else:
                adata.obs[col] = adata.obs[col].astype('category')

    adata.write_h5ad(os.path.join(save_dir, 'cistopic_objects/all_sample_cistopic_anndata.h5ad'), compression='gzip')   

    del adata


# %%
# 
#     
### **Topic modeling**

# if use_mallet:
#     target_dir = os.path.join(os.getcwd())
#     os.chdir(target_dir)
#     os.system('wget https://github.com/mimno/Mallet/releases/download/v202108/Mallet-202108-bin.tar.gz')       
#     os.system('tar -xf Mallet-202108-bin.tar.gz')
#     os.chdir('../')

os.makedirs(os.path.join(save_dir, f'models/{cell_type_column}_models_500_iter_LDA.pkl'), exist_ok=True)

if use_mallet:

    from pycisTopic.lda_models import run_cgs_models_mallet

    mallet_path=os.path.join(os.getcwd(), "Mallet-202108/bin/mallet")

    os.makedirs(os.path.join(tmp_dir, 'ray_spill', 'mallet'), exist_ok=True)
    os.environ['MALLET_MEMORY'] = '200G'

    os.environ.update({'MALLET_HOME': r'/work_bgfs/t/tadeoye/scRNAseq_AD_meta_analysis/scripts/Mallet-202108/'})

    models=run_cgs_models_mallet(
                                cistopic_obj,
                                n_topics=[48],
                                n_cpu=40,
                                n_iter=500,
                                random_state=555,
                                alpha=50,
                                alpha_by_topic=True,
                                eta=0.1,
                                eta_by_topic=False,
                                tmp_path=os.path.join(tmp_dir, 'ray_spill', 'mallet'),
                                save_path=os.path.join(save_dir, f'models/{cell_type_column}_models_500_iter_LDA.pkl'),
                                mallet_path=mallet_path,
                                reuse_corpus=True,
                                run_models_in_parallel=False
                            )
    
else:
    from pycisTopic.lda_models import run_cgs_models

    models=run_cgs_models(cistopic_obj,
                        n_topics=[2, 4, 10, 16, 24, 32, 48],
                        n_cpu=40,
                        n_iter=500,
                        random_state=555,
                        alpha=50,
                        alpha_by_topic=True,
                        eta=0.1,
                        eta_by_topic=False,
                        save_path=os.path.join(save_dir, f'models/{cell_type_column}_models_500_iter_LDA.pkl'),
                    )