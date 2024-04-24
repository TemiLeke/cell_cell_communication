# %%
import os
import sys
import pickle
import warnings
import numpy as np
import pandas as pd
import pandas as pd
import scanpy as sc
import pyranges as pr
import matplotlib.pyplot as plt
from matplotlib.pyplot import rcParams
from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import *
from pycisTopic.clust_vis import *
from pycisTopic.topic_binarization import *
from pycisTopic.diff_features import *
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk


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
load_cell_data = True
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
save_dir = f'/Volumes/Seagate/SEA-AD/{region_name}/ATACseq/results'
tmp_dir = f'/Volumes/Seagate/SEA-AD/{region_name}/ATACseq/temp_files'

if not os.path.exists(save_dir):
    os.makedirs(save_dir)


fragments_dir = f'/Volumes/Seagate/SEA-AD/{region_name}/ATACseq/fragments/'
fragments_dict = {fragment.split('-')[0]: fragments_dir + fragment for fragment in os.listdir(fragments_dir) if fragment.endswith('.tsv.gz')}
fragments_metadata = pd.read_csv(fragments_dir+f'{save_prefix}_fragments_metadata.csv')
fragments_metadata['fileID'] = fragments_metadata['fileID'].astype(str)

# %% [markdown]
# ## **Generate pseudobulk ATAC-seq profiles, call peaks and generate a consensus peak set**
# 
# We will use the cell type labels from the scRNA-seq side of the data to generate pseudobulk profiles per cell type. These pseudobulk profiles will be used to call peaks which will be merged in one consensus peak set.
# 
# The advantage of calling peaks per cell type (instead of calling peaks on the whole bulk profile at once) is that, for rare cell types, there is a risk that the ATAC signal of these cells might be confused for noise by the peak caller when viewed in the whole bulk profile. Calling peaks per cell type helps resolving these rare peaks and thus increases the resolution of the data.
# 
# We will first load the cell type annotation in the scRNA-seq data

# %% [markdown]
# ### **Get Cell Type Annotations**
# 

# %%
if load_cell_data:
    cell_data = pd.read_csv(f'../data/raw/{save_prefix}/{save_prefix}_cell_data.csv', dtype=str)
    cell_data.set_index('cell_barcode', inplace=True)
else:
    
    obs_data = {}

    for cell_group in subclass[save_prefix].keys():
        adata = sc.read_h5ad(f'../data/raw/{save_prefix}/anndata/atac/{cell_group}_raw_anndata.h5ad')
        adata.obs_names_make_unique()
        adata.var_names_make_unique()
        obs_data[cell_group] = adata.obs

    cell_data = pd.concat([df for df in obs_data.values()], ignore_index=False)
    cell_data = pd.merge(cell_data, fragments_metadata, left_on='Donor ID', right_on='individualID', how='outer')

    # adata = sc.read_h5ad('/Volumes/Seagate/SEA-AD/MTG/ATACseq/anndata/SEAAD_MTG_ATACseq_all-nuclei.2024-02-13.h5ad')
    # cell_data = pd.merge(adata.obs, fragments_metadata, left_on='Donor ID', right_on='individualID', how='outer')

    del(obs_data)
    del (adata)

    cell_data = cell_data[~cell_data['sample_id'].isna()]
    cell_data['cell_barcode'] = cell_data['sample_id'].apply(lambda x: "-".join([x.split('-')[0], '1']))
    cell_data['cell_barcode'] = cell_data['cell_barcode'] + '.' + cell_data['fileID']
    cell_data.set_index('cell_barcode', inplace=True)

    cell_data['sample_id'] = cell_data['fileID']
    cell_data['celltype'] = cell_data[cell_type_column] # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.

    cell_data = cell_data[cell_data['Neurotypical reference'] == 'False']
    cell_data = cell_data[cell_data['Severely Affected Donor'] == 'N']
    cell_data = cell_data[~cell_data['fileID'].isna()]
    cell_data = cell_data[~cell_data['celltype'].isna()]

    cell_data.to_csv(f'../data/raw/{save_prefix}/{save_prefix}_cell_data.csv')

# %%
for key in list(fragments_dict.keys()):
    if key not in cell_data['sample_id'].unique():
        del fragments_dict[key]

# %% [markdown]
# ### **Get Chromosome Sizes**

# %%
# Get chromosome sizes (for hg38 here)
target_url='http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'
chromsizes=pd.read_csv(target_url, sep='\t', header=None)
chromsizes.columns=['Chromosome', 'End']
chromsizes['Start']=[0]*chromsizes.shape[0]
chromsizes=chromsizes.loc[:,['Chromosome', 'Start', 'End']]
# Exceptionally in this case, to agree with CellRangerARC annotations
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].replace('v', '.') for x in range(len(chromsizes['Chromosome']))]
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].split('_')[1] if len(chromsizes['Chromosome'][x].split('_')) > 1 else chromsizes['Chromosome'][x] for x in range(len(chromsizes['Chromosome']))]
chromsizes=pr.PyRanges(chromsizes)

# %% [markdown]
# ### **Generate Pseudobulk Profile**
# 
# Next, we will generate the pseudobulk profiles. 
# 
# This will generate two sets of files:
# 
# 1. pseudobulk_bed_files: pseudobulk profiles stored in [bed format](https://genome.ucsc.edu/FAQ/FAQformat.html).
# 2. pseudobulk_bw_files: pseudobulk profiles stored in [BigWig format](https://genome.ucsc.edu/goldenpath/help/bigWig.html).
# 
# The BigWig files are useful for visualization in [IGV](https://software.broadinstitute.org/software/igv/) or [UCSC genome browser](https://genome.ucsc.edu/).
# 

# %%
bw_paths, bed_paths = export_pseudobulk(input_data = cell_data,
                                        variable = 'celltype',                                                                     # variable by which to generate pseubulk profiles, in this case we want pseudobulks per celltype
                                        sample_id_col = 'sample_id',
                                        chromsizes = chromsizes,
                                        bed_path = os.path.join(save_dir, 'consensus_peak_calling/pseudobulk_bed_files/'),  # specify where pseudobulk_bed_files should be stored
                                        bigwig_path = os.path.join(save_dir, 'consensus_peak_calling/pseudobulk_bw_files/'), # specify where pseudobulk_bw_files should be stored
                                        path_to_fragments = fragments_dict,                                                        # location of fragment fiels
                                        n_cpu = 2,                                                                                 # specify the number of cores to use, we use ray for multi processing
                                        normalize_bigwig = True,
                                        #remove_duplicates = True,
                                        temp_dir = os.path.join(tmp_dir, 'ray_spill'),
                                        split_pattern = '.')

print(bw_paths)
print(bed_paths)