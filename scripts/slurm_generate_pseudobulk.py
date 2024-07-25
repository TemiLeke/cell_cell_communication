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

    for key in list(fragments_dict.keys()):
        if key not in cell_data['sample_id'].unique():
            del fragments_dict[key]

    if exclude_cells != None:
        cells_to_include = [ctype for ctype in list(itertools.chain(*list(subclass[save_prefix].values()))) if ctype not in exclude_cells]
        cell_data = cell_data[cell_data['celltype'].isin(cells_to_include)]
    
    if include_cells != None:
        cells_to_include = [ctype for ctype in list(itertools.chain(*list(subclass[save_prefix].values()))) if ctype in include_cells]
        cell_data = cell_data[cell_data['celltype'].isin(cells_to_include)]

else:
    
    adata = sc.read_h5ad(f'../data/SEA-AD/{region_name}/ATACseq/anndata/SEAAD_MTG_ATACseq_final-nuclei.2024-02-13.h5ad')
    # adata = sc.read_h5ad('../data/SEA-AD/{region_name}/ATACseq/anndata/SEAAD_MTG_ATACseq_all-nuclei.2024-02-13.h5ad')
    cell_data = adata.obs.copy()
    del (adata)

    # cell_data = cell_data[~cell_data['Supertype confidence'].isna()]
    cell_data = cell_data[cell_data['Neurotypical reference'] == 'False']
    cell_data = cell_data[cell_data['Severely Affected Donor'] == 'N']
    cell_data = cell_data[~cell_data['sample_id'].isna()]
    cell_data['celltype'] = cell_data[cell_type_column] # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.
    cell_data = cell_data[~cell_data['celltype'].isna()]


    cell_data['cell_barcode'] = cell_data['sample_id'].apply(lambda x: "-".join([x.split('-')[0], '1']))
    cell_data['sample_id'] = cell_data['sample_id'].apply(lambda x: x.split('-')[-1]).astype(str)
    cell_data['cell_barcode'] = cell_data['cell_barcode'] + '.' + cell_data['sample_id']

    files = list(set(cell_data['sample_id'].unique()) & set(fragments_metadata['fileID'].unique()))
    cell_data = cell_data[cell_data['sample_id'].isin(files)]

    cell_data.set_index('cell_barcode', inplace=True)

    for key in list(fragments_dict.keys()):
        if key not in cell_data['sample_id'].unique():
            del fragments_dict[key]

    cell_data.to_csv(f'../data/raw/{save_prefix}/{save_prefix}_cell_data.csv')

    if exclude_cells != None:
        cells_to_include = [ctype for ctype in list(itertools.chain(*list(subclass[save_prefix].values()))) if ctype not in exclude_cells]
        cell_data = cell_data[cell_data['celltype'].isin(cells_to_include)]
    
    if include_cells != None:
        cells_to_include = [ctype for ctype in list(itertools.chain(*list(subclass[save_prefix].values()))) if ctype in include_cells]
        cell_data = cell_data[cell_data['celltype'].isin(cells_to_include)]

print(f'Total number of cells: {len(cell_data)}')

cell_data[cell_type_column].value_counts()

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

if load_pseudobulk_profiles:
    
    log.info(f"loading bigwig and bed files in {os.path.join(save_dir, 'consensus_peak_calling/')}")

    bed_paths = pickle.load(open(os.path.join(save_dir, 'consensus_peak_calling/pseudobulk_bed_files/bed_paths.pkl'), 'rb'))
    bw_paths =  pickle.load(open(os.path.join(save_dir, 'consensus_peak_calling/pseudobulk_bw_files/bw_paths.pkl'), 'rb'))

else:
    bw_paths, bed_paths = export_pseudobulk(input_data = cell_data,
                                            variable = 'celltype',                                                                     # variable by which to generate pseubulk profiles, in this case we want pseudobulks per celltype
                                            sample_id_col = 'sample_id',
                                            chromsizes = chromsizes,
                                            bed_path = os.path.join(save_dir, 'consensus_peak_calling/pseudobulk_bed_files/'),  # specify where pseudobulk_bed_files should be stored
                                            bigwig_path = os.path.join(save_dir, 'consensus_peak_calling/pseudobulk_bw_files/'), # specify where pseudobulk_bw_files should be stored
                                            path_to_fragments = fragments_dict,                                                        # location of fragment fiels
                                            n_cpu = 1,                                                                                 # specify the number of cores to use, we use ray for multi processing
                                            normalize_bigwig = True,
                                            #remove_duplicates = True,
                                            temp_dir = os.path.join(tmp_dir, 'ray_spill'),
                                            split_pattern = '.',
                                            stage = stage,
                                            verbose = True)
    
    for key in list(bed_paths.keys()):
        bed_paths[key.replace(' ', '_').replace('/', '_')] = bed_paths[key]
        del bed_paths[key]
        bw_paths[key.replace(' ', '_').replace('/', '_')] = bw_paths[key]
        del bw_paths[key]
        
    pickle.dump(bed_paths, 
            open(os.path.join(save_dir, 'consensus_peak_calling/pseudobulk_bed_files/bed_paths.pkl'), 'wb'))
    pickle.dump(bw_paths,
            open(os.path.join(save_dir, 'consensus_peak_calling/pseudobulk_bw_files/bw_paths.pkl'), 'wb'))


# %% [markdown]

# ### **Call peaks per pseudobulk profile**

macs_path='/work_bgfs/t/tadeoye/scRNAseq_AD_meta_analysis/scripts/scenic_env/bin/macs2'
# # tmp_dir = f'/Volumes/Seagate/temp_files's
# # Run peak calling
# narrow_peaks_dict = peak_calling(macs_path,
#                                 bed_paths,
#                                 os.path.join(save_dir, 'consensus_peak_calling/MACS/'),
#                                 genome_size='hs',
#                                 n_cpu=1,
#                                 input_format='BEDPE',
#                                 shift=73, 
#                                 ext_size=146,
#                                 keep_dup = 'all',
#                                 q_value = 0.05,
#                                 _temp_dir = os.path.join(tmp_dir, 'ray_spill'))