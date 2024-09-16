# %%
import os
import sys
import warnings
import subprocess


#supress warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

_stderr = sys.stderr                                                         
null = open(os.devnull,'wb')

# %%
os.getcwd()

# %%
save_prefix = 'seaad_mtg' # this takes the format '{StudyName}_{ThreeLetterAccronymForBrainRegion}'

cell_type_column = 'Subclass' # 'Supertype (non-expanded)', 'Subclass'

region_name = save_prefix.split('_')[-1].upper()
save_dir = f'/media/tadeoye/Volume1/SEA-AD/{region_name}/ATACseq/results'
tmp_dir = f'/media/tadeoye/Volume1/SEA-AD/{region_name}/ATACseq/temp_files'

if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# %% [markdown]
# # **Creating a custom cistarget database**
# 
# Here we will create the custom cistarget database using consensus peaks.
# 
# This involves precomputed scores for all the motifs in our motif collection on a predefined set of regions
# 
# The developers already provided precomputed databases for [human](https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/), [mouse](https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/) and [fly](https://resources.aertslab.org/cistarget/databases/drosophila_melanogaster/dm6/flybase_r6.02/mc_v10_clust/region_based/). These databases are computed on regulatory regions spanning the genome. For the best result the authors recommend generating a custom database, given that it is highly likely that the precomputed databases don't cover all the regions in our consensus peak set.

# %% [markdown]
# ## **Download `create_cistarget_database`**
# 
# We will start by downloading and installing the `create_cistarget_database` repository.

# %%
target_dir = os.path.join(os.getcwd(), 'functions')
working_dir = os.getcwd()
os.chdir(target_dir)
os.system('git clone https://github.com/aertslab/create_cisTarget_databases')       
os.chdir(working_dir)

# %% [markdown]
# ## **Download `cluster-buster`**
# 
# [Cluster-buster](https://github.com/weng-lab/cluster-buster) will be used to score the regions using the motif collection. We used the precompiled binary of cluster buster at **https://github.com/weng-lab/cluster-buster**
# 

# %%
target_dir = os.path.join(os.getcwd(), 'functions')
working_dir = os.getcwd()
os.chdir(target_dir )
os.system('wget https://resources.aertslab.org/cistarget/programs/cbust -O cbust')       
os.system('chmod a+x cbust')
os.chdir(working_dir)

# %% [markdown]
# ## **Download motif collection**
# 
# Next, we will download the motif collection.

# %%
os.makedirs(f"{save_dir}/motif_collection", exist_ok=True)
os.system(f'wget -O {save_dir}/motif_collection/v10nr_clust_public.zip https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip')

# %%
os.chdir(f"{save_dir}/motif_collection")
os.system('yes Y | unzip -q v10nr_clust_public.zip')       
os.chdir(working_dir)

# %% [markdown]
# These are the motif-to-TF annotations for:
# 
# - **`Chicken`**: motifs-v10-nr.chicken-m0.00001-o0.0.tbl
# - **`fly`**: motifs-v10-nr.flybase-m0.00001-o0.0.tbl
# - **`human`**: motifs-v10-nr.hgnc-m0.00001-o0.0.tbl
# - **`mouse`**: motifs-v10-nr.mgi-m0.00001-o0.0.tbl

# %%
os.system(f'ls {save_dir}/motif_collection/v10nr_clust_public/snapshots/')

# %% [markdown]
# Here are some example motifs, they are stored in cb format.

# %%
os.system(f'ls -l {save_dir}/motif_collection/v10nr_clust_public/singletons | head')

# %%
os.system(f'cat {save_dir}/motif_collection/v10nr_clust_public/singletons/bergman__Adf1.cb')

# %% [markdown]
# ## **Prepare fasta from consensus regions**
# 
# Next we will get sequences for all the consensus peaks. We will also add 1kb of background padding, this will be used as a background sequence for cluster-buster. It is completely optional to add this padding, the authors have noticed that it does not affect the analyses a lot.

# %%
os.makedirs(f"{save_dir}/fasta", exist_ok=True)

target_url='https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips'

REGION_BED = f"{save_dir}/consensus_peak_calling/consensus_regions_modified.bed"
GENOME_FASTA = f"{save_dir}/fasta/hg38.fa"
CHROMSIZES = f"{save_dir}/fasta/hg38.chrom.sizes"
DATABASE_PREFIX = f"{save_prefix}_1kb_bg_with_mask"
SCRIPT_DIR = f"{working_dir}/functions/create_cisTarget_databases"
SAVE_PREFIX = save_prefix

# %% [markdown]
# We can run this in python using

# %%

os.system(f'wget {target_url}/hg38.fa.gz -O {GENOME_FASTA}.gz')
os.system(f'yes Y | gunzip -c {GENOME_FASTA}.gz > {GENOME_FASTA}') 
os.system(f'wget {target_url}/hg38.chrom.sizes -O {CHROMSIZES}')


shell_command = f"""
{SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \\
    {GENOME_FASTA} \\
    {CHROMSIZES} \\
    {REGION_BED} \\
    {save_dir}/fasta/hg38.{SAVE_PREFIX}.with_1kb_bg_padding.fa \\
    1000 \\
    yes
"""

subprocess.run(shell_command, shell=True, check=True)

# %%
os.system(f"head -n 2 {save_dir}/fasta/hg38.{SAVE_PREFIX}.with_1kb_bg_padding.fa")

# %% [markdown]
# ## Create cistarget databases
# 
# Now we can create the ranking and score database. This step will take some time so we recommend to run it as a job (i.e. not in jupyter notebooks).

# %%
os.system(f'ls {save_dir}/motif_collection/v10nr_clust_public/singletons > {save_dir}/motif_collection/motifs.txt')

# %%
os.system(f'export PATH=$PATH:/{working_dir}/functions/cbust')

# %%
OUT_DIR=f"{save_dir}/motif_collection"
CBDIR=f"{working_dir}/functions/cbust"
FASTA_FILE=f"{save_dir}/fasta/hg38.{SAVE_PREFIX}.with_1kb_bg_padding.fa"
MOTIF_LIST=f"{OUT_DIR}/motifs.txt"
MOTIF_DIR = f"{save_dir}/motif_collection/v10nr_clust_public/singletons"

cmd = [
    f"{SCRIPT_DIR}/create_cistarget_motif_databases.py",
    "-f", FASTA_FILE,
    "-c", CBDIR,
    "-M", MOTIF_DIR,
    "-m", MOTIF_LIST,
    "-o", f"{OUT_DIR}/{DATABASE_PREFIX}",
    "--bgpadding", "1000",
    "-t", "40"
]

try:
    result = subprocess.run(cmd, check=True, text=True, capture_output=True)
    print("Command executed successfully")
    print("Output:", result.stdout)
except subprocess.CalledProcessError as e:
    print("An error occurred while running the command")
    print("Error output:", e.stderr)
    print("Return code:", e.returncode)
    print("Command output:", e.stdout)
except FileNotFoundError:
    print(f"The script {cmd[0]} was not found. Please check the path.")

# If you want to see the command that was run:
print("Command:", " ".join(cmd))

# %% [markdown]
# ## **Get Moftif Enrichment**

# %% [markdown]
# We next use this custom cistarget database to estimate motif enrichment in `/scripts/eGRN_motif_enrichment.ipynb`

# %% [markdown]
# 


