#!/bin/bash

# Uncomment the following line if needed for your cluster
# module load cluster/wice/bigmem

# Set variables
save_dir="/work_bgfs/t/tadeoye/scRNAseq_AD_meta_analysis/data/SEA-AD/MTG/ATACseq/results"  # Replace ${region_name} with actual value or set it before this line
save_prefix="seaad_mtg"
working_dir="/work_bgfs/t/tadeoye/scRNAseq_AD_meta_analysis/scripts"

# Set URLs and file paths
target_url="https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips"
REGION_BED="${save_dir}/consensus_peak_calling/consensus_regions.bed"
GENOME_FASTA="${save_dir}/fasta/hg38.fa"
CHROMSIZES="${save_dir}/fasta/hg38.chrom.sizes"
DATABASE_PREFIX="${save_prefix}_1kb_bg_with_mask"
SCRIPT_DIR="${working_dir}/create_cisTarget_databases"
SAVE_PREFIX="${save_prefix}"

# Creating a custom cistarget database
git clone https://github.com/aertslab/create_cisTarget_databases

# Download `cluster-buster`
cd functions/
wget https://resources.aertslab.org/cistarget/programs/cbust -O cbust
chmod a+x cbust
cd ..

# Download motif collection
mkdir -p "${save_dir}/motif_collection"
wget -O "${save_dir}/motif_collection/v10nr_clust_public.zip" https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip
cd "${save_dir}/motif_collection"
yes Y | unzip -q v10nr_clust_public.zip
cd "${working_dir}"

# Create directory
mkdir -p "${save_dir}/fasta"

# Download and extract genome fasta
wget "${target_url}/hg38.fa.gz" -O "${GENOME_FASTA}.gz"
yes Y | gunzip -c "${GENOME_FASTA}.gz" > "${GENOME_FASTA}"

# Download chrom sizes
wget "${target_url}/hg38.chrom.sizes" -O "${CHROMSIZES}"

# Load BEDTools module
# module load apps/bedtools/2.30.0

export PATH="/work_bgfs/t/tadeoye/scRNAseq_AD_meta_analysis/scripts/bedtools/bin:$PATH"
chmod a+x /work_bgfs/t/tadeoye/scRNAseq_AD_meta_analysis/scripts/bedtools/bin/bedtools

# Run the create_fasta_with_padded_bg_from_bed.sh script
"${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh" \
    "${GENOME_FASTA}" \
    "${CHROMSIZES}" \
    "${REGION_BED}" \
    "hg38.${SAVE_PREFIX}.with_1kb_bg_padding.fa" \
    1000 \
    yes