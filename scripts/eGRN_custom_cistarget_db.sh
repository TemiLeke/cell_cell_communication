#!/bin/bash

set -e

# Suppress warnings
export PYTHONWARNINGS="ignore:FutureWarning"

# Set variables
SAVE_PREFIX="seaad_mtg"
CELL_TYPE_COLUMN="Subclass"
REGION_NAME=$(echo "${SAVE_PREFIX}" | awk -F'_' '{print toupper($NF)}')
SAVE_DIR="/media/tadeoye/Volume1/SEA-AD/${REGION_NAME}/ATACseq/results"
TMP_DIR="/media/tadeoye/Volume1/SEA-AD/${REGION_NAME}/ATACseq/temp_files"

# Create directories
mkdir -p "${SAVE_DIR}"
mkdir -p "${TMP_DIR}"

# Clone create_cisTarget_databases repository
TARGET_DIR="${PWD}/functions"
mkdir -p "${TARGET_DIR}"
cd "${TARGET_DIR}"
git clone https://github.com/aertslab/create_cisTarget_databases
cd -

# Download cluster-buster
cd "${TARGET_DIR}"
wget https://resources.aertslab.org/cistarget/programs/cbust -O cbust
chmod a+x cbust
cd -

# Download motif collection
mkdir -p "${SAVE_DIR}/motif_collection"
wget -O "${SAVE_DIR}/motif_collection/v10nr_clust_public.zip" https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip
cd "${SAVE_DIR}/motif_collection"
unzip -q v10nr_clust_public.zip
cd -

# Prepare fasta from consensus regions
mkdir -p "${SAVE_DIR}/fasta"
TARGET_URL="https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips"
REGION_BED="${SAVE_DIR}/consensus_peak_calling/consensus_regions_modified.bed"
GENOME_FASTA="${SAVE_DIR}/fasta/hg38.fa"
CHROMSIZES="${SAVE_DIR}/fasta/hg38.chrom.sizes"
DATABASE_PREFIX="${SAVE_PREFIX}_1kb_bg_with_mask"
SCRIPT_DIR="${PWD}/functions/create_cisTarget_databases"

wget "${TARGET_URL}/hg38.fa.gz" -O "${GENOME_FASTA}.gz"
gunzip -c "${GENOME_FASTA}.gz" > "${GENOME_FASTA}"
wget "${TARGET_URL}/hg38.chrom.sizes" -O "${CHROMSIZES}"

"${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh" \
    "${GENOME_FASTA}" \
    "${CHROMSIZES}" \
    "${REGION_BED}" \
    "${SAVE_DIR}/fasta/hg38.${SAVE_PREFIX}.with_1kb_bg_padding.fa" \
    1000 \
    yes

# Create cistarget databases
ls "${SAVE_DIR}/motif_collection/v10nr_clust_public/singletons" > "${SAVE_DIR}/motif_collection/motifs.txt"

export PATH="${PATH}:${PWD}/functions/cbust"

OUT_DIR="${SAVE_DIR}/motif_collection"
CBDIR="${PWD}/functions/cbust"
FASTA_FILE="${SAVE_DIR}/fasta/hg38.${SAVE_PREFIX}.with_1kb_bg_padding.fa"
MOTIF_LIST="${OUT_DIR}/motifs.txt"
MOTIF_DIR="${SAVE_DIR}/motif_collection/v10nr_clust_public/singletons"

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f "${FASTA_FILE}" \
    -c "${CBDIR}" \
    -M "${MOTIF_DIR}" \
    -m "${MOTIF_LIST}" \
    -o "${OUT_DIR}/${DATABASE_PREFIX}" \
    --bgpadding 1000 \
    -t 40

echo "Script completed successfully"