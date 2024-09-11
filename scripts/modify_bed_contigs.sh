#!/bin/bash

# Set variables
save_dir="/media/tadeoye/Volume1/SEA-AD/MTG/ATACseq/results"  # Replace ${region_name} with actual value or set it before this line
save_prefix="seaad_mtg"
working_dir="/home/tadeoye/Documents/research_codes/cell_cell_communication/scripts"

# Set input files
INPUT_BED="${save_dir}/consensus_peak_calling/consensus_regions.bed"
FASTA_FILE="${save_dir}/fasta/hg38.fa"
OUTPUT_BED="${save_dir}/consensus_peak_calling/consensus_regions_modified.bed"
TEMP_BED=$(mktemp)

# Function to get the full contig name from FASTA
get_full_contig_name() {
    local short_name=$1
    grep "^>" "$FASTA_FILE" | grep -i "$short_name" | sed 's/^>//' | head -n 1
}

# Process the BED file
while IFS=$'\t' read -r chrom start end rest; do
    if [[ $chrom == GL* || $chrom == KI* ]]; then
        full_name=$(get_full_contig_name "$chrom")
        if [[ -n $full_name ]]; then
            echo -e "${full_name}\t${start}\t${end}\t${rest}"
        else
            echo "Warning: No match found for $chrom" >&2
            echo -e "${chrom}\t${start}\t${end}\t${rest}"
        fi
    else
        echo -e "${chrom}\t${start}\t${end}\t${rest}"
    fi
done < "$INPUT_BED" > "$TEMP_BED"

# Sort the BED file
sort -V -k1,1 -k2,2n "$TEMP_BED" > "$OUTPUT_BED"

# Remove temporary file
rm "$TEMP_BED"

echo "Modified and sorted BED file saved as $OUTPUT_BED"