#!/bin/bash

# Bash script implementing a pipeline to:
# - Filter snpeff-consistency result by mode and effect using 'filter-consistency.py';
# - Convert filtered consistency table to BED format with 'consistency-to-bed.py'
# - Filter in SNPs with 'vcftools' and the bed file from above;
# - Rename newly generated vcfs and clean up

# This is a suggestion in how to run all these steps at once.
# If you have one vcf for each chromosome, you can loop through them
# or run them in parallel, by running this script in diffrent terminal
# windows (or using GNU parallel).

# Stop on errors
set -e

# Parameters
FILTERCONSISTENCY="filter-consistency.py"
MODE="rule"
EFFECTS='NON_SYNONYMOUS_CODING SYNONYMOUS_CODING INTRON'
CONSISTENCYTOBED="consistency-to-bed.py"
FILENAMEPREFIX="test_ann_dataset"
INPUTDIR="../testdata/snpeff-consistency"
OUTPUTDIR="../testdata/snpeff-consistency/filtered-consistency"
VCFDIR="../testdata"

# Check if the scripts and other programs exists

## Check filter-consistency.py
if [ ! -f "$FILTERCONSISTENCY" ]; then
    echo "Error: Script not found at $FILTERCONSISTENCY"
    exit 1
fi

## Check consistency-to-bed.py
if [ ! -f "$CONSISTENCYTOBED" ]; then
    echo "Error: Script not found at $CONSISTENCYTOBED"
    exit 1
fi

## Check vcftools
if ! command -v vcftools &> /dev/null; then
    echo "Error: vcftools is required but not found in PATH"
    echo "Please install vcftools (e.g., 'brew install vcftools' on macOS)"
    exit 1
fi

# Check if input file exists
if [ ! -f "${INPUTDIR}/${FILENAMEPREFIX}_consistency_${MODE}.txt" ]; then
    echo "Warning: Input file ${FILENAMEPREFIX}_consistency_${MODE}.txt not found, skipping"
    exit 1
fi

# Check if output directory exists, create if it doesn't
if [ -d "$OUTPUTDIR" ]; then
    echo "Warning: Output directory already exists: $OUTPUTDIR"
else
    echo "Creating output directory: $OUTPUTDIR"
    mkdir -p "$OUTPUTDIR"
fi

# Filter snpeff-consistency by mode and effect 
python $FILTERCONSISTENCY ${INPUTDIR}/${FILENAMEPREFIX}_consistency_${MODE}.txt \
                          -o ${OUTPUTDIR}/${FILENAMEPREFIX}_consistency_${MODE}_fltreffects.txt \
                          --mode ${MODE} \
                          --effects $EFFECTS

echo "Done filtering..."

# Convert filtered consistency table to BED format
python $CONSISTENCYTOBED -o ${OUTPUTDIR}/${FILENAMEPREFIX}_consistency_${MODE}_fltreffects.bed \
                          ${OUTPUTDIR}/${FILENAMEPREFIX}_consistency_${MODE}_fltreffects.txt

echo "Done converting..."

# Filter in the initial vcf file for SNPs in the BED file
vcftools --vcf ${VCFDIR}/${FILENAMEPREFIX}.vcf \
         --bed ${OUTPUTDIR}/${FILENAMEPREFIX}_consistency_${MODE}_fltreffects.bed \
         --recode --recode-INFO-all \
         --out ${OUTPUTDIR}/${FILENAMEPREFIX}_consistency_${MODE}_fltreffects

echo "Done with vcftools..."

# Rename the vcf file to remove 'encode' from the filename
mv ${OUTPUTDIR}/${FILENAMEPREFIX}_consistency_${MODE}_fltreffects.recode.vcf ${OUTPUTDIR}/${FILENAMEPREFIX}_consistency_${MODE}_fltreffects.vcf

echo "Processing finished!"