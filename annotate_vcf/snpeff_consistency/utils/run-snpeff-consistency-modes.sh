#!/bin/bash

# Bash script to run the snpeff-consistency program.
# This is a suggestion in how to run all modes at once.
# If you have one vcf for each chromosome, you can loop through them
# or run them in parallel, by running this script in diffrent terminal
# windows (or using GNU parallel).

# Stop on errors
set -e

# Parameters
SNPEFFCONSISTENCY="../snpeff-consistency.py"
FILENAMEPREFIX="test_ann_dataset"
INPUTDIR="../testdata"
OUTPUTDIR="../testdata/snpeff-consistency"
THRESHOLD=0.75
DISTANCE=100

# Function to check if output files exist for a specific mode
check_output_files() {
    local mode=$1
    local base_output="${OUTPUTDIR}/${FILENAMEPREFIX}_consistency_${mode}"
    
    if [ -f "${base_output}.txt" ] || [ -f "${base_output}_stats.txt" ]; then
        echo "Error: Output files already exist for ${mode} mode:"
        [ -f "${base_output}.txt" ] && echo "- ${base_output}.txt"
        [ -f "${base_output}_stats.txt" ] && echo "- ${base_output}_stats.txt"
        exit 1
    fi
}

# Check if script exists
if [ ! -f "$SNPEFFCONSISTENCY" ]; then
    echo "Error: Script not found at $SNPEFFCONSISTENCY"
    exit 1
fi

# Check if input file exists
if [ ! -f "${INPUTDIR}/${FILENAMEPREFIX}.vcf" ]; then
    echo "Warning: Input file ${FILENAMEPREFIX}.vcf not found, skipping"
    exit 1
fi

# Check if output directory exists, create if it doesn't
if [ -d "$OUTPUTDIR" ]; then
    echo "Warning: Output directory already exists: $OUTPUTDIR"
else
    echo "Creating output directory: $OUTPUTDIR"
    mkdir -p "$OUTPUTDIR"
fi

# Check for existing output files before running each mode
check_output_files "strict"
check_output_files "rule"
check_output_files "specific"

# Program execution for each of the available modes
echo "Running the strict mode..."
python $SNPEFFCONSISTENCY -o ${OUTPUTDIR}/${FILENAMEPREFIX}_consistency_strict --mode strict --stats --codon_stats ${INPUTDIR}/${FILENAMEPREFIX}.vcf
echo "Strict mode finished."

echo "Running the rule mode with threshold ${THRESHOLD} and distance ${DISTANCE}..."
python $SNPEFFCONSISTENCY -o ${OUTPUTDIR}/${FILENAMEPREFIX}_consistency_rule --mode rule -t $THRESHOLD -d $DISTANCE --stats --codon_stats ${INPUTDIR}/${FILENAMEPREFIX}.vcf
echo "Rule mode finished."

echo "Running the specific mode..."
python $SNPEFFCONSISTENCY -o ${OUTPUTDIR}/${FILENAMEPREFIX}_consistency_specific --mode specific --stats --codon_stats ${INPUTDIR}/${FILENAMEPREFIX}.vcf
echo "Specific mode finished."
    
echo "Processing complete!"