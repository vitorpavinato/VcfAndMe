#!/bin/bash

# Bash script to run the snpeff-consistency script.
# This is a suggestions how to run all modes at once.
# If you have one vcf for each chromosome, you can loop through them
# or run them in parallel.

# Stop on errors
set -e

# Parameters
FILENAME="test_ann_dataset"
SNPEFFCONSISTENCY=snpeff-consistency.py

# Check if script exists
if [ ! -f "$SNPEFFCONSISTENCY" ]; then
    echo "Error: Script not found at $SNPEFFCONSISTENCY"
    exit 1
fi


THRESHOLD=0.75
DISTANCE=100


# Check if input file exists
if [ ! -f "testdata/${FILENAME}.vcf" ]; then
    echo "Warning: Input file ${FILENAME}.vcf not found, skipping"
    continue
fi

python $SNPEFFCONSISTENCY -o testdata/${FILENAME}_consistency_strict --mode strict --stats testdata/${FILENAME}.vcf

python $SNPEFFCONSISTENCY -o testdata/${FILENAME}_consistency_rule --mode rule -t $THRESHOLD -d $DISTANCE --stats testdata/${FILENAME}.vcf

python $SNPEFFCONSISTENCY -o testdata/${FILENAME}_consistency_specific --mode specific --stats testdata/${FILENAME}.vcf

    

echo "Processing complete!"