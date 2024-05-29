# How to annotate a vcf and convert it to a tsv table

Here I show how to run the annotation pipeline on a vcf file (that was filtered out for tri-allelic SNPs).
```zsh
python pipeline_to_annotate_vcf.py -i examples/biallelic/example_remade_rooted_lifted_filtered.vcf -d Drosophila_melanogaster -b examples/intervals/short_introns.bed -o examples/snpeff -s /Users/tur92196/local/sift4g/BDGP6.83 -f examples/sift4g
```

To run this pipeline like above, you should have installed [SNPEff](https://pcingola.github.io/SnpEff/) and [SIFT4g](https://sift.bii.a-star.edu.sg/sift4g/). SIFT4 is optional if you don't provide a PATH to a database. But if you provided, you should provid the PATH for an output folder. Another optional argument is the PATH for the file containing intervals you want to include in SNPEff annotation. It should be a BED-like file containing somehow "custom" annotations. The pipeline looks for the presence of a file PATH and when triggered, it implements SNPEff `-interval` argument.

Make sure to change the necessary SNPEff and SIFT4g PATHs in the `annotate_vcf/config.ini` file.

The other parameters are self explained when you have the pipeline help message:
```zsh
usage: python pipeline_to_annotate_vcf.py [-h] -i INPUT_FILE -d SNPEFF_DATABASE [-b SNPEFF_INTERVAL_FILE] -o SNPEFF_OUTPUT_FOLDER [-s SIFT4G_DATABASE]
                                          [-f SIFT4G_OUTPUT_FOLDER]

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE         input vcf file not annotated (default: None)
  -d SNPEFF_DATABASE    SNPEff database (default: None)
  -b SNPEFF_INTERVAL_FILE
                        SNPEff interval BED file (default: None)
  -o SNPEFF_OUTPUT_FOLDER
                        output folder for SNPEff annotations (default: None)
  -s SIFT4G_DATABASE    SIFT4G database (default: None)
  -f SIFT4G_OUTPUT_FOLDER
                        Output folder for SIFT4G annotations (default: None)
```

An attentive reader should have noticed that to run the pipeline above I provided in the command-line a file containing intervals of short-introns. You can provide any BED-like file with custom annotations, it doesn't need to be short-introns. Becaue I included an example with short-introns, for completeness I provide here a script that I developed to create the BED file containing short-introns. It can be used as a templete for any interval present in a BED file derived from a GFF you want to retain or filter out. Here is the command:
```zsh
python get_short_introns_from_bed.py -i  examples/intervals/introns.bed -o examples/intervals/short_introns.bed -c chr2L -s 86 -t 8
```

Here is the complete list of arguments:
```zsh
usage: python get_short_introns_from_bed.py [-h] -i INPUTFILE -o OUTPUTFILE -c CHROM_LIST [CHROM_LIST ...] [-s SHORT_INTRON_SIZE] [-t TRAILLING_SIZE]

options:
  -h, --help            show this help message and exit
  -i INPUTFILE          Input BED file name (default: None)
  -o OUTPUTFILE         Output BED file name (default: None)
  -c CHROM_LIST [CHROM_LIST ...]
                        List of chroms to process (default: None)
  -s SHORT_INTRON_SIZE  Short intron size (default: 86)
  -t TRAILLING_SIZE     Trailling size (default: 8)
```

After running the pipeline, you can convert the annotated vcf to a tsv table:
```zsh
python vcf_to_tsv.py -i examples/sift4g/example_remade_rooted_lifted_filtered_ann_simplified_SIFTpredictions.vcf -o examples/tables/example_remade_rooted_lifted_filtered_ann_table_snpeff_sift4g.vcf -r PATH/TO/REFERENCE -s PATH/TO/SAMTOOLS -f 3 -c short_introns.bed -n SI -e
```

Here is the complete list of arguments:
```zsh
usage: python vcf_to_table.py [-h] -i INPUTFILE [-o OUTPUTFILE] -r REFERENCE -s SAMTOOLS_PATH [-f NFLANKINBPS] [-c CUSTOM_EFFECT_NAME] [-n NEW_CUSTOM_EFFECT_NAME] [-e]

options:
  -h, --help            show this help message and exit
  -i INPUTFILE          Input vcf file name (default: None)
  -o OUTPUTFILE         Output tsv file name (default: None)
  -r REFERENCE          Path to the reference genome of the vcf file (default: None)
  -s SAMTOOLS_PATH      Path to the samtools (default: None)
  -f NFLANKINBPS        Number of bases flanking each targeted SNP (default: 3)
  -c CUSTOM_EFFECT_NAME
                        Custom effect name (default: None)
  -n NEW_CUSTOM_EFFECT_NAME
                        New custom effect name (default: None)
  -e                    Input vcf with SIFT4G annotations (default: False)
```

### Issues:
- `simplify_snpeff.py` only takes the first term from SNPEff annotation. For bi-allelic SNPs, the first term should be the important one (with the annotation of the most impacted change); but for tri-allelic there are the issue with the most important annotation per allele, and which allele is the most (or least) common one. Right now, the script picks the first term, regardless of all these issues for tri-allelic(s). ADVICE: remove tri-allelics before running this pipeline.
- `vcf_to_tsv.py` also had to deal with tri-allelic SNPs in SIFT4g annotations. I did a quick fix to baypass it. It basica disregard the second annotated term, basically discarding the functional annotation of the second alternative allele. Because of it, the logic would be to remove tri-allelic SNPs.
- `get_reversed_complementary_strand()` also doesn't know how to deal with tri-allelics.