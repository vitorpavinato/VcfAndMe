# Annotate a vcf

Here we have a pipeline that takes the file produced by [remake_vcf](https://github.com/vitorpavinato/VcfAndMe/tree/main/remake_vcf) pipeline and annotate it with [SNPEff](https://pcingola.github.io/SnpEff/) and [SIFT4g](https://sift.bii.a-star.edu.sg/sift4g/). 

Here is how you can run the pipeline in the script `pipeline_to_annotate_vcf.py`.
(Note: VCF files from `remake_vcf` should be filter out for >2 alleles SNPs)

```zsh
python pipeline_to_annotate_vcf.py -i examples/biallelic/example_remade_rooted_lifted_filtered.vcf -d Drosophila_melanogaster -b examples/intervals/short_introns.bed -o examples/snpeff -s /Users/tur92196/local/sift4g/BDGP6.83 -f examples/sift4g
```

To run above pipeline, you should have installed SNPEff and SIFT4g. SIFT4 is optional if you don't provide a PATH to a database. But if you provided, you should provid the PATH for an output folder. Another optional argument is the PATH for the file containing intervals you want to include in SNPEff annotation. It should be a BED-like file containing somehow "custom" annotations. The pipeline looks for the presence of a file PATH and when triggered, it implements SNPEff `-interval` argument.

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

An attentive reader should have noticed that to run the pipeline above you need to provide in the command-line a file containing intervals of short-introns. You can provide any BED-like file with custom annotations, it doesn't need to be short-introns. For completeness we provide here a script that we developed for creating a BED file containing short-introns intervals. It can be used as a templete for any interval present in a BED file derived from a GFF you want to retain or filter out. Here is the command:
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

The pipeline first annotates the input VCF with SNPEff and produce a VCF for each chromosome named like this `*_remade_rooted_lifted_filtered_ann.vc`. If no interval file is passed, the pipeline only runs SNPEff using the *Drosophila melanogaster* genome build 6 derived database. \[OPTIONAL\]  If an interval file is provided, it takes the annotated VCFs and run SNPEff again with an interval file, but do not use the *D. melanogaster* database. It also creates annother annotated VCF with the `_simplified.vcf` suffix. This VCF contains the same variants, but only the first effect of each variant annotated is taken. The "simplified" VCF files can be annotated with SIFT4g \[OPTIONAL\]. The simplified VCFs annotated with SIFT4g can also be converted into a table \[OPTIONAL\].

Here is command used to annotate with SNPEff:
```zsh
java -jar snpeff ann -c snpeff_config -classic
          -noStats -csvStats csv_stats_file -v database
          input_file > output_file
        
```

\[OPTIONAL\] After running the pipeline, you can convert the annotated VCF to a TSV table:
```zsh
python vcf_to_tsv.py -i examples/sift4g/example_remade_rooted_lifted_filtered_ann_simplified_SIFTpredictions.vcf -o examples/tables/example_remade_rooted_lifted_filtered_ann_table_snpeff_sift4g.vcf -r PATH/TO/REFERENCE -s PATH/TO/SAMTOOLS -f 3 -c short_introns.bed -n SI -e
```

Here is the complete list of arguments:
```zsh
usage: python vcf_to_tsv.py [-h] -i INPUTFILE [-o OUTPUTFILE] -r REFERENCE -s SAMTOOLS_PATH [-f NFLANKINBPS] [-c CUSTOM_EFFECT_NAME] [-n NEW_CUSTOM_EFFECT_NAME] [-e]

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
- `simplify_snpeff.py` only takes the first term from SNPEff annotation. For bi-allelic SNPs, the first term should be the important one (with the annotation of the most impacted change); but for tri-allelic there are the issue with the most important annotation per allele, and which allele is the most (or least) common one. Right now, the script picks the first term, regardless of all these issues for tri-allelic(s). ADVICE: remove tri-allelics before running this pipeline. In [annotate_vcf/snpeff_consistency](https://github.com/vitorpavinato/VcfAndMe/tree/main/annotate_vcf/snpeff_consistency) there is framework to check and obtain SNP annotations based on different consisteny rules (strictly one term, a dominant term etc) that can be modified to deal with different annotations for each of the different alleles for a tri-allelic SNP. I didnd't tried it yet. If you found this framework useful for this, fell free to create a pull request with these additions.
- `vcf_to_tsv.py` also had to deal with tri-allelic SNPs in SIFT4g annotations. I did a quick fix to baypass it. It basically disregard the second annotated term, discarding the functional annotation of the second alternative allele. Because of it, the logic would be to remove tri-allelic SNPs (this was removed to force the use of biallelic-only vcf).
- `get_reversed_complementary_strand()` also doesn't know how to deal with tri-allelics, so keep only biallelic SNPs.
