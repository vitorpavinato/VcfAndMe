# Extra annotations:
This script helps you process additional annotations to a table file. By doing it, you can easely merge the additional annotation to an already produced `.TSV` with `vcf_to_tsv.py`. To keep the things consistent, it only takes the first item of the SNPEff annotation. This is the same behaviour of the annotation pipeline, where `simplify_snpeff.py` takes the first item, and in the `vcf_to_tsv.py`, where it takes the first item of an optional SIFT4 annotation

You should run SNPEff as below, with the flag `-interval` where you parse a `.BED` file and with the flag `-noGenome` to avoid using database-based annotations. You still need to parse the database `-v` is a required flag.

Here is an example, where I created a `.BED` file from UCSC containing only exon-intron junctions:
```zsh
java -jar snpEff.jar ann -classic -noStats -noGenome -interval dm6_eijunctions.bed -v Drosophila_melanogaster rooted_lifted_filtered.vcf > snpeff/extra_annotations/exon-introns_junctions/rooted_lifted_filtered_eijann.vcf
```

By running like this, SNPEff you annotate any SNPs that overlaps any interval in the `.BED` file as: `EFF=CUSTOM[dm6_eijunctions](MODIFIER...`. See that with brackets, SNPEff adds the name of the file used with the `-interval` flag. If you run `vcf_ann_to_table.py` parsing only the input and ouput file, the value inside the brackets will be the annotation in the output file. You can modify it by parsing any string in `-n` (`annotation_name`).

```zsh
usage: python vcf_ann_to_table.py [-h] -i INPUTFILE [-o OUTPUTFILE] [-n ANNOTATION_NAME]

options:
  -h, --help          show this help message and exit
  -i INPUTFILE        Input vcf file name (default: None)
  -o OUTPUTFILE       Output tsv file name (default: None)
  -n ANNOTATION_NAME  New custom annotation name (default: None)
```
