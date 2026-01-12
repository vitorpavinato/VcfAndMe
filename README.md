## VcfAndMe

Here, we have some scripts that we used to convert [DGN](https://www.johnpool.net/genomes.html) whole-genome sequence aligmnets to VCFs. The workflow presented in `remake_vcf` was motivated by the fact that *Drosophila melanogaster* data on DGN weren't ready to use as polymorhism-containing files (aka VCFs). The below workflow shows how we converted multiple sequence alignment (MSA) files to VCFs, and how we rooted the variants in those VCFs by parsimony – using "ancestral states" of alleles found in th e*Drosophila simulans* genome.


In `annotate_vcf` we provide some tools for variant annotation. We developed a pipeline for annotating variants using [SNPEff](https://pcingola.github.io/SnpEff/) and [SIFT4g](https://sift.bii.a-star.edu.sg/sift4g/), and a pipeline to make sense of the multiple effects SNPEff determines to a variant (`snpeff_consistency`). In addition, it has a script to convert a VCF file into a TSV table (be aware that the variant effect from the VCF that goes into the table is the first effect — some work should be done to use the results from `snpeff_consistency` in the script that convert the VCF into a table). 


Table of Content:

- [Remake vcf files](https://github.com/vitorpavinato/VcfAndMe/tree/main/remake_vcf)
- [Annotate SNPs](https://github.com/vitorpavinato/VcfAndMe/tree/main/annotate_vcf)
    - [snpeff_consistency](https://github.com/vitorpavinato/VcfAndMe/tree/main/annotate_vcf/snpeff_consistency)