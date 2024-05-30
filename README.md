## VcfAndMe

Here I have developed some tools to work with `vcf` files. I started working on the ideas presented here in a project I need to deal with polymorphism data of two *Drosophila melanogaster* populations. The tools present in `remake_vcf` were motivated by the fact that Dmel data on [DGN](https://www.johnpool.net/genomes.html), besides being a great database that makes available a ton of population genomics data, all of that aligned and process with the same bioinformatics tools, weren't ready to use as vcfs files. I provide the steps and some tools to have these vcfs from multiple sequence alignment (MSA) files. I also provide a tool for rooting a vcf, by parsimony, using the alignment between Dmel and *Drosophila simulans*. 

In `annotate_vcf` I tried to go a step further providing tools and ideas for variant annotation. Basically I developed a pipeline for annotation using [SNPEff](https://pcingola.github.io/SnpEff/) and [SIFT4g](https://sift.bii.a-star.edu.sg/sift4g/), and a tool to convert an annotated vcf to a table.

Table of Content:

- [Remake vcf files](#Remake-vcf-files)
- [Annotate SNPs](#Annotate-SNPs)


### Remake vcf files
Follow the steps on the `remake_vcf` page to make DGN files usable (might seems restricted to Dmel world, but if you have individual genomes in MSA files you should be able to use most of the tools).

### Annotate SNPs
Tools for variant annotation and vcf to tsv conversion (should work with any species)