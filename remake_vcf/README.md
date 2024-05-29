# Remake vcf file

Some ideas to make DGN data usable for vcf-based tools.

I recommend you structure your working directory as something like this:
```bash
.    # My root folder is named DGN (Here I also have data from other populations)
├── dpgp3
│   ├── dmel_data # This is a copy of the GitHub repository
│   ├── dpgp3.txt
│   ├── masked
│   │   ├── Chr2L
│   │   ├── Chr2R
│   │   ├── Chr3L
│   │   ├── Chr3R
│   │   ├── fastas
│   │   └── vcfs
│   │       ├── liftedover
│   │       ├── remade
│   │       └── rooted
│   └── originals # Is where the original data from DGN are stored
├── data
├── reference # Reference genome and associated files
├── packages # Where I placed the masking packages
└── simulans_sequences # This is for rooting (optional showed below)
```

### Download DGN data
First download the consensus sequence FASTA from [DGN](https://www.johnpool.net/genomes.html) for the target population (sample). Go to the bottom of the page where you see a bunch of links like \<name\>_SEQ. These are the links to a compressed folder that contains `.seq` files with consensus FASTA sequences for each individual organized for each chromosome, for each collection of samples. Here I am showing as an example the analysis of data from DPGP3 that contains a collection of 197 Zambia genomess

```zsh
# From the root folder of your project:
mkdir data
cd data
wget http://pooldata.genetics.wisc.edu/dpgp3_sequences.tar.bz2

md5sum-lite dpgp3_sequences.tar.bz2 906d282740a56e5273a4cbc5abfb61f9
```

After checking the md5sum and decompressing the file, you should place the folder containing the chromosomes FASTAs inside the folder I named `originals`.

```zsh
cd ../
mkdir -p dpgp3/originals
mv data/dpgp3_sequences/dpgp3_Chr* dpgp3/originals
```

### Run DGN masking scripts
Second, download the [masking package](http://johnpool.net/masking.zip) also from DGN. It contains a set of files and scripts used to mask problematic sites previously identified containing traces of identity-by-descent (IBD) or admixture. Copy both scripts and the interval files (two `.csv` files) inside each chromosome folder (that should now be inside `originals`), and run the two `Perl`` scripts. Here I am showing how to run them inside chromosome 2L folder (I am showing only the analysis of Chr2L, for the rest of this document).

You can either delete the compressed files or keep them. You can create a folder in the root of your project called `packages` to save the scrips, CSVs and the original `.zip` files in there.

Download the masking package:
```zsh
mkdir packages
cd packages
wget http://johnpool.net/masking.zip
```

Then copy these files to each of the folder containing sequences for each chromosome (I am oly showing for the chromosome 2L here).

```zsh
cd ..
cp -b packages/masking/*.pl dpgp3/originals/Chr2L
cp -b packages/masking/*.csv dpgp3/originals/Chr2L
cd dpgp3/originals/Chr2L

# Run the IBD masking script first:
perl ibd_mask_seq.pl

# Wait until it finish, then run the Admixture masking script:
perl admixture_mask_seq.pl
```

You need to repeat it for all folders (I know, it is tedious, I stopped at chr3R after running the other three major chromosomes - maybe you can open four terminal windows and run then one script at time, but in four different chromosomes...)

If the individual has IBD or Admixture traces, an `N` will replace the nucleotide on the problematic position. For each chromosome, the result is a set of individuals genomes (masked or not). The unmasked files will be moved to the folder corresponding to the script you ran. After running the two script you should have the masked files and two folders: `unmasked_admixture` and `unmasked_ibd` and the masked files, one for each individual. Move the masked files to the corresponding folder named with the chromosome name inside the folder `masked`.

```zsh
mv *_Chr2L.seq ../../masked/Chr2L/
cd ../..masked/Chr2L
```

### Add sample names to each sequence
These `.seq` files don't have header (I don't know why?). Copy the name of each individual and place it as a header (using bash please with something like ```ls -1 *.seq | rev | cut -c11- | rev > ../../dpgp3.txt```). You can keep this files at the `dpgp3` directory of your project, as it is shown in the above directory structure. I prepared a script that takes the names of each individual placed in a one-column `.txt` file (as the one created above), opens the corresponding `.seq`, and adds the corresponding name of the sample as the header:

```zsh
# You should run it for each chromosome file. Inside the folder for the chromosome, type:
python ../../dmel_data/remake_vcf/add_head_to_fasta_files.py -s _Chr2L.seq -n ../../dpgp3.txt
```
Where `-s` means suffix and `-n` means the file with the sample names.

### Run snp-site
Then concatenate all `.seq` file by chromosome using `cat`. I did something like this. Inside the folder with all sequences from the chromosome 2L, type:

```zsh
mkdir ../fastas
cat *.seq > ../fastas/ZI_Chr2L.fasta
cd ../fastas
```

Then, run [snp-site](https://sanger-pathogens.github.io/snp-sites/) to convert the *multi-fasta alignment* to a somehow descent and usable `.vcf` file. Here is how the command looks like (make sure that `snp-sites` is in your PATH, either absolute or when a conda environment is active):

```zsh
# Run snp-sites in each chromosome alignment
snp-sites -v -o ZI_Chr2L.vcf ZI_Chr2L.fasta

# Move the VCFs to another folder
mkdir ../vcfs
mv ZI_Chr2L.vcf ../vcfs
```

A note here: If you are trying to install snp-site from conda on a Apple computer with Mx chips, run this command before (with the target conda environemtn activated):
```zsh
conda config --env --set subdir osx-64
```

This will create a `.vcf` from the `.fasta` file. Here is important to have all sequences with their corresponding header. Usually it is the sample name. The VCF will have then the same ID as the FASTA for each sample. 

Snp-site has no way to know which base in the multi-fasta alignment is in the reference base, so I prepared a script to deal with this (see [remake_vcf](https://github.com/vitorpavinato/dmel_data/tree/main/remake_vcf) folder). The script also identify which genotype number corresponds to the missing character `*` and change the number to `./.`, effectively making the masking done before usable. And it re-sort the order of alternative alleles and make genotype codes accordingly to the re-ordering. The most important ordering is the true reference and the alternatives, so the allele that is the reference in the genome (dm3 or the flybase release 5) will be `0` and reference genotypes will be `0/0`. 

Here is how you run `remake_vcf`. There are two examples below: 1- for the entire 2L vcf, another with a sample of SNP from chromosome 2L (called example). For the rest of this tutorial I will be using this vcf.

```zsh
cd ..
mkdir -p vcfs/remade # to store the fixed vcf files
cd vcfs

# An example on how to run on a vcf derived from the aligments of chromosome 2L samples
python ../../dmel_data/remake_vcf/remake_vcf.py -i ZI_Chr2L.vcf -r ../../../reference/dm3.fa -c chr2L -l 23011544 -o remade/ZI_Chr2L_remade.vcf -s /Users/tur92196/local/samtools1_8/bin/samtools

# To run the examples provided with the GitHub repository, you need to move insed the copy of the repository you downloaded.
cd ../../dmel_data/remake_vcf

python remake_vcf.py -i example/example.vcf -r ../../../reference/dm3.fa -c chr2L -l 23011544 -o example/example_remade.vcf -s /Users/tur92196/local/samtools1_8/bin/samtools
```

Where `-i`is the imput vcf you whant to fix, `-r` is the PATH for the reference genome the original DGN data was aligned to (dm3 also named as dmelr5), the `-c` chromosome name; `-l` is the length of the chromosome (you can put any number here, it is just to annoate the vcf with it); `-o`is the output and `-s` is the PATH for samtools as it is required (later I will probably remove this from the CLI interface, so make sure samtools is in your PATH).

I advise to use the reference genomes (and annotations) stored at [UCSC genome browse](https://genome.ucsc.edu/cgi-bin/hgGateway) because the chromosomes are named with the same naming system as in the chain file. The chain files will be used in later steps to liftover the SNPs from dm3 (dmel r.5.x) to dm6 (dmel r.6.x).

### Rooting or annotating the SNPs with the ancestral state
Now that we have a `better DGN vcf files`, we can determine the ancestral state of each SNP using a simple parsimony strategy. It is optional. If you don't need to know the AAs, you can jump to the next section. The idea here is to use the reference genome sequencing of one of *D. melanogaster* sister species to find the state before the two species splited from a common ancestor. I use *D. simulans* because, conviently, there is a version of this species genome that was aligned to *D. melanogaster* genome. This whole genome aligment underwent some manipulation to make sure the syntenic regions between the two species was present along with masked regions, in a way that the new genome of *D. simulans* has the same number of bp as in *D. melanogaster* genome. This make any sequence comparison straigthford. This alignment was done on dm3, so were are good to go.

```zsh
python root_vcf_by_parsimony.py -i example/example_remade.vcf -o example/example_remade_rooted.vcf -r ../../../simulans_sequences/dsim2_as_dmelr5.fasta -s /Users/tur92196/local/samtools1_8/bin/samtools
```

### Liftover VCF
Now is time to lift the positions over to the newest (at least the newest at UCSC genome browser) Dmel genome (release 6). We are going to use GATK `Picard` [LiftoverVcf](https://gatk.broadinstitute.org/hc/en-us/articles/360037060932-LiftoverVcf-Picard). The input VCF here can be the one from above (rooted) or the one not rooted.

First create variables to store the PATH of GATK and the JAVA required to run the program.

```zsh
JAVA='/Library/Java/JavaVirtualMachines/temurin-17.jdk/Contents/Home/bin/java'
GATK='/Users/tur92196/local/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar'
```

To run GATK `LiftoverVcf`, you should have:
- your `.vcf` you want to lift to the newest genome;
- the reference sequence of the target genome (the one you are lifting the SNPs to);
- the associated `.dict` file (see how to get one below);
- and the chain files, downloaded from UCSC (found [here](https://hgdownload.soe.ucsc.edu/goldenPath/dm3/liftOver/)), containing the chains from the oldest (here dm3) to the newest target genome (dm6).

Go inside the folder you have the genome file of the target genome. There we will run, to obtain the `.dict` file. You can find the [documentation here](https://gatk.broadinstitute.org/hc/en-us/articles/360037422891-CreateSequenceDictionary-Picard-):

```bash
cd ../../../reference
$JAVA -jar $GATK CreateSequenceDictionary \\ 
    -R dm6.fa \\ 
    -O dm6.fa.dict
```

With everything in place, then type, to have lift SNPs from one genome to the target genome:
```bash
cd ../dpgp3/dmel_data/remake_vcf/

$JAVA -jar $GATK LiftoverVcf \\
-I example/example_remade_rooted.vcf \\
-O example/example_remade_rooted_lifted.vcf \\
-C ../../../reference/dm3ToDm6.over.chain \\
--REJECT example/example_remade_rooted_rejected.vcf \\
-R ../../../reference/dm6.fa \\
--WARN_ON_MISSING_CONTIG true \\
--RECOVER_SWAPPED_REF_ALT true &
```

### Keep only bi-allelic SNPs

This is mandatory step: the [annotation pipeline](https://github.com/vitorpavinato/dmel_data/tree/main/annotate_vcf) only deals with bi-allelic SNPs.

Now that everthing is in place, we are going to filter out SNPs with more than one alternative allele to retain only bi-allelic. You can also filter out a lot of missing data (more than 50%). We are going to use `vcftools` so make sure it is in your PATH for this task (this is optional depending on your project context):

```zsh
vcftools --vcf example/example_remade_rooted_lifted.vcf \\
--out example/example_remade_rooted_lifted_fltr \\
--min-alleles 2 \\
--max-alleles 2 \\
--max-missing 0.5 \\
--recode \\
--recode-INFO-all
```

Remove the `.recode.` part of the file name to kept file naming simple.
```zsh
mv example/example_remade_rooted_lifted_fltr.recode.vcf example/example_remade_rooted_lifted_fltr.vcf
```