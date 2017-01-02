Output
======

The output structure from each Tychus module is described below.

Alignment Module
----------------

Quality Filtered Read Pairs
``````````````````````

Quality filtered reads from each pair of FASTQ inputs are provided as output in the *PreProcessing* directory. Each pair of FASTQ pairs is prefixed with the dataset ID from which it came from and suffixed with a *1P* and *2P* extension for each forward and reverse read in the pair.

Alignment Files
```````````````

FASTQ inputs that pass quality filtering are then aligned to four reference databases: a reference genome, an antimicrobial resistance database (AMR), a virulence factor database, and a plasmid database. Outputs produced from each of these steps are the corresponding BAM formatted alignment and index files. Each file is prefixed with the dataset ID from which it came from and suffixed with a *.bam* and *.bai* file extension. These files can be found in the *GenomeAlignment*, *AMRAlignment*, *VirulenceAlignment*, and *PlasmidAlignment* directories. For more information about the contents and structure of the alignment file, please see this helpful wiki.

Consensus Files
```````````````

Each *BAM* file produced from alignment against the reference genome in the previous step is then provided as input to *Freebayes*, a haplotype-based variant caller. The variants produced are then used to create a FASTA formatted consensus sequence with *Bcftools*. Each consensus sequence is prefixed with the dataset ID from the previous step and suffixed with a *.fa* extension. These files can be found in the *Consensus* directory.

SNPs and Phylogenies
````````````````````

The FASTA formatted consensus sequences and reference genome are then provided as input to a program called kSNP3, which is able to identify SNPs and build phylogenies. The SNPs produced are output in a variety of formats: (1) SNP allele fasta alignment, (2) SNP matrix, and (3) SNPs with allele positions corresponding to where they were found in the genome [citation]. Newick formatted phylogenies are also produced. The phylogenies are built using a variety of different algorithms including: (1) parsimony, (2) maximum likelihood, and (3) neighbor joining methods.

Annotated Phylogenies
`````````````````````

Each Newick formatted phylogeny is then provided as input to Figtree, a program for visualizing text-based phylogenies. The annotated phylogenies output from Figtree will come in a format specified by you, the user (JPEG|PDF|PNG|SVG). These files can be found in the *PhylogeneticTreeImages* directory.
