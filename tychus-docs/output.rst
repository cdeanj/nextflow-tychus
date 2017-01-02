Output
======

The output structure from each Tychus module is described below.

Alignment Module
----------------

Quality Filtered Read Pairs
```````````````````````````

**Description**: User-input reads are filtered and trimmed with Trimmomatic

**Directory**: PreProcessing/

**File Names**: {dataset_id}_1P.fastq, {dataset_id}_2P.fastq

.. code-block:: console
   :linenos:

   SFBRL043-M3237-14-001_S7_L001_1P.fastq
   SFBRL043-M3237-14-001_S7_L001_2P.fastq


Quality filtered reads from each pair of FASTQ inputs are provided as output in the **PreProcessing/** directory. Each pair of FASTQ pairs is prefixed with the ``dataset ID`` from which it came from and suffixed with a ``1P`` and ``2P`` extension for each forward and reverse read in the pair.

Alignment Files
```````````````

FASTQ inputs that pass quality filtering are then aligned to four reference databases: a reference genome, an antimicrobial resistance database (AMR), a virulence factor database, and a plasmid database. Outputs produced from each of these steps are the corresponding BAM formatted alignment and index files. Each file is prefixed with the ``dataset ID`` from which it came from and suffixed with a ``.bam`` and ``.bai`` file extension. These files can be found in the **GenomeAlignment/**, **AMRAlignment/**, **VirulenceAlignment/**, and **PlasmidAlignment/** directories. For more information about the Sequence Alignment/Map Format, please see this helpful `page <https://samtools.github.io/hts-specs/SAMv1.pdf>`_.

Consensus Files
```````````````

Each ``BAM`` file produced from the alignment against the reference genome in the previous step is then provided as input to `Freebayes <https://github.com/ekg/freebayes>`_, a haplotype-based variant caller. The variants produced are then incorporated into the reference genmone, creating a FASTA formatted consensus sequence with `Bcftools <https://samtools.github.io/bcftools/bcftools.html>`_. Each consensus sequence is prefixed with the ``dataset ID`` from the previous step and suffixed with a ``.fa`` extension. These files can be found in the **Consensus/** directory.

SNPs and Phylogenies
````````````````````

The FASTA formatted consensus sequences and reference genome are then provided as input to a program called `kSNP3 <https://sourceforge.net/projects/ksnp/>`_, which identifies SNPs and builds phylogenies. The SNPs produced are provided as output in a variety of formats: (1) SNP allele fasta alignment, (2) SNP matrix, and (3) SNPs with allele positions corresponding to where they were found in the reference genome [citation]. `Newick <https://en.wikipedia.org/wiki/Newick_format>`_ formatted phylogenies are also produced using a variety of different algorithms including by: (1) parsimony, (2) maximum likelihood, and (3) neighbor joining methods.

Annotated Phylogenies
`````````````````````

Each Newick formatted phylogeny is then provided as input to `Figtree <http://tree.bio.ed.ac.uk/software/figtree/>`_, a program for visualizing text-based phylogenies. The annotated phylogenies output from Figtree will come in a format specified by you, the user (JPEG|PDF|PNG|SVG). Each figure produced is prefixed with the basename of the file used to build the phylogeny and suffixed with the user-specified image format. These files can be found in the **PhylogeneticTreeImages/** directory.
