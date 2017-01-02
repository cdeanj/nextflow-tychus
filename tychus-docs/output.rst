Output
======

The output structure from each Tychus module is described below.

Alignment Module
----------------

Quality Filtered Read Pairs
```````````````````````````

* **Description**: User-input reads are filtered and trimmed with Trimmomatic
* **Directory**: PreProcessing/
* **File Names**: {dataset_id}_1P.fastq, {dataset_id}_2P.fastq

.. code-block:: console

   SFBRL043-M3237-14-001_S7_L001_1P.fastq
   SFBRL043-M3237-14-001_S7_L001_2P.fastq


Alignment Files
```````````````

* **Description**: FASTQ inputs that pass quality filtering are then aligned to four reference databases: a reference genome, an antimicrobial resistance database (AMR), a virulence factor database, and a plasmid database. Outputs produced from each of these steps are the corresponding BAM formatted alignment and index files.
* **Directory**: GenomeAlignment/, AMRAlignment/, VirulenceAlignment/, and PlasmidAlignment/.
* **File Names**: {dataset_id}_alignment.bam, {dataset_id}_alignment.bai

.. code-block:: console

   SFBRL043-M3237-14-001_S7_L001_alignment.bam
   SFBRL043-M3237-14-001_S7_L001_alignment.bai


* **More Info**: For more information about the Sequence Alignment/Map Format, please see this helpful `page <https://samtools.github.io/hts-specs/SAMv1.pdf>`_.

Consensus Files
```````````````

* **Description**: Each ``BAM`` file produced from the alignment against the reference genome in the previous step is then provided as input to `Freebayes <https://github.com/ekg/freebayes>`_, a haplotype-based variant caller. The variants produced are then incorporated into the reference genmone, creating a FASTA formatted consensus sequence with `Bcftools <https://samtools.github.io/bcftools/bcftools.html>`_.
* **Directory**: Consensus/
* **File Names**: {dataset_id}_consensus.fa

.. code-block:: console

   SFBRL043-M3237-14-001_S7_L001_consensus.fa


SNPs and Phylogenies
````````````````````

* **Description**: The FASTA formatted consensus sequences and reference genome are then provided as input to a program called `kSNP3 <https://sourceforge.net/projects/ksnp/>`_, which identifies SNPs and constructs `Newick <https://en.wikipedia.org/wiki/Newick_format>`_ formatted phylogenies.
* **Directory**: Polymorphisms/, Trees/
* **File Names**: 


Annotated Phylogenies
`````````````````````

* **Description**: Each Newick formatted phylogeny is then provided as input to `Figtree <http://tree.bio.ed.ac.uk/software/figtree/>`_, a program for visualizing text-based phylogenies. The annotated phylogenies output from Figtree will come in a format specified by you, the user (JPEG|PDF|PNG|SVG). 
* **Directory**: PhylogeneticTreeImages/
* **File Names**: {basename}*.svg

.. code-block:: console

   tree_AlleleCounts.core.NodeLabel.svg		tree_AlleleCounts.core.svg
   tree_AlleleCounts.majority0.75.NodeLabel.svg	tree_AlleleCounts.majority0.75.svg
   tree_AlleleCounts.NJ.NodeLabel.svg		tree_AlleleCounts.NJ.svg
   tree_AlleleCounts.parsimony.NodeLabel.svg	tree_AlleleCounts.parsimony.svg
   tree.core.svg				tree.majority0.75.svg
   tree.NJ.svg					tree.parsimony.svg
   tree_tipAlleleCounts.core.svg		tree_tipAlleleCounts.majority0.75.svg
   tree_tipAlleleCounts.NJ.svg			tree_tipAlleleCounts.parsimony.svg
