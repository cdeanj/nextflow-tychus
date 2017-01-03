Output
======

The output structure and files produced from each Tychus module are described below.

Alignment Module
----------------

Quality Filtering
`````````````````

* **Description**: User-input reads are filtered and trimmed with a program called `Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_.
* **Directory**: PreProcessing/
* **File Names**: {dataset_id}_1P.fastq, {dataset_id}_2P.fastq

.. code-block:: console

   SFBRL043-M3237-14-001_S7_L001_1P.fastq
   SFBRL043-M3237-14-001_S7_L001_2P.fastq


Sequence Alignment
``````````````````

* **Description**: FASTQ inputs that pass quality filtering are then aligned to four reference databases: a reference genome, an antimicrobial resistance database (AMR), a virulence factor database, and a plasmid database. Outputs produced from each of these steps are the corresponding BAM formatted alignment and index files.
* **Directory**: GenomeAlignment/, AMRAlignment/, VirulenceAlignment/, and PlasmidAlignment/.
* **File Names**: {dataset_id}_alignment.bam, {dataset_id}_alignment.bai

.. code-block:: console

   SFBRL043-M3237-14-001_S7_L001_alignment.bam
   SFBRL043-M3237-14-001_S7_L001_alignment.bai


* **More Info**: For more information about the Sequence Alignment/Map Format, please see this helpful `page <https://samtools.github.io/hts-specs/SAMv1.pdf>`_.

Consensus Calling
`````````````````

* **Description**: Each BAM file produced from the alignment against the reference genome in the previous step is then provided as input to `FreeBayes <https://github.com/ekg/freebayes>`_, a haplotype-based variant caller. The variants produced are then incorporated into the reference genome, creating a FASTA formatted consensus sequence with `BCFtools <https://samtools.github.io/bcftools/bcftools.html>`_.
* **Directory**: Consensus/
* **File Names**: {dataset_id}_consensus.fa

.. code-block:: console

   SFBRL043-M3237-14-001_S7_L001_consensus.fa


SNPs and Phylogenies
````````````````````

* **Description**: The FASTA formatted consensus sequences and reference genome are then provided as input to a program called `kSNP3 <https://sourceforge.net/projects/ksnp/>`_, which identifies SNPs and constructs `Newick <https://en.wikipedia.org/wiki/Newick_format>`_ formatted phylogenies.
* **Directory**: Polymorphisms/, Trees/
* **File Names**: The files produced by this process are abundant and will not be covered in detail. These include SNP matrices, SNP counts shared by each produced genome in the ``consensus step``, phylogenetic trees, core SNPs, majority SNPs, andhomoplasy groups.
* **More Info:** Users are encouraged to read ``The Output Files`` section of the kSNP3 `user manual <http://download2.nust.na/pub4/sourceforge/k/ks/ksnp/kSNP3.02%20User%20Guide.pdf>`_ to learn more about the output files produced by kSNP3.


Phylogeny Annotation
````````````````````

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


Assembly Module
----------------

Quality Filtering
`````````````````

* **Description**: User-input reads are filtered and trimmed with Trimmomatic
* **Directory**: PreProcessing/
* **File Names**: {dataset_id}_1P.fastq, {dataset_id}_2P.fastq

.. code-block:: console

   SFBRL043-M3237-14-001_S7_L001_1P.fastq
   SFBRL043-M3237-14-001_S7_L001_2P.fastq

Assembly Contigs
````````````````

* **Description**: FASTQ inputs that pass quality filtering are then used as input to four *de novo* genome assemblers (Abyss, IDBA, SPades, and Velvet), which are used to build genome assemblies from short-read sequence data.
* **Directory**: AbyssContigs/, IDBAContigs/, SPadesContigs/, VelvetContigs/
* **File Names**: {dataset_id}_{assembler-name}-contigs.fa

.. code-block:: console

   SFBRL043-M3237-14-001_S7_L001_abyss-contigs.fa
   SFBRL043-M3237-14-001_S7_L001_idba-contigs.fa
   SFBRL043-M3237-14-001_S7_L001_spades-contigs.fa
   SFBRL043-M3237-14-001_S7_L001_velvet-contigs.fa


Contig Integration
``````````````````

* **Description**: Contigs produced from each of the four genome assemblers are then used as input to a program called `CISA <http://sb.nhri.org.tw/CISA/en/CISA;jsessionid=A1D6759153E59AA27F20E0EB4E537E66>`_, which produces a kind of ``super assembly`` of higher contiguity and accuracy.
* **Directory**: IntegratedContigs/
* **File Names**: {dataset_id}_master_integrated_contigs.fa

.. code-block:: console

   SFBRL043-M3237-14-001_S7_L001_master_integrated_contigs.fa


Contig Annotation
`````````````````

* **Description**: The integrated contigs from the previous step are used as input to `Prokka <https://github.com/tseemann/prokka>`_, a prokaryotic genome annotation tool used to identify genomic features of interest.
* **Directory**: AnnotatedContigs/
* **File Names**: {dataset_id}*

.. code-block:: console

   SFBRL043-M3237-14-001_S7_L001.err		SFBRL043-M3237-14-001_S7_L001.fna  
   SFBRL043-M3237-14-001_S7_L001.gff		SFBRL043-M3237-14-001_S7_L001.tbl
   SFBRL043-M3237-14-001_S7_L001.faa		SFBRL043-M3237-14-001_S7_L001.fsa
   SFBRL043-M3237-14-001_S7_L001.log		SFBRL043-M3237-14-001_S7_L001.txt
   SFBRL043-M3237-14-001_S7_L001.ffn		SFBRL043-M3237-14-001_S7_L001.gbk
   SFBRL043-M3237-14-001_S7_L001.sqn

* **More Info**: For more information about each of the output files produced from Prokka, please see their output files description `page <https://github.com/tseemann/prokka#output-files>`_.

Assembly Evaluation
````````````````````

* **Description**: Assemblies produced from each assembler (including the ``super assembly``) are then evaluated using a genome evaluation tool called `QUAST <https://github.com/ablab/quast>`_. The reports produced can be used to evaluate each assembly based on a variety of metrics such as contig length, number of contigs, and N50. They can also be used to come up with your own assembly score function if you're into that sort of thing.
* **Directory**: AssemblyReport/
* **File Names**: {dataset_id}*

.. code-block:: console

   SFBRL043-M3237-14-001_S7_L001_quast.log		SFBRL043-M3237-14-001_S7_L001_report.html 
   SFBRL043-M3237-14-001_S7_L001_report.tex		SFBRL043-M3237-14-001_S7_L001_report.tsv 
   SFBRL043-M3237-14-001_S7_L001_report.txt		SFBRL043-M3237-14-001_S7_L001_transposed_report.tex 
   SFBRL043-M3237-14-001_S7_L001_transposed_report.tsv	SFBRL043-M3237-14-001_S7_L001_transposed_report.txt

* **More info**: For more information about how to interpret the files produced by QUAST, please see the QUAST output `page <http://quast.bioinf.spbau.ru/manual.html#sec3>`_.
