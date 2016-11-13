Output Files
============

Assembly Module
---------------

============================ ===============================================================
File Name                    Description
============================ ===============================================================
abyss-contigs.fa            
velvet-contigs.fa
spades-contigs.fa
idba-contigs.fa
master_integrated_contigs.fa
============================ ===============================================================

Dependencies
------------
Tychus utilizes a number of open source projects, which are all resolved by Docker:

* `Nextflow <https://www.nextflow.io>`_ Workflow framework
* `Docker <https://www.docker.com/what-docker>`_ Software containerization platform
* `Trimmomatic <http://bioinformatics.oxfordjournals.org/content/early/2014/04/01/bioinformatics.btu170>`_ Read trimmer and quality control
* `Bowtie2 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3322381/>`_ Short-read sequence aligner
* `Samtools <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2723002/>`_ SAM/BAM alignment processor
* `Freebayes <https://arxiv.org/abs/1207.3907>`_ Probabilistic variant caller
* `Prokka <https://www.ncbi.nlm.nih.gov/pubmed/24642063>`_ Prokaryotic genome annotation tool
* `CoverageSampler <https://github.com/cdeanj/coverage_sampler>`_ Resistome analyzer
* `kSNP3 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3857212/>`_ Phylogenetic analysis
* `KmerGenie <https://arxiv.org/pdf/1304.5665.pdf>`_ Optimal kmer selection for building De-Bruijn graphs
* `Abyss <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2694472/>`_ *De novo* sequence assembler for short-paired reads
* `SPades <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3342519/>`_ Assembler for single-celled bacterial genomes
* `IDBA-UD <http://i.cs.hku.hk/~chin/paper/idba_ud-revised-latest.pdf>`_ Genome assembler for short reads
* `Velvet <http://genome.cshlp.org/content/genome/18/5/821.full.html>`_ *De novo* short-read assembler
* `CISA <http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0060843q>`_ Contig integrator
