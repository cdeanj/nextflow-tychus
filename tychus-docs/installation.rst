Installation
============

Tychus can be run on all major operating systems, including ``MAC OS X``, ``Linux 64-bit``, and ``Windows``. If running on Windows, you will need access to a virtual machine, such as `Virtual Box <https://www.virtualbox.org>`_.

To run Tychus, you will need `Nextflow <https://www.nextflow.io>`_, a DSL workflow framework for running computational pipelines, and `Docker <https://www.docker.com>`_, a software containerization platform that resolves dependency nightmares for complex workflows such as Tychus. See below for installation instructions.

Nextflow Installation
---------------------
.. code-block:: console
   :linenos:

   curl -fsSL get.nextflow.io | bash
   ./nextflow

Docker Installation
-------------------
.. code-block:: console
   :linenos:

   sudo apt-get update
   sudo apt-get install docker-compose
   sudo groupadd docker
   sudo usermod -aG docker $USER

Pull Images
-----------
Before we pull, let's understand what Docker images and containers are. You can think of a Docker image as a class, similar to what you see in programming languages like Java, C++, or Python. Classes have attributes. In our case, these attributes will be our dependencies (i.e., IDBA, KmerGenie, Prokka, etc..). Furthermore, we know that an instance of a class is an object. For our case, the object is the running container. To sum up, an image is a description of a class and a container is an instance of an object, or an instance of an image.

If you wish to run the ``alignment`` pipeline you will need the latest ``tychus-alignment`` Docker image.

.. code-block:: console
   :linenos:

   docker pull abdolab/tychus-alignment

If you wish to run the ``assembly`` pipeline, you will need the latest ``tychus-assembly`` Docker image.

.. code-block:: console
   :linenos:

   docker pull abdolab/tychus-assembly

======================== =============== =============== ================= =============
Repository               Tag             Image OS        Image Size        Download Time
======================== =============== =============== ================= =============
abdolab/tychus-alignment Latest          Ubuntu 15.10    1.436 GB          4 minutes
abdolab/tychus-assembly  Latest          Ubuntu 15.10    2.270 GB          7 minutes
======================== =============== =============== ================= =============

Dependencies
------------
Tychus utilizes a number of open source projects, which are all resolved by Docker:

* `Nextflow <https://www.nextflow.io>`_ Workflow framework
* `Docker <https://www.docker.com/what-docker>`_ Software containerization platform
* Trimmomatic `[Paper] <http://bioinformatics.oxfordjournals.org/content/early/2014/04/01/bioinformatics.btu170>`_ `[Code] <https://github.com/timflutre/trimmomatic>`_ Read trimmer and quality control
* Bowtie2 `[Paper] <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3322381/>`_ `[Code] <https://github.com/BenLangmead/bowtie2>`_ Short-read sequence aligner
* Samtools `[Paper] <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2723002/>`_ `[Code]<https://github.com/samtools/samtools>`_ SAM/BAM alignment processor
* Freebayes `[Paper] <https://arxiv.org/abs/1207.3907>`_ `[Code] <https://github.com/ekg/freebayes>`_ Probabilistic variant caller
* `Prokka `[Paper] <https://www.ncbi.nlm.nih.gov/pubmed/24642063>`_ `[Code] <https://github.com/tseemann/prokka>`_ Prokaryotic genome annotation tool
* CoverageSampler `[Code] <https://github.com/cdeanj/coverage_sampler>`_ Resistome analyzer
* kSNP3 `[Paper] <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3857212/>`_ `[Code] <https://sourceforge.net/projects/ksnp/>`_ Phylogenetic analysis
* KmerGenie `[Paper] <https://arxiv.org/pdf/1304.5665.pdf>`_ `[Code]<http://kmergenie.bx.psu.edu/>`_ Optimal kmer selection for building De-Bruijn graphs
* Abyss `[Paper] <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2694472/>`_ `[Code] <https://github.com/bcgsc/abyss>`_ *De novo* sequence assembler for short-paired reads
* SPades `[Paper] <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3342519/>`_ `[Code] <https://github.com/ablab/spades>`_ Assembler for single-celled bacterial genomes
* IDBA-UD `[Paper] <http://i.cs.hku.hk/~chin/paper/idba_ud-revised-latest.pdf>`_ `[Code] <https://github.com/loneknightpy/idba>`_ Genome assembler for short reads
* Velvet `[Paper] <http://genome.cshlp.org/content/genome/18/5/821.full.html>`_ `[Code] <https://github.com/dzerbino/velvet>`_ *De novo* short-read assembler
* CISA `[Paper] <http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0060843q`_ `[Code] <http://sb.nhri.org.tw/CISA/en/CISA;jsessionid=125169F363E3D18705C397E7C6F68C8E>`_ Contig integrator
