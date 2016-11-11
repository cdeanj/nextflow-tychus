Installation
============

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
If you wish to run the alignment pipeline you will need to obtain the latest tychus-alignment Docker image

.. code-block:: console
   :linenos:

   docker pull chrisd/tychus-alignment

If you wish to run the assembly pipeline, you will need to obtain the latest tychus-assembly Docker image

.. code-block:: console
   :linenos:

   docker pull chrisd/tychus-assembly

======================= =============== 
Image Name              Image Size     
======================= =============== 
chrisd/tychus-alignment 1.779 GB
chrisd/tychus-assembly  1.000 GB
======================= ===============

Dependencies
------------
Tychus utilizes a number of open source projects, which are all resolved by Docker:

* `Nextflow <https://www.nextflow.io>`_ Workflow framework
* `Docker <https://www.docker.com/what-docker>`_ Software containerization platform
* `Trimmomatic <https://github.com/timflutre/trimmomatic>`_ Read trimmer and quality control
* `Bowtie2 <https://github.com/BenLangmead/bowtie2>`_ Short-read sequence aligner
* `Samtools <https://github.com/samtools/samtools>`_ SAM/BAM alignment processor
* `Freebayes <https://github.com/ekg/freebayes>`_ Probabilistic variant caller
* `Prokka <https://github.com/tseemann/prokka>`_ Prokaryotic genome annotation tool
* `CoverageSampler <https://github.com/cdeanj/coverage_sampler>`_ Resistome analyzer
* `kSNP3 <https://sourceforge.net/projects/ksnp/>`_ Phylogenetic analysis
* `KmerGenie <http://kmergenie.bx.psu.edu/>`_ Optimal kmer selection for build De-Bruijn graphs
* `Abyss <https://github.com/bcgsc/abyss>`_ *De novo* sequence assembler for short-paired reads
* `SPades <http://spades.bioinf.spbau.ru/release3.9.0/manual.html>`_ Assembler for single-celled bacterial genomes
* `IDBA-UD <https://github.com/loneknightpy/idba>`_ Genome assembler for short reads
* `Velvet <https://github.com/dzerbino/velvet>`_ *De novo* short-read assembler
* `CISA <http://sb.nhri.org.tw/CISA/en/CISA;jsessionid=125169F363E3D18705C397E7C6F68C8E>`_ Contig integrator
