Installation
============

Tychus can be run on all major operating systems, including ``MAC OS X``, ``Linux 64-bit``, and ``Windows``. If running on Windows, you will need access to a virtual machine, such as `Virtual Box <https://www.virtualbox.org>`_.

To run Tychus, you will need `Nextflow <https://www.nextflow.io>`_, a DSL workflow framework for running computational pipelines, and `Docker <https://www.docker.com>`_, a software containerization platform that resolves dependency nightmares for complex workflows such as Tychus. See below for installation instructions.

Clone Repository
----------------
.. code-block:: console
   :linenos:

   git clone https://github.com/cdeanj/nextflow-tychus.git


Nextflow Installation
---------------------
.. code-block:: console
   :linenos:

   cd nextflow-tychus
   curl -fsSL get.nextflow.io | bash
   ./nextflow


If you are installing nextflow behind a proxy server, you may encounter the following error:

.. code-block:: console
   :linenos:

   Unable to initialize nextflow environment


In this case, open a terminal and type the following commands:

.. code-block:: console
   :linenos:

   wget -O nextflow http://www.nextflow.io/releases/v0.20.1/nextflow-0.20.1-all
   chmod u+x nextflow
   ./nextflow


Docker Installation
-------------------
Docker can be installed on both Linux and Mac operating systems. See below for installation instructions.

Linux
`````
.. code-block:: console
   :linenos:

   # Update package manager
   sudo apt-get update
   # Add the GPG key for the official Docker repository
   sudo apt-key adv --keyserver hkp://p80.pool.sks-keyservers.net:80 --recv-keys 58118E89F3A912897C070ADBF76221572C52609D
   # Add the Docker repository to APT sources
   sudo apt-add-repository 'deb https://apt.dockerproject.org/repo ubuntu-xenial main'
   # Update package manager
   sudo apt-get update
   # Install Docker
   sudo apt-get install -y docker-engine
   # Run docker command to make sure installation went as planned
   docker

MAC OS X
````````
.. code-block:: console
   :linenos:

    # Install Homebrew
    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    # Install Cask
    brew install caskroom/cask/brew-cask
    # Install docker toolbox
    brew cask install docker-toolbox
    # create the docker machine
    docker-machine create --driver "virtualbox" myBoxName
    # start the docker machine
    docker-machine start myBoxName
    # this command allows the docker commands to be used in the terminal
    eval "$(docker-machine env myBoxName)"
    # at this point can run any "docker" or "docker-compose" commands you want
    docker-compose up

Pull Images
-----------
To run the ``alignment`` pipeline you will need the latest ``tychus-alignment`` Docker image.

.. code-block:: console
   :linenos:

   docker pull abdolab/tychus-alignment

To run the ``assembly`` pipeline, you will need the latest ``tychus-assembly`` Docker image.

.. code-block:: console
   :linenos:

   docker pull abdolab/tychus-assembly

======================== =============== =============== ================= =============
Repository               Tag             Image OS        Image Size        Download Time
======================== =============== =============== ================= =============
abdolab/tychus-alignment Latest          Ubuntu 16.04    3.5 GB            4 minutes
abdolab/tychus-assembly  Latest          Ubuntu 16.04    3.4 GB            7 minutes
======================== =============== =============== ================= =============

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
* `kSNP3 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3857212/>`_ Phylogenetic analysis
* `KmerGenie <https://arxiv.org/pdf/1304.5665.pdf>`_ Optimal kmer selection for building De-Bruijn graphs
* `Abyss <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2694472/>`_ *De novo* sequence assembler for short-paired reads
* `SPades <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3342519/>`_ Assembler for single-celled bacterial genomes
* `IDBA-UD <http://i.cs.hku.hk/~chin/paper/idba_ud-revised-latest.pdf>`_ Genome assembler for short reads
* `Velvet <http://genome.cshlp.org/content/genome/18/5/821.full.html>`_ *De novo* short-read assembler
* `CISA <http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0060843q>`_ Contig integrator
* `FigTree <http://tree.bio.ed.ac.uk/software/figtree/>`_ Newick to image converter
