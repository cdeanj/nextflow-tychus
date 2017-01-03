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
Pull the latest ``tychus-alignment`` Docker image.

.. code-block:: console
   :linenos:

   docker pull abdolab/tychus-alignment

Pull the latest ``tychus-assembly`` Docker image.

.. code-block:: console
   :linenos:

   docker pull abdolab/tychus-assembly

======================== ====== =============== ================= ================
Repository               Tag    Image OS        Image Size        Download Time
======================== ====== =============== ================= ================
abdolab/tychus-alignment Latest Ubuntu 16.04    3.5 GB            ~ 5 - 15 minutes
abdolab/tychus-assembly  Latest Ubuntu 16.04    3.1 GB            ~ 5 - 15 minutes
======================== ====== =============== ================= ================

