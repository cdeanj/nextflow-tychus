Usage
=====

Assembly Parameters
-------------------

========== ==================== ================================================================================================================== ==================================================
**Option** **Parameter**        **Description**                                                                                                    **Default behavior**
-forward   forward read file(s) Read forward FASTQ files from the specified directory. These files must have an R1 substring within the file name. Program will use tutorial data.
-reverse   reverse read file(s) Read reverse FASTQ files from the specified directory. These files must have an R2 substring within the file name. Program will use tutorial data.
-threads   integer              The number of threads to use for each process.                                                                     The program will use a single thread.
-output    dir                  Create a folder called dir and write output files to it.                                                           The program will create a directory called output.
========== ==================== ================================================================================================================== ==================================================

Alignment Parameters
--------------------




The ``forward`` parameter allows you to specify the location of the forward read pairs

.. code-block:: console

   -forward /path/to/your/fastq/*.fq

The ``reverse`` parameter allows you to specify the location of the reverse read pairs

.. code-block:: console

   -reverse /path/to/your/fastq/*.fq

The ``amr_db`` parameter allows you to specify the location of your fasta formatted antimicrobial resistance database

.. code-block:: console

   -amr_db /path/to/your/resistance/database.fa

The ``vf_db`` parameter allows you to specify the location of your fasta formatted virulence factor database

.. code-block:: console

   -vf_db /path/to/your/virulence/database.fa

The ``plasmid_db`` parameter allows you to specify the location of your fasta formatted plasmid database

.. code-block:: console

   -plasmid_db /path/to/your/plasmid/database.fa

The ``threads`` parameter can be used to specify the number of threads to use for each process (default: 10)

.. code-block:: console

   -threads 10

The ``output`` parameter can be used to specify the location for writing output files. Note, the output directory does not need to exist prior to pipeline execution

.. code-block:: console

   -output /path/to/where/you/want/your/output/
