Usage
=====

Assembly Parameters
-------------------

The assembly module can run reference-based or *de novo* assemblies. If you want to run de novo assemblies, simply omit the reference parameter. If you wish to run a reference-based assembly, specify the location of your fasta formatted reference database.

The ``forward`` parameter allows you to specify the location of the forward read pairs

.. code-block:: console

   -forward /path/to/your/forward/fastq/*.fq

The ``reverse`` parameter allows you to specify the location of the reverse read pairs

.. code-block:: console

   -reverse /path/to/your/reverse/fastq/*.fq

The ``reference`` parameter allow you to specify the location of fasta formatted reference database

.. code-block:: console

   -reference /path/to/your/reference/database.fa

The ``threads`` parameter can be used to specify the number of threads to use for each process (default: 10)

.. code-block:: console

   -threads 10

The ``output`` parameter can be used to specify the location for writing output files. Location for output files to be written to. Note, the output directory does not need to exist prior to pipeline execution

.. code-block:: console

   -output /path/to/where/you/want/your/output/

Alignment Parameters
--------------------

It is not required to have access to a resistance, virulence, or plasmid database to run the alignment module. These databases have been custom curated and provided in the tutorial directory of the Tychus repository. If you have access to your own databases, simply specify the appropriate parameter and location of each database on your file system.

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
