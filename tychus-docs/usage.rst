Usage
=====

Assembly Parameters
-------------------

Location of the forward read pairs

.. code-block:: console

   -forward /path/to/your/forward/fastq/*.fq

Location of the reverse read pairs

.. code-block:: console

   -reverse /path/to/your/reverse/fastq/*.fq

Location of fasta formatted reference database

.. code-block:: console

   -reference /path/to/your/reference/database.fa

Number of threads to use (default: 10)

.. code-block:: console

   -threads 10

Location for output files to be written to. Note, the output directory does not need to exist prior to pipeline execution

.. code-block:: console

   -output /path/to/where/you/want/your/output/

Alignment Parameters
--------------------

Location of the forward read pairs

.. code-block:: console

   -forward /path/to/your/fastq/*.fq

Location of the reverse read pairs

.. code-block:: console

   -reverse /path/to/your/fastq/*.fq

Location of fasta formatted antimicrobial resistance database

.. code-block:: console

   -amr_db /path/to/your/resistance/database.fa

Location of fasta formatted virulence factor database

.. code-block:: console

   -vf_db /path/to/your/virulence/database.fa

Location of fasta formatted plasmid database

.. code-block:: console

   -plasmid_db /path/to/your/plasmid/database.fa

Number of threads to use (default: 10)

.. code-block:: console

   -threads 10

Location for output files to be written to. Note, the output directory does not need to exist prior to pipeline execution

.. code-block:: console

   -output /path/to/where/you/want/your/output/
