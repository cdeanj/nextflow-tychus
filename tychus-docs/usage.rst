Usage
=====

Assembly Parameters
-------------------

-forward

Location of the forward read pairs

.. code-block:: console

   -forward /path/to/your/fastq/*.fq

-reverse

Location of the reverse read pairs

.. code-block:: console

   -reverse /path/to/your/fastq/*.fq

-amr_db

Location of fasta formatted antimicrobial resistance database

.. code-block:: console

   -amr_db /path/to/your/resistance/database.fa

-vf_db

Location of fasta formatted virulence factor database

.. code-block:: console

   -vf_db /path/to/your/virulence/database.fa

-plasmid_db

Location of fasta formatted plasmid database

.. code-block:: console

   -plasmid_db /path/to/your/plasmid/database.fa

-threads

Number of threads to use (default: 10)

.. code-block:: console

   -threads 10

-output

Location for output files to be written to. Note, the output directory does not need to exist prior to pipeline execution

.. code-block:: console

   -output /path/to/where/you/want/your/output/

Alignment Parameters
--------------------
