Usage
=====

Assembly Parameters
-------------------
R = Required Parameters

NR = Not Required Parameters

Location of the forward read pairs (R)

.. code-block:: console

   -forward /path/to/your/fastq/*.fq

Location of the reverse read pairs (R)

.. code-block:: console

   -reverse /path/to/your/fastq/*.fq

Location of fasta formatted antimicrobial resistance database (NR)

.. code-block:: console

   -amr_db /path/to/your/resistance/database.fa

Location of fasta formatted virulence factor database (NR)

.. code-block:: console

   -vf_db /path/to/your/virulence/database.fa

Location of fasta formatted plasmid database (NR)

.. code-block:: console

   -plasmid_db /path/to/your/plasmid/database.fa

Number of threads to use (default: 10) (NR)

.. code-block:: console

   -threads 10

Location for output files to be written to. Note, the output directory does not need to exist prior to pipeline execution (R)

.. code-block:: console

   -output /path/to/where/you/want/your/output/

Alignment Parameters
--------------------
