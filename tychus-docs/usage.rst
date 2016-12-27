Usage
=====

The alignment module requires five inputs: a FASTA formatted reference, virulence, plasmid, and AMR database, as well as a pair of FASTQ files. If you don't have access to any of these, you can simply use the reference databases and sequence files provided in the ``tutorial`` directory of the github repository. For example:

.. code-block:: console
   :linenos:

   ./nextflow run alignment.nf -profile alignment --with-docker --output dir

This will run the alignment pipeline with the appropriate Dockerfile and default sequence datasets and write output files to the ``dir`` directory. Users may run any number of FASTQ files in parallel by including an appropriate commandline wildcard with the ``--read_pairs`` option. For example:

.. code-block:: console
   :linenos:

   ./nextflow run alignment.nf -profile alignment --with-docker --read_pairs=/raw_sequence_data/_R{1,2}_001.fastq

The assembly module requires a single input: a pair of FASTQ files. As mentioned above, any number of FASTQ files can be run in parallel by specifying an appropriate commandline wildcard. For example:

.. code-block:: console
   :linenos:

   ./nextflow run assembly.nf -profile assembly --with-docker --read_pairs=/raw_sequence_data/_R{1,2}_001.fastq

For more information about available parameters and options, you can ask for help:

.. code-block:: console
   :linenos:

   ./nextflow run assembly.nf --help


Alignment Module Examples
-------------------------

Assembly Module Examples
------------------------
