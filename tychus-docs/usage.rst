Usage
=====

The alignment module requires five inputs: a FASTA formatted reference, virulence, plasmid, and AMR database, as well as a pair of FASTQ files. If you don't have access to any of these, you can simply use the reference databases and sequence files provided in the ``tutorial`` directory of the github repository. For example:

.. code-block:: console
   :linenos:

   ./nextflow run alignment.nf -profile alignment --with-docker --threads 10 --output dir

This will run the alignment pipeline with the appropriate Dockerfile and default sequence datasets and write output files to the ``dir`` directory.

The assembly module requires a single input: a pair of FASTQ files. To run the tutorial data, you can simply type the following command into a terminal:

.. code-block:: console
   :linenos:

   ./nextflow run assembly.nf -profile assembly --with-docker --threads 10 --output dir

For more information about available parameters and options, you can ask for help:

.. code-block:: console
   :linenos:

   ./nextflow run [nextflow-script-name.nf] --help


Alignment Module Examples
-------------------------

Assembly Module Examples
------------------------
