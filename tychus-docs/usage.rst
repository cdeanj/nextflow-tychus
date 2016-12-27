Usage
=====

Tychus can be utilized as an alignment-based or assembly-based pipeline. The alignment module requires five inputs: a FASTA formatted reference, virulence, plasmid, and AMR database, as well as a pair of FASTQ files. If you don't have access to any of these, you can simply use the reference databases and sequence files provided in the ``tutorial`` directory of the github repository. For example:

.. code-block:: console
   :linenos:

   ./nextflow run alignment.nf -profile alignment --with-docker --output dir

This will run the alignment pipeline with the appropriate Dockerfile and default sequence datasets and write output files to the ``dir`` directory. Users may run any number of FASTQ files in parallel by including an appropriate commandline wildcard with the ``--read_pairs`` option. For example:

.. code-block:: console
   :linenos:

   ./nextflow run alignment.nf -profile alignment --with-docker --read_pairs=/raw_sequence_data/\*_R{1,2}_001.fastq

Assembly Parameters
-------------------

============ ==================== ================================================================================================================== ==================================================
**Option**   **Parameter**        **Description**                                                                                                    **Default behavior**
--read_pairs read pair files      Directory of paired (forward and reverse) FASTQ files.                                                             Program will use tutorial data.
--threads    integer              The number of threads to use for each process.                                                                     Program will use a single thread.
--output     dir                  Create a folder called dir and write output files to it.                                                           Program will create a directory called output.
============ ==================== ================================================================================================================== ==================================================

Alignment Parameters
--------------------

============ =========================== ================================================================================================================== ==============================================
**Option**   **Parameter**               **Description**                                                                                                    **Default behavior**
--read_pairs read pair files             Directory of paired (forward and reverse) FASTQ files.                                                             Program will use tutorial data.
--ref_db     primary reference database  Name of the FASTA formatted reference database.                                                                    Program will use tutorial reference database.
--amr_db     primary resistance database Name of the FASTA formatted resistance database.                                                                   Program will use tutorial resistance database.
--vf_db      primary virulence database  Name of the FASTA formatted virulence database.                                                                    Program will use tutorial virulence database.
--plasmid_db primary plasmid database    Name of the FASTA formatted plasmid database.                                                                      Program will use tutorial plasmid database.
--threads    integer                     Number of threads to use for each process.                                                                         Program will use a single thread.
--output     dir                         Create a folder called dir and write output files to it.                                                           Program will create a directory called output.
============ =========================== ================================================================================================================== ==============================================
