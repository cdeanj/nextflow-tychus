Example Usage
=============

Below we provide common example usages for the Tychus pipeline.

Alignment Module
````````````````

Assembly Module
```````````````

Run a directory of FASTQ files with an R1 and R2 strand specifier.

.. code-block:: console
   :linenos:

   ./nextflow run assembly.nf -profile assembly --with-docker --read_pairs tutorial/raw_sequence_data/\*_R{1,2}_001.fastq
