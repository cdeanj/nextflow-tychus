Example Usage
=============

Below we provide common example usages for the Tychus pipeline.

Alignment Module
````````````````

Run a directory of FASTQ files with an R1 and R2 strand specifier.

.. code-block:: console
   :linenos:

   ./nextflow run assembly.nf \
   -profile assembly \
   --with-docker \
   --read_pairs tutorial/raw_sequence_data/_R{1,2}_001.fastq

.. note ::

    It is not necessary to run these commands as multi-line arguments. You can instead enter them on a single line as shown below.

.. code-block:: console
   :linenos:

   ./nextflow run assembly.nf -profile assembly --with-docker --read_pairs tutorial/raw_sequence_data/_R{1,2}_001.fastq


Perform quality filtering for each pair of FASTQ files, removing leading and trailing bases below an average quality of 3, averaging across 5 bases with a minimum average quality of 15, and dropping reads below 36 base pairs.

.. code-block:: console
   :linenos:

   ./nextflow run assembly.nf \
   -profile assembly \
   --with-docker \
   --read_pairs tutorial/raw_sequence_data/_R{1,2}_001.fastq \
   --leading 3 \
   --trailing 3 \
   --slidingwindow 4:15 \
   --minlen 36


Assembly Module
```````````````

Run a directory of FASTQ files with an R1 and R2 strand specifier.

.. code-block:: console
   :linenos:

   ./nextflow run assembly.nf -profile assembly --with-docker --read_pairs tutorial/raw_sequence_data/_R{1,2}_001.fastq

Perform quality filtering for each pair of FASTQ files, removing leading and trailing bases below an average quality of 3, averaging across 5 bases that must have an average quality of 15, and dropping reads below 36 base pairs.

.. code-block:: console
   :linenos:

   ./nextflow run assembly.nf \
   -profile assembly \
   --with-docker \
   --read_pairs tutorial/raw_sequence_data/_R{1,2}_001.fastq \
   --leading 3 \
   --trailing 3 \
   --slidingwindow 4:15 \
   --minlen 36


