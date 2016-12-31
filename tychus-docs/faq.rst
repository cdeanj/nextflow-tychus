Frequently Asked Questions
==========================

Q: Can you explain how to process multiple FASTQ files again?

A: Sure! Often times in bioinformatics, you want to work with multiple files at once. To do this we can take advantage of command line globs. A glob is simply a way of grouping file names based on a particular search pattern. For example, suppose we want to provide two files as input (**SRR532663_1.fastq** and **SRR532663_2.fastq**) which are stored in a directory called **raw/**. To do this we could write:

.. code-block:: console
   :linenos:

   --read_pairs = raw/_{1,2}.fastq

The ***** is a wildcard character that can recognize any string of characters. In this case, we want to recognize any characters ending with **_1** or **_2** that occur before the **.fastq** suffix and the curly brackets allow us to match strings with either a **1** or **2** character.

Q: Can I use compressed FASTQ data as input to each module of the Tychus pipeline?

A: Currently, no. These files must be uncompressed.

Q: Can I use the Tychus assembly module to run my Eukaryotic sequence data?

A: No. Many of the tools used for assembly and annotation are specifically designed to process Prokaryotic genomes.

Q: After running the Tychus pipeline I see that a *work* directory is created. What is that and can I delete it?

A: The *work* directory is created by Nextflow to store all of the intermediate outputs produced by each process as the pipeline runs. Many of these outputs shouldn't be useful to you as a user and may be safely deleted to clean up disk space.

Q: How long will the pipeline take to complete on a large number of FASTQ datasets?

A: Running ~200 Listeria *monocytogenes* samples took about 8 hours to complete using a 64-core server.
