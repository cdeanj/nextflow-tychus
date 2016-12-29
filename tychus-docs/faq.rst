Frequently Asked Questions
==========================

Q: Can you explain how to process multiple FASTQ files again?

A: Sure!

Q: Can I use the Tychus assembly module to run my Eukaryotic sequence data?

A: No. Many of the tools used for assembly and annotation are specifically designed to process Prokaryotic genomes.

Q: After running the Tychus pipeline I see that a *work* directory is created. What is that and can I delete it?

A: The *work* directory is created by Nextflow to store all of the intermediate outputs produced by each process as the pipeline runs. Many of these outputs shouldn't be useful to you as a user and may be safely deleted to clean up disk space.

Q: Can I use compressed FASTQ data as input to each module of the Tychus pipeline?

A: Currently, no. These files must be uncompressed.

Q: How long will the pipeline take to complete on a large number of FASTQ datasets?

A: Running ~200 Listeria *monocytogenes* samples took about 8 hours to complete using a 64-core server.
