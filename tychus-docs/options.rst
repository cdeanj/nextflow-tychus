Options
=======

There are a number of available options that can be utilized throughout the Tychus pipeline. Below we describe the function and use of each parameter.

Alignment Module
----------------

General Options
```````````````

1. **help** - Display help text
 - Default is to not display help text.
 - Optional parameter.

.. code-block:: console

   ./nextflow run alignment.nf --help


2. **read_pairs** - Directory of FASTQ formatted sequence data.

 - To run data in parallel, the FASTQ file names *must* have a strand identifier such as R1 or R2.
 - Required parameter.

.. code-block:: console

   ./nextflow run alignment.nf --read_pairs /tutorial/raw_sequence_data/_R{1,2}.fastq


3. **genome** - Location of the FASTA formatted reference database.

 - Default reference is a Listeria *monocytogenes* database which can be found in the tutorial/genome_reference/ directory.
 - Required parameter.

.. code-block:: console

   ./nextflow run alignment.nf --genome /tutorial/genome_reference/listeriadb.fa


4. **amr_db** - Location of the FASTA formatted antimicrobial resistance database.

 - Default AMR database is the newly published MEGARes database, which can be found in the tutorial/amr_reference/ directory.
 - Required parameter.

.. code-block:: console

   ./nextflow run alignment.nf --amr_db /tutorial/amr_reference/megaresdb.fa


5. **vf_db** - Location of the FASTA formatted virulence factor database.

 - Defaults to a custom curated virulence factor database, which can be found in the tutorial/virulence_reference/ directory.
 - Required parameter.

.. code-block:: console

   ./nextflow run alignment.nf --vf_db /tutorial/virulence_reference/virulencedb.fa


6. **plasmid_db** - Location of the FASTA formatted plasmid database.

 - Defaults to a custom curated plasmid database, which can be found in the tutorial/plasmid_reference/ directory.
 - Required parameter.

.. code-block:: console

   ./nextflow run alignment.nf --plasmid_db /tutorial/plasmid_reference/plasmiddb.fa


7. **threads** - The number of threads to use for each process.

 - Any number of threads can be used.
 - Defaults to 1.
 - Optional parameter.

.. code-block:: console

   ./nextflow run alignment.nf --threads 8


8. **out_dir** - Name of the directory to write output files to.

 - Default is to publish results to the *tychus_alignment_output/* directory.
 - Optional parameter.

.. code-block:: console

   ./nextflow run alignment.nf --output tychus_alignment_output


QC Options
``````````

1. **leading** - Remove leading low quality or N bases.

 - Default is to remove leading low quality or N bases below quality 3.
 - Optional parameter.

2. **trailing** - Remove trailing low quality or N bases.

 - Default is to remove trailing low quality or N bases below quality 3.
 - Optional parameter.

3. **slidingwindow** - Scan read with a sliding window.

 - Default is to scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15.
 - Optional parameter.

4. **minlen** - Name of the directory to write output files to.

 - Defaults to removing reads which are less than 36 bases long.
 - Optional parameter.

Phylogeny Options
`````````````````

1. **ML** - Calculates an ML tree.

 - Default is to calculate an ML tree.
 - Optional parameter.

2. **NJ** - Calculates an NJ tree.

 - Default is to not calculate an NJ tree.
 - Optional parameter.

3. **min_frac** - Calculates a tree based on only SNP loci occurring in at least this fraction of genomes.

 - Default is to calculate a tree based on SNP loci occurring in atleast 0.75 of genomes.
 - Optional parameter.

Assembly Module
----------------

General Options
```````````````

1. **help** - Display help text
 - Default is to not display help text.
 - Optional parameter.

2. **read_pairs** - Directory of FASTQ formatted sequence data.

 - To run data in parallel, the FASTQ file names *must* have a strand identifier such as R1 or R2.
 - Required parameter.

3. **threads** - The number of threads to use for each process.

 - Any number of threads can be used.
 - Defaults to 1.
 - Optional parameter.

4. **out_dir** - Name of the directory to write output files to.

 - Default is to publish results to the *tychus_assembly_output/* directory.
 - Optional parameter.

QC Options
``````````

1. **leading** - Remove leading low quality or N bases.

 - Default is to remove leading low quality or N bases below quality 3.
 - Optional parameter.

2. **trailing** - Remove trailing low quality or N bases.

 - Default is to remove trailing low quality or N bases below quality 3.
 - Optional parameter.

3. **slidingwindow** - Scan read with a sliding window.

 - Default is to scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15.
 - Optional parameter.

4. **minlen** - Name of the directory to write output files to.

 - Defaults to removing reads which are less than 36 bases long.
 - Optional parameter.

Annotation Options
``````````````````

1. **genus** - Genus name (first character of genus name must be capitalized).

 - Must be used with *species* parameter.
 - Optional parameter.

2. **species** - Species name (first character of species name must be lowercase).

 - Must be used with *genus* parameter.
 - Optional parameter.
