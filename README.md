<!---
(Modelled after Will-Rowe's beautiful SEAR README https://github.com/will-rowe/SEAR/blob/master/README.md)
-->

Tychus: A tool to characterize the bacterial genome.
====================================================

Table of Contents
-----------------

[Overview](https://github.com/cdeanj/nextflow-tychus#overview)

[Requirements](https://github.com/cdeanj/nextflow-tychus#requirements)

[Quickstart](https://github.com/cdeanj/nextflow-tychus#quickstart)

[Options](https://github.com/cdeanj/nextflow-tychus#options)

[Example Usage](https://github.com/cdeanj/nextflow-tychus#example-usage)

[Results](https://github.com/cdeanj/nextflow-tychus#results)

[Dependencies](https://github.com/cdeanj/nextflow-tychus#dependencies)

[Contact](https://github.com/cdeanj/nextflow-tychus#contact)

-----------------

Overview
========
Tychus is a tool that allows researchers to perform massively parallel sequence data analysis with the goal of producing a high confidence and comprehensive description of the bacterial genome. Key features of the Tychus pipeline include the assembly, annotation, and phylogenetic inference of large numbers of WGS isolates in parallel using open-source bioinformatics tools and virtualization technology.

-------

Requirements
============
  - Java 7+
  - Docker

To check your Java version, type the following command into a terminal:
```
$ java -version
```
-----------

Quickstart
==========
Install Nextflow
----------------
Open a terminal and type the following commands:
```
$ curl -fsSL get.nextflow.io | bash
$ ./nextflow
```

If installing Nextflow behind a proxy server, you may encounter the following error message:
```
$ Unable to initialize nextflow environment
```
In this case, you can type the following commands to obtain the Nextflow executable.
```
$ wget -O nextflow http://www.nextflow.io/releases/v0.20.1/nextflow-0.20.1-all
$ chmod u+x nextflow
$ ./nextflow
```

It may be beneficial for you to move the Nextflow executable to somewhere you have global execute permissions, such as `/usr/local/bin`. If you don't have the appropriate permissions to do so, you can add the following line to your `.bashrc` file, which can be found inside your `home` directory.
```
export PATH=$PATH:/path/to/your/nextflow/executable
```

Install Tychus Pipeline
-----------------------
The Tychus pipeline can be pulled and installed from Github with the following command:
```
$ git clone https://github.com/cdeanj/nextflow-tychus.git
$ cd nextflow-tychus/
```

----------

Pipeline Options
================
To view available pipeline options for each of the Tychus modules, you can type the following command into a terminal:
```
./nextflow run <module-name> --help
```
Alignment Module
----------------
```
./nextflow run alignment.nf --help

N E X T F L O W  ~  version 0.23.0
Launching `alignment.nf` [tender_wing] - revision: aa90f777d3

Tychus - Alignment Pipeline

Usage: 
    nextflow run alignment.nf -profile alignment [options]

General Options: 
    --read_pairs      DIR		Directory of paired FASTQ files
    --genome          FILE		Path to the FASTA formatted reference database
    --amr_db          FILE		Path to the FASTA formatted resistance database
    --vf_db           FILE		Path to the FASTA formatted virulence database
    --plasmid_db      FILE		Path to the FASTA formatted plasmid database
    --threads         INT		Number of threads to use for each process
    --out_dir         DIR		Directory to write output files to

Trimmomatic Options: 
    --leading         INT		Remove leading low quality or N bases
    --trailing        INT		Remove trailing low quality or N bases
    --slidingwindow   INT		Scan read with a sliding window
    --minlen          INT		Drop reads below INT bases long

kSNP Options: 
    --ML              BOOL		Estimate maximum likelihood tree
    --NJ              BOOL		Estimate neighbor joining tree
    --min_frac        DECIMAL	Minimum fraction of genomes with locus
```

Assembly Module
---------------
```
./nextflow run assembly.nf --help

N E X T F L O W  ~  version 0.23.0
Launching `assembly.nf` [sleepy_bohr] - revision: 05adc382a5

Tychus - Assembly Pipeline

Usage: 
    nextflow run assembly.nf -profile assembly [options]

General Options: 
    --read_pairs      DIR		Directory of paired FASTQ files
    --threads         INT       Number of threads to use for each process
    --output          DIR       Directory to write output files to

Trimmomatic Options: 
    --leading         INT		Remove leading low quality or N bases
    --trailing        INT		Remove trailing low quality or N bases
    --slidingwindow   STR		Scan read with a sliding window
    --minlen          INT		Drop reads below INT bases long
    
Prokka Options:
    --genus           STR		Target genus
    --species         STR		Target species
```

----------------

Example Usage
=============
FASTQ Input
-----------
The most useful command for both modules will be to read in your sequence data. With Nextflow, we can specify a command line glob to provide a directory of FASTQ files as input. Doing so will allow Nextflow to process data in parallel, using multiple processors. For example, a typical command may look like the following:

```
./nextflow run <module-name> -profile <profile-name> --read_pairs "tutorial/raw_sequence_data/*_R{1,2}_001.fastq.gz"
```

Here, we are using the `*` wildcard to grab all files within the `tutorial/raw_sequence_data/` directory. The `{1,2}` wildcards allows us to further group the files based on the presence of an `_R1` or `_R2` substring. What is returned is essentially a sorted list of files that Nextflow can group together and process appropriately.

Trimmomatic Operations
----------------------
You may want to use your own trimming operation instead of the defaults providided by each module. To change them you can enter the following command:
```
./nextflow run <module-name> -profile <profile-name> --read_pairs "tutorial/raw_sequence_data/*_R{1,2}_001.fastq.gz" --leading 5 --trailing 5 --slidingwindow 5:16 --minlen 45
```

kSNP Operations
---------------
By default, maximum likelihood (ML) trees are computed with kSNP. Although this is the `recommended` tree format to produce, you can specify the neighbor joining (NJ) method by including the `--NJ` option. Furthermore, you can enter a decimal number between 0 and 1 specifying the fraction of loci that must be present in all genomes to be included in the resulting SNP phylogeny.
```
./nextflow run alignment.nf -profile alignment --read_pairs "tutorial/raw_sequence_data/*_R{1,2}_001.fastq.gz" --NJ --min_frac 0.85
```

Figtree Options
---------------
By deafult the SNP phylogenies produced by kSNP and written to a Newick formatted `.tre` file. Figtree is used to produce phylogenies in the image format of your choosing. By default, SNP phylognies are annotated and saved as scalable vector graphic [SVG](https://en.wikipedia.org/wiki/Scalable_Vector_Graphics) images. To change this, simply specify an alternative image format (JPEG,PDF,PNG,SVG).
```
./nextflow run alignment.nf -profile alignment --read_pairs "tutorial/raw_sequence_data/*_R{1,2}_001.fastq.gz" --JPEG
```

Prokka Options
--------------
We allow users to annotate contigs using specific BLAST databases. To do this, you must specify both the `genus` and `species` parameters.
```
./nextflow run assembly.nf -profile assembly --read_pairs "tutorial/raw_sequence_data/*_R{1,2}_001.fastq.gz" --genus Listeria --species monocytogenes
```

Database Options
----------------
If you would like to specify an alternative `reference`, `virulence`, `plasmid` or `resistance` database than the ones already provided, you can do that as well.
```
./nextflow run alignment.nf -profile alignment --read_pairs "tutorial/raw_sequence_data/*_R{1,2}_001.fastq.gz" --ref_db "path/to/your/reference/db/ref.fa" --vf_db "path/to/your/virulence/db/vf.fa" --plasmid path/to/your/plasmid/db/plasmid.fa --amr_db path/to/your/resistance/db/resistance.fa
```

Other Options
-------------
By now you should be familiar with how to specify the various options provided by each module. Here are some more options. The `threads` parameter allows you to control how many threads each process will use. By default, this value is set to 1. The `output` directory allows you to specify where outputs will be stored.
```
./nextflow run <module-name> -profile <profile-name> --read_pairs "tutorial/raw_sequence_data/*_R{1,2}_001.fastq.gz" --threads 4 --output dir
```

----------------

Results
=======
Alignment Module
----------------
Directory | Description
--------- | -----------
Alignment | `Contains all BAM formatted alignment files produced by the alignment of reads against the user-input reference, plasmid, resistance, and virulence databases.`
Consensus| `Contains all FASTA formatted consensus sequences produced by the VCF formatted SNPs called by FreeBayes.`
PreProcessing | `Contains all FASTQ formatted trimmed sequence files produced by Trimmomatic.`
Resistome | `Contains all TSV formatted resistome files.`
SNPsAndPhylogenies | `Contains all SNPs and Newick formatted Phylogenies produced by kSNP3. The SNP files can be found in the SNPs/ subdirectory. The Newick formatted phylogenies can be found in the Trees/ directory. The Newick formatted image files can be found in the TreeImages/ directory.`

Assembly Module
---------------
Directory | Description
--------- | -----------
AbyssContigs | `Contains all FASTA formatted contigs produced by the Abyss assembler.`
IDBAContigs| `Contains all FASTA formatted contigs produced by the IDBA-UD assembler.`
SPAdesContigs | `Contains all FASTA formatted contigs produced by the SPAdes assembler.`
VelvetContigs | `Contains all FASTA formatted contigs produced by the Velvet assembler.`
IntegratedContigs | `Contains all super assembly contigs produced by the CISA contig integrator.`
AnnotatedContigs | `Contains all annotation files produced by Prokka.`
AssemblyReport | `Contains all assembly evaulation files produced by QUAST.`
PreProcessing | `Contains all FASTQ formatted trimmed sequence files produced by Trimmomatic.`

-------

Dependencies
============

Tychus utilizes a number of open-source bioinformatics tools to run. Please click on the tool names below to learn more about each tool.

Software | Function
--------- | --------
[Abyss](https://github.com/bcgsc/abyss) | `Used to produce assembly contigs.`
[BCFtools](https://github.com/samtools/bcftools) | `Used to generate consensus sequences from VCF formatted SNPs.`
[Bowtie2](https://github.com/BenLangmead/bowtie2) | `Used to align short fragments of DNA to a reference genome.`
[CISA](http://sb.nhri.org.tw/CISA/en/CISA) | `Used to integrate assembly contigs into a super assembly.`
[Docker](https://www.docker.com/what-docker) | `Software containerization platform used to resolve the dependencies listed here.`
[Figtree](http://tree.bio.ed.ac.uk/software/figtree/) | `Used to create images from Newick formatted phylogenies.`
[IDBA-UD](https://github.com/loneknightpy/idba) | `Used to produce assembly contigs.`
[KmerGenie](http://kmergenie.bx.psu.edu) | `Used for optimizing chosen values of k (kmer) for non-iterative genome assemblers.`
[kSNP3](https://sourceforge.net/projects/ksnp/files/) | `Used to generate SNPs and SNP phylogenies.`
[Nextflow](https://www.nextflow.io) | `Used as the backend framework for the Tychus pipeline.`
[Prokka](https://github.com/tseemann/prokka) | `Used to identify genomic features of interest.`
[QUAST](https://github.com/ablab/quast) | `Used for the evaluation and interrogation of assembly contigs.`
[Samtools](https://github.com/samtools/samtools) | `Used for manipulating SAM/BAM formatted alignment files.`
[SPAdes](https://github.com/ablab/spades) | `Used to produce assembly contigs.`
[Trimmomatic](https://github.com/timflutre/trimmomatic) | `Used for the removal of adapter sequences and low quality base pairs.`
[Velvet](https://github.com/dzerbino/velvet) | `Used to produce assembly contigs.`
[VelvetOptimiser](https://github.com/tseemann/VelvetOptimiser) | `Used to optimize parameter values for the Velvet assembler.`

------------

Contact
=======
Questions, bugs, or feature requests should be directed to Chris Dean at cdean11 at rams dot colostate dot edu. Alternatively, you can [Submit an Issue](https://github.com/cdeanj/nextflow-tychus/issues) on Github.
