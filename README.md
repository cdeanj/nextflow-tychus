# What is Tychus?
That's a good question

## Installation
### [Nextflow](http://www.nextflow.io/) Installation
```bash
$ curl -fsSL get.nextflow.io | bash 
```

## Running the Pipeline
```bash
$ ./nextflow run tychus.nf -with-docker
```

## Parameters

##### `--forward`
Location of the forward read pairs

##### `--reverse`
Location of the reverse read pairs

##### `--amr_db`
Location of the fasta formatted antimicrobial resistance database

##### `--vf_db`
Location of the fasta formatted virulence factor database

##### `--plasmid_db`
Location of the fasta formatted plasmid database

##### `--threads`
Number of threads to use for each process (where applicable)

##### `--output`
Location for output files to be written to. Note, directory does not need to exist prior to pipeline execution


## Dependencies

Tychus utilizes a number of open source projects to run:

* [Docker](https://www.docker.com/what-docker) - Software containerization platform
* [Trimmomatic](https://github.com/timflutre/trimmomatic) - Read trimmer
* [Bowtie2](https://github.com/BenLangmead/bowtie2) - Sequence aligner
* [Samtools](https://github.com/samtools/samtools) - Alignment processor
* [Freebayes](https://github.com/ekg/freebayes) - Variant caller
* [Coverage Sampler](https://github.com/cdeanj/coverage_sampler) - Rarefaction calculator
* [kSNP3](https://sourceforge.net/projects/ksnp/) KMER selection

