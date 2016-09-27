[![nextflow](https://img.shields.io/badge/nextflow-≥0.14.3-brightgreen.svg)](http://nextflow.io)

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
```bash
$ --forward /path/to/your/fastq/*.fq
```

##### `--reverse`
Location of the reverse read pairs
```bash
$ --reverse /path/to/your/fastq/*.fq
```

##### `--amr_db`
Location of the fasta formatted antimicrobial resistance database
```bash
$ --amr_db /path/to/your/resistance/database/*.fa
```

##### `--vf_db`
Location of the fasta formatted virulence factor database
```bash
$ --vf_db /path/to/your/virulence/database/*.fa
```

##### `--plasmid_db`
Location of the fasta formatted plasmid database
```bash
$ --plasmid_db /path/to/your/plasmid/database/*.fa
```

##### `--threads`
Number of threads to use for each process (where applicable)
```bash
$ --threads 10
```

##### `--output`
Location for output files to be written to. Note, the output directory does not need to exist prior to pipeline execution
```bash
$ --output /path/to/where/you/want/your/output/files
```

## Dependencies

Tychus utilizes a number of open source projects to run:

* [Nextflow](https://www.nextflow.io/index.html) - Workflow framework
* [Docker](https://www.docker.com/what-docker) - Software containerization platform
* [Trimmomatic](https://github.com/timflutre/trimmomatic) - Read trimmer
* [Bowtie2](https://github.com/BenLangmead/bowtie2) - Sequence aligner
* [Samtools](https://github.com/samtools/samtools) - Alignment processor
* [Freebayes](https://github.com/ekg/freebayes) - Variant caller
* [Coverage Sampler](https://github.com/cdeanj/coverage_sampler) - Rarefaction calculator
* [kSNP3](https://sourceforge.net/projects/ksnp/) - KMER selection

