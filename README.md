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

## Dependencies

Tychus utilizes a number of open source projects to run:

* [Docker](https://www.docker.com/what-docker) - Software containerization platform
* [Bowtie2](https://github.com/BenLangmead/bowtie2) - Sequence aligner
* [Samtools](https://github.com/samtools/samtools) - Alignment processor
* [Freebayes](https://github.com/ekg/freebayes) - Variant caller
* [Coverage Sampler](https://github.com/cdeanj/coverage_sampler) - Rarefaction calculator
* [kSNP3](https://sourceforge.net/projects/ksnp/) KMER selection

