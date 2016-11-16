#!/usr/bin/env nextflow

params.pair1 = "/home/chris_dean/nextflow/tychus/tutorial/reads/SFBRL001-M3237-15-001_S1_L001_R1_001.fastq"
params.pair2 = "/home/chris_dean/nextflow/tychus/tutorial/reads/SFBRL001-M3237-15-001_S1_L001_R2_001.fastq"
params.threads = 1
params.outdir = "${PWD}/results"
results_path = './results'

threads = params.threads

forward_reads = Channel
		.fromPath(params.pair1)
		.map { path -> [ path.toString().replace('_R1', '_RX'), path ] }

reverse_reads = Channel
		.fromPath(params.pair2)
		.map { path -> [ path.toString().replace('_R2', '_RX'), path ] }

params.read_pairs = forward_reads
	     .phase(reverse_reads)
             .map { pair1, pair2 -> [ extractSampleName(pair1[1]), pair1[1], pair2[1] ] }

params.read_pairs.into {
	abyss_read_pairs;
	velvet_read_pairs;
	spades_read_pairs;
	idba_read_pairs;
	kmer_genie_read_pairs
}

process run_kmer_genie {
	input:
	set dataset_id, file(forward), file(reverse) from kmer_genie_read_pairs

	output:
	set dataset_id, file("${dataset_id}_best-k.txt") into (best_abyss_kmer_results, best_velvet_kmer_results)

	"""
	echo $forward > ${dataset_id}_read_pair_list.txt
	echo $reverse >> ${dataset_id}_read_pair_list.txt
	kmergenie "${dataset_id}_read_pair_list.txt" -t ${threads} | tail -n 1 | awk '{print \$3}' > ${dataset_id}_best-k.txt
	"""
}


// Thanks Paolo
best_abyss_kmer_results.map{ id, file -> [id, file.text.trim()] } .set { best_abyss_k } 
best_velvet_kmer_results.map{ id, file -> [id, file.text.trim()] } .set { best_velvet_k }

process run_abyss_assembly {
        input:
        set dataset_id, file(forward), file(reverse) from abyss_read_pairs
	val best from best_abyss_k

        output:
        set dataset_id, file("${dataset_id}_abyss-contigs.fa") into abyss_assembly_results

        """
        abyss-pe k=${best[1]} name=abyss np=${threads} in='${forward} ${reverse}'
	mv abyss-contigs.fa ${dataset_id}_abyss-contigs.fa
        """
}

process run_velvet_assembly {
	input:
	set dataset_id, file(forward), file(reverse) from velvet_read_pairs
	val best from best_velvet_k

	output:
	set dataset_id, file("${dataset_id}_velvet-contigs.fa") into velvet_assembly_results
	
	"""
	velveth auto ${best[1]} -fastq -shortPaired ${forward} -fastq -shortPaired2 ${reverse}
	velvetg auto -exp_cov auto
	mv auto/contigs.fa ${dataset_id}_velvet-contigs.fa
	"""
}


process run_spades_assembly {
	input:
	set dataset_id, file(forward), file(reverse) from spades_read_pairs

	output:
	set dataset_id, file("${dataset_id}_spades-contigs.fa") into spades_assembly_results

	"""
	spades.py --pe1-1 ${forward} --pe1-2 ${reverse} -t ${threads} -o spades_output
	mv spades_output/contigs.fasta ${dataset_id}_spades-contigs.fa
	"""
}

process run_idba_assembly {
	input:
	set dataset_id, file(forward), file(reverse) from idba_read_pairs

	output:
	set dataset_id, file("${dataset_id}_idba-contigs.fa") into idba_assembly_results

	"""
	fq2fa --merge --filter ${forward} ${reverse} ${dataset_id}_idba-paired-contigs.fa
	idba_ud -r ${dataset_id}_idba-paired-contigs.fa --num_threads ${threads} -o idba_output
	mv idba_output/contig.fa ${dataset_id}_idba-contigs.fa	
	"""
}

process run_cisa_contig_integrator {
	publishDir "$results_path/assembly_contigs"

	input:
	set dataset_id, spades_contigs from spades_assembly_results
	set dataset_id, idba_contigs from idba_assembly_results
	set dataset_id, abyss_contigs from abyss_assembly_results
	set dataset_id, velvet_contigs from velvet_assembly_results

	output:
	set dataset_id, file("${dataset_id}_master_contigs.fa") into cisa_merged_contigs
	set dataset_id, file("${dataset_id}_master_integrated_contigs.fa") into cisa_integrated_contigs

	shell:
	'''
	#!/bin/sh
	echo count=4 > Merge.config
	echo data=!{spades_contigs},title=SPades >> Merge.config
	echo data=!{idba_contigs},title=IDBA >> Merge.config
	echo data=!{abyss_contigs},title=Abyss >> Merge.config
	echo data=!{velvet_contigs},title=Velvet >> Merge.config
	echo Master_file=!{dataset_id}_master_contigs.fa >> Merge.config
	Merge.py Merge.config
	echo genome=`grep 'Whole Genome' 'Merge_info' | cut -d ':' -f2 | sort -rn | head -n 1 | tr -d [:space:]` > CISA.config
	echo infile=!{dataset_id}_master_contigs.fa >> CISA.config
	echo outfile=!{dataset_id}_master_integrated_contigs.fa >> CISA.config
	echo nucmer=`which nucmer` >> CISA.config
	echo R2_Gap=0.95 >> CISA.config
	echo CISA=/home/chris_dean/nextflow/assembly_pipeline/assemblers/CISA1.2 >> CISA.config
	echo makeblastdb=`which makeblastdb` >> CISA.config
	echo blastn=`which blastn` >> CISA.config
	CISA.py CISA.config
	'''
}

process run_prokka_annotation {
	publishDir "$results_path/annotations"

	input:
	set dataset_id, file(cisa_contigs) from cisa_integrated_contigs
	
	output:
	file("${dataset_id}.*") into prokka_annotations

	"""
	prokka ${cisa_contigs} --prefix ${dataset_id} --outdir annotations
	mv annotations/* .
	"""
}

def extractSampleName(s) {
	ret = s =~ /\/(.+)_R/;
	basepath = ~/.+\//
  	return ret[0][1] - basepath;
}
