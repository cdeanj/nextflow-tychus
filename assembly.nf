#!/usr/bin/env nextflow

params.read_pairs = "/home/chris_dean/nextflow/assembly_pipeline/tutorial/reads/test/test-data/*_R{1,2}_001.fastq"
params.threads = 1
results_path = './results'

threads = params.threads

Channel
	.fromFilePairs(params.read_pairs, flat: true)
	.into { abyss_read_pairs; velvet_read_pairs; spades_read_pairs; idba_read_pairs; kmer_genie_read_pairs }

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

best_abyss_kmer_results.map{ id, file -> [id, file.text.trim()] } .set { best_abyss_k } 
best_velvet_kmer_results.map{ id, file -> [id, file.text.trim()] } .set { best_velvet_k }

process run_abyss_assembly {
	maxForks 1

        input:
        set dataset_id, file(forward), file(reverse) from abyss_read_pairs
	val best from best_abyss_k

        output:
        set dataset_id, file("${dataset_id}_abyss-contigs.fa") into abyss_assembly_results

	"""
        abyss-pe k=${best[1]} name=abyss j=${threads} in='${forward} ${reverse}'
        mv abyss-contigs.fa ${dataset_id}_abyss-contigs.fa
	"""
}

process run_velvet_assembly {
	maxForks 1

	input:
	set dataset_id, file(forward), file(reverse) from velvet_read_pairs
	val best from best_velvet_k

	output:
	set dataset_id, file("${dataset_id}_velvet-contigs.fa") into velvet_assembly_results

	"""
	velveth auto ${best[1]} -fastq -shortPaired ${forward} -fastq -shortPaired2 ${reverse}
	velvetg auto
	mv auto/contigs.fa ${dataset_id}_velvet-contigs.fa
	"""
}

process run_spades_assembly {
	maxForks 1

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
	maxForks 1

	input:
	set dataset_id, file(forward), file(reverse) from idba_read_pairs

	output:
	set dataset_id, file("${dataset_id}_idba-contigs.fa") into idba_assembly_results

	"""
	fq2fa --merge --filter ${forward} ${reverse} ${dataset_id}_idba-paired-contigs.fa
	idba_ud -r ${dataset_id}_idba-paired-contigs.fa --num_threads ${threads} -o ${dataset_id}_idba_output
	mv ${dataset_id}_idba_output/contig.fa ${dataset_id}_idba-contigs.fa	
	"""
}

process run_cisa_contig_integrator {
	maxForks 1

	publishDir "$results_path/assembly_contigs"

	input:
	set dataset_id, file(spades_contigs) from spades_assembly_results
	set dataset_id, file(idba_contigs) from idba_assembly_results
	set dataset_id, file(abyss_contigs) from abyss_assembly_results
	set dataset_id, file(velvet_contigs) from velvet_assembly_results

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
	echo CISA=/CISA1.2 >> CISA.config
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
