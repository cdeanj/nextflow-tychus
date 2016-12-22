#!/usr/bin/env nextflow

/*
 * Copyright (c) 2017 Chris Dean

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

// General configuration variables
params.help = false
params.read_pairs = "tutorial/raw_sequence_data/*_R{1,2}_001.fastq"
params.threads = 1

threads = params.threads

if(params.help) {
	log.info ''
	log.info 'Tychus - Assembly Pipeline'
	log.info ''
	log.info 'Usage: '
	log.info '    nextflow run assembly.nf -profile assembly -with-docker [options]'
	log.info ''
	log.info 'General Options: '
	log.info '    --read_pairs      DIR		Directory of paired FASTQ files'
	log.info '    --threads         INT             Number of threads to use for each process'
	log.info '    --output          DIR             Directory to write output files to'
	log.info 'Abyss Options: '
	log.info 'Velvet Options: '
	log.info 'SPades Options: '
	log.info 'IDBA-UD Options: '
	log.info 'Prokka Options: '
	log.info '
	return
}

Channel
	.fromFilePairs(params.read_pairs, flat: true)
	.into { trimmomatic_read_pairs }

process run_trimmomatic {
        input:
        set dataset_id, file(forward), file(reverse) from trimmomatic_read_pairs

        output:
        set dataset_id, file("${dataset_id}_1P.fastq"), file("${dataset_id}_2P.fastq") into (abyss_read_pairs, velvet_read_pairs, spades_read_pairs, idba_read_pairs, kmer_genie_read_pairs)

        """
        java -jar /Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads ${threads} $forward $reverse -baseout ${dataset_id} ILLUMINACLIP:Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10:3:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        mv ${dataset_id}_1P ${dataset_id}_1P.fastq
        mv ${dataset_id}_2P ${dataset_id}_2P.fastq
        """
}

process run_kmer_genie {
	input:
	set dataset_id, file(forward), file(reverse) from kmer_genie_read_pairs

	output:
	file("${dataset_id}_best-k.txt") into (best_abyss_kmer_results, best_velvet_kmer_results)
	set dataset_id, file("${dataset_id}_forward_kg.fq"), file("${dataset_id}_reverse_kg.fq") into (abyss_kg_pairs, velvet_kg_pairs)

	"""
	echo $forward > ${dataset_id}_read_pair_list.txt
	echo $reverse >> ${dataset_id}_read_pair_list.txt
	kmergenie "${dataset_id}_read_pair_list.txt" -t ${threads} | tail -n 1 | awk '{print \$3}' > ${dataset_id}_best-k.txt
	cp $forward "${dataset_id}_forward_kg.fq"
	cp $reverse "${dataset_id}_reverse_kg.fq"
	"""
}

process run_abyss_assembly {
        input:
        set dataset_id, file(forward), file(reverse) from abyss_kg_pairs
	val best from best_abyss_kmer_results

	output:
        set dataset_id, file("${dataset_id}_abyss-contigs.fa") into abyss_assembly_results

	shell:
	'''
	#!/bin/sh
	best_kmer=`cat !{best}`
	abyss-pe k=$best_kmer name=abyss j=!{threads} in='!{forward} !{reverse}'
        mv abyss-contigs.fa !{dataset_id}_abyss-contigs.fa
	'''
}

process run_velvet_assembly {
	input:
	set dataset_id, file(forward), file(reverse) from velvet_kg_pairs
	val best from best_velvet_kmer_results

	output:
	set dataset_id, file("${dataset_id}_velvet-contigs.fa") into velvet_assembly_results

	shell:
	'''
	#!/bin/sh
	best_kmer=`cat !{best}`
	velveth auto $best_kmer -fastq -shortPaired !{forward} -fastq -shortPaired2 !{reverse}
	velvetg auto
	mv auto/contigs.fa !{dataset_id}_velvet-contigs.fa
	'''
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
	idba_ud -r ${dataset_id}_idba-paired-contigs.fa --num_threads ${threads} -o ${dataset_id}_idba_output
	mv ${dataset_id}_idba_output/contig.fa ${dataset_id}_idba-contigs.fa	
	"""
}

process run_cisa_contig_integrator {
	publishDir "${params.out_dir}/assembly_contigs"

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
	publishDir "${params.out_dir}/annotations"

	input:
	set dataset_id, file(cisa_contigs) from cisa_integrated_contigs
	
	output:
	file("${dataset_id}.*") into prokka_annotations

	"""
	prokka ${cisa_contigs} --prefix ${dataset_id} --cpus ${threads} --outdir annotations
	mv annotations/* .
	"""
}
