#!/usr/bin/env nextflow

//General configuation variables
params.work_dir = "$PWD/temporary_files"
params.read_pairs = "tutorial/raw_sequence_data/*_R{1,2}_001.fastq"
params.genome = "tutorial/genome_reference/listeriadb.fa"
params.amr_db = "tutorial/amr_reference/megaresdb.fa"
params.vf_db = "tutorial/virulence_reference/virulencedb.fa"
params.plasmid_db = "tutorial/plasmid_reference/plasmiddb.fa"
params.threads = 1

genome = file(params.genome)
amr_db = file(params.amr_db)
vf_db = file(params.vf_db)
plasmid_db = file(params.plasmid_db)
threads = params.threads

// Trimmomatic configuration variables
params.leading = 3
params.trailing = 3
params.slidingwindow = "4:15"
params.minlen = 36

leading = params.leading
trailing = params.trailing
slidingwindow = params.slidingwindow
minlen = params.minlen

// kSNP3 configuration variables
params.ML = true
params.NJ = false
params.min_frac = 0.75

ML = params.ML
NJ = params.NJ
min_frac = params.min_frac


Channel
        .fromFilePairs(params.read_pairs, flat: true)
        .into { trimmomatic_read_pairs }

process build_genome_index {
	input:
	file genome

	output:
	file 'genome.index*' into genome_index

	"""
	bowtie2-build $genome genome.index
	"""
}

process build_amr_index {
	input:
        file amr_db

        output:
        file 'amr.index*' into amr_index

        """
        bowtie2-build $amr_db amr.index
	"""
}

process build_vf_index {
	input:
        file vf_db

        output:
        file 'vf.index*' into vf_index

        """
        bowtie2-build $vf_db vf.index
	"""
}

/*process build_plasmid_index {
	input:
        file plasmid_db

        output:
        file 'plasmid.index*' into plasmid_index

        """
        bowtie2-build $plasmid_db plasmid.index
	"""
}*/

process run_trimmomatic {
	input:
        set dataset_id, file(forward), file(reverse) from trimmomatic_read_pairs

        output:
        set dataset_id, file("${dataset_id}_1P.fastq"), file("${dataset_id}_2P.fastq") into (amr_read_pairs, plasmid_read_pairs, vf_read_pairs, genome_read_pairs)

        """
        java -jar /Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads ${threads} $forward $reverse -baseout ${dataset_id} ILLUMINACLIP:Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10:3:TRUE LEADING:${leading} TRAILING:${trailing} SLIDINGWINDOW:${slidingwindow} MINLEN:${minlen}
        mv ${dataset_id}_1P ${dataset_id}_1P.fastq
        mv ${dataset_id}_2P ${dataset_id}_2P.fastq
        """
}

process bowtie2_genome_alignment {
	input:
	set dataset_id, file(forward), file(reverse) from genome_read_pairs
	file index from genome_index.first()

	output:
	set dataset_id, file("${dataset_id}_genome_alignment.sam") into genome_sam_files
	set dataset_id, file("${dataset_id}_genome_alignment.bam") into genome_bam_files
        set dataset_id, file("${dataset_id}_genome_alignment.bai") into genome_index_files

	"""
	bowtie2 -p ${threads} -x genome.index -1 $forward -2 $reverse -S ${dataset_id}_genome_alignment.sam
	samtools view -bS ${dataset_id}_genome_alignment.sam | samtools sort -@ ${threads} -o ${dataset_id}_genome_alignment.bam
        samtools index ${dataset_id}_genome_alignment.bam ${dataset_id}_genome_alignment.bai
	"""
}

process bowtie2_amr_alignment {
	input:
	set dataset_id, file(forward), file(reverse) from amr_read_pairs
	file index from amr_index.first()

	output:
	set dataset_id, file("${dataset_id}_amr_alignment.sam") into amr_sam_files

	"""
	bowtie2 -p ${threads} -x amr.index -1 $forward -2 $reverse -S ${dataset_id}_amr_alignment.sam
	"""
}

process bowtie2_vfdb_alignment {
        input:
        set dataset_id, file(forward), file(reverse) from vf_read_pairs
        file index from vf_index.first()

        output:
        set dataset_id, file("${dataset_id}_vf_alignment.sam") into vf_sam_files

        """
        bowtie2 -p ${threads} -x vf.index -1 $forward -2 $reverse -S ${dataset_id}_vf_alignment.sam
        """
}

/*process bowtie2_plasmid_alignment {
	input:
	set dataset_id, file(forward), file(reverse) from plasmid_read_pairs
	file index from plasmid_index.first()

	output:
	set dataset_id, file("${dataset_id}_plasmid_alignment.sam") into plasmid_sam_files

	"""
	bowtie2 -p ${threads} -x plasmid.index -1 $forward -2 $reverse -S ${dataset_id}_plasmid_alignment.sam
	"""
}*/

process freebayes_snp_caller {
	storeDir 'temporary_files'

	input:
	set dataset_id, file(bam) from genome_bam_files
	set dataset_id, file(bai) from genome_index_files
	file genome

	output:
	file("${dataset_id}_consensus.fa") into consensus_files
	file("${dataset_id}_in_list.txt") into ksnp3_configuration

	"""
	freebayes -p 1 -f ${genome} $bam | bgzip -c > ${dataset_id}_genome_variants.vcf.gz
	tabix ${dataset_id}_genome_variants.vcf.gz
	cat $genome | bcftools consensus ${dataset_id}_genome_variants.vcf.gz > ${dataset_id}_consensus.fa
	echo -e "$params.work_dir/${dataset_id}_consensus.fa\t$dataset_id" >> ${dataset_id}_in_list.txt
	"""
}

process prepare_ksnp3_configuration {
	storeDir 'temporary_files'

	input:
	file genome

	output:
	file("genome_path.txt") into kchooser_configuration
	file("$genome") into kchooser_genome

	shell:
	'''
	#!/bin/sh
	genome_prefix=`echo !{genome} | cut -f1 -d '.'`
	genome_fp=`readlink !{genome}`
	echo "${genome_fp}\t${genome_prefix}" > genome_path.txt
	'''
}

process run_ksnp3 {
	input:
	file kchooser_config from kchooser_configuration
	file ksnp3_config from ksnp3_configuration.toList()

	shell:
	'''
	#!/bin/sh
	/usr/local/kSNP3/MakeFasta !{kchooser_config} MF.fasta > /dev/null
	/usr/local/kSNP3/Kchooser MF.fasta > /dev/null
	optimum_k=`grep "The optimum" Kchooser.report | tr -dc '0-9'`
	cat !{kchooser_config} > in_list
	cat !{ksnp3_config} >> in_list
	if [ !{ML} && !{NJ} ]
	then
		/usr/local/kSNP3/kSNP3 -in in_list -outdir kSNP3_results -k ${optimum_k} -ML -NJ -core -min_frac !{min_frac} >> /dev/null
	elif [ !{ML} ]
	then
		/usr/local/kSNP3/kSNP3 -in in_list -outdir kSNP3_results -k ${optimum_k} -ML -core -min_frac !{min_frac} >> /dev/null
	else
		/usr/local/kSNP3/kSNP3 -in in_list -outdir kSNP3_results -k ${optimum_k} -NJ -core -min_frac !{min_frac} >> /dev/null
	fi
	'''
}
