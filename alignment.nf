#!/usr/bin/env nextflow

params.genome = "/home/chris_dean/nextflow/assembly_pipeline/tutorial/listeriadb/listeriadb.fa"
params.amr_db = "/home/chris_dean/nextflow/assembly_pipeline/tutorial/amrdb/amrdb.fa"
params.vf_db = "/home/chris_dean/nextflow/assembly_pipeline/tutorial/vfdb/vfdb.fa"
params.plasmid_db = "/home/chris_dean/nextflow/assembly_pipeline/tutorial/plasmiddb/plasmiddb.fa"
params.threads = 1
params.keep = false
params.outdir = "${PWD}/alignment_results"

genome = file(params.genome)
amr_db = file(params.amr_db)
vf_db = file(params.vf_db)
plasmid_db = file(params.plasmid_db)
threads = params.threads

Channel
        .fromFilePairs("/home/chris_dean/nextflow/assembly_pipeline/tutorial/reads/test/test-data/*_R{1,2}_001.fastq", flat: true)
        .into { genome_read_pairs; amr_read_pairs; vf_read_pairs; plasmid_read_pairs }

process build_genome_index {
	input:
	file genome

	output:
	file 'genome.index*' into genome_index

	"""
	bowtie2-build --threads ${threads} $genome genome.index
	"""
}

process build_amr_index {
	input:
        file amr_db

        output:
        file 'amr.index*' into amr_index

        """
        bowtie2-build --threads ${threads} $amr_db amr.index
	"""
}

process build_vf_index {
	input:
        file vf_db

        output:
        file 'vf.index*' into vf_index

        """
        bowtie2-build --threads ${threads} $vf_db vf.index
	"""
}

/*process build_plasmid_index {
	input:
        file plasmid_db

        output:
        file 'plasmid.index*' into plasmid_index

        """
        bowtie2-build --threads ${threads} $plasmid_db plasmid.index
	"""
}*/

process bowtie2_genome_alignment {
	maxForks 1

	input:
	set dataset_id, file(forward), file(reverse) from genome_read_pairs
	file index from genome_index.first()

	output:
	set dataset_id, file("${dataset_id}_genome_alignment.sam") into genome_sam_files
	set dataset_id, file("${dataset_id}_genome_alignment.bam") into genome_bam_files
        set dataset_id, file("${dataset_id}_genome_alignment.bai") into genome_index_files

	"""
	bowtie2 --threads ${threads} -x genome.index -1 $forward -2 $reverse -S ${dataset_id}_genome_alignment.sam
	samtools view -bS ${dataset_id}_genome_alignment.sam | samtools sort -@ ${threads} - ${dataset_id}_genome_alignment
        samtools index ${dataset_id}_genome_alignment.bam ${dataset_id}_genome_alignment.bai
	"""
}

process bowtie2_amr_alignment {
	maxForks 1

	input:
	set dataset_id, file(forward), file(reverse) from amr_read_pairs
	file index from amr_index.first()

	output:
	set dataset_id, file("${dataset_id}_amr_alignment.sam") into amr_sam_files

	"""
	bowtie2 --threads ${threads} -x amr.index -1 $forward -2 $reverse -S ${dataset_id}_amr_alignment.sam
	"""
}

process bowtie2_vfdb_alignment {
	maxForks 1

        input:
        set dataset_id, file(forward), file(reverse) from vf_read_pairs
        file index from vf_index.first()

        output:
        set dataset_id, file("${dataset_id}_vf_alignment.sam") into vf_sam_files

        """
        bowtie2 --threads ${threads} -x vf.index -1 $forward -2 $reverse -S ${dataset_id}_vf_alignment.sam
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
	input:
	set dataset_id, file(bam) from genome_bam_files
	set dataset_id, file(bai) from genome_index_files
	file genome

	output:
	file("${dataset_id}_genome_variants.vcf.gz") into compressed_variants
	file("${dataset_id}_genome_variants.vcf.gz.tbi") into indexed_variants

	"""
	freebayes -p 1 -f ${genome} $bam | bgzip -c > ${dataset_id}_genome_variants.vcf.gz
	tabix ${dataset_id}_genome_variants.vcf.gz
	"""
}

process bcftools_isec_and_consensus {
	input:
	file genome
	file(vcf) from compressed_variants.toList()
	file(idx) from indexed_variants.toList()

	output:
	file('consensus.fa') into consensus_sequence

	"""
	bcftools isec -p dir --nfiles=3 -w1 $vcf
	mv dir/0000.vcf isect.vcf
	bgzip isect.vcf
	tabix isect.vcf.gz
	cat $genome | bcftools consensus isect.vcf.gz > consensus.fa
	"""
}

/*process genome_coverage_sampler {
	input:
	set dataset_id, file(genome_sam_alignment) from genome_sam_files
	file genome

	output:
	set dataset_id, file('coverage_sampler_genome.tab') into genome_csa_files

	"""
	csa -ref_fp $genome -sam_fp $genome_sam_alignment -min 100 -max 100 -skip 5 -t 80 -samples 1 -out_fp coverage_sampler_genome.tab
	"""
}*/

/*process amr_coverage_sampler {
	input:
	set dataset_id, file(amr_sam_alignment) from amr_sam_files
	file amrdb from amr_db

	output:
	set dataset_id, file('coverage_sampler_amr.tab') into amr_csa_files

	"""
	csa -ref_fp $amrdb -sam_fp $amr_sam_alignment -min 100 -max 100 -skip 5 -t 80 -samples 1 -out_fp coverage_sampler_amr.tab
	"""
}*/

/*process vf_coverage_sampler { 
	input:
	set dataset_id, file(vf_sam_alignment) from vf_sam_files
	file vfdb from vf_db

	output:
	set dataset_id, file('coverage_sampler_vf.tab') into vf_csa_files

	"""
	csa -ref_fp $vfdb -sam_fp $vf_sam_alignment -min 100 -max 100 -skip 5 -t 80 -samples 1 -out_fp coverage_sampler_vf.tab
	"""
}*/

def extractSampleName(s) {
        ret = s =~ /\/(.+)_R/;
        basepath = ~/.+\//
        return ret[0][1] - basepath;
}
