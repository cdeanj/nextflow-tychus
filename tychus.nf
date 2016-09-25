#!/usr/bin/env nextflow

params.pair1 = "/home/chris_dean/nextflow/tychus/tutorial/listeria_reads/*_R1*.fastq"
params.pair2 = "/home/chris_dean/nextflow/tychus/tutorial/listeria_reads/*_R2*.fastq"
params.amr_db = "/home/chris_dean/nextflow/tychus/tutorial/amrdb/amrdb.fa"
params.vf_db = "/home/chris_dean/nextflow/tychus/tutorial/vfdb/vfdb.fa"
params.plasmid_db = "/home/chris_dean/nextflow/tychus/tutorial/plasmiddb/plasmiddb.fa"
params.threads = 1

log.info "T Y C H U S - NF ~ version 1.0.0"
log.info "================================"
log.info "AMR Database       : ${params.amr_db}"
log.info "Virulence Database : ${params.vf_db}"
log.info "Plasmid Database   : ${params.plasmid_db}"
log.info "Forward Reads      : ${params.pair1}"
log.info "Reverse Reads      : ${params.pair2}"
log.info "Number of Threads  : ${params.threads}"
log.info "\n"

amr_db = file(params.amr_db)
vf_db = file(params.vf_db)
plasmid_db = file(params.plasmid_db)
threads = params.threads

if(!amr_db.exists()) {
        exit 1, "Unable to find genome file: {params.amr_db}"
}
if(!vf_db.exists()) {
        exit 1, "Unable to find virulence file: {params.vf_db}"
}
if(!plasmid_db.exists()) {
        exit 1, "Unable to find plasmid file: {params.plasmid_db}"
}

forward_reads = Channel
		.fromPath(params.pair1)
		.map { path -> [ path.toString().replace('_R1', '_RX'), path ] }

reverse_reads = Channel
		.fromPath(params.pair2)
		.map { path -> [ path.toString().replace('_R2', '_RX'), path ] }

params.read_pairs = forward_reads
	     .phase(reverse_reads)
             .map { pair1, pair2 -> [ pathToDatasetID(pair1[1]), pair1[1], pair2[1] ] }

params.read_pairs.into {
	read_files_amr;
	read_files_vf;
	read_files_plasmid
}

process build_amr_index {
	input:
	file amr_db

	output:
	file 'amr.index*' into amr_index

	"""
	bowtie2-build ${amr_db} amr.index
	"""
}

process build_vfdb_index {
        input:
        file vf_db

        output:
        file 'vf.index*' into vf_index

        """
        bowtie2-build ${vf_db} vf.index
        """
}

process build_plasmiddb_index {
	input:
	file plasmid_db

	output:
	file 'plasmid.index*' into plasmid_index

	"""
	bowtie2-build ${plasmid_db} plasmid.index
	"""
}

process bowtie_amr_alignment {
	input:
	set dataset_id, file(forward), file(reverse) from read_files_amr
	file index from amr_index.first()

	output:
	set dataset_id, file('amr_alignment.sam') into amr_sam_files

	"""
	bowtie2 -p ${threads} -x amr.index -1 $forward -2 $reverse -S amr_alignment.sam
	"""
}

process bowtie2_vfdb_alignment {
        input:
        set dataset_id, file(forward), file(reverse) from read_files_vf
        file index from vf_index.first()

        output:
        set dataset_id, file('vf_alignment.sam') into vf_sam_files

        """
        bowtie2 -p ${threads} -x vf.index -1 $forward -2 $reverse -S vf_alignment.sam
        """
}

process bowtie2_plasmid_alignment {
	input:
	set dataset_id, file(forward), file(reverse) from read_files_plasmid
	file index from plasmid_index.first()

	output:
	set dataset_id, file('plasmid_alignment.sam') into plasmid_sam_files

	"""
	bowtie2 -p ${threads} -x plasmid.index -1 $forward -2 $reverse -S plasmid_alignment.sam
	"""
}

process amr_coverage_sampler {
	input:
	set dataset_id, file(amr_sam_alignment) from amr_sam_files
	file amrdb from amr_db

	output:
	set dataset_id, file('coverage_sampler_amr.tab') into amr_csa_files

	"""
	csa -ref_fp $amrdb -sam_fp $amr_sam_alignment -min 100 -max 100 -skip 5 -t 80 -samples 1 -out_fp coverage_sampler_amr.tab
	"""
}

process vf_coverage_sampler {
	input:
	set dataset_id, file(vf_sam_alignment) from vf_sam_files
	file vfdb from vf_db

	output:
	set dataset_id, file('coverage_sampler_vf.tab') into vf_csa_files

	"""
	csa -ref_fp $vfdb -sam_fp $vf_sam_alignment -min 100 -max 100 -skip 5 -t 80 -samples 1 -out_fp coverage_sampler_vf.tab
	"""
}

def pathToDatasetID(path) {
  	return path.getParent().toString();
}
