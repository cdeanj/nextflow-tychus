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

//General configuation variables
params.help = ""
params.pwd = "$PWD"
params.output = "tychus_alignment_output"
params.work_dir = "$PWD/temporary_files"
params.read_pairs = "tutorial/raw_sequence_data/*_R{1,2}_001.fastq"
params.genome = "tutorial/genome_reference/listeriadb.fa"
params.amr_db = "tutorial/amr_reference/megaresdb.fa"
params.vf_db = "tutorial/virulence_reference/virulencedb.fa"
params.plasmid_db = "tutorial/plasmid_reference/plasmiddb.fa"
params.out_dir = params.pwd + "/" + params.output
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
params.ML = ""
params.NJ = ""
params.min_frac = 0.75

ML = params.ML
NJ = params.NJ
min_frac = params.min_frac

// Figtree configuration variables
params.PNG = ""
params.JPEG = ""
params.PDF = ""
params.SVG = ""

PNG = params.PNG
JPEG = params.JPEG
PDF = params.PDF
SVG = params.SVG

// Display help message
if(params.help) {
	log.info ''
	log.info 'Tychus - Alignment Pipeline'
	log.info ''
	log.info 'Usage: '
	log.info '    nextflow run alignment.nf -profile alignment [options]'
	log.info ''
	log.info 'General Options: '
	log.info '    --read_pairs      DIR		Directory of paired FASTQ files'
	log.info '    --genome          FILE		Path to the FASTA formatted reference database'
	log.info '    --amr_db          FILE		Path to the FASTA formatted resistance database'
	log.info '    --vf_db           FILE		Path to the FASTA formatted virulence database'
	log.info '    --plasmid_db      FILE		Path to the FASTA formatted plasmid database'
	log.info '    --threads         INT		Number of threads to use for each process'
	log.info '    --out_dir         DIR		Directory to write output files to'
	log.info ''
	log.info 'Trimmomatic Options: '
	log.info '    --leading         INT		Remove leading low quality or N bases'
	log.info '    --trailing        INT		Remove trailing low quality or N bases'
	log.info '    --slidingwindow   INT		Scan read with a sliding window'
	log.info '    --minlen          INT		Drop reads below INT bases long'
	log.info ''
	log.info 'kSNP Options: '
	log.info '    --ML              BOOL		Estimate maximum likelihood tree'
	log.info '    --NJ              BOOL		Estimate neighbor joining tree'
	log.info '    --min_frac        DECIMAL		Minimum fraction of genomes with locus'
	log.info ''
	log.info 'Figtree Options: '
	log.info '    --JPEG            BOOL		Convert newick tree to annotated JPEG'
	log.info '    --PDF             BOOL		Convert newick tree to annotated PDF'
	log.info '    --PNG             BOOL		Convert newick tree to annotated PNG'
	log.info '    --SVG             BOOL		Convert newick tree to annotated SVG'
	log.info ''
	return
}


// Let's group the read pairs and place them into a channel
Channel
        .fromFilePairs(params.read_pairs, flat: true)
        .into { trimmomatic_read_pairs }

process BuildGenomeIndex {
	tag { "${genome.baseName}" }

	input:
	file genome

	output:
	file 'genome.index*' into genome_index

	"""
	bowtie2-build $genome genome.index
	"""
}

process BuildAMRIndex {
	tag { "${amr_db.baseName}" }

	input:
        file amr_db

        output:
        file 'amr.index*' into amr_index

        """
        bowtie2-build $amr_db amr.index
	"""
}

process BuildVFIndex {
	tag { "${vf_db.baseName}" }

	input:
        file vf_db

        output:
        file 'vf.index*' into vf_index

        """
        bowtie2-build $vf_db vf.index
	"""
}

/*process BuildPlasmidIndex {
	tag { "${plasmid_db.baseName}" }

	input:
        file plasmid_db

        output:
        file 'plasmid.index*' into plasmid_index

        """
        bowtie2-build $plasmid_db plasmid.index
	"""
}*/

process RunQC {
	publishDir "${params.out_dir}/Preprocessing", mode: "copy"

	tag { dataset_id }

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

process GenomeAlignment {
	publishDir "${params.out_dir}/Genome_Alignment", mode: "copy"

	tag { dataset_id }

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

process AMRAlignment {
	publishDir "${params.out_dir}/AMR_Alignment", mode: "copy"

	tag { dataset_id }

	input:
	set dataset_id, file(forward), file(reverse) from amr_read_pairs
	file index from amr_index.first()

	output:
	set dataset_id, file("${dataset_id}_amr_alignment.sam") into amr_sam_files

	"""
	bowtie2 -p ${threads} -x amr.index -1 $forward -2 $reverse -S ${dataset_id}_amr_alignment.sam
	"""
}

process VFAlignment {
	publishDir "${params.out_dir}/Virulence_Alignment", mode: "copy"

	tag { dataset_id }

        input:
        set dataset_id, file(forward), file(reverse) from vf_read_pairs
        file index from vf_index.first()

        output:
        set dataset_id, file("${dataset_id}_vf_alignment.sam") into vf_sam_files

        """
        bowtie2 -p ${threads} -x vf.index -1 $forward -2 $reverse -S ${dataset_id}_vf_alignment.sam
        """
}

/*process PlasmidAlignment {
	publishDir "${params.out_dir}/Plasmid_Alignment", mode: "copy"

	tag { dataset_id }

	input:
	set dataset_id, file(forward), file(reverse) from plasmid_read_pairs
	file index from plasmid_index.first()

	output:
	set dataset_id, file("${dataset_id}_plasmid_alignment.sam") into plasmid_sam_files

	"""
	bowtie2 -p ${threads} -x plasmid.index -1 $forward -2 $reverse -S ${dataset_id}_plasmid_alignment.sam
	"""
}*/

process BuildConesnsusSequence {
	tag { dataset_id }

	publishDir "${params.out_dir}/Consensus", mode: "copy"

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
	echo -e "$params.out_dir/Consensus/${dataset_id}_consensus.fa\t$dataset_id" >> ${dataset_id}_in_list.txt
	"""
}

process PreparePhylogeneticAnalysis {
	tag { "genome_configuration" }

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

process BuildPhylogenies {
	publishDir "${params.out_dir}/Trees", mode: "copy"

	tag { "configuration_files" }

	input:
	file kchooser_config from kchooser_configuration
	file ksnp3_config from ksnp3_configuration.toList()

	output:
	file("trees/*.tre") into phylogenetic_trees
	file("polymorphisms/*") into polymorphisms

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
	mkdir trees
	mkdir polymorphisms
	mv kSNP3_results/*.tre trees
	mv kSNP3_results/* polymorphisms
	rm -rf !{params.work_dir}
	'''
}

process ConvertNewickToPDF {
	publishDir "${params.out_dir}/Phylogenetic_Tree_Images", mode: "copy"

	input:
	file tree from phylogenetic_trees.flatten()

	output:
	file "${base}*"

	script:
	base = tree.baseName

	shell:
        '''
        #!/bin/sh
        if [ !{PDF} ]
        then
                java -jar /figtree/lib/figtree.jar -graphic PDF !{tree} !{base}.pdf
        elif [ !{PNG} ]
        then
                java -jar /figtree/lib/figtree.jar -graphic PNG !{tree} !{base}.png
        elif [ !{JPEG} ]
        then
                java -jar /figtree/lib/figtree.jar -graphic JPEG !{tree} !{base}.jpg
        else
                java -jar /figtree/lib/figtree.jar -graphic SVG !{tree} !{base}.svg
        fi
        '''
}

// Display information about the completed run
// See https://www.nextflow.io/docs/latest/metadata.html for more
// information about available onComplete options
workflow.onComplete {
	log.info "Nextflow Version:	$workflow.nextflow.version"
  	log.info "Command Line:		$workflow.commandLine"
	log.info "Container:		$workflow.container"
	log.info "Duration:		$workflow.duration"
	log.info "Output Directory:	$params.out_dir"
}
