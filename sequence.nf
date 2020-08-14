#!/usr/bin/env nextflow 


// To see all intermediate files, uncomment the publishDir lines 
// in each process


// set to corresponding directory paths before using
work_dir    = "/gpfs/loomis/project/sherman/hl723/nf"
data_dir    = "${work_dir}/sfBacteria/fastq"
ref_dir     = "${work_dir}/sfBacteria/reference"
ref         = "${ref_dir}/PA14_genome.fa"
ref_name    = "PA14_genome"


// initial channel that emits fastq read pairs
// the sampleId emitted and carried on by this channel is actually 
// the largest matching group of the basename and may not be as 
// you may expect
Channel.fromFilePairs("${data_dir}/*_L???_R{1,2}_???.fastq").set{fastq_ch}


// index reference before proceeding with the pipeline
// rm commands are to clean the directory of other extraneous files 
// before indexing
process prepareRef {
	label "grace" 

	publishDir "${ref_dir}"

	output:
	val true into indexedref_ch

	"""
	module load BWA/0.7.17-foss-2018a
	module load GATK/4.1.8.1-Java-1.8
	module load SAMtools/1.7-foss-2018a

	rm -f ${ref_dir}/${ref_name}.fai
	rm -f ${ref_dir}/${ref_name}.d*
	rm -f ${ref_dir}/${ref_name}.fa.*

	bwa index ${ref} 
	samtools faidx ${ref}
	gatk CreateSequenceDictionary -R ${ref}
	"""
}


// convert pairs of reads into sam files
process reads_to_sam {
	label "grace" 

	// publishDir "${work_dir}/sam/"

	input: 
	set sampleId, file(reads) from fastq_ch
	val x from indexedref_ch

	output:
	set sampleId, file("${sampleId}.sam") into sam_ch

	"""
	module load BWA/0.7.17-foss-2018a

	bwa mem ${ref} $reads > ${sampleId}.sam
	""" 
}


// convert each sam to bam file 
process sam_to_bam {
	label "grace" 

	// publishDir "${work_dir}/bam/"

	input: 
	set sampleId, file(input_sam) from sam_ch

	output:
	set sampleId, file("${sampleId}.bam") into bam_ch

	"""
	module load SAMtools/1.7-foss-2018a

	samtools view -bu ${input_sam} | samtools sort -o ${sampleId}.bam
	""" 
}


// get picard.jar with wget
process get_picard {
	label "grace"

	output:
	file "picard.jar" into picard_ch

	"""
	wget https://github.com/broadinstitute/picard/releases/download/2.23.3/picard.jar
	"""
}


// mark duplicates for the bam files
process mark_duplicates {
	label "grace" 

	// publishDir "${work_dir}/bam_dup/"

	input:
	file "picard.jar" from picard_ch
	set sampleId, file(input_bam) from bam_ch

	output:
	set sampleId, file("${sampleId}_marked.bam") into dup_ch

	"""
	java -jar picard.jar MarkDuplicates \
      I=${input_bam} \
      O=${sampleId}_marked.bam \
      M=${sampleId}_marked_dup_metrics.txt
    """

}


// add read groups with the tags specified
// may want to change for different projects
process add_read_groups {
	label "grace" 

	// publishDir "${work_dir}/bam_dup_rg/"

	input:
	file "picard.jar" from picard_ch
	set sampleId, file(input_bam) from dup_ch

	output: 
	set sampleId, file("${sampleId}_final.bam") into bam_rg_ch

	"""
	java -jar picard.jar AddOrReplaceReadGroups \
      I=${input_bam} \
      O=${sampleId}_final.bam \
      RGID=${sampleId} \
      RGSM=${sampleId} \
      RGPL=ILLUMINA \
      RGLB=${sampleId} \
      RGPU=FLOWCELL1.LANE1.${sampleId}
    """
}


// indexes each bam file and calls gatk HaplotypeCaller
// the publishDir command converts the filename to the actual sample name
process haploTypeCaller {
	label "grace" 

	// publishDir "${work_dir}/g.vcf/", saveAs: {filename -> "${filename[0..-15]}.g.vcf.gz"}

	input:
	set sampleId, file(input_bam) from bam_rg_ch

	output:
	file("${sampleId}.g.vcf.gz") into gvcf_ch

	"""
	module load GATK/4.1.8.1-Java-1.8
	module load SAMtools/1.7-foss-2018a

	samtools index ${input_bam} 

	gatk --java-options "-Xmx4g" HaplotypeCaller  \
   		-R ${ref} \
   		-I ${input_bam}  \
   		-O ${sampleId}.g.vcf.gz \
   		-ERC GVCF
	"""
}


// combine all the gvcf files together
// and converts the output channel to a string --variant [filename]
// so that it can be directly used in the command after joining each argument
// the .tbi file is needed in the next process
process combineGVCFs {
	label "grace_more_time" 

	// publishDir "${work_dir}/combined_g.vcf/"
	
	input:
	val vcfs from gvcf_ch.map{file -> "--variant ${file}"}.collect()

	output:
	file("cohort.g.vcf.gz") into combined_vcfs_ch
	file("*.tbi") into vcfs_index_ch

	
	"""
	module load GATK/4.1.8.1-Java-1.8

	gatk CombineGVCFs \
   		-R ${ref} \
   		${vcfs.join(" ")} \
	    -O cohort.g.vcf.gz
	"""
}


// takes the combined file and genotypes it
// stores the final file in the vcf.gz folder in your work directory
// *** BE SURE TO CLEAN THE vcf.gz DIRECTORY FOR THIS PROCESS TO RUN ***
process genotypeGVCFs {
	label "grace_more_time" 

	storeDir "${work_dir}/vcf.gz/"
	
	input:
	file(cohort) from combined_vcfs_ch
	file(index) from vcfs_index_ch

	output:
	file("final.vcf.gz") into genotype_vcfs_ch

	"""
	module load GATK/4.1.8.1-Java-1.8

	gatk --java-options "-Xmx4g" GenotypeGVCFs \
   		-R ${ref} \
   		-V ${cohort} \
   		-O final.vcf.gz
	"""
}