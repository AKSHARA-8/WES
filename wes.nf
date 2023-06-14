#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
BED file: ${params.bedfile}.bed
Sequences in:${params.sequences}
"""
process trimming_trimmomatic {
 
        input:
			val Sample
        output:
                        tuple val (Sample), file("*1P.fq.gz"), file("*2P.fq.gz")
        script:
        """
        ${params.trimmomatic_path}trimmomatic PE \
        ${params.sequences}/${Sample}_*R1_*.fastq.gz ${params.sequences}/${Sample}_*R2_*.fastq.gz \
        -baseout ${Sample}.fq.gz \
        ILLUMINACLIP:${params.adaptors}:2:30:10:2:keepBothReads \
        LEADING:3 SLIDINGWINDOW:4:15 MINLEN:40
        sleep 5s
        """
}
process pair_assembly_pear {
	
        input:
			tuple val (Sample), file(paired_forward), file(paired_reverse)

        output:
			tuple val (Sample), file("*.assembled.fastq")

        script:
        """
        ${params.pear_path} -f ${paired_forward} -r ${paired_reverse} -o ${Sample} -n 53 -j 25
	
        """
}
process mapping_reads {
	input:
			tuple val (Sample), file(pairAssembled)

	output: 
		        tuple val (Sample), file ("*.sam")
	
	script:
        """
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200" -M -t 20 ${params.genome} ${pairAssembled} > ${Sample}.sam 
        """     
}
process sam_conversion{
	input:
			tuple val (Sample), file(samfile)
	
	output:
			tuple val (Sample), file("*.sorted.bam"), file("*.sorted.bam.bai")

	script:
	"""
	${params.samtools} view -bT ${params.genome} ${samfile} > ${Sample}.bam 
	${params.samtools} sort ${Sample}.bam > ${Sample}.sorted.bam
        ${params.samtools} index ${Sample}.sorted.bam > ${Sample}.sorted.bam.bai

	"""
}
workflow WES {

    Channel
                .fromPath(params.input)
                .splitCsv(header:false)
                .flatten()
                .map{ it }
                .set { samples_ch }

        main:
        trimming_trimmomatic(samples_ch) | pair_assembly_pear | mapping_reads | sam_conversion

}

