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

process sam_conversion {
	input:
		tuple val (Sample), file (samfile)
	output:
		tuple val (Sample), file ("*.sorted.bam"), file("*.sorted.bam.bai")
	script:
	"""
	${params.samtools} view -bT ${params.genome} ${samfile} > ${Sample}.bam 
	${params.samtools} sort ${Sample}.bam > ${Sample}.sorted.bam
        ${params.samtools} index ${Sample}.sorted.bam > ${Sample}.sorted.bam.bai
	"""
}

process mark_duplicates {
	input:
		tuple val (Sample), file (sorted_bam), file (sorted_bam_index)
	output:
		tuple val (Sample), file ("*.bam"), file ("*.bam.bai"),  file ("*.txt")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} MarkDuplicates \
	I=${sorted_bam} \
	O=${Sample}_sorted_marked.bam \
	M=${Sample}_picard.info.txt
	REMOVE_DUPLICATES=true 
	${params.samtools} index ${Sample}_sorted_marked.bam > ${Sample}_sorted_marked.bam.bai
	"""
}

process RealignerTargetCreator {
	input:
		tuple val (Sample), file (bamFile), file (bamBai), file (pi_card_info)
	output:
		tuple val (Sample), file ("*.intervals")
	script:
	"""
	 ${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T RealignerTargetCreator -R ${params.genome} -nt 10 -I ${bamFile} --known ${params.site1} -o ${Sample}_target.intervals
	"""
}

process IndelRealigner {
	input:
		tuple val (Sample), file (target_interval), file (bamFile), file (bamBai), file (pi_card_info)
	output:
		tuple val (Sample), file ("*.realigned.bam")
	script:
	"""
	echo ${Sample} ${target_interval} ${bamFile}
        ${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T IndelRealigner -R ${params.genome} -I ${bamFile} -known ${params.site1} --targetIntervals ${target_interval} -o ${Sample}.realigned.bam
	"""
}

process BaseRecalibrator {
        input:
                tuple val (Sample), file (realignedBam)
        output:
                tuple val(Sample), file ("*.recal_data.table")
        script:
        """
        ${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T BaseRecalibrator -R ${params.genome} -I ${realignedBam} -knownSites ${params.site2} -knownSites ${params.site3} -maxCycle 600 -o ${Sample}.recal_data.table
        """
}

process PrintReads {
	input:
                tuple val (Sample), file (realigned_Bam), file (recal_data_table)
        output:
                tuple val (Sample), file ("*.aligned.recalibrated.bam")
        script:
        """
        ${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T PrintReads -R ${params.genome} -I ${realigned_Bam} --BQSR ${recal_data_table} -o ${Sample}.aligned.recalibrated.bam
	"""
}

process generatefinalbam {
	input:
		tuple val (Sample), file (aligned_recalibrated_bam)
	output:
                tuple val(Sample), file ("*.final.bam"), file ("*.final.bam.bai")

        script:
	"""
	${params.java_path}/java -Xmx16G -jar ${params.abra2_path}/abra2-2.23.jar --in ${aligned_recalibrated_bam} --out ${Sample}.abra.bam --ref ${params.genome} --threads 8 --tmpdir ./ > abra.log
	${params.samtools} sort ${Sample}.abra.bam > ${Sample}.final.bam
	${params.samtools} index ${Sample}.final.bam > ${Sample}.final.bam.bai
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
		trimming_trimmomatic(samples_ch) | pair_assembly_pear | mapping_reads | sam_conversion | mark_duplicates
		RealignerTargetCreator(mark_duplicates.out)
		IndelRealigner(RealignerTargetCreator.out.join(mark_duplicates.out)) | BaseRecalibrator
		PrintReads(IndelRealigner.out.join(BaseRecalibrator.out)) | generatefinalbam
}
