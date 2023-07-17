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
	${params.sequences}/${Sample}_*R1_*.fastq.gz ${params.sequences}/${Sample}_*R2_*.fastq.gz -baseout ${Sample}.fq.gz \
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
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final.bam*'
	input:
		tuple val (Sample), file (aligned_recalibrated_bam)
	output:
		tuple val(Sample), file ("*.final.bam"), file ("*.final.bam.bai"), file ("*.old_final.bam"), file ("*.old_final.bam.bai")
	script:
	"""
	${params.java_path}/java -Xmx16G -jar ${params.abra2_path}/abra2-2.23.jar --in ${aligned_recalibrated_bam} --out ${Sample}.abra.bam --ref ${params.genome} --threads 8 --tmpdir ./ > abra.log
	${params.samtools} sort ${aligned_recalibrated_bam} > ${Sample}.old_final.bam
	${params.samtools} index ${Sample}.old_final.bam > ${Sample}.old_final.bam.bai
	${params.samtools} sort ${Sample}.abra.bam > ${Sample}.final.bam
	${params.samtools} index ${Sample}.final.bam > ${Sample}.final.bam.bai
	"""

}

process coverage_mosdepth {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*cov*'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*")
	script:
	"""
	${params.mosdepth_script} ${finalBam} ${Sample}_cov
	"""
}

process freebayes {
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.freebayes.vcf")
	script:
	"""
	${params.freebayes_path} -f ${params.genome} -b ${finalBam} > ${Sample}.freebayes.vcf
	"""
}

process haplotypecaller {
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.haplotypecaller.vcf")
	script:
	"""
	${params.java_path}/java -Xmx10G -jar ${params.GATK38_path} -T HaplotypeCaller -R ${params.genome} -I ${finalBam} -o ${Sample}.haplotypecaller.vcf --dbsnp ${params.site2} 
	"""
}

process strelka {
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.strelka.vcf")
	script:
	"""
	${params.strelka_path}/configureStrelkaGermlineWorkflow.py --bam ${finalBam} --referenceFasta ${params.genome} --targeted --exome --runDir ./
	./runWorkflow.py -m local -j 20
	gunzip -f ./results/variants/variants.vcf.gz
	mv ./results/variants/variants.vcf ${Sample}.strelka.vcf
	"""
}

process platypus {
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file("*.platypus.vcf")
	script:
	"""
	 python2.7 ${params.platypus_path} callVariants --bamFiles=${finalBam[0]} --refFile=${params.genome} --output=${Sample}.platypus.vcf --nCPU=15 --minFlank=0 --filterDuplicates=0 --maxVariants=6 --minReads=6
	"""
}

process vardict {
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.vardict.vcf")
	script:
	"""
	VarDict -G ${params.genome} -f 0.03 -N ${Sample} -b ${finalBam} -c 1 -S 2 -E 3 -g 4 ${params.bedfile} | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.03 > ${Sample}.vardict.vcf
	"""
}

process varscan {
	input:
		 tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.varscan.vcf")
	script:
	"""
	${params.samtools} mpileup -f ${params.genome} ${finalBam} > ${Sample}.mpileup
	${params.java_path}/java -jar ${params.varscan_path} mpileup2snp ${Sample}.mpileup --min-coverage 10 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.003 --p-value 1e-4 --output-vcf 1 > ${Sample}.varscan_snp.vcf
	${params.java_path}/java -jar ${params.varscan_path} mpileup2indel ${Sample}.mpileup --min-coverage 10 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.003 --p-value 1e-4 --output-vcf 1 > ${Sample}.varscan_indel.vcf
	bgzip -c ${Sample}.varscan_snp.vcf > ${Sample}.varscan_snp.vcf.gz
	bgzip -c ${Sample}.varscan_indel.vcf > ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_snp.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} concat -a ${Sample}.varscan_snp.vcf.gz ${Sample}.varscan_indel.vcf.gz -o ${Sample}.varscan.vcf
	"""
}

process pindel {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_pindel_SI.vcf'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_vep_delheaders.txt'
	input:
		tuple val (Sample), file (finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file("*")
	script:
	"""
	export BAM_2_PINDEL_ADAPT=${params.pindel}/Adaptor.pm
	sh ${params.pindel_config_script} -s ${Sample} -b ${finalBam} -c config.txt
	${params.pindel}pindel -f ${params.genome} -i config.txt -c chr13 -o ${Sample}_pindel
	${params.pindel}pindel2vcf -r ${params.genome} -P ${Sample}_pindel -R hg19 -d 07102019 -v ${Sample}_pindel_SI.vcf

	${params.vep_script_path} ${Sample}_pindel_SI.vcf ${Sample}
	#${params.vep_extract_path} ${Sample}.final.concat_append.csv ${PWD}/Final_Output/${Sample}/${Sample}_vep_delheaders.txt > ${Sample}.vep
	"""
}

process lofreq {
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file("*.lofreq.filtered.vcf")
	script:
	"""
	${params.lofreq_path} viterbi -f ${params.genome} -o ${Sample}.lofreq.pre.bam ${oldfinalBam}
	${params.samtools} sort ${Sample}.lofreq.pre.bam > ${Sample}.lofreq.bam
	${params.lofreq_path} call -b dynamic -C 50 -a 0.00005 -q 30 -Q 30 -m 50 -f ${params.genome} -l ${params.bedfile}.bed -o ${Sample}.lofreq.vcf ${Sample}.lofreq.bam
	${params.lofreq_path} filter -a 0.005 -i ${Sample}.lofreq.vcf -o ${Sample}.lofreq.filtered.vcf
	"""
}

process somaticSeq_run {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.somaticseq.vcf'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_vep_delheaders.txt'
	input:
		tuple val (Sample), file (lofreqVcf), file (varscanVcf), file (platypusVcf), file (strelkaVcf), file (haplotypecallerVcf), file (freebayesVcf), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.somaticseq.vcf")
	script:
	"""
	${params.vcf_sorter_path} ${freebayesVcf} ${Sample}.freebayes.sorted.vcf
	${params.vcf_sorter_path} ${platypusVcf} ${Sample}.platypus.sorted.vcf
	${params.vcf_sorter_path} ${haplotypecallerVcf} ${Sample}.haplotypecaller.sorted.vcf

	python3 ${params.splitvcf_path} -infile ${Sample}.platypus.sorted.vcf -snv ${Sample}_platypus_cnvs.vcf -indel ${Sample}_platypus_indels.vcf
	python3 ${params.splitvcf_path} -infile ${Sample}.freebayes.sorted.vcf -snv ${Sample}_freebayes_cnvs.vcf -indel ${Sample}_freebayes_indels.vcf
	python3 ${params.splitvcf_path} -infile ${Sample}.haplotypecaller.sorted.vcf -snv ${Sample}_haplotypecaller_cnvs.vcf -indel ${Sample}_haplotypecaller_indels.vcf

	${params.vcf_sorter_path} ${Sample}_platypus_cnvs.vcf ${Sample}_platypus_cnvs_sort.vcf
	${params.vcf_sorter_path} ${Sample}_platypus_indels.vcf ${Sample}_platypus_indels_sort.vcf
	${params.vcf_sorter_path} ${Sample}_freebayes_cnvs.vcf ${Sample}_freebayes_cnvs_sort.vcf
	${params.vcf_sorter_path} ${Sample}_freebayes_indels.vcf ${Sample}_freebayes_indels_sort.vcf
	${params.vcf_sorter_path} ${Sample}_haplotypecaller_cnvs.vcf ${Sample}_haplotypecaller_cnvs_sort.vcf
	${params.vcf_sorter_path} ${Sample}_haplotypecaller_indels.vcf ${Sample}_haplotypecaller_indels_sort.vcf
	
	somaticseq_parallel.py --output-directory ${Sample}.somaticseq --genome-reference ${params.genome} --inclusion-region ${params.bedfile}.bed --threads 25 --algorithm xgboost  --dbsnp-vcf  /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file ${finalBam} --varscan-vcf ${varscanVcf} --lofreq-vcf ${lofreqVcf} --strelka-vcf ${strelkaVcf} --sample-name ${Sample} --arbitrary-snvs ${Sample}_freebayes_cnvs_sort.vcf ${Sample}_platypus_cnvs_sort.vcf ${Sample}_haplotypecaller_cnvs_sort.vcf --arbitrary-indels ${Sample}_freebayes_indels_sort.vcf ${Sample}_platypus_indels_sort.vcf ${Sample}_haplotypecaller_indels_sort.vcf

	${params.vcf_sorter_path} ${Sample}.somaticseq/Consensus.sSNV.vcf ${Sample}.somaticseq/somaticseq_snv.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_snv.vcf > ${Sample}.somaticseq/somaticseq_snv.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_snv.vcf.gz

	${params.vcf_sorter_path} ${Sample}.somaticseq/Consensus.sINDEL.vcf ${Sample}.somaticseq/somaticseq_indel.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_indel.vcf > ${Sample}.somaticseq/somaticseq_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_indel.vcf.gz

	${params.bcftools_path} concat -a ${Sample}.somaticseq/somaticseq_snv.vcf.gz ${Sample}.somaticseq/somaticseq_indel.vcf.gz -o ${Sample}.somaticseq.vcf

	sed -i 's/##INFO=<ID=VLK012,Number=6,Type=Integer,Description="Calling decision of the 6 algorithms: VarScan2, LoFreq, Strelka, SnvCaller_0, SnvCaller_1, SnvCaller_2">/##INFO=<ID=VLSFPH,Number=6,Type=String,Description="Calling decision of the 6 algorithms:  VarScan2, LoFreq, Strelka, Freebayes, Platypus, Haplotypecaller">/g' ${Sample}.somaticseq.vcf

	sed -i 's/VLK012/VLSFPH/g' ${Sample}.somaticseq.vcf
	
	#adding vep
	${params.vep_script_path} ${Sample}.somaticseq.vcf ${Sample}
	"""
}

process cava {
	input:
		 tuple val(Sample), file (somaticVcf)
	output:
		tuple val(Sample), file ("*")
	script:
	"""
	${params.cava_path}/cava -c ${params.cava_path}/config_v2.txt -t 10 -i ${Sample}.somaticseq.vcf -o ${Sample}.somaticseq
	python3 ${params.cava_script_path} ${Sample}.somaticseq.txt ${Sample}.cava.csv
	"""
}

process merge_csv {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.xlsx'
	input:
		tuple val (Sample), file (cava_csv)
	output:
		val Sample
	script:
	"""
	python3 ${params.pharma_marker_script} ${Sample} $PWD/Final_Output/${Sample}/${Sample}_vep_delheaders.txt ${params.pharma_input_xlxs} ${PWD}/${Sample}/${Sample}_pharma.csv
	python3 ${params.merge_csvs_script} ${Sample} ${cava_csv} ${PWD}/Final_Output/${Sample}/${Sample}.xlsx  $PWD/Final_Output/${Sample}/${Sample}_cov.mosdepth.summary.txt $PWD/Final_Output/${Sample}/${Sample}_cov.regions.bed $PWD/Final_Output/${Sample}/${Sample}_median50 $PWD/Final_Output/${Sample}/${Sample}_pindel.vep.txt $PWD/Final_Output/${Sample}/${Sample}_somaticseq.vep.txt
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
		coverage_mosdepth(generatefinalbam.out)
		freebayes(generatefinalbam.out)
		haplotypecaller(generatefinalbam.out)
		strelka(generatefinalbam.out)
		platypus(generatefinalbam.out)
		//vardict(generatefinalbam.out)
		varscan(generatefinalbam.out)
		lofreq(generatefinalbam.out)
		pindel(generatefinalbam.out)
		somaticSeq_run(lofreq.out.join(varscan.out.join(platypus.out.join(strelka.out.join(haplotypecaller.out.join(freebayes.out.join(generatefinalbam.out)))))))
		cava(somaticSeq_run.out)
		merge_csv(cava.out)
}
