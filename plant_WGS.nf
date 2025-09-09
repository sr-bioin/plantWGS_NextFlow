#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Plant WGS Analysis Pipeline (Steps 1–13)
 */

params.reads    = "/HiFi/m64069_200705_101404.Q20.fastq"
params.hic_r1   = "/HiC/raw_data/S25_L001_R1_001.fastq.gz"
params.hic_r2   = "/HiC/raw_data/S25_L001_R2_001.fastq.gz"
params.bam      = "input_reads.bam"
params.genome   = "assembly.fasta"
params.outdir   = "results"

/****************************
 * Step 1–4: QC
 ****************************/

process FASTQC {
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
	input: path reads from params.reads
    
	output: path "*.zip", path "*.html"
    
	script:
    """
	module load fastqc
	
    fastqc ${reads}
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
	input: path reports from FASTQC.out.collect()
    
	output: path "multiqc_report.html"
    
	script:
    """
    multiqc .
    """
}

process LONGQC {
    publishDir "${params.outdir}/longqc", mode: 'copy'
    
	input: path bam from params.bam
    
	output: path "out_dir"
	
	// activate conda environment
	conda '/home/user/.conda/envs/longQC'
    
	script:
	
    """
    longQC.py sampleqc -x pb-sequel -o out_dir ${bam}
    """
}

process NANOPLOT {
    publishDir "${params.outdir}/nanoplot", mode: 'copy'
    
	input: path reads from params.reads
    
	output: path "_HiFi"
    
	script:
    """
    NanoPlot -t 2 --fastq ${reads} --maxlength 40000 --plots kde --legacy kde -o _HiFi
    """
}

/****************************
 * Step 5: Assembly
 ****************************/

process HIFIASM {
    publishDir "${params.outdir}/assembly", mode: 'copy'
    
	input: path reads from params.reads
           path hic_r1 from params.hic_r1
           path hic_r2 from params.hic_r2
    
	output: path "*.asm*"
    
	script:
    """
    hifiasm -o asm --primary --n-perturb 20000 --f-perturb 0.15 --seed 11 \
        -l3 --n-weight 6 -s 0.55 -k 60 --h1 ${hic_r1} --h2 ${hic_r2} ${reads}
    """
}

/****************************
 * Step 6: QUAST
 ****************************/

process QUAST {
    publishDir "${params.outdir}/quast", mode: 'copy'
    
	input: path assemblies from HIFIASM.out.collect()
    
	output: path "quast_results"
    
	script:
    """
	module load quast/5.3.0
	
    quast.py ${assemblies.join(" ")} -o quast_results
    """
}

/****************************
 * Step 7: Merqury
 ****************************/

process MERYL {
    publishDir "${params.outdir}/meryl", mode: 'copy'
    
	input: path reads from params.reads
    
	output: path "reads.meryl"
    
	script:
    """
    meryl count k=21 ${reads} output reads.meryl
    """
}

process MERQURY {
    publishDir "${params.outdir}/merqury", mode: 'copy'
    input: path assemblies from HIFIASM.out.collect()
           path meryl_db from MERYL.out
    
	output: path "merqury_results"
    
	conda activate merqury
	
	script:
    """
    merqury.sh ${meryl_db} ${assemblies.join(" ")} merqury_results
    """
}

/****************************
 * Step 8: BUSCO
 ****************************/

process BUSCO {
    publishDir "${params.outdir}/busco", mode: 'copy'
    
	input: path genome from HIFIASM.out.collect().first()
    
	output: path "busco_results"
    
	script:
    """
	module load busco/5.4.3
	
    busco -i ${genome} -l embryophyta_odb10 -m genome -o busco_results
    """
}

/****************************
 * Step 9: Repeat Annotation
 ****************************/

process REPEATMODELER {
    publishDir "${params.outdir}/repeat/repeatmodeler", mode: 'copy'

    input:
    path genome from HIFIASM.out.collect().first()

    output:
    path "consensi.fa.classified"

    script:
    """
	module load repeatmodeler/2.0.4
	
    BuildDatabase -name genome_db ${genome}
    RepeatModeler -database genome_db -pa 8 -LTRStruct
    cp genome_db-families.fa consensi.fa.classified
    """
}

process REPEATMASKER {
    publishDir "${params.outdir}/repeat/repeatmasker", mode: 'copy'

    input:
    path genome from HIFIASM.out.collect().first()
    path repeatlib from REPEATMODELER.out

    output:
    path "repeatmasker.out"
    path "genome.masked.fasta"

    script:
    """
	module load repeatmasker/4.1.5
	
    RepeatMasker -pa 8 -lib ${repeatlib} -xsmall -gff ${genome}
    mv ${genome}.masked genome.masked.fasta
    """
}

process LONGREPEATS {
    publishDir "${params.outdir}/repeat/longrepeats", mode: 'copy'

    input:
    path genome from HIFIASM.out.collect().first()

    output:
    path "longrepeats.gff3"

    script:
    """
    # Long terminal repeat (LTR) detection
    LTR_retriever -genome ${genome} -threads 8 -out longrepeats
    mv longrepeats.gff3 longrepeats.gff3
    """
}

process EDTA {
    publishDir "${params.outdir}/repeat/edta", mode: 'copy'

    input:
    path genome from HIFIASM.out.collect().first()

    output:
    path "genome.EDTA.fa"
    path "EDTA.anno"
	
	conda '/home/user/.conda/envs/EDTA'
	
    script:
    """
    EDTA.pl --genome ${genome} --species others --step all --threads 8 --anno 1 --overwrite 1
    mv ${genome%.fa}.mod.EDTA.fa genome.EDTA.fa
    mv ${genome%.fa}.EDTA.anno EDTA.anno
    """
}

process COMBINE_REPEATS {
    publishDir "${params.outdir}/repeat/final", mode: 'copy'

    input:
    path masked from REPEATMASKER.out.collect()
    path edta from EDTA.out.collect()
    path ltrs from LONGREPEATS.out.collect()

    output:
    path "final_repeat_masked_genome.fa"
    path "combined_repeats.gff3"

    script:
    """
    # Combine masked genomes: prefer EDTA as primary masked version
    cp ${edta.find { it.name.endsWith('.fa') }} final_repeat_masked_genome.fa

    # Merge annotations (GFF3 from RepeatMasker + EDTA + LongRepeats)
    cat ${masked.find { it.name.endsWith('.out') }} ${edta.find { it.name.endsWith('.anno') }} ${ltrs} > combined_repeats.gff3
    """
}


/****************************
 * Step 10: MAKER annotation
 ****************************/

process MAKER {
    publishDir "${params.outdir}/maker", mode: 'copy'
    
	input: path masked_genome from REPEATMASKER.out
    
	output: path "maker_out"
    
	script:
    """
	
    maker -genome ${masked_genome} -base maker_out
    """
}

/****************************
 * Step 11: Functional annotation
 ****************************/

process INTERPROSCAN {
    publishDir "${params.outdir}/interproscan", mode: 'copy'
    
	input: path proteins from MAKER.out
    
	output: path "interproscan_out"
    
	script:
    """
	module load interproscan/5.63-95.0
	
    interproscan.sh -i ${proteins} -f TSV -o interproscan_out
    """
}

/****************************
 * Step 12: Gene model QC
 ****************************/

process ANNOT_STATS {
    publishDir "${params.outdir}/annotation_stats", mode: 'copy'
    
	input: path gff from MAKER.out
    
	output: path "annotation_stats.txt"
    
	script:
    """
    gffstats ${gff} > annotation_stats.txt
    """
}

/****************************
 * Step 13: Final summary
 ****************************/

process FINAL_REPORT {
    publishDir "${params.outdir}/final_report", mode: 'copy'
    
	input: path reports from (MULTIQC.out.collect() +
                              QUAST.out.collect() +
                              BUSCO.out.collect() +
                              ANNOT_STATS.out.collect())
    
	output: path "final_summary.html"
    
	script:
    """
    multiqc . -o final_summary
    mv final_summary/multiqc_report.html final_summary.html
    """
}

/****************************
 * Workflow
 ****************************/

workflow {
    reads   = file(params.reads)
    bam     = file(params.bam)
    hic_r1  = file(params.hic_r1)
    hic_r2  = file(params.hic_r2)

    // QC
    fq = FASTQC(reads)
    MULTIQC(fq)
    LONGQC(bam)
    NANOPLOT(reads)

    // Assembly
    asm = HIFIASM(reads, hic_r1, hic_r2)
    QUAST(asm)

    // Merqury
    mdb = MERYL(reads)
    MERQURY(asm, mdb)

    // Assembly QC
    BUSCO(asm)

	// Repeat annotation
	rep_lib  = REPEATMODELER(asm)
	rep_mask = REPEATMASKER(asm, rep_lib)
	ltr_out  = LONGREPEATS(asm)
	edta_out = EDTA(asm)
	repeat_final = COMBINE_REPEATS(rep_mask, edta_out, ltr_out)

// Annotation uses final masked genome
	maker_out = MAKER(repeat_final)

    // Final integration
    FINAL_REPORT()
}
