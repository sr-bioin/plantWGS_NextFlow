<h2>Nextflow pipeline for whole genome sequence analysis of Plant </h2>

<h3>Software used in the pipeline.</h3>

**FastQC** is widely used and is robust, efficient, and versatile quality control software for a varied range of raw genetic data. It outputs a quality report which can be viewed to give an indication on how good the respective reads are. For more information, please visit the FastQC website https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.

**MultiQC** is used create a single report with interactive plots for multiple bioinformatics analyses across many samples. It reports can describe multiple analysis steps and large numbers of samples within a single plot, and multiple analysis tools making it ideal for routine fast quality control. For more information, please visit the MultiQC website https://github.com/MultiQC/MultiQC.

**NanoQC** is a plotting tool for long read sequencing data and alignments. For more information, please check its website: https://biocontainers.pro/tools/nanoplot and https://github.com/wdecoster/NanoPlot. <br>
<img width="400" height="400" alt="LengthvsQualityScatterPlot_legacy_kde" src="https://github.com/user-attachments/assets/1bcb8557-c3c7-451b-a7db-7ea96fab3aa5" />

**Hifiasm** 
is a fast haplotype-resolved de novo assembler initially designed for PacBio HiFi reads. It produces arguably the best single-sample telomere-to-telomere assemblies combing HiFi, ultralong and Hi-C reads, and it is one of the best haplotype-resolved assemblers for the trio-binning assembly given parental short reads. For more information, please visit https://github.com/chhylp123/hifiasm. <br> 
  &ensp; **Completeness:** <br>
  &emsp; &emsp; **Merqury** is a tool for validating genomes assembled using long read sequencing. The algorithm uses k-mer to evaluate base accuracy and <br> &emsp; &emsp; completeness of a genome by comparing the de novo assembled genome with high accuracy reads that are not used in the assembly. <br>
  &emsp; &emsp; The program does not require another reference assembly for its evaluation. It is able to evaluate the haplotype-specific accuracy, <br>
  &emsp; &emsp; completeness, phase block continuity and switch errors in trio binned assemblies. https://github.com/marbl/merqury <br>
  &emsp; &emsp; **QUAST** (QUality ASsessment Tool) is a widely used software for evaluating and comparing genome assemblies in bioinformatics, <br>
  &emsp; &emsp; providing metrics on assembly quality, such as contiguity and completeness. The tool can be used with or without a reference genome <br>
  &emsp; &emsp; and generates detailed reports, tables, and plots for analysis.

**Maker2** 
is a genome annotation and data management tool designed for second-generation genome projects. It is a multi-threaded, parallelized application that can process second-generation datasets of virtually any size. It can produce accurate annotations for novel genomes where training-data are limited, of low quality or even non-existent. . For more information, please visit the website https://github.com/Yandell-Lab/maker.

**QUAST** 
is a quality assessment tool for evaluating and comparing genome assemblies by computing various metrics and works both with and without reference genomes. It produces many reports, summary tables and plots to help scientists in their research and in their publications. For more information, please visit the Quast website https://github.com/ablab/quast.

**BUSCO**
 (Benchmarking Universal Single-Copy Orthologs) provides measures for quantitative assessment of genome assembly, gene set, and transcriptome completeness based on evolutionarily informed expectations of gene content from near-universal single-copy orthologs, please visit https://busco.ezlab.org/.


