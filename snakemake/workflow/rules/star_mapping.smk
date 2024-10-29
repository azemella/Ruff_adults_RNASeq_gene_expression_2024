# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Calidris pugnax RNA-Seq reads were mapped using STAR.
# In step 1, the Ruff reference genome was indexed.
# In step 2, Independent RNA-Seq libraries were mapped on the reference genome using STAR's 2nd pass mapping strategy.
# In step 3, Independent RNA-Seq libraries were mapped on the reference genome using STAR's 2nd pass mapping strategy.
# In step 4, Faeder RNA-Seq libraries were mapped on the reference genome using STAR's 2nd pass mapping strategy.
# The final outputs of this pipeline consist of alignment BAM files and gene count files for Independents, Satellites, and Faeders.

# Step 1
rule STAR_index_with_annotations:
    input:
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
        annotation = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.gtf",
    output:
        directory("star_ncbi/GenomeDir"),
    params:
        extra1 = "--runMode genomeGenerate",
        extra2 = "--genomeDir star_ncbi/GenomeDir",
        extra3 = "--genomeFastaFiles genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
        extra4 = "--sjdbGTFfile genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.gtf",
        extra5 = "--sjdbOverhang 99",
        extra6 = "--runThreadN 10",
    threads: 10
    conda:
        "../envs/star.yaml",
    shell:
        "STAR {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5} {params.extra6}"

# Step 2
rule STAR_mapping_1Pass_Independent:
    input:
        fastq1 = "trimgalore/{ind_sample}_R1_val_1.fq.gz",
        fastq2 = "trimgalore/{ind_sample}_R2_val_2.fq.gz",
        index = "star_ncbi/GenomeDir",
    output:
        junctions = "star_ncbi/junctions/{ind_sample}_SJ.out.tab",
    params:
        extra1 = "--runMode alignReads",
        extra2 = "--readFilesCommand zcat",
        extra3 = "--outSAMtype None",
        extra4 = "--runThreadN 3",
        extra5 = "--outFileNamePrefix star_ncbi/junctions/{ind_sample}_",
    threads: 3
    conda:
        "../envs/star.yaml",
    shell:
        "STAR --genomeDir {input.index} --readFilesIn {input.fastq1} {input.fastq2} "
        " {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5}"

rule junctions_filtering_Independent:
    input:
        SJ = expand("star_ncbi/junctions/{ind_sample}_SJ.out.tab", ind_sample = config["independent_samples"]),
    output:
        filteredSJ = "star_ncbi/junctions/SJ_ind_merged_filtered_junctions.tab",
    threads: 1
    shell:
        "cat {input.SJ} | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > {output.filteredSJ}"

rule STAR_mapping_2Pass_Independent:
    input:
        fastq1 = "trimgalore/{ind_sample}_R1_val_1.fq.gz",
        fastq2 = "trimgalore/{ind_sample}_R2_val_2.fq.gz",
        index = "star_ncbi/GenomeDir",
        junctions = "star_ncbi/junctions/SJ_ind_merged_filtered_junctions.tab",
    output:
        aligns = "star_ncbi/mapped/{ind_sample}_Aligned.out.bam",
        transbam = "star_ncbi/mapped/{ind_sample}_ReadsPerGene.out.tab",
    params:
        extra1 = "--runMode alignReads",
        extra2 = "--quantMode GeneCounts",
        extra3 = "--readFilesCommand zcat",
        extra4 = "--outSAMtype BAM Unsorted",
        extra5 = "--outSAMmapqUnique 60",
        extra6 = "--outSAMattributes NH HI AS nM NM MD jM jI XS MC ch",
        extra7 = "--runThreadN 3",
        extra8 = "--outSAMattrRGline ID:{ind_sample}",
        extra9 = "--outFileNamePrefix star_ncbi/mapped/{ind_sample}_",
    threads: 3
    conda:
        "../envs/star.yaml",
    shell:
        "STAR --genomeDir {input.index} --readFilesIn {input.fastq1} {input.fastq2} --sjdbFileChrStartEnd {input.junctions}"
        " {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5} {params.extra6}"
        " {params.extra7} {params.extra8} {params.extra9}"

# Step 3
rule STAR_mapping_1Pass_Satellite:
    input:
        fastq1 = "trimgalore/{sat_sample}_R1_val_1.fq.gz",
        fastq2 = "trimgalore/{sat_sample}_R2_val_2.fq.gz",
        index = "star_ncbi/GenomeDir",
   output:
        junctions = "star_ncbi/junctions/{sat_sample}_SJ.out.tab",
    params:
        extra1 = "--runMode alignReads",
        extra2 = "--readFilesCommand zcat",
        extra3 = "--outSAMtype None",
        extra4 = "--runThreadN 3",
        extra5 = "--outFileNamePrefix star_ncbi/junctions/{sat_sample}_",
    threads: 3
    conda:
        "../envs/star.yaml",
    shell:
        "STAR --genomeDir {input.index} --readFilesIn {input.fastq1} {input.fastq2} "
        " {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5}"

rule junctions_filtering_Satellite:
    input:
        SJ = expand("star_ncbi/junctions/{sat_sample}_SJ.out.tab", sat_sample = config["satellite_samples"]),
    output:
        filteredSJ = "star_ncbi/junctions/SJ_sat_merged_filtered_junctions.tab",
    threads: 1
    shell:
        "cat {input.SJ} | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > {output.filteredSJ}"

rule STAR_mapping_2Pass_Satellite:
    input:
        fastq1 = "trimgalore/{sat_sample}_R1_val_1.fq.gz",
        fastq2 = "trimgalore/{sat_sample}_R2_val_2.fq.gz",
        index = "star_ncbi/GenomeDir",
        junctions = "star_ncbi/junctions/SJ_sat_merged_filtered_junctions.tab",
    output:
        aligns = "star_ncbi/mapped/{sat_sample}_Aligned.out.bam",
       transbam = "star_ncbi/mapped/{sat_sample}_ReadsPerGene.out.tab",
    params:
        extra1 = "--runMode alignReads",
        extra2 = "--quantMode GeneCounts",
        extra3 = "--readFilesCommand zcat",
        extra4 = "--outSAMtype BAM Unsorted",
        extra5 = "--outSAMmapqUnique 60",
        extra6 = "--outSAMattributes NH HI AS nM NM MD jM jI XS MC ch",
        extra7 = "--runThreadN 3",
        extra8 = "--outSAMattrRGline ID:{sat_sample}",
        extra9 = "--outFileNamePrefix star_ncbi/mapped/{sat_sample}_",
    threads: 3
    conda:
        "../envs/star.yaml",
    shell:
        "STAR --genomeDir {input.index} --readFilesIn {input.fastq1} {input.fastq2} --sjdbFileChrStartEnd {input.junctions}"
        " {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5} {params.extra6}"
        " {params.extra7} {params.extra8} {params.extra9}"

# Step 4
rule STAR_mapping_1Pass_Faeder:
    input:
        fastq1 = "trimgalore/{fae_sample}_R1_val_1.fq.gz",
        fastq2 = "trimgalore/{fae_sample}_R2_val_2.fq.gz",
        index = "star_ncbi/GenomeDir",
    output:
        junctions = "star_ncbi/junctions/{fae_sample}_SJ.out.tab",
    params:
        extra1 = "--runMode alignReads",
        extra2 = "--readFilesCommand zcat",
        extra3 = "--outSAMtype None",
        extra4 = "--runThreadN 3",
        extra5 = "--outFileNamePrefix star_ncbi/junctions/{fae_sample}_",
    threads: 3
    conda:
        "../envs/star.yaml",
    shell:
        "STAR --genomeDir {input.index} --readFilesIn {input.fastq1} {input.fastq2} "
        " {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5}"

rule junctions_filtering_Faeder:
    input:
        SJ = expand("star_ncbi/junctions/{fae_sample}_SJ.out.tab", fae_sample = config["faeder_samples"]),
    output:
        filteredSJ = "star_ncbi/junctions/SJ_fae_merged_filtered_junctions.tab",
    threads: 1
    shell:
        "cat {input.SJ} | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > {output.filteredSJ}"

rule STAR_mapping_2Pass_Faeder:
    input:
        fastq1 = "trimgalore/{fae_sample}_R1_val_1.fq.gz",
        fastq2 = "trimgalore/{fae_sample}_R2_val_2.fq.gz",
        index = "star_ncbi/GenomeDir",
        junctions = "star_ncbi/junctions/SJ_fae_merged_filtered_junctions.tab",
    output:
        aligns = "star_ncbi/mapped/{fae_sample}_Aligned.out.bam",
        transbam = "star_ncbi/mapped/{fae_sample}_ReadsPerGene.out.tab",
    params:
        extra1 = "--runMode alignReads",
        extra2 = "--quantMode GeneCounts",
       extra3 = "--readFilesCommand zcat",
        extra4 = "--outSAMtype BAM Unsorted",
       extra5 = "--outSAMmapqUnique 60",
       extra6 = "--outSAMattributes NH HI AS nM NM MD jM jI XS MC ch",
        extra7 = "--runThreadN 3",
        extra8 = "--outSAMattrRGline ID:{fae_sample}",
        extra9 = "--outFileNamePrefix star_ncbi/mapped/{fae_sample}_",
    threads: 3
    conda:
        "../envs/star.yaml",
    shell:
        "STAR --genomeDir {input.index} --readFilesIn {input.fastq1} {input.fastq2} --sjdbFileChrStartEnd {input.junctions}"
        " {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5} {params.extra6}"
        " {params.extra7} {params.extra8} {params.extra9}"
