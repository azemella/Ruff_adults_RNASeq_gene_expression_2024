# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Faeder and Satellite RNA-Seq reads were mapped on each inversion morph-specific consensus genome using STARconsensus mapping.
# In step 1, Faeder RNA-Seq libraries were mapped on the new Faeder's consensus genome using STAR 2nd pass mapping strategy.
# In step 2, Satellite RNA-Seq libraries were mapped on the new Satellite's consensus genome using STAR 2nd pass mapping strategy.
# The final output of this file consists of alignment BAM files and gene count files for Faeders and Satellites.

# Step 1
rule STARconsensus_genome_index_Faeder:
    input:
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
        annotation = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.gtf",
    output:
        directory("star_consensus/faeder/GenomeDir"),
    params:
        extra1 = "--runMode genomeGenerate",
        extra2 = "--genomeDir star_consensus/faeder/GenomeDir",
        extra3 = "--genomeFastaFiles genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
        extra4 = "--sjdbGTFfile genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.gtf",
        extra5 = "--sjdbOverhang 99",
        extra6 = "--runThreadN 20",
        extra7 = "--genomeTransformType Haploid",
        extra8 = "--genomeTransformVCF bcftools/consensus_inversion_hetSNPs_and_hetINDELs_merged_faeder.vcf",
    threads: 20
    conda:
        "../envs/star.yaml",
    shell:
        "STAR {params.extra1} {params.extra2} {params.extra3} {params.extra4}"
        " {params.extra5} {params.extra6} {params.extra7} {params.extra8}"

rule STARconsensus_mapping_1Pass_Faeder:
    input:
        fastq1 = "trimgalore/{fae_sample}_R1_val_1.fq.gz",
        fastq2 = "trimgalore/{fae_sample}_R2_val_2.fq.gz",
        index = "star_consensus/faeder/GenomeDir",
    output:
        junctions = "star_consensus/faeder/junctions/{fae_sample}_SJ.out.tab",
    params:
        extra1 = "--runMode alignReads",
        extra2 = "--readFilesCommand zcat",
        extra3 = "--outSAMtype None",
        extra4 = "--runThreadN 5",
        extra5 = "--outFileNamePrefix star_consensus/faeder/junctions/{fae_sample}_",
    threads: 5
    conda:
        "../envs/star.yaml",
    shell:
        "STAR --genomeDir {input.index} --readFilesIn {input.fastq1} {input.fastq2} "
        " {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5}"

rule STARconsensus_junctions_filtering_Faeder:
    input:
        SJ = expand("star_consensus/faeder/junctions/{fae_sample}_SJ.out.tab", fae_sample = config["faeder_samples"]),
    output:
        filteredSJ = "star_consensus/faeder/junctions/SJ.faeder.merged.filtered.tab",
    threads: 1
    shell:
        "cat {input.SJ} | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > {output.filteredSJ}"

rule STARconsensus_mapping_2Pass_Faeder:
    input:
        fastq1 = "trimgalore/{fae_sample}_R1_val_1.fq.gz",
        fastq2 = "trimgalore/{fae_sample}_R2_val_2.fq.gz",
        index = "star_consensus/faeder/GenomeDir",
        junctions = "star_consensus/faeder/junctions/SJ.faeder.merged.filtered.tab",
    output:
        aligns = "star_consensus/faeder/mapped/{fae_sample}_Aligned.out.bam",
        transbam = "star_consensus/faeder/mapped/{fae_sample}_ReadsPerGene.out.tab",
    params:
        extra1 = "--runMode alignReads",
        extra2 = "--quantMode GeneCounts",
        extra3 = "--readFilesCommand zcat",
        extra4 = "--outSAMtype BAM Unsorted",
        extra5 = "--outSAMmapqUnique 60",
        extra6 = "--outSAMattributes NH HI AS nM NM MD jM jI XS MC ch",
        extra7 = "--runThreadN 5",
        extra8 = "--outSAMattrRGline ID:{fae_sample}",
        extra9 = "--outFileNamePrefix star_consensus/faeder/mapped/{fae_sample}_",
        extra10 = "--genomeTransformOutput SAM SJ",
    threads: 5
    conda:
        "../envs/star.yaml",
    shell:
        "STAR --genomeDir {input.index} --readFilesIn {input.fastq1} {input.fastq2} --sjdbFileChrStartEnd {input.junctions}"
        " {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5} {params.extra6}"
        " {params.extra7} {params.extra8} {params.extra9} {params.extra10}"

# Step 2
rule STARconsensus_genome_index_Satellite_:
    input:
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
        annotation = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.gtf",
    output:
        directory("star_consensus/satellite/GenomeDir"),
    params:
        extra1 = "--runMode genomeGenerate",
        extra2 = "--genomeDir star_consensus/satellite/GenomeDir",
        extra3 = "--genomeFastaFiles genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
        extra4 = "--sjdbGTFfile genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.gtf",
        extra5 = "--sjdbOverhang 99",
        extra6 = "--runThreadN 20",
        extra7 = "--genomeTransformType Haploid",
        extra8 = "--genomeTransformVCF bcftools/consensus_inversion_hetSNPs_and_hetINDELs_merged_satellite.vcf",
    threads: 20
    conda:
        "../envs/star.yaml",
    shell:
        "STAR {params.extra1} {params.extra2} {params.extra3} {params.extra4}"
        " {params.extra5} {params.extra6} {params.extra7} {params.extra8}"

rule STARconsensus_mapping_1Pass_Satellite:
    input:
        fastq1 = "trimgalore/{sat_sample}_R1_val_1.fq.gz",
        fastq2 = "trimgalore/{sat_sample}_R2_val_2.fq.gz",
        index = "star_consensus/satellite/GenomeDir",
    output:
        junctions = "star_consensus/satellite/junctions/{sat_sample}_SJ.out.tab",
    params:
        extra1 = "--runMode alignReads",
        extra2 = "--readFilesCommand zcat",
        extra3 = "--outSAMtype None",
        extra4 = "--runThreadN 5",
        extra5 = "--outFileNamePrefix star_consensus/satellite/junctions/{sat_sample}_",
    threads: 5
    conda:
        "../envs/star.yaml",
    shell:
        "STAR --genomeDir {input.index} --readFilesIn {input.fastq1} {input.fastq2} "
        " {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5}"

rule STARconsensus_junctions_filtering_Satellite:
    input:
        SJ = expand("star_consensus/satellite/junctions/{sat_sample}_SJ.out.tab", sat_sample = config["satellite_samples"]),
    output:
        filteredSJ = "star_consensus/satellite/junctions/SJ.satellite.merged.filtered.tab",
    threads: 1
    shell:
        "cat {input.SJ} | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > {output.filteredSJ}"

rule STARconsensus_mapping_2Pass_Satellite:
    input:
        fastq1 = "trimgalore/{sat_sample}_R1_val_1.fq.gz",
        fastq2 = "trimgalore/{sat_sample}_R2_val_2.fq.gz",
        index = "star_consensus/satellite/GenomeDir",
        junctions = "star_consensus/satellite/junctions/SJ.satellite.merged.filtered.tab",
    output:
        aligns = "star_consensus/satellite/mapped/{sat_sample}_Aligned.out.bam",
        transbam = "star_consensus/satellite/mapped/{sat_sample}_ReadsPerGene.out.tab",
    params:
        extra1 = "--runMode alignReads",
        extra2 = "--quantMode GeneCounts",
        extra3 = "--readFilesCommand zcat",
        extra4 = "--outSAMtype BAM Unsorted",
        extra5 = "--outSAMmapqUnique 60",
        extra6 = "--outSAMattributes NH HI AS nM NM MD jM jI XS MC ch",
        extra7 = "--runThreadN 5",
        extra8 = "--outSAMattrRGline ID:{sat_sample}",
        extra9 = "--outFileNamePrefix star_consensus/satellite/mapped/{sat_sample}_",
        extra10 = "--genomeTransformOutput SAM SJ",
    threads: 3
    conda:
        "../envs/star.yaml",
    shell:
        "STAR --genomeDir {input.index} --readFilesIn {input.fastq1} {input.fastq2} --sjdbFileChrStartEnd {input.junctions}"
        " {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5} {params.extra6}"
        " {params.extra7} {params.extra8} {params.extra9} {params.extra10}"