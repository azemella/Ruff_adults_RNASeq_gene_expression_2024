# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# GATK RNA-Seq short variant discovery (SNPs + Indels) workflow
# In step 1, all the reads in each input library were assigned a new read-group name.
# In step 2, sequencing duplicates were identified and removed.
# In step 3, the input reads were sorted into coordinate-ordered reads for the next step.
# In step 4, reads containing Ns in their cigar string were split.
# In step 5, germline SNPs and INDELs were called via local re-assembly of haplotypes.
# In step 6, a list of input files for GATK CombineGVCFs was created.
# In step 7, joint genotyping was performed on the single multi-sample GVCF created by CombineGVCFs.
# The final output of this pipeline is a VCF in which all samples have been jointly genotyped.

# Step 1
rule AddOrReplaceReadGroups:
    input:
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
        bam = "star_ncbi/mapped/{sample}_Aligned.out.bam",
    output:
        out = temp("gatk4/addorreplacereadgroups/{sample}_addorreplacereadgroups.bam"),
    params:
        java = "--java-options -Xmx50G",
        extra0 = "--TMP_DIR gatk4/addorreplacereadgroups",
        extra1 = "--SORT_ORDER queryname",
        extra2 = "-LB 1",
        extra3 = "-PL illumina",
        extra4 = "-PU 1",
        extra5 = "-SM {sample}",
        extra6 = "--VALIDATION_STRINGENCY LENIENT",
    threads: 10
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} AddOrReplaceReadGroups -I {input.bam} -O {output.out} -R {input.genome}"
        " {params.extra0} {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5} {params.extra6}"

# Step 2
rule MarkDuplicates:
    input:
        bam = "gatk4/addorreplacereadgroups/{sample}_addorreplacereadgroups.bam",
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
    output:
        out = temp("gatk4/markduplicates/{sample}_markdups.bam"),
        txt = "gatk4/markduplicates/{sample}_marked_dup_metrics.txt",
    params:
        java = "--java-options -Xmx50G",
        extra0 = "--TMP_DIR gatk4/markduplicates",
        extra1 = "--REMOVE_SEQUENCING_DUPLICATES true",
        extra2 = "--TAGGING_POLICY OpticalOnly",
        extra3 = "--CREATE_INDEX true",
    threads: 10
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} MarkDuplicates -I {input.bam} -O {output.out} -R {input.genome} -M {output.txt}"
        " {params.extra0} {params.extra1} {params.extra2} {params.extra3}"

# Step 3
rule SortSam:
    input:
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
        bam = "gatk4/markduplicates/{sample}_markdups.bam",
    output:
        out = temp("gatk4/sortsam/{sample}_coordinate_sorted.bam"),
    params:
        java = "--java-options -Xmx50G",
        extra0 = "--TMP_DIR gatk4/sortsam",
        extra1 = "--SORT_ORDER coordinate",
    threads: 10
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} SortSam -I {input.bam} -O {output.out} -R {input.genome}"
        " {params.extra0} {params.extra1}"

# Step 4
rule SplitNCigarReads:
    input:
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
        bam = "gatk4/sortsam/{sample}_coordinate_sorted.bam",
    output:
        splitncigar = temp("gatk4/splitncigarreads/{sample}_splitncigar.bam"),
    params:
        java = "--java-options -Xmx50G",
        extra1 = "--tmp-dir gatk4/splitncigarreads",
    threads: 10
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} SplitNCigarReads -R {input.genome} -I {input.bam} -O {output.splitncigar} {params.extra1}"

# Step 5
rule HaplotypeCaller:
    input:
        bam = "gatk4/splitncigarreads/{sample}_splitncigar.bam",
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
    output:
        variants = "gatk4/haplotypecaller/{sample}_haplotypecaller.g.vcf.gz",
    params:
        java = "--java-options -Xmx50G",
        extra0 = "--tmp-dir gatk4/haplotypecaller",
        extra1 = "-ERC GVCF",
        extra2 = "-G StandardAnnotation",
        extra3 = "-G AS_StandardAnnotation",
    threads: 10
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} HaplotypeCaller -R {input.genome} -I {input.bam} -O {output.variants}"
        " {params.extra0} {params.extra1} {params.extra2} {params.extra3}"

# Step 6
rule CreateList:
    input:
        variants = expand("gatk4/haplotypecaller/{sample}_haplotypecaller.g.vcf.gz", sample = config["samples"]),
    output:
        list = "config/haplotypecaller_vcfs.list",
    threads: 4
    shell:
        "ls {input.variants} > {output.list}"

# Step 7
rule CombineGVCFs:
    input:
        variants = "config/haplosamples1.list",
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
    output:
        combine = "gatk4/combinegvcfs/combinegvcfs.vcf",
    params:
        java = "--java-options -Xmx50G",
        extra1 = "--tmp-dir gatk4/combinegvcfs",
        extra2 = "-G StandardAnnotation",
        extra3 = "-G AS_StandardAnnotation",
    threads: 50
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} CombineGVCFs -R {input.genome} -V {input.variants} -O {output.combine}"
        " {params.extra1} {params.extra2} {params.extra3}"

# Step 8
rule GenotypeGVCFs:
    input:
        combine = "gatk4/combinegvcfs/combinegvcfs.vcf",
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
    output:
        genovcf = "gatk4/genotypegvcfs/genotypegvcfs.vcf",
    params:
        java = "--java-options -Xmx50G",
        extra1 = "--tmp-dir gatk4/genotypegvcfs",
        extra2 = "-G StandardAnnotation",
        extra3 = "-G AS_StandardAnnotation",
    threads: 50
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} GenotypeGVCFs -V {input.combine} -R {input.genome} -O {output.genovcf} {params.extra1} {params.extra2} {params.extra3}"
