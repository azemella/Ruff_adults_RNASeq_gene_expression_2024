# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Read counts per HetSNPs allele were calculated using the GATK ASEReadCounter workflow.
# In step 1, all the Satellite reads were assigned a new read-group name.
# In step 2, all the Faeder reads were assigned a new read-group name.
# In step 3, sequencing duplicates in Satellite samples were identified and removed.
# In step 4, sequencing duplicates in Faeder samples were identified and removed.
# In step 5, a table of filtered base counts at each Bi-allelic HetSNPs was generated for Satellite samples.
# In step 6, a table of filtered base counts at each Bi-allelic HetSNPs was generated for Faeder samples.
# The final output of this file consists of tables of allele counts at each Bi-allelic HetSNPs for downstream Allele-Specific Expression (ASE) analysis in R.

# Step 1
rule AddOrReplaceReadGroups_Satellite:
    input:
        bam = "star_consensus/satellite/mapped/{sat_sample}_Aligned.out.bam",
    output:
        out = "gatk4/asereadcounter/satellite/{sat_sample}_Aligned_sorted_RG.bam",
    params:
        java = "--java-options -Xmx20G",
        extra0 = "--TMP_DIR gatk4/addorreplacereadgroups",
        extra1 = "--SORT_ORDER coordinate",
        extra2 = "-LB 1",
        extra3 = "-PL Illumina",
        extra4 = "-PU 1",
        extra5 = "-SM {sat_sample}",
        extra6 = "--VALIDATION_STRINGENCY LENIENT",
        extra7 = "--CREATE_INDEX true",
    threads: 5
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} AddOrReplaceReadGroups -I {input.bam} -O {output.out} {params.extra0} {params.extra1}"
        " {params.extra2} {params.extra3} {params.extra4} {params.extra5} {params.extra6} {params.extra7}"

# Step 2
rule AddOrReplaceReadGroups_Faeder:
    input:
        bam = "star_consensus/faeder/mapped/{fae_sample}_Aligned.out.bam",
    output:
        out = "gatk4/asereadcounter/faeder/{fae_sample}_Aligned_sorted_RG.bam",
    params:
        java = "--java-options -Xmx20G",
        extra0 = "--TMP_DIR gatk4/addorreplacereadgroups",
        extra1 = "--SORT_ORDER coordinate",
        extra2 = "-LB 1",
        extra3 = "-PL Illumina",
        extra4 = "-PU 1",
        extra5 = "-SM {fae_sample}",
        extra6 = "--VALIDATION_STRINGENCY LENIENT",
        extra7 = "--CREATE_INDEX true",
    threads: 5
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} AddOrReplaceReadGroups -I {input.bam} -O {output.out} {params.extra0} {params.extra1}"
        " {params.extra2} {params.extra3} {params.extra4} {params.extra5} {params.extra6} {params.extra7}"

# Step 3
rule MarkDuplicates_Satellite:
    input:
        bam = "gatk4/asereadcounter/satellite/{sat_sample}_Aligned_sorted_RG.bam",
    output:
        out = "gatk4/asereadcounter/satellite/{sat_sample}_markdups.bam",
        txt = "gatk4/asereadcounter/{sat_sample}_marked_dup_metrics.txt",
    params:
        java = "--java-options -Xmx20G",
        extra0 = "--TMP_DIR gatk4/markduplicates",
        extra1 = "--REMOVE_SEQUENCING_DUPLICATES true",
        extra2 = "--TAGGING_POLICY OpticalOnly",
        extra3 = "--CREATE_INDEX true",
    threads: 5
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} MarkDuplicates -I {input.bam} -O {output.out} -M {output.txt}"
        " {params.extra0} {params.extra1} {params.extra2} {params.extra3}"

# Step 4
rule MarkDuplicates_Faeder:
    input:
        bam = "gatk4/asereadcounter/faeder/{fae_sample}_Aligned_sorted_RG.bam",
    output:
        out = "gatk4/asereadcounter/faeder/{fae_sample}_markdups.bam",
        txt = "gatk4/asereadcounter/{fae_sample}_marked_dup_metrics.txt",
    params:
        java = "--java-options -Xmx20G",
        extra0 = "--TMP_DIR gatk4/markduplicates",
        extra1 = "--REMOVE_SEQUENCING_DUPLICATES true",
        extra2 = "--TAGGING_POLICY OpticalOnly",
        extra3 = "--CREATE_INDEX true",
    threads: 5
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} MarkDuplicates -I {input.bam} -O {output.out} -M {output.txt}"
        " {params.extra0} {params.extra1} {params.extra2} {params.extra3}"

# Step 5
rule ASEReadCounter_Satellite:
    input:
        vcf = "gatk4/ruff_adults_variants/consensus_inversion_hetSNPs_satellite.vcf",
        bam = "gatk4/asereadcounter/satellite/{sat_sample}_markdups.bam",
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
    output:
        aser = "gatk4/asereadcounter/satellite/{sat_sample}_ASEReadCounter.table",
    params:
        java = "--java-options '-Xmx20G'",
        extra0 = "--min-mapping-quality '50'",
        extra1 = "--min-base-quality '25'",
        extra2 = "--tmp-dir gatk4/asereadcounter/satellite",
        extra3 = "--intervals 'config/inversion_genes_coordinates.list'",
    threads: 5
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} ASEReadCounter -V {input.vcf} -I {input.bam} -R {input.genome} -O {output.aser}"
        " {params.extra0} {params.extra1} {params.extra2} {params.extra3}"

# Step 6
rule ASEReadCounter_Faeder:
    input:
        vcf = "gatk4/ruff_adults_variants/consensus_inversion_hetSNPs_faeder.vcf",
        bam = "gatk4/asereadcounter/faeder/{fae_sample}_markdups.bam",
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
    output:
        aser = "gatk4/asereadcounter/faeder/{fae_sample}_ASEReadCounter.table",
    params:
        java = "--java-options '-Xmx20G'",
        extra0 = "--min-mapping-quality '50'",
        extra1 = "--min-base-quality '25'",
        extra2 = "--tmp-dir gatk4/asereadcounter/faeder",
        extra3 = "--intervals 'config/inversion_genes_coordinates.list'",
    threads: 5
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} ASEReadCounter -V {input.vcf} -I {input.bam} -R {input.genome} -O {output.aser}"
        " {params.extra0} {params.extra1} {params.extra2} {params.extra3}"
