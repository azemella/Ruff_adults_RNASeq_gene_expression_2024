# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Faeder and Satellite RNA-Seq reads were mapped on each inversion morph-specific consensus genome using STARconsensus mapping.
# In step 1, Faeder consensus INDELs VCF was compressed.
# In step 2, Faeder consensus INDELs VCF was indexed.
# In step 3, Faeder consensus SNPs VCF was compressed.
# In step 4, Faeder consensus SNPs VCF was indexed.
# In step 5, Faeder INDELs and SNPs were combined into a single VCF.
# In step 6, Satellite consensus INDELs VCF was compressed.
# In step 7, Satellite consensus INDELs VCF was indexed.
# In step 8, Satellite consensus SNPs VCF was compressed.
# In step 9, Satellite consensus SNPs VCF was indexed.
# In step 10, Satellite INDELs and SNPs were combined into a single VCF.
# The final output of this file consists of two VCFs of SNPs and INDELs for Faeders and Satellites.

# Step 1
rule bcftools_bgzip_inversion_hetINDELs_Faeder:
    input:
        indels = "gatk4/ruff_adults_variants/consensus_inversion_hetINDELs_faeder.vcf",
    output:
        bgzip = "bcftools/consensus_inversion_hetINDELs_faeder.vcf.gz",
    threads: 1
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bgzip -c {input.indels} > {output.bgzip}"

# Step 2
rule bcftools_tabix_inversion_hetINDELs_Faeder:
    input:
        bgzip = "bcftools/consensus_inversion_hetINDELs_faeder.vcf.gz",
    output:
        tabix = "bcftools/consensus_inversion_hetINDELs_faeder.vcf.gz.tbi",
    threads: 1
    priority: 1000
    conda:
        "../envs/bcftools.yaml",
    shell:
        "tabix -p vcf {input.bgzip}"

# Step 3
rule bcftools_bgzip_inversion_hetSNPs_Faeder:
    input:
        snps = "gatk4/ruff_adults_variants/consensus_inversion_hetSNPs_faeder.vcf",
    output:
        bgzip = "bcftools/consensus_inversion_hetSNPs_faeder.vcf.gz",
    threads: 1
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bgzip -c {input.snps} > {output.bgzip}"

# Step 4
rule bcftools_tabix_inversion_hetSNPs_Faeder:
    input:
        bgzip = "bcftools/consensus_inversion_hetSNPs_faeder.vcf.gz",
    output:
        tabix = "bcftools/consensus_inversion_hetSNPs_faeder.vcf.gz.tbi",
    threads: 1
    priority: 1000
    conda:
        "../envs/bcftools.yaml",
    shell:
        "tabix -p vcf {input.bgzip}"

# Step 5
rule bcftools_concat_hetVARIANTS_Faeder:
    input:
        snps = "bcftools/consensus_inversion_hetSNPs_faeder.vcf.gz",
        indels = "bcftools/consensus_inversion_hetINDELs_faeder.vcf.gz",
    output:
        merged = "bcftools/consensus_inversion_hetSNPs_and_hetINDELs_merged_faeder.vcf",
    params:
        extra = "--allow-overlaps",
    threads: 1
    priority: 500
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bcftools concat {params.extra} {input.snps} {input.indels} -o {output.merged}"

# Step 6
rule bcftools_bgzip_inversion_hetINDELs_Satellite:
    input:
        indels = "gatk4/ruff_adults_variants/consensus_inversion_hetINDELs_satellite.vcf",
    output:
        bgzip = "bcftools/consensus_inversion_hetINDELs_satellite.vcf.gz",
    threads: 1
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bgzip -c {input.indels} > {output.bgzip}"

# Step 7
rule bcftools_tabix_inversion_hetINDELs_Satellite:
    input:
        bgzip = "bcftools/consensus_inversion_hetINDELs_satellite.vcf.gz",
    output:
        tabix = "bcftools/consensus_inversion_hetINDELs_satellite.vcf.gz.tbi",
    threads: 1
    priority: 1000
    conda:
        "../envs/bcftools.yaml",
    shell:
        "tabix -p vcf {input.bgzip}"

# Step 8
rule bcftools_bgzip_inversion_hetSNPs_Satellite:
    input:
        snps = "gatk4/ruff_adults_variants/consensus_inversion_hetSNPs_satellite.vcf",
    output:
        bgzip = "bcftools/consensus_inversion_hetSNPs_satellite.vcf.gz",
    threads: 1
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bgzip -c {input.snps} > {output.bgzip}"

# Step 9
rule bcftools_tabix_inversion_hetSNPs_Satellite:
    input:
        bgzip = "bcftools/consensus_inversion_hetSNPs_satellite.vcf.gz",
    output:
        tabix = "bcftools/consensus_inversion_hetSNPs_satellite.vcf.gz.tbi",
    threads: 1
    priority: 1000
    conda:
        "../envs/bcftools.yaml",
    shell:
        "tabix -p vcf {input.bgzip}"

# Step 10
rule bcftools_concat_hetVARIANTS_Satellite:
    input:
        snps = "bcftools/consensus_inversion_hetSNPs_satellite.vcf.gz",
        indels = "bcftools/consensus_inversion_hetINDELs_satellite.vcf.gz",
    output:
        merged = "bcftools/consensus_inversion_hetSNPs_and_hetINDELs_merged_satellite.vcf",
    params:
        extra = "--allow-overlaps",
    threads: 1
    priority: 500
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bcftools concat {params.extra} {input.snps} {input.indels} -o {output.merged}"