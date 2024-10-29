# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)
        
# Heterozygous INDELs were extracted for the Satellite inversion morph.
# In step 1, all the bi-allelic Het (GT 0/1) and HomVar INDELs (GT 1/1) were extracted from the Satellite population.
# In step 2, recommended GATK filters were applied, and then HomRef INDELs (GT 0/0) were removed from the Satellite VCF.
# In step 3, the HomVar INDELs (GT 1/1) were set to Heterozygous INDELs (GT 0/1) for the Satellite VCF.
# In step 4, Het INDELs found in at least 5 Satellite samples were kept, and the flagged genotypes were transformed with a null genotype call.
# In step 5, the Satellite hetINDELs VCF file was compressed.
# In step 6, the Satellite hetINDELs VCF file was indexed.
# In step 7, overlapping INDELs between the Satellite and Independent populations were removed.
# In step 8, unnecessary output files were deleted.
# The final output of this file is a VCF of consensus Heterozygous Satellite INDELs.

# NOTE1:
# Homozygous variants for the alternate allele (GT 1/1) have been changed to heterozygous variants (GT 0/1).
# This was only possible due to the Faeder and Satellite morphs being heterozygous for the autosomal inversion.

# NOTE2:
# This file contains fewer steps than the Faeder one because Independent INDELs files had been prepared in the other one.

# Step 1
rule SelectVariants_inversion_HetINDELs_and_HomVarINDELs_Satellite:
    input:
        vcf = "gatk4/genotypegvcfs/genotypegvcfs.vcf",
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
    output:
        snps = "gatk4/ruff_adults_variants/selectvariants_inversion_hetINDELs_and_homvarINDELs_satellite.vcf",
    params:
        java = "--java-options '-Xmx10G'",
        extra0 = "--select-type-to-include 'INDEL'",
        extra1 = "--tmp-dir 'gatk4/selectvariants'",
        extra2 = "--sample-name 'config/samples_satellite.args'",
        extra3 = "--intervals 'config/inversion_genes_coordinates.list'",
        extra4 = "--restrict-alleles-to 'BIALLELIC'",
        extra5 = "-select 'vc.getHetCount() >= 1 || vc.getHomVarCount() >= 1'",
    threads: 10
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} SelectVariants -R {input.genome} -V {input.vcf} -O {output.snps} {params.extra0}"
        " {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5}"

# Step 2
rule VariantFiltration_inversion_hetINDELs_and_homvarINDELs_Satellite:
    input:
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
        vcf = "gatk4/ruff_adults_variants/selectvariants_inversion_hetINDELs_and_homvarINDELs_satellite.vcf",
    output:
        snps = "gatk4/ruff_adults_variants/variantfiltration_inversion_hetINDELs_and_homvarINDELs_satellite.vcf",
    params:
        java = "--java-options '-Xmx10G'",
        # GATK recommended parameters
        extra0 = "-filter 'QD < 2.0' --filter-name 'QD2'",
        extra1 = "-filter 'QUAL < 30.0' --filter-name 'QUAL30'",
        extra2 = "-filter 'SOR > 3.0' --filter-name 'SOR3'",
        extra3 = "-filter 'FS > 60.0' --filter-name 'FS60'",
        extra4 = "-filter 'MQ < 40.0' --filter-name 'MQ40'",
        extra5 = "-filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5'",
        extra6 = "-filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8'",
        # Extra parameters
        extra7 = "-filter 'DP < 3.0' --filter-name 'DP3'",
        extra8 = "--genotype-filter-expression 'isHomRef == 1'",
        extra9 = "--genotype-filter-name 'isHomRefFilter'",
    threads: 10
    conda:
        "../envs/gatk4.yaml", 
    shell:
        "gatk {params.java} VariantFiltration -V {input.vcf} -R {input.genome} -O {output.snps} {params.extra0} {params.extra1} {params.extra2}"
        " {params.extra3} {params.extra4} {params.extra5} {params.extra6} {params.extra7} {params.extra8} {params.extra9}"

# Step 3
rule sed_HomVarINDELs_to_HetINDELs_transformation_Satellite:
    input:
        vcf = "gatk4/ruff_adults_variants/variantfiltration_inversion_hetINDELs_and_homvarINDELs_satellite.vcf",
    output:
        out = "gatk4/ruff_adults_variants/transformed_inversion_hetINDELs_satellite.vcf",
    threads: 10
    shell:
        "sed 's/1\/1:/0\/1:/g; s/1|1:/0|1:/g' {input.vcf} > {output.out}"

# Step 4
rule SelectVariants_transformed_inversion_HetINDELs_Satellite:
    input:
        vcf = "gatk4/ruff_adults_variants/transformed_inversion_hetINDELs_satellite.vcf",
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
    output:
        snps = "gatk4/ruff_adults_variants/selectvariants_transformed_inversion_hetINDELs_satellite.vcf",
    params:
        java = "--java-options '-Xmx10G'",
        extra1 = "--tmp-dir 'gatk4/selectvariants'",
        extra2 = "-select 'vc.getHetCount() > 5'",
        extra3 = "--set-filtered-gt-to-nocall",
    threads: 10
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} SelectVariants -R {input.genome} -V {input.vcf} -O {output.snps} {params.extra1} {params.extra2} {params.extra3}"

# Step 5
rule bcftools_bgzip_inversion_HetINDELs_Satellite:
    input:
        snps = "gatk4/ruff_adults_variants/selectvariants_transformed_inversion_hetINDELs_satellite.vcf",
    output:
        bgzip = "gatk4/ruff_adults_variants/selectvariants_transformed_inversion_hetINDELs_satellite.vcf.gz",
    threads: 10
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bgzip -c {input.snps} > {output.bgzip}"

# Step 6
rule bcftools_tabix_inversion_HetINDELs_Satellite:
    input:
        bgzip = "gatk4/ruff_adults_variants/selectvariants_transformed_inversion_hetINDELs_satellite.vcf.gz",
    output:
        tabix = "gatk4/ruff_adults_variants/selectvariants_transformed_inversion_hetINDELs_satellite.vcf.gz.tbi",
    threads: 10
    priority: 100
    conda:
        "../envs/bcftools.yaml",
    shell:
        "tabix -p vcf {input.bgzip}"

# Step 7
rule bcftools_intersect_consensus_inversion_hetINDELs_Satellite:
    input:
        vcf_sat = "gatk4/ruff_adults_variants/selectvariants_transformed_inversion_hetINDELs_satellite.vcf.gz",
        vcf_ind = "gatk4/ruff_adults_variants/selectvariants_inversion_INDELs_independent_filtered.vcf.gz",
        dir = "gatk4/ruff_adults_variants",
    output:
        out = "gatk4/ruff_adults_variants/0000.vcf",
        readme = "gatk4/ruff_adults_variants/README.txt",
        sites = "gatk4/ruff_adults_variants/sites.txt",
    params:
        extra1 = "--complement",
    threads: 10
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bcftools isec {params.extra1} -p {input.dir} {input.vcf_sat} {input.vcf_ind}"

# Step 8
rule rename_output_consensus_inversion_hetINDELs_Satellite:
    input:
        vcf = "gatk4/ruff_adults_variants/0000.vcf",
        readme = "gatk4/ruff_adults_variants/README.txt",
        sites = "gatk4/ruff_adults_variants/sites.txt",
    output:
        out = "gatk4/ruff_adults_variants/consensus_inversion_hetINDELs_satellite.vcf",
    threads: 1
    shell:
        "mv {input.vcf} {output.out} | rm {input.sites} {input.readme}"