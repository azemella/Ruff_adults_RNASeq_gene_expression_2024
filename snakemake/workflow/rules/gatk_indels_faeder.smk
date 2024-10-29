# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Heterozygous INDELs were extracted for the Faeder inversion morph.
# In step 1, all the bi-allelic Het (GT 0/1) and HomVar INDELs (GT 1/1) were extracted from the Faeder population.
# In step 2, recommended GATK filters were applied, and then HomRef INDELs (GT 0/0) were removed from the Faeder VCF.
# In step 3, the HomVar INDELs (GT 1/1) were set to Heterozygous INDELs (GT 0/1) for the Faeder VCF.
# In step 4, Het INDELs found in more than 5 Faeder samples were kept, and the flagged genotypes were transformed with a null genotype call.
# In steps 5-8, the same approach was employed to retain INDELs from the Independent population
# In steps 9 and 10, bcftools was used to compress and tabix-index the Faeder and Independent variant calling files.
# In step 11, overlapping INDELs between the Faeder and Independent populations were filtered out.
# In step 12, I deleted unnecessary output files and renamed the consensus variant calling file.
# The final output of this file is a variant calling file of consensus Het Faeder INDELs.

# NOTE:
# Homozygous variants for the alternate allele (GT 1/1) have been manually changed to heterozygous variants in the 4th step.
# This was only possible due to the Faeder and Satellite morphs being heterozygous obliged for the autosomal inversion.
# The output file has been named "transformed_inversion_hetINDELs_faeder.vcf".

# Step 1
rule SelectVariants_inversion_HetINDELs_and_HomVarINDELs_Faeder:
    input:
        vcf = "gatk4/genotypegvcfs/genotypegvcfs.vcf",
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
    output:
        indels = "gatk4/ruff_adults_variants/selectvariants_inversion_hetINDELs_and_homvarINDELs_faeder.vcf",
    params:
        java = "--java-options '-Xmx10G'",
        extra0 = "--select-type-to-include 'INDEL'",
        extra1 = "--tmp-dir 'gatk4/ruff_adults_variants'",
        extra2 = "--sample-name 'config/samples_faeder.args'",
        extra3 = "-select 'vc.getHetCount() >= 1 || vc.getHomVarCount() >= 1'",
        extra4 = "--restrict-alleles-to 'BIALLELIC'",
        extra5 = "--intervals 'config/inversion_genes_coordinates.list'",
    threads: 10
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} SelectVariants -R {input.genome} -V {input.vcf} -O {output.indels}"
        " {params.extra0} {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5}"

# Step 2
rule VariantFiltration_HetINDELs_and_HomVarINDELs_Faeder:
    input:
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
        vcf = "gatk4/ruff_adults_variants/selectvariants_inversion_hetINDELs_and_homvarINDELs_faeder.vcf",
    output:
        indels = "gatk4/ruff_adults_variants/variantfiltration_inversion_hetINDELs_and_homvarINDELs_faeder.vcf",
    params:
        java = "--java-options '-Xmx10G'",
        # GATK recommended parameters
        extra0 = "-filter 'QD < 2.0' --filter-name 'QD2'",
        extra1 = "-filter 'QUAL < 30.0' --filter-name 'QUAL30'",
        extra2 = "-filter 'FS > 200.0' --filter-name 'FS200'",
        extra3 = "-filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20'",
        # Extra parameters
        extra4 = "-filter 'DP < 3.0' --filter-name 'DP3'",
        extra5 = "--genotype-filter-expression 'isHomRef == 1'",
        extra6 = "--genotype-filter-name 'isHomRefFilter'",
    threads: 10
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} VariantFiltration -V {input.vcf} -R {input.genome} -O {output.indels} {params.extra0} {params.extra1}"
        " {params.extra2} {params.extra3} {params.extra4} {params.extra5} {params.extra6}"

# Step 3
rule sed_HomVarINDELs_to_HetINDELs_transformation_Faeder:
    input:
        vcf = "gatk4/ruff_adults_variants/variantfiltration_inversion_hetINDELs_and_homvarINDELs_faeder.vcf",
    output:
        out = "gatk4/ruff_adults_variants/transformed_inversion_hetINDELs_faeder.vcf",
    threads: 10
    shell:
        "sed 's/1\/1:/0\/1:/g; s/1|1:/0|1:/g' {input.vcf} > {output.out}"

# Step 4
rule SelectVariants_transformed_inversion_HetINDELs_Faeder:
    input:
        vcf = "gatk4/ruff_adults_variants/transformed_inversion_hetINDELs_faeder.vcf",
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
    output:
        indels = "gatk4/ruff_adults_variants/selectvariants_transformed_inversion_hetINDELs_faeder.vcf",
    params:
        java = "--java-options '-Xmx10G'",
        extra1 = "--tmp-dir 'gatk4/selectvariants'",
        extra2 = "-select 'vc.getHetCount() > 5'",
        extra3 = "--set-filtered-gt-to-nocall 'true'",
    threads: 10
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} SelectVariants -R {input.genome} -V {input.vcf} -O {output.indels} {params.extra1} {params.extra2} {params.extra3}"

# Step 5
rule SelectVariants_inversion_INDELs_Independent:
    input:
        vcf = "gatk4/genotypegvcfs/genotypegvcfs.vcf",
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
    output:
        indels = "gatk4/ruff_adults_variants/selectvariants_inversion_INDELs_independent.vcf",
    params:
        java = "--java-options '-Xmx10G'",
        extra0 = "--select-type-to-include 'INDEL'",
        extra1 = "--tmp-dir 'gatk4/selectvariants'",
        extra2 = "--sample-name 'config/samples_independent.args'",
        extra3 = "--intervals 'config/inversion_genes_coordinates.list'",
        extra4 = "--restrict-alleles-to 'BIALLELIC'",
        extra5 = "-select 'vc.getHetCount() >= 1 || vc.getHomVarCount() >= 1'",
    threads: 10
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} SelectVariants -R {input.genome} -V {input.vcf} -O {output.indels} {params.extra0}"
        " {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5}"

# Step 6
rule VariantFiltration_inversion_INDELs_Independent:
    input:
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
        vcf = "gatk4/ruff_adults_variants/selectvariants_inversion_INDELs_independent.vcf",
    output:
        indels = "gatk4/ruff_adults_variants/variantfiltration_inversion_INDELs_independent.vcf",
    params:
        java = "--java-options '-Xmx10G'",
        # GATK recommended parameters
        extra0 = "-filter 'QD < 2.0' --filter-name 'QD2'",
        extra1 = "-filter 'QUAL < 30.0' --filter-name 'QUAL30'",
        extra2 = "-filter 'FS > 200.0' --filter-name 'FS200'",
        extra3 = "-filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20'",
        # Extra parameters
        extra4 = "-filter 'DP < 3.0' --filter-name 'DP3'",
        extra5 = "--genotype-filter-expression 'isHomRef == 1'",
        extra6 = "--genotype-filter-name 'isHomRefFilter'",
    threads: 10
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} VariantFiltration -V {input.vcf} -R {input.genome} -O {output.indels} {params.extra0} {params.extra1}"
        " {params.extra2} {params.extra3} {params.extra4} {params.extra5} {params.extra6}"

# Step 7
rule sed_HomVarINDELs_to_HetINDELs_transformation_Independent:
    input:
        vcf = "gatk4/ruff_adults_variants/variantfiltration_inversion_INDELs_independent.vcf",
    output:
        out = "gatk4/ruff_adults_variants/transformed_inversion_INDELs_independent.vcf",
    threads: 10
    shell:
        "sed 's/1\/1:/0\/1:/g; s/1|1:/0|1:/g' {input.vcf} > {output.out}"

# Step 8
rule SelectVariants_inversion_INDELs_Independent_transformed:
    input:
        indels = "gatk4/ruff_adults_variants/transformed_inversion_INDELs_independent.vcf",
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
    output:
        filtered_indels = "gatk4/ruff_adults_variants/selectvariants_inversion_INDELs_independent_filtered.vcf",
    params:
        java = "--java-options '-Xmx10G'",
        extra1 = "--tmp-dir 'gatk4/selectvariants'",
        extra2 = "-select 'vc.getHetCount() >= 3'",
        extra3 = "--set-filtered-gt-to-nocall",
    threads: 10
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk {params.java} SelectVariants -V {input.indels} -R {input.genome} -O {output.filtered_indels} {params.extra1} {params.extra2} {params.extra3}"

# Step 9
rule bcftools_bgzip_inversion_HetINDELs_Faeder:
    input:
        indels = "gatk4/ruff_adults_variants/selectvariants_transformed_inversion_hetINDELs_faeder.vcf",
    output:
        bgzip = "gatk4/ruff_adults_variants/selectvariants_transformed_inversion_hetINDELs_faeder.vcf.gz",
    threads: 10
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bgzip -c {input.indels} > {output.bgzip}"

rule bcftools_bgzip_inversion_INDELs_Independent:
    input:
        indels = "gatk4/ruff_adults_variants/selectvariants_inversion_INDELs_independent_filtered.vcf",
    output:
        bgzip = "gatk4/ruff_adults_variants/selectvariants_inversion_INDELs_independent_filtered.vcf.gz",
    threads: 10
    conda:
        "../envs/bcftools.yaml",
    shell:
        "bgzip -c {input.indels} > {output.bgzip}"

# Step 10        
rule bcftools_tabix_inversion_HetINDELs_Faeder:
    input:
        bgzip = "gatk4/ruff_adults_variants/selectvariants_transformed_inversion_hetINDELs_faeder.vcf.gz",
    output:
        tabix = "gatk4/ruff_adults_variants/selectvariants_transformed_inversion_hetINDELs_faeder.vcf.gz.tbi",
    threads: 10
    priority: 100
    conda:
        "../envs/bcftools.yaml",
    shell:
        "tabix -p vcf {input.bgzip}"

rule bcftools_tabix_inversion_INDELs_Independent:
    input:
        bgzip = "gatk4/ruff_adults_variants/selectvariants_inversion_INDELs_independent_filtered.vcf.gz",
    output:
        tabix = "gatk4/ruff_adults_variants/selectvariants_inversion_INDELs_independent_filtered.vcf.gz.tbi",
    threads: 10
    priority: 100
    conda:
        "../envs/bcftools.yaml",
    shell:
        "tabix -p vcf {input.bgzip}"

# Step 11
rule bcftools_intersect_consensus_inversion_hetINDELs_Faeder:
    input:
        vcf_fae = "gatk4/ruff_adults_variants/selectvariants_transformed_inversion_hetINDELs_faeder.vcf.gz",
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
        "bcftools isec {params.extra1} -p {input.dir} {input.vcf_fae} {input.vcf_ind}"

# Step 12
rule rename_output_consensus_inversion_hetINDELs_Faeder:
    input:
        vcf = "gatk4/ruff_adults_variants/0000.vcf",
        readme = "gatk4/ruff_adults_variants/README.txt",
        sites = "gatk4/ruff_adults_variants/sites.txt",
    output:
        out = "gatk4/ruff_adults_variants/consensus_inversion_hetINDELs_faeder.vcf",
    threads: 1
    shell:
        "mv {input.vcf} {output.out} | rm {input.sites} {input.readme}"