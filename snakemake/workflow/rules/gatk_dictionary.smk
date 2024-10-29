# Author: Alex Zemella
# Copyright: Copyright 2023, Alex Zemella
# Email: alex.zemella@bi.mpg.de
# License: Max Planck Institute for Biological Intelligence (MPIBI)

# Rule to create a sequence dictionary file for the Ruff's genome.

rule CreateSequenceDictionary:
    input:
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
    output:
        dict = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.dict",
    threads: 10
    conda:
        "../envs/gatk4.yaml",
    shell:
        "gatk CreateSequenceDictionary -R {input.genome} -O {output.dict}"
