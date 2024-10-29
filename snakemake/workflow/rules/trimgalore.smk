# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# RNA-Seq libraries quality trimming and adapter sequences removal using Trim Galore!

rule trim_galore_pe:
    input:
        ["/group/ag_nowick/data/Alex_RNASeq_Ruff/RNASeq_Adults/RNASeq_Blood_Adults/{sample}_R1.fastq.gz", 
        "/group/ag_nowick/data/Alex_RNASeq_Ruff/RNASeq_Adults/RNASeq_Blood_Adults/{sample}_R2.fastq.gz"],
    output:
        "trimgalore/{sample}_R1_val_1.fq.gz",
        "trimgalore/{sample}_R1.fastq.gz_trimming_report.txt",
        "trimgalore/{sample}_R2_val_2.fq.gz",
        "trimgalore/{sample}_R2.fastq.gz_trimming_report.txt",
    params:
        extra = "--illumina -q 30 --length 50",
    threads: 4
    log:
        "logs/trim_galore/{sample}.log",
    wrapper:
        "v1.7.0/bio/trim_galore/pe"
