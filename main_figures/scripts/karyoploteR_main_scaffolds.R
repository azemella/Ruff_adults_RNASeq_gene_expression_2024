# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

### This script was used to generate the left-side illustrations in Fig.1C-E
### Run it three times after changing lines #56, #57 and #69

# Set working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Load required R libraries
library("DESeq2")
library("BiocParallel")
library("AnnotationDbi")
library("GenomicFeatures")
library("karyoploteR")
library("BSgenome.Cpugnax.NCBI.ASM143184v1")
library("plyranges")

# Set up parallel computing
register(MulticoreParam(4))

# Read metadata
samples <- read.csv(file.path(dir, "metadata/metadata.csv"))
samples$Phenotype <- factor(samplesPhenotype)
samples$Tissue <- factor(samples$Tissue)
rownames(samples) <- samples$Sample_ID

# Read count matrix
count_matrix <- read.csv("DESeq2_gene_expression_analysis/STAR_gene_counts.csv", row.names = 1)
count_matrix <- count_matrix[, rownames(samples)]

# Define experimental design
design <- ~ Tissue + Phenotype
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = samples, design = design)

# Filter low-count genes
ind_sample_count <- sum(samples$Phenotype == "Independent")
sat_sample_count <- sum(samples$Phenotype == "Satellite")
fae_sample_count <- sum(samples$Phenotype == "Faeder")
smallestGroupSize <- min(ind_sample_count, sat_sample_count, fae_sample_count)
keep <- rowSums(counts(dds) >= 5) >= smallestGroupSize
dds <- dds[keep,]

# Perform differential expression analysis
dds <- DESeq(dds, parallel = TRUE)
dds <- dds[which(mcols(dds)$betaConv),]

# Load genomic features
gff_file <- file.path(dir, "metadata/NCBI_RefSeq_annotations.gtf")
txdb <- makeTxDbFromGFF(gff_file, format = "gtf")
cp.genes <- genes(txdb)
cp.genes

# Analyze differential expression results
res <- results(dds, contrast = c("Phenotype", "Faeder", "Satellite"), alpha = 0.05, parallel = TRUE)
res <- lfcShrink(dds, contrast = c("Phenotype", "Faeder", "Satellite"), res = res, type = "ashr", parallel = TRUE)

# Prepare results for plotting on karyotype
cp.genes_res <- genes(txdb)
mcols(cp.genes_res) <- res[names(cp.genes), c("log2FoldChange", "pvalue", "padj")]
cp.genes_names <- join_overlap_left(cp.genes, cp.genes_res)

# Read Calidris pugnax BSgenome and gene locations
genome <- getBSgenome(genome = "BSgenome.Cpugnax.NCBI.ASM143184v1", masked = FALSE)
custom.cytobands <- read.delim("others/Cpugnax_genome_cytobands.txt", header = T, sep = "\t")

# Set up SVG output
svg("main_figures/plots/karyoploteR_main_scaffolds_faeder_vs_satellite.svg", 
    width = 13, height = 6)

# Define parameters for plotting
filtered.cp.genes <- cp.genes_names[!is.na(cp.genes_names$padj)]
log.pval <- -log10(filtered.cp.genes$padj)
mcols(filtered.cp.genes)$log.pval <- log.pval
sign.genes <- filtered.cp.genes[filtered.cp.genes$padj < 0.05,]
sign.genes <- unique(sign.genes)
sign.genes_df <- as.data.frame(sign.genes)
sign.genes_df$log.pval <- ifelse(sign.genes_df$log.pval > 30, 30, sign.genes_df$log.pval)
sign.genes_df <- sign.genes_df %>% mutate(color = ifelse(log2FoldChange > 0, adjustcolor("red", alpha.f = 0.5), adjustcolor("blue", alpha.f = 0.5)))
fc.ymax <- 8
fc.ymin <- -fc.ymax
cex.val <- sqrt(sign.genes_df$log.pval)/1.5
points.top <- 0.4
sign.col <- sign.genes_df$color
regs <- toGRanges("NW_015090842.1:5801914-10348254")

# Set karyoploteR plot parameters
pp <- getDefaultPlotParams(plot.type = 2)
pp$leftmargin <- 0.065

# Plot scaffolds of similar length as the one carrying the autosomal inversion for comparison
kp <- plotKaryotype(genome = genome, cytobands = custom.cytobands, 
                    chromosomes = c("NW_015090839.1",
                                    "NW_015094507.1",
                                    "NW_015090893.1",
                                    "NW_015090842.1",
                                    "NW_015090812.1",
                                    "NW_015090995.1",
                                    "NW_015090990.1",
                                    "NW_015090810.1"), plot.type = 2, labels.plotter = NULL, plot.params = pp)

kpAddChromosomeNames(kp, chr.names = c(" s56", "s219", " s27", " s28", " s18", " s38", "s103", " s13"), cex = 1.7)
#Data panel 1
kpPlotRegions(kp, data = regs, r0 = 0, r1 = 0.75, col = lighter("purple", amount = 200), border = lighter("purple", amount = 100), lwd = 0.75)
kpPoints(kp, data = sign.genes, y = sign.genes$log2FoldChange, cex = cex.val, ymax = fc.ymax, ymin = -fc.ymax, 
         r0 = points.top, r1 = points.top, col = sign.col)
kp <- kpPlotDensity(kp, data = cp.genes, window.size = 10e4, data.panel = 2)

# Close SVG device
dev.off()
