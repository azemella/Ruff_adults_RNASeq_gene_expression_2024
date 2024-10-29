# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

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
samples$Phenotype <- factor(samples$Phenotype)
samples$Tissue <- factor(samples$Tissue)
rownames(samples) <- samples$Sample_ID

# Read count matrix
count_matrix <- read.csv("DESeq2_gene_expression_analysis/STAR_gene_counts.csv", row.names = 1)
count_matrix <- count_matrix[, rownames(samples)]

# Define experimental design
design <- ~ Tissue + Phenotype
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = samples, design = design)

# Filter low-count genes
ind_sample_count <- sum(samples$phenotype == "Independent")
sat_sample_count <- sum(samples$phenotype == "Satellite")
fae_sample_count <- sum(samples$phenotype == "Faeder")
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
res <- results(dds, contrast = c("Phenotype", "Satellite", "Independent"), alpha = 0.05, parallel = TRUE)
res <- lfcShrink(dds, contrast = c("Phenotype", "Satellite", "Independent"), res = res, type = "ashr", parallel = TRUE)

# Prepare results for plotting on karyotype
cp.genes_res <- genes(txdb)
mcols(cp.genes_res) <- res[names(cp.genes), c("log2FoldChange", "pvalue", "padj")]
cp.genes_names <- join_overlap_left(cp.genes, cp.genes_res)

# Read Calidris pugnax BSgenome and gene locations
genome <- getBSgenome(genome = "data/BSgenome.Cpugnax.NCBI.ASM143184v1", masked = FALSE)
custom.cytobands <- read.delim("others/Cpugnax_genome_cytobands.txt", header = T, sep = "\t")

# Set up SVG output
svg("main_figures/plots/karyoploteR_inversion_satellite_vs_independent_zoom.svg", 
    width = 28, height = 14)

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
cex.val <- sqrt(sign.genes_df$log.pval)*1.7
points.top <- 0.8
sign.col <- sign.genes_df$color
top.genes <- sign.genes[sign.genes$gene_id %in% c("TUBB3", "CENPN", "DPEP1", "HSD17B2", "BCO1")]
points.top <- 0.8
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2

# Set karyoploteR plot parameters
pp <- getDefaultPlotParams(plot.type = 4)
pp$leftmargin <- 0.085
pp$data1outmargin <- 0
pp$ideogramheight <- 2.5
pp$topmargin <- 0.01
pp$rightmargin <- 0.01
pp$data1inmargin <- 2.5
pp$data1height <- 100
pp$bottommargin <- 5
regs1 <- toGRanges("NW_015090842.1:5801914-7300084")
regs2 <- toGRanges("NW_015090842.1:7316310-8917817")
regs3 <- toGRanges("NW_015090842.1:8934521-9102274")
regs4 <- toGRanges("NW_015090842.1:9118156-9838999")
regs5 <- toGRanges("NW_015090842.1:9855539-10348254")
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2

# Plot inversion region
kp <- plotKaryotype(genome = genome, cytobands = custom.cytobands, chromosomes = "NW_015090842.1", 
                    plot.type = 4, plot.params = pp, labels.plotter = NULL, 
                    zoom = toGRanges(data.frame("NW_015090842.1", 5650000, 10790000)))

kpPlotRegions(kp, data = regs1, r0 = 0, r1 = 0.80, col = lighter("#d9d9d9", amount = 20))
kpPlotRegions(kp, data = regs2, r0 = 0, r1 = 0.80, col = lighter("#bdbdbd", amount = 20))
kpPlotRegions(kp, data = regs3, r0 = 0, r1 = 0.80, col = lighter("#d9d9d9", amount = 20))
kpPlotRegions(kp, data = regs4, r0 = 0, r1 = 0.80, col = lighter("#bdbdbd", amount = 20))
kpPlotRegions(kp, data = regs5, r0 = 0, r1 = 0.80, col = lighter("#d9d9d9", amount = 20))

kpAbline(kp, chr = "NW_015090842.1", v = c(5801914, 7300084), col=lighter("#d9d9d9", amount = 0), 
         lwd=1.5, lty=3, r1=0.801)
kpAbline(kp, chr = "NW_015090842.1", v = c(7316310, 8917817), col=lighter("#bdbdbd", amount = 0), 
         lwd=1.5, lty=3, r1=0.801)
kpAbline(kp, chr = "NW_015090842.1", v = c(8934521, 9102274), col=lighter("#d9d9d9", amount = 0), 
         lwd=1.5, lty=3, r1=0.801)
kpAbline(kp, chr = "NW_015090842.1", v = c(9118156, 9838999), col=lighter("#bdbdbd", amount = 0), 
         lwd=1.5, lty=3, r1=0.801)
kpAbline(kp, chr = "NW_015090842.1", v = c(9855539, 10348254), col=lighter("#d9d9d9", amount = 0), 
         lwd=1.5, lty=3, r1=0.801)

kpPoints(kp, data = sign.genes, y = sign.genes$log2FoldChange, cex = cex.val, ymax = fc.ymax, ymin = -fc.ymax, 
         r1 = points.top, col = sign.col)
kpAxis(kp, ymax = fc.ymax, ymin = fc.ymin, r1 = points.top, cex = 3.75, numticks = 5, lwd = 3.5)
kpAddLabels(kp, labels = expression('Log'[2]*'FC Satellite vs. Independent  '), srt = 90, pos = 1, label.margin = 0.07, 
            ymax = fc.ymax, ymin = fc.ymin, r0 = 0.15, r1 = points.top, cex = 5.5)
kpSegments(kp, chr = as.character(seqnames(top.genes)), x0 = gene.mean, x1 = gene.mean, y0 = top.genes$log2FoldChange, 
           y1 = 6, ymax = fc.ymax, ymin = fc.ymin, r1 = points.top, col = darker("#777777", amount = 25), lwd = 2.5)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0 = 0.6, r1 = 0.87, label.dist = 0.05, label.margin = 2, max.iter = 1000,
              label.color = "#000000", line.color = darker("#777777", amount = 25), marker.parts = c(0.7,0.3,0), srt = 30, cex = 4.75, font = 3, lwd = 2.5)

# Close SVG device
dev.off()
