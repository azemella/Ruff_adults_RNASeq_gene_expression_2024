# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")

# Load required libraries
library("AnnotationForge")

# Read input files
cSym <- read.csv(file.path(dir, "csym.csv"))
cTx <- read.csv(file.path(dir, "ctx.csv"))
cGOs <- read.csv(file.path(dir, "cgo.csv"))
cType <- read.csv(file.path(dir, "ctype.csv"))
cEns <- read.csv(file.path(dir, "cens.csv"))
cKegg <- read.csv(file.path(dir, "ckegg.csv"))
cAltSym <- read.csv(file.path(dir, "caltsym.csv"))

# Create custom ruff Org.Db object
makeOrgPackage(gene_info = cSym, 
               tx_info = cTx, 
               go = cGOs,
               gene_type = cType,
               ensembl_info = cEns,
               kegg_info = cKegg,
               alt_symbols = cAltSym,
               version = "0.1.1",
               maintainer = "Alex Zemella <azemella@bi.mpg.de>",
               author = "Alex Zemella <azemella@bi.mpg.de>",
               outputDir = ".",
               tax_id = "198806",
               genus = "Calidris",
               species = "pugnax",
               goTable = "go")

# Install object in R as a package
install.packages("./org.Cpugnax.eg.db", repos = NULL)
