# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set the working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")

# Load required libraries
library("tidyverse")
library("PMCMRplus")
library("stats")

# Read hormone data
data <- read.csv(file = "data/hormones_and_behaviors_data/Testes_and_blood_hormone_data.csv")

# Subset data by morph
ind_data <- subset.data.frame(data, morph=="Ind", select = c(gon_A4, gon_T, blo_A4, blo_T))
sat_data <- subset.data.frame(data, morph=="Sat", select = c(gon_A4, gon_T, blo_A4, blo_T))
fae_data <- subset.data.frame(data, morph=="Fae", select = c(gon_A4, gon_T, blo_A4, blo_T))

options(digits = 10)

# Create correlation matrices with p-values
# Independent
calculate_spearman <- function(ind_data) {
  n <- ncol(ind_data)
  cor_matrix <- matrix(NA, nrow = n, ncol = n)
  p_values <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      cor_result <- cor.test(ind_data[[i]], ind_data[[j]], method = "spearman", use = "pairwise.complete.obs")
      cor_matrix[i, j] <- cor_matrix[j, i] <- cor_result$estimate
      p_values[i, j] <- p_values[j, i] <- cor_result$p.value
    }
  }
  list(cor_matrix = cor_matrix, p_values = p_values)
}

# Spearman correlation coefficients and p-values (Independent)
spearman_results <- calculate_spearman(ind_data)

# Satellite
calculate_spearman <- function(sat_data) {
  n <- ncol(sat_data)
  cor_matrix <- matrix(NA, nrow = n, ncol = n)
  p_values <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      cor_result <- cor.test(sat_data[[i]], sat_data[[j]], method = "spearman", use = "pairwise.complete.obs")
      cor_matrix[i, j] <- cor_matrix[j, i] <- cor_result$estimate
      p_values[i, j] <- p_values[j, i] <- cor_result$p.value
    }
  }
  list(cor_matrix = cor_matrix, p_values = p_values)
}

# Spearman correlation coefficients and p-values (Satellite)
spearman_results <- calculate_spearman(sat_data)

# Faeder
calculate_spearman <- function(fae_data) {
  n <- ncol(fae_data)
  cor_matrix <- matrix(NA, nrow = n, ncol = n)
  p_values <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      cor_result <- cor.test(fae_data[[i]], fae_data[[j]], method = "spearman", use = "pairwise.complete.obs")
      cor_matrix[i, j] <- cor_matrix[j, i] <- cor_result$estimate
      p_values[i, j] <- p_values[j, i] <- cor_result$p.value
    }
  }
  list(cor_matrix = cor_matrix, p_values = p_values)
}

# Spearman correlation coefficients and p-values (Faeder)
spearman_results <- calculate_spearman(fae_data)

# Get exact p-values (Independent)
cor.test(ind_data$gon_A4, ind_data$gon_T, method = "spearman", use = "pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  ind_data$gon_A4 and ind_data$gon_T
# S = 56, p-value < 2.2204e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.931372549 

cor.test(ind_data$gon_A4, ind_data$blo_A4, method = "spearman", use = "pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  ind_data$gon_A4 and ind_data$blo_A4
# S = 32, p-value = 0.008235571
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.8060606061 

cor.test(ind_data$gon_A4, ind_data$blo_T, method = "spearman", use = "pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  ind_data$gon_A4 and ind_data$blo_T
# S = 116, p-value < 2.2204e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.8578431373 

cor.test(ind_data$gon_T, ind_data$blo_A4, method = "spearman", use = "pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  ind_data$gon_T and ind_data$blo_A4
# S = 12, p-value = 0.0001301624
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.9272727273 

cor.test(ind_data$gon_T, ind_data$blo_T, method = "spearman", use = "pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  ind_data$gon_T and ind_data$blo_T
# S = 54, p-value < 2.2204e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.9338235294 

cor.test(ind_data$blo_A4, ind_data$blo_T, method = "spearman", use = "pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  ind_data$blo_A4 and ind_data$blo_T
# S = 4, p-value < 2.2204e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.9757575758 

# Function to correct pvalues
bh_correction <- function(p_values) {
  n <- length(p_values)
  ranks <- rank(p_values)
  adj_p_values <- p_values * (n / ranks)
  adj_p_values[adj_p_values > 1] <- 1
  return(adj_p_values)
}

# Pairwise comparisons with p-value adjustment (Independent)
p_values <- c(2.22e-16, 0.008, 2.22e-16, 0.00013, 2.22e-16, 2.22e-16)
adjusted_p_values <- bh_correction(p_values)
# 2.22e-16 8.00e-03 2.22e-16 1.30e-04 2.22e-16 2.22e-16
# 5.328e-16 8.000e-03 5.328e-16 1.560e-04 5.328e-16 5.328e-16

# Get exact p-values (Satellite)
cor.test(sat_data$gon_A4, sat_data$gon_T, method = "spearman", use = "pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  sat_data$gon_A4 and sat_data$gon_T
# S = 20, p-value = 0.001977059
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.8787878788 

cor.test(sat_data$gon_A4, sat_data$blo_A4, method = "spearman", use = "pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  sat_data$gon_A4 and sat_data$blo_A4
# S = 24, p-value = 0.05758929
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.7142857143 

cor.test(sat_data$gon_A4, sat_data$blo_T, method = "spearman", use = "pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  sat_data$gon_A4 and sat_data$blo_T
# S = 112, p-value = 0.367684
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.3212121212 

cor.test(sat_data$gon_T, sat_data$blo_A4, method = "spearman", use = "pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  sat_data$gon_T and sat_data$blo_A4
# S = 42, p-value = 0.2161706
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.5 

cor.test(sat_data$gon_T, sat_data$blo_T, method = "spearman", use = "pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  sat_data$gon_T and sat_data$blo_T
# S = 94, p-value = 0.2180285
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.4303030303 

cor.test(sat_data$blo_A4, sat_data$blo_T, method = "spearman", use = "pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  sat_data$blo_A4 and sat_data$blo_T
# S = 58, p-value = 0.4618056
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.3095238095 

# Pairwise comparisons with p-value adjustment (Satellite)
p_values <- c(0.0019, 0.0576, 0.368, 0.216, 0.218, 0.462)
adjusted_p_values <- bh_correction(p_values)
# 2.22e-16 8.00e-03 2.22e-16 1.30e-04 2.22e-16 2.22e-16
# 5.328e-16 8.000e-03 5.328e-16 1.560e-04 5.328e-16 5.328e-16

# Get exact p-values (Faeder)
cor.test(fae_data$gon_A4, fae_data$gon_T, method = "spearman", use = "pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  fae_data$gon_A4 and fae_data$gon_T
# S = 12, p-value = 0.002028219
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.9 

cor.test(fae_data$gon_A4, fae_data$blo_A4, method = "spearman", use = "pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  fae_data$gon_A4 and fae_data$blo_A4
# S = 12, p-value = 0.04801587
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.7857142857 

cor.test(fae_data$gon_A4, fae_data$blo_T, method = "spearman", use = "pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  fae_data$gon_A4 and fae_data$blo_T
# S = 62, p-value = 0.1937996
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.4833333333 

cor.test(fae_data$gon_T, fae_data$blo_A4, method = "spearman", use = "pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  fae_data$gon_T and fae_data$blo_A4
# S = 6, p-value = 0.01230159
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.8928571429 

cor.test(fae_data$gon_T, fae_data$blo_T, method = "spearman", use = "pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  fae_data$gon_T and fae_data$blo_T
# S = 42, p-value = 0.06656195
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.65 

cor.test(fae_data$blo_A4, fae_data$blo_T, method = "spearman", use = "pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  fae_data$blo_A4 and fae_data$blo_T
# S = 16, p-value = 0.08809524
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.7142857143 

# Pairwise comparisons with p-value adjustment (Faeder)
p_values <- c(0.00203, 0.048, 0.194, 0.012, 0.06656, 0.088)
adjusted_p_values <- bh_correction(p_values)
# 0.00203 0.04800 0.19400 0.01200 0.06656 0.08800
# 0.01218 0.09600 0.19400 0.03600 0.09984 0.10560