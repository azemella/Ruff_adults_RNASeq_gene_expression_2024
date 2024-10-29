library("rstudioapi")

setwd("/home/alex-zemella/paml-4.10.7/paml-tutorial-main/positive-selection/")

lnL_vals <- read.table(file = "HSD17B2_positive_selection_analysis/Branch/Inversion_morphs_clade/9_LRT/lnL_inversion_morphs_clade_branch_mods.txt", 
                      sep = " ", stringsAsFactors = FALSE, header = FALSE)

rownames( lnL_vals ) <- c("BS_morphs", "BS-w1-morphs")

# Compute the LRT statistic.
diff_BSvsBSw1 <- 2*(lnL_vals[1,] - lnL_vals[2,])
# diff =  8.345988
pchisq(diff_BSvsBSw1, df = 1, lower.tail = F)
# p-val = 0.003865399

Chisq.crit1 <- 3.84
Chisq.crit2 <- 5.99

# Visualization results
curve( dchisq( x, df = 1 ), from = 0, to =  7 )
abline( v = c( Chisq.crit1, Chisq.crit2, diff_BSvsBSw1 ), col = c( "darkgray", "brown", "red" ) )
coords_dev     <- c( 1.35, 1.3 )
coords_pval    <- c( 1.25, 1.2 )
coords_alphac  <- c( 1.28, 1.0 )
coords_alphac2 <- c( 1.28, 0.9 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = 1.6604', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["1,0.05"]^"2", "= 3.84", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_alphac2[1], y = coords_alphac2[2],
      labels = expression( atop( paste( chi["1,0.01"]^"2", "= 5.99", sep = " " ) ) ),
      cex = 1.2, col = "brown" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 0.1975493' ) ),
      cex = 1.2, col = "black" )
title( expression( 'Faeder morph: branch model A vs. branch model A with '*omega*'=1' ) )
