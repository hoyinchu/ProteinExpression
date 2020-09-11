library(WGCNA)
library(data.table)
pearson_corr_df <- read.csv('./pairwise_csv/proteomeHD_pearson_corr.csv')
pearson_corr_mat <- as.matrix(pearson_corr_df)
colnames(pearson_corr_mat) <- colnames(pearson_corr_df)
row.names(pearson_corr_mat) <- colnames(pearson_corr_df)
pearson_corr_mat[is.nan(pearson_corr_mat)] <- 0
pearson_corr_dist_mat <- pearson_corr_mat^2 # so it ranges between 0 and 1
pearson_corr_dist_mat_adj <- sigmoidAdjacencyFunction(pearson_corr_dist_mat)

## Get the Topological Overlap Matrix
pearson_tom_sim <- TOMsimilarity(pearson_corr_dist_mat_adj, TOMDenom = "mean")
colnames(pearson_tom_sim) <- colnames(pearson_corr_df)
row.names(pearson_tom_sim) <- colnames(pearson_corr_df)

## Test if network is scale-free
connectivity <- colSums(pearson_tom_sim, na.rm = TRUE)
scaleFreePlot(connectivity)

## Turn similarity matrices into long-format, remove duplicates and merge into final table
pearson_sim_dt <- as.data.table( reshape2::melt(pearson_corr_dist_mat)) 
pearson_sim_dt <- pearson_sim_dt[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), pearson_sim = value ) ]
pearson_sim_dt <- pearson_sim_dt[ Protein_1 > Protein_2 ]                  # Removes duplicate pairs (incl. self-comparisons)

pearson_adj_mat_dt <- as.data.table( reshape2::melt( pearson_corr_dist_mat_adj )) 
pearson_adj_mat_dt <- pearson_adj_mat_dt[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), pearson_adj = value ) ]
pearson_adj_mat_dt <- pearson_adj_mat_dt[ Protein_1 > Protein_2 ]                # Removes duplicate pairs (incl. self-comparisons)

pearson_tom_dt <- as.data.table( reshape2::melt( pearson_tom_sim )) 
pearson_tom_dt <- pearson_tom_dt[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), pearson_tom = value ) ]
pearson_tom_dt <- pearson_tom_dt[ Protein_1 > Protein_2 ]                  # Removes duplicate pairs (incl. self-comparisons)

pearson_dt <- merge( pearson_sim_dt, pearson_adj_mat_dt, by = c("Protein_1", "Protein_2"))
pearson_dt <- merge( pearson_dt, pearson_tom_dt,      by = c("Protein_1", "Protein_2"))

## Write out the combined result file
fwrite(pearson_dt, "./tom_optimized_csv/pearson_tom.csv")

## Write out simplified result file containing only the final scores
pearson_tom_dt_final <- pearson_dt[, .(Protein_1, Protein_2, pearson_tom_score = pearson_tom)]
fwrite(pearson_tom_dt_final, "./tom_optimized_csv/pearson_tom_scores.csv")
