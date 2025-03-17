# Script: Network metrics function
# Author: Iago Ferreiro-Arias
# Date: 15th July, 2024

# Define a function to retrieve network metrics from bipartite,
# perform bootstrap and extract CI values

library(dplyr)
library(bipartite)

network_metrics <- function(matrix, index, frac, resamp, null_model_method, N_permutations) {
  # Observed values
  Observed <- networklevel(matrix, index = index)
  
  boot_metrics <- NULL
  set.seed(123)
  
  # Bootstrap loop
  for (i in 1:resamp) {
    mat <- sample_frac(as.data.frame(matrix), frac, replace = TRUE) %>% as.matrix() # Sampling
    metrics_mat <- networklevel(mat, index = index) # Bootstrapped values
    boot_metrics <- rbind(boot_metrics, metrics_mat) # Bind results
  }
  
  boot_metrics <- as.data.frame(boot_metrics)
  colnames(boot_metrics) <- index # Set column names
  
  # Calculate mean, randomized SE, and CI levels
  boot_mean <- colMeans(boot_metrics)
  rnd.se <- apply(boot_metrics, 2, sd) / sqrt(resamp)
  ci_lower <- Observed - (1.962 * rnd.se)
  ci_upper <- Observed + (1.962 * rnd.se)
  
  # Null model significance
  nulls <- nullmodel(matrix, N = N_permutations, method = null_model_method)
  null_metrics <- sapply(1:N_permutations, function(i) {
    networklevel(nulls[[i]], index = index)
  })
  
  null_metrics <- as.data.frame(t(null_metrics))
  colnames(null_metrics) <- index
  
  null_mean <- colMeans(null_metrics)
  null_sd <- apply(null_metrics, 2, sd)
  z_scores <- (Observed - null_mean) / null_sd
  p_values <- 2 * pnorm(-abs(z_scores))
  
  # Store results in a single data.frame with the specified column order
  results <- data.frame(
    Metric = index,
    Observed = as.vector(Observed),
    Expected = null_mean, # Add expected values from chance (i.e, mean values from null models)
    Z_Score = z_scores,
    P_value = p_values,
    Boot_Mean = boot_mean,
    Boot_SE = rnd.se,
    CI_lower = ci_lower,
    CI_upper = ci_upper
  )
  
  rownames(results) <- NULL
  return(results)
}

# Example usage of the function
# Suppose that overhunted is the example matrix

#frac <- 0.9  # Sampling fraction for bootstrap
#resamp <- 100  # Number of resamples for bootstrap
#null_model_method <- 3  # Null model method (i.e., vaznull, see bipartite function)
#N_permutations <- 100  # Number of permutations for the null model

#result <- network_metrics(overhunted, frac = 0.9, 
#                          resamp = 100, null_model_method = 3, 
#                          N_permutations = 100,
#                          index =c("weighted connectance", "weighted NODF", "modularity"))

#print(result)
