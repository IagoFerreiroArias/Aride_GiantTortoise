library(igraph)
library(bipartite)

# computeModules from bipartite return inconsistent number of modules, including its
# composition. Hence, we can not trust. Greedy algorithm for community detection
#returned similar values of modularity while keeping constant the number of modules
# as well as it composition. We will create a function to work with igraph objects
# applied to our case

# We retrieve the formula to compute z and c values from Olensen et al 2006.
# z = (Kis- mean(ks))/sd(ks),
# c = 1 - Σ(Kit/Ki)^2

# kis is number of links of i to other species in its own module "s",
# mean(ks) and sd(ks) are average and standard deviation of within-module k of all species in s
# kit is number of links from i to species in module t (including i's own module)
# ki is degree of species i

# Function to calculate c (among module connectivity) and z (within module connectivity) values
zc_igraph <- function(graph, membership_vector) {
  
  degree_vector <- degree(graph) #extract degree (n of links) from each node
  z_values <- numeric(vcount(graph)) #create vector to store z values
  c_values <- numeric(vcount(graph)) #create vector to store c values
  
  # Iterate over each node to calculate z and c values
  for (i in V(graph)) {
    
    node <- i #actual node
    module <- membership_vector[node] #extract membership of each node
    same_module_nodes <- V(graph)[membership_vector == module] #nodes within actual module
    
    kis <- sum(neighbors(graph, node) %in% same_module_nodes) #degree within module (kis)
    k_mean <- mean(degree(graph, v = same_module_nodes)) #mean of kis
    k_sd <- sd(degree(graph, v = same_module_nodes)) #sd of jis
    
    # Calculate within module connectivity: z-values. If sd = 0, return Inf, so set 0.
    if (k_sd > 0) {z_values[node] <- (kis - k_mean) / k_sd} else {z_values[node] <- 0}
    
    # Calculate among module connectivity: c values.
    total_degree <- degree_vector[node]
    sum_term <- 0
    
    for (mod in unique(membership_vector)) {
      nodes_in_module <- V(graph)[membership_vector == mod]
      kit <- sum(neighbors(graph, node) %in% nodes_in_module)
      sum_term <- sum_term + (kit / total_degree) ^ 2
    }
    c_values[node] <- 1 - sum_term
  }
  
  # Return a data.frame with z and c values 
  return(data.frame(name = V(graph)$name, z = z_values, c = c_values))
}

# First, calculate modularity and extract membership of each node 
# greddy <- cluster_fast_greedy(full_network)
# membership_vector <- membership(greddy)

#Then, apply function
# result <- zc_igraph(full_network, membership_vector)

