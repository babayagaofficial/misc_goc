# Load necessary libraries
library(phylofactor)
library(ape)  # For phylogenetic tree handling
library(ggtree)
library(ggplot2)

# Example input (replace with your own data)
tree <- read.tree("/home/daria/Documents/projects/ABC/host_presence/host_tree/pruned_95.nwk")  # Load your tree
print(is.recursive(tree))
all_traits <- read.csv("/home/daria/Documents/projects/ABC/host_presence/presence_per_host.tsv", row.names = 1, sep = "\t")  # Load your traits matrix
traits <- all_traits[, "community_0_subcommunity_503"]

names(traits) <- rownames(all_traits)  # Ensure species names are set
#traits <- as.data.frame(traits)

# Perform phylogenetic factorization
phylofactor_result <- twoSampleFactor(traits,tree, method = "Fisher", nfactors = 50)

gtree <- pf.tree(phylofactor_result)
gtree$ggplot
ggsave("plot.png")

