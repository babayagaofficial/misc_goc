library(ape)

clades <- c("st73_cl733_community_0_subcommunity_493", "st73_cl568_community_0_subcommunity_503", "st95_cl378_community_0_subcommunity_503", "st73_cl429_community_0_subcommunity_501", "st73_cl703_community_0_subcommunity_493", "st131_cl494_community_0_subcommunity_490", "st95_cl400_community_0_subcommunity_503", "st69_cl493_community_0_subcommunity_487", "st73_cl597_community_0_subcommunity_501", "st95_cl461_community_0_subcommunity_503", "st131_cl416_community_0_subcommunity_503", "st95_cl605_community_0_subcommunity_498", "st73_cl613_community_0_subcommunity_501", "st95_cl389_community_0_subcommunity_502", "st131_cl386_community_0_subcommunity_485", "st69_cl461_community_0_subcommunity_489", "st131_cl406_community_0_subcommunity_496", "st73_cl455_community_0_subcommunity_494", "st69_cl375_community_0_subcommunity_499")

for (clade in clades){
tree <- read.tree(paste0("/home/daria/Documents/projects/ABC/pangraph_workflow/fitch/relabelled_trees/", clade, ".nw"))  

png(file=paste0("tree_vis/og/", clade, ".png"))

plot(tree, show.tip.label = TRUE)

dev.off()

edge_table <- read.csv(paste0("/home/daria/Documents/projects/ABC/pangraph_workflow/spp_dcj_sol/res_adj/",clade,"_edists.txt"), sep=" ", header=FALSE, col.names=c("child","parent","dist"))

# Step 1: Create a lookup table from node names to internal node numbers
node_names <- c(tree$tip.label, tree$node.label)  # all node labels
node_ids <- 1:(length(tree$tip.label) + tree$Nnode)
name_to_id <- setNames(node_ids, node_names)

# Step 2: Map parent and child names in your table to node numbers
edge_table$parent_id <- name_to_id[as.character(edge_table[[2]])]
edge_table$child_id  <- name_to_id[as.character(edge_table[[1]])]

# Step 3: Match to tree$edge
# For each row in tree$edge, find the matching row in edge_table
edge_indices <- match(
  paste(tree$edge[, 1], tree$edge[, 2], sep = "-"),
  paste(edge_table$parent_id, edge_table$child_id, sep = "-")
)

# Step 4: Assign distances
tree$edge.length <- edge_table[[3]][edge_indices]

png(file=paste0("tree_vis/dcj/",clade,".png"))

plot(tree, show.tip.label = TRUE)

edgelabels(tree$edge.length)

dev.off()
}