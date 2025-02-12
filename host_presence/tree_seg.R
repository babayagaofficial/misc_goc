library(treeSeg)
library(ape)
library(stepR)
#library(treeio)

args <- commandArgs()
print(args)
subcom <- args[2]
st <-  args[3]

path <- "/hps/nobackup/iqbal/daria/game_of_clones/treeseg"
out_dir <- paste0(path, "/st", st, "/", subcom)
dir.create(paste0(path, "/st", st))
dir.create(out_dir)

# Input
tree <- read.tree(paste0("host_tree/pruned_dated_", st, ".nwk"))  # Load your tree
print(is.recursive(tree))
all_traits <- read.csv("presence_per_host.tsv", row.names = 1, sep = "\t")  # Load your traits matrix
traits <- all_traits[, subcom]

names(traits) <- rownames(all_traits)  # Ensure species names are set

traits <- traits[tree$tip.label]

seg <- treeSeg(traits, tree, alpha = 0.05)

#store output

output <- data.frame(Genome_ID=tree$tip.label,mLP=seg$mlP,lower_conf=seg$confBandP[,1], upper_conf=seg$confBandP[,2])
write.csv(x = output, file = paste0(out_dir, "/probabilities.csv"))

dir.create(paste0(out_dir, "/clades"))
for(i in seq_along(seg$mlAN)){
  node <- seg$mlAN[i]
  members <- extract.clade(tree, node)$tip.label
  writeLines(members, paste0(out_dir, "/clades/", node, ".txt"))
}

#plot

pdf(file=paste0(out_dir,"/tree.pdf"))

n <- length(tree$tip.label)
lwdEdge <- rep(1.5, dim(tree$edge)[1])
layout(mat=matrix(1:3,ncol = 1),heights = rep(c(3,0.1, 1),6))
par(mar = c(0.1,4,0.1,0.1))

# plot tree
plot(tree, type = "phylogram", show.tip.label = F , use.edge.length = F, 
     node.pos = 1, direction = "downwards", show.node.label = F, edge.width = lwdEdge) 

# plot maximum likelihood estimate for the nodes where the distribution of phenotype changes.
for(i in seq_along(seg$mlAN)){
  nodelabels("", seg$mlAN[i], pch = 18, col = "red", frame = "none",cex = 2)
}
text(180,-8, 'community_0_subcommunity_503')

# plot tip phenos
plot(1:n, rep(0,n), axes = F, col = c("lightgray", "black")[traits+1], pch = "|", cex = 2,xlab = '',ylab = '') 

# plot the maximum likelihood estimates for the phenotype distributions for each segment.
par(mar = c(2,4,0.1,0.1))
plot(1:n, seg$mlP,type = "l", lwd = 1, axes = F,ylim = c(-0.1,1.1),col='red',xlab = '',ylab = '')

# plot the confidence interval for the phenotype distributions.
polygon(c(1:n, rev(1:n)), c(seg$confBandP[,2], rev(seg$confBandP[,1])),
        col = rgb(255, 165, 0, alpha = 100, maxColorValue = 255), border = NA)

#axis(1, at = c(1,50,100,150,200), cex = 2)
axis(2, at = c(0,0.5,1),labels = c('0','0.5','1') )
mtext(side=2,'Prob.',line=2.5,cex=0.8)
points(1:n, rep(0,n), ylim = c(-0.1,1.1), type = "l", lwd = 1, axes = F)

dev.off()