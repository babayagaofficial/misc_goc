library(treeSeg)
library(ape)
library(stepR)
library(hash)

colours <- hash()
colours[["community_0_subcommunity_503"]] = "#ff28ba"
colours[["community_0_subcommunity_496"]] = "#1c007d"
colours[["community_0_subcommunity_502"]] = '#ba10c2'
colours[["community_0_subcommunity_500"]] = '#00650c'
colours[["community_0_subcommunity_493"]] = '#49caff'
colours[["community_0_subcommunity_494"]] = '#ae0069'
colours[["community_0_subcommunity_501"]] = '#35926d'
colours[["community_0_subcommunity_499"]] = '#20aa00'
colours[["community_0_subcommunity_498"]] = '#f7b69a'
colours[["community_0_subcommunity_485"]] = '#ffaeeb'
colours[["community_0_subcommunity_490"]] = '#1c5951'
colours[["community_0_subcommunity_487"]] = '#416dff'
colours[["community_0_subcommunity_489"]] = '#590000'
colours[["community_0_subcommunity_486"]] = '#ca2800'
colours[["community_0_subcommunity_474"]] = '#aeff0c'
colours[["community_0_subcommunity_453"]] = '#ff316d'
colours[["community_0_subcommunity_481"]] = '#510039'
colours[["community_0_subcommunity_470"]] = '#0096a6'
colours[["community_0_subcommunity_478"]] = '#65008e'
colours[["community_0_subcommunity_495"]] = '#0431ff'
colours[["community_0_subcommunity_482"]] = '#31e7ce'

pdf(file="tree.pdf")

# Define layout matrix for 2x2 panels, each with 2 rows (tree + trait)
layout_matrix <- rbind(
  matrix(c(
    1, 2,
    5, 6,
    3, 4,
    7, 8
  ), ncol = 2, byrow = TRUE),
  c(9, 9)  # <- Legend panel (row 5)
)


# Each entry is a plot slot (top: tree, bottom: trait), arranged row-wise
layout(layout_matrix, heights = rep(c(4, 1), 2))  # top plots taller than bottom


for (st in c(131,95,73,69)){
if (st==131) {subcomms <- c("community_0_subcommunity_485","community_0_subcommunity_490", "community_0_subcommunity_496", "community_0_subcommunity_503", "community_0_subcommunity_453", "community_0_subcommunity_470", "community_0_subcommunity_482" )
                location_1 <- c(1,1)
                location_2 <- c(3,1)}
else if (st==95) {subcomms <- c("community_0_subcommunity_498","community_0_subcommunity_500","community_0_subcommunity_501","community_0_subcommunity_502","community_0_subcommunity_503", "community_0_subcommunity_495")
                location_1 <- c(1,2)
                location_2 <- c(3,2)}
else if (st==73) {subcomms <- c("community_0_subcommunity_493","community_0_subcommunity_494","community_0_subcommunity_501","community_0_subcommunity_503","community_0_subcommunity_478")
                location_1 <- c(2,1)
                location_2 <- c(4,1)}
else if (st==69) {subcomms <- c("community_0_subcommunity_487","community_0_subcommunity_489","community_0_subcommunity_499", "community_0_subcommunity_474", "community_0_subcommunity_486")
                location_1 <- c(2,2)
                location_2 <- c(4,2)}


tree <- read.tree(paste0("/home/daria/Documents/projects/ABC/host_presence/host_tree/ST", st, "_100000000_1000_arc_BD.tre"))  # Load your tree
all_traits <- read.csv("/home/daria/Documents/projects/ABC/host_presence/presence_per_host.tsv", row.names = 1, sep = "\t")  # Load your traits matrix
dtip_dates <- read.csv("/home/daria/Documents/projects/ABC/goc/1-s2.0-S2666524721000318-mmc2.csv", sep = ",")
dtip_dates <- dtip_dates[c("lane","year")]
dtip_dates <- dtip_dates[!duplicated(dtip_dates),]
tip_dates <- as.vector(dtip_dates[, "year"])
names(tip_dates) <- dtip_dates$lane
tip_dates <- tip_dates[tree$tip.label]
traits <- hash()
for(i in seq_along(subcomms)){
  traits[[subcomms[i]]] <- as.vector(all_traits[, subcomms[i]])
  names(traits[[subcomms[i]]]) <- rownames(all_traits)  # Ensure species names are set
  traits[[subcomms[i]]] <- traits[[subcomms[i]]][tree$tip.label]
}

seg <- hash()
mlAN <- hash()
path <- paste0("/home/daria/Documents/projects/ABC/host_presence/phylofactor/st", st, "/")
for(i in seq_along(subcomms)){
  seg[[subcomms[i]]] <- read.csv(paste0(path, subcomms[i], "/rates.csv"), row.names = 1, header=TRUE)
  mlAN[[subcomms[i]]] <- as.integer(tools::file_path_sans_ext(list.files(path = paste0(path, subcomms[i],"/correct_clades"))))
}


# Distance from root to each tip (in time units)
depths <- node.depth.edgelength(tree)
tip_depths <- depths[1:length(tree$tip.label)]

# Match tree tips to their sample dates
tip_dates_ordered <- tip_dates[tree$tip.label]

print(tip_dates_ordered)

# Estimate root year: max tip date - branch length to tip
estimated_root_year <- max(tip_dates_ordered) - max(tip_depths)
print(estimated_root_year)

#plot

n <- length(tree$tip.label)
lwdEdge <- rep(0.8, dim(tree$edge)[1])
par(mfg=location_1, mar = c(0.1,4,1,0.1))

# plot tree
plot(tree, type = "phylogram", show.tip.label = F , node.pos = 1, direction = "downwards", show.node.label = F, edge.color = rep("#515151", dim(tree$edge)[1]),  edge.width = lwdEdge, main=st) 

# Get the full vertical range
tree_height <- max(node.depth.edgelength(tree))

# Build calendar year labels
start_year <- ceiling(estimated_root_year)
end_year <- floor(max(tip_dates_ordered))
print(start_year)
print(end_year)
year_labels <- seq(start_year, end_year, by = 20)
year_labels <- append(year_labels,end_year)
ticks <- end_year - year_labels

# Add vertical axis (left side) with calendar years
axis(2, at = ticks, labels = year_labels, las = 1, cex.axis=0.6)
mtext("Year", side = 2, line = 3,cex=0.6)

# plot maximum likelihood estimate for the nodes where the distribution of phenotype changes.
for(j in seq_along(subcomms)){
  nodelabels(text="", node = mlAN[[subcomms[j]]], pch = 18, col = colours[[subcomms[j]]], frame = "none",cex = 1.5)
}

# plot the maximum likelihood estimates for the phenotype distributions for each segment.
par(mfg=location_2, mar = c(2,4,0.1,0.1))
plot(1:n, seg[[subcomms[1]]]$rate,type = "l", lwd = 1, axes = F,ylim = c(-0.1,1.1),xlab = '',ylab = '')
for(j in seq_along(subcomms)){
  frame <- seg[[subcomms[j]]]
  lines(1:n, frame$rate, col=colours[[subcomms[j]]])
}


#axis(1, at = c(1,50,100,150,200), cex = 2)
axis(2, at = c(0,0.5,1),labels = c('0','0.5','1'), las = 1, cex.axis=0.6 )
mtext(side=2,'Proportion',line=2.5,cex=0.4)
points(1:n, rep(0,n), ylim = c(-0.1,1.1), type = "l", lwd = 1, axes = F)


}

names <- c("c_0_sc_503","c_0_sc_496","c_0_sc_502","c_0_sc_500","c_0_sc_493","c_0_sc_494","c_0_sc_501","c_0_sc_499","c_0_sc_498","c_0_sc_485","c_0_sc_490","c_0_sc_487","c_0_sc_489","c_0_sc_486","c_0_sc_474","c_0_sc_453","c_0_sc_481","c_0_sc_470","c_0_sc_478","c_0_sc_495","c_0_sc_482")
colours_vec <- c("#ff28ba","#1c007d",'#ba10c2','#00650c','#49caff','#ae0069','#35926d','#20aa00','#f7b69a','#ffaeeb','#1c5951','#416dff','#590000','#ca2800','#aeff0c','#ff316d','#510039','#0096a6','#65008e','#0431ff','#31e7ce')
par(mar = c(0.1, 1, 0.1, 1))
plot.new()
legend("center", legend = names, fill = colours_vec, ncol=7, cex = 1, bty = "n")

dev.off()