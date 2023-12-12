source("nn_ws_cluster.R")

pts <- read.table("set_dom_pp.txt")

## std is the standard deviation of a normal distribution
## used as a kernel to define the local density at points
std <- 0.05
nn <- define.neighbors( pts, std=0.05 )
## nn is a list containing neighbor information (in $nbor)
## and potentially scaled points in $pts
## $nbor is a matrix of three columns giving
## edge (nbor), edge length (dist) and point density
## (dens) for each point.

## nn.clust follows edges that are no longer than max.d
## in order to connect members of clusters; It returns
## a numeric vector giving the cluster membership
max.d <- 0.25
nn.clust <- define.clusters( nn$nbor, max.d=max.d )

## draw.nn draws the points with arrows indicating edges
## and their directions.
png("clustered_points.png", width=1024, height=1024)
draw.nn(pts, nn$nbor, max.d=max.d, clust=nn.clust, cex=0.75)
dev.off()
