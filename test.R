dyn.load("src/nearest_neighbor.so")
pts <- read.table("set_dom_pp.txt")


source("nn_ws_cluster.R")

max.d=0.25
nn <- define.neighbors( pts, std=0.05 )
nn.clust <- define.clusters( nn$nbor, max.d=max.d )
draw.nn(pts, nn$nbor, max.d=max.d, clust=nn.clust, cex=0.75)

pts.sc <- scale(pts)

pts.nn <- .Call("nearest_neighbor", t(pts.sc), 0.1)

pts.clusters <- .Call("cluster_by_path", pts.nn, 0.25)

i <- which(pts.nn[,1] > 0 & pts.nn[,2] < 0.25)
j <- pts.nn[i,1]
plot(pts, type='n')
arrows( pts[i,1], pts[i,2], pts[j,1], pts[j,2], lwd=1, length=0.1, angle=20 )
points( pts, cex=0.1, col='red')

## we can see the biggest clusters from:
sort(table(pts.clusters))
b <- pts.clusters == 6 ## 5905 points
points( pts[b,], col='green', pch=19, cex=0.5 )

b <- pts.clusters == 4 ## 3684 points
points( pts[b,], col='blue', pch=19, cex=0.5 )

b <- pts.clusters == 7 ## 1504 points
points( pts[b,], col='cyan', pch=19, cex=0.5 )

b <- pts.clusters == 5 ## 959 points
points( pts[b,], col='purple', pch=19, cex=0.5 )

b <- pts.clusters == 1 ## 896 points
points( pts[b,], col='brown', pch=19, cex=0.5 )

b <- pts.clusters == 9 ## 382 points
points( pts[b,], col='yellow', pch=19, cex=0.5 )

b <- pts.clusters == 8 ## 354 points
points( pts[b,], col='dark blue', pch=19, cex=0.5 )

b <- pts.clusters == 3 ## 224 points
points( pts[b,], col='dark blue', pch=19, cex=0.5 )

b <- pts.clusters == 2 ## 146 points
points( pts[b,], col='dark green', pch=19, cex=0.5 )


## or we can try to be a bit more systematic:
cluster.order <- as.integer(names(sort(table(pts.clusters), decreasing=TRUE)))
cluster.cols <- hsv( h=abs(sin(1:length(cluster.order)/0.75)), s=(length(cluster.order):1) / length(cluster.order),
                    v=0.4 + 0.5 * abs(cos(1:length(cluster.order))) )

plot(pts.sc, asp=1, type='n')
i <- which(pts.nn[,'nbor'] > 0 & pts.nn[,'dist'] < 0.25)
j <- pts.nn[i,'nbor']
arrows( pts.sc[i,1], pts.sc[i,2], pts.sc[j,1], pts.sc[j,2], lwd=1, length=0.1, angle=20 )
points(pts.sc, asp=1, col=cluster.cols[pts.clusters], cex=0.5, pch=19)


## We can then try what happens with a smaller standard deviation for the kernel:
## don't understand why this is so slow.
pts.nn.2 <- .Call("nearest_neighbor", t(pts.sc), 0.05)
pts.clusters.2 <- .Call("cluster_by_path", pts.nn.2, 0.25)

cluster.order <- as.integer(names(sort(table(pts.clusters.2), decreasing=TRUE)))
cluster.cols <- hsv( h=abs(sin(1:length(cluster.order)/0.75)), s=(length(cluster.order):1) / length(cluster.order),
                    v=0.4 + 0.5 * abs(cos(1:length(cluster.order))) )

plot(pts.sc, asp=1, type='n')
i <- which(pts.nn.2[,'nbor'] > 0 & pts.nn.2[,'dist'] < 0.25)
j <- pts.nn.2[i,'nbor']
arrows( pts.sc[i,1], pts.sc[i,2], pts.sc[j,1], pts.sc[j,2], lwd=1, length=0.1, angle=20 )
points(pts.sc, asp=1, col=cluster.cols[pts.clusters.2], cex=0.5, pch=19)



pdf("points_arrows.pdf", width=10, height=10)
i <- which(pts.nn[,2] < 0.25 & pts.nn[,1] > 0)
j <- pts.nn[i,1]
plot(pts.sc, type='n', asp=1)
arrows( pts.sc[i,1], pts.sc[i,2], pts.sc[j,1], pts.sc[j,2], lwd=0.25, length=0.025, angle=10 )
points( pts.sc, cex=0.001, col='red', pch=19)
dev.off()

b <- pts.nn[,2] >= 0
plot(log10(1/pts.nn[b,2]), pts.nn[b,3], cex=0.5, xlab="1 / distance", ylab="density")

## x and y should be matrices of the same dimensions
edist <- function(x, y){
    sqrt( rowSums((x-y)^2) )
}

## this confirms that this is OK
dd <- edist( pts.sc, pts.sc[ pts.nn[,1], ] )
plot(dd, pts.nn[,2])
