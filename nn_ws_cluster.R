dyn.load( paste(dirname(sys.frame(1)$ofile), "src/nearest_neighbor.so", sep="/") )

## points: a matrix of points. The points should be in rows
##         with one column per dimension. If data is n-dimensional then it may
##         make sense to perform a PCA on it first.
##
## std: the standard deviation of the normal distribution used to define
##      a kernel that is used to calculated local densities at each point
define.neighbors <- function(points, std, scale=TRUE){
    if(!is.numeric(std))
        stop("std should be a numeric value")
    if(!is.double(std))
        std <- as.double(std)
    points <- as.matrix(points)
    if(!is.matrix(points) || !is.numeric(points))
        stop("points should be a numeric matrix")
    if(!is.double(points)) ## there really should be a better way
        points <- matrix(as.double(points), nrow=nrow(points), ncol=ncol(points),
                         dimnames=list(rownames(points), colnames(points)))
    if(scale)
        points <- scale(points)
    nn <- .Call("nearest_neighbor", t(points), std)
    list(pts=points, nbor=nn)
}

define.clusters <- function(nbors, max.d){
    if(!is.matrix(nbors) || !is.double(nbors))
        stop("nbors should be a matrix of real doubles")
    if(any( colnames(nbors) != c("nbor", "dist", "dens")))
        stop("unexpected column names")
    .Call("cluster_by_path", nbors, max.d)
}

cluster.cols <- function(clusters){
    clust.counts <- sort(table(clusters), decreasing=TRUE)
    clust.o <- as.integer(names(clust.counts))
    cols <- hsv(h=abs(sin(1:length(clust.o)/0.75)),
                s=length(clust.o):1/length(clust.o),
                v=0.4 + 0.5 * abs(cos(1:length(clust.o))))
    cols
}

draw.nn <- function(pts, pts.nn, max.d=0.25, dims=1:2, ar.lwd=1, ar.angle=10,
                    ar.length=0.1, pch=19, cex=0.5, clust=NULL, asp=NA, draw.ar=TRUE,
                    ...){
    pts <- pts[,dims]
    plot(pts, asp=asp, type='n', ...)
    i <- which(pts.nn[,'nbor'] > 0 & pts.nn[,'dist'] < max.d)
    j <- pts.nn[i,'nbor']
    if(draw.ar)
        arrows( pts[i,1], pts[i,2], pts[j,1], pts[j,2], lwd=ar.lwd, length=ar.length, angle=ar.angle )
    if(is.null(clust))
        points(pts, cex=cex, pch=pch)
    if(!is.null(clust)){
        clust.cols <- cluster.cols(clust)
        points(pts, cex=cex, pch=pch, col=clust.cols[ clust ])
        invisible( clust.cols )
    }
}

    
