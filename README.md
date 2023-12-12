# Nearest neighbor watershed clustering (nn_ws_cluster)

This repository provides functions that allow points in N-dimensional
space to be clustered into groups based on a nearest neighbour path
through the point set. Edges connect points to their nearest neighbor
that has a higher localised point density. The density at each point
is defined using a normal kernel defined by its standard deviation and
is in effect a Gaussian blur of all the points.

This means that points are naturally connected to their nearest high
density centers in a way that is analogous to watershed
segmentation. The advantage over other forms of clustering is that
there is no limitation to the shapes of clusters; they can arise from
any form of a connected space, and points that are very close to each
other can still be divided into separate clusters.

## Implementation

The primary functions that define the graph and use it do cluster
points are defined in `src/nearest_neighbor.c` and are compiled to a
shared object file using `R CMD SHLIB nearest_neighbor.c`. The
resulting shared object file is loaded by wrapper functions defined in
`nn_ws_cluster.R`.

Note that the function appears to run somewhat slower than I
expected. I have for this reason included a `Makevars` file that sets
the optimisation level for the compilation. However, it does not seem
to make that much difference.

## Wrapper functions

`nn_ws_cluster.R` defines two wrapper functions that call functions
defined in `src/nearest_neighbor.c`. These are:

### `define_neighbors(points, std, scale=TRUE)`

This takes three arguments:

1. points: a numeric matrix where each row defines one point in an
   n-dimensional space (one column per dimension).
2. std: the standard deviation of the kernel used to define local
   density. Smaller values result in smaller and more arbitrarily
   shaped clusters.
3. scale: A logical indicating whether the data should be scaled.
   Note that the standard deviation should fit the scale of the
   data used. If `scale` is `TRUE`, then columns of the data will
   be scaled using the built in `scale` function (this sets the mean
   to 0 and the standard devation to 1).

The function returns a list containing the points used (`pts`) and a
matrix giving edge information (`nbor`) in a named list. `$nbor` is
a matrix with three columns:

1. `nbor`: the index of the points pointed to from each point. Note that
   there will always be one point which doesn't point at any other points.
   This will be the one with the highest density in the set.
2. `dist`: the lengths of the edges.
3. `dens`: the density at each point.

### `define.clusters(nbors max.d)`

This function takes the `$nbors` matrix produced by `define.neighbors`
and a maximum distance. Edges are followed from nodes as long as their
lengths are smaller than or equal to `max.d`. The function returns a
numeric vector giving the resulting cluster of each point. Note that
edges that are distant from other points will end up in a large number
of singleton clusters.

## Plotting

`draw.nn(pts, pts.nn, max.d=0.25, dims=1:2, ar.lwd=1, ar.angle=10, ar.length=0.1, pch=19, cex=0.5, clust=NULL, asp=NA, ...)`

Will draw the points in their first two dimensions by default with
arrows connecting the points as defined in the `$nbor` matrix. The points
will be coloured by their cluster identity if a `clust` vector is passed
to the function.

![Example of proteins clustered by the scores of alignments to two different
descriptions of the same protein motif (SET-domain).](clustered_points.png)

## Bugs

There is an unfortunate problem if two points or more coexist at exactly the same point
in space. Since they will have exactly the same density it is not possible to add
them to the same cluster. The simple solution to allow addition of points that have
larger than or equal density doesn't work as that immediately leads to circularity
at these points. That means that a solution would need to keep track of all points that
have been added to a cluster. Unfortunately I have not come up with a good solution to
this problem. It is probably OK to simply add a little noise to the system or to detect
points that have 0 distance between them and perturbing one of them a little bit. But
that is not ideal either. 