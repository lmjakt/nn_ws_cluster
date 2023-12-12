#include <R.h>
#include <Rinternals.h>

// Given a matrix of numeric (double) values:
// For each point find the closest point and its index
// and its distance. Note that both values will needt
// to be expressed as doubles.

// This is simply the Euclidean distance
double dist(double *r1, double *r2, size_t n){
  double d = 0;
  for(size_t i=0; i < n; ++i)
    d += (r1[i] - r2[i]) * (r1[i] - r2[i]);
  return(sqrt(d));
}

// A pure C-function..
// points should be in columns; there are n columns and dim_n rows.
// neighbors should be a matrix of n rows and two columns.
// this will do the stupid way that takes twice as much time to
// simplify the code.
void nearest_neighbors(double *points, size_t n, size_t dim_n, double *neighbors, double std){
  // shorthand for writing the data
  double *ni = neighbors;
  double *nd = neighbors + n;
  double *dens = neighbors + 2 * n;
  memset( neighbors, 0, sizeof(double) * n * dim_n );
  // We will use a normal distribution to get a kernel density for each point
  // The equation for the normal distribution is:
  // (1/std * sqrt(2 * pi)) * exp( -0.5 * (d/std)^2 )
  // the first part of this is constant
  double c = 1 / (2 * sqrt( 2 * M_PI) );
  double var = std * std;
  for(size_t i=0; i < n; ++i){
    double *pt1 = points + dim_n * i;
    dens[i] = 0;
    for(size_t j=0; j < n; ++j){
      double d = dist( pt1, points + dim_n * j, dim_n );
      dens[i] += ( c * exp( -0.5 * ( (d * d) / var )) );
    }
  }
  for(size_t i=0; i < n; ++i){
    double *pt1 = points + dim_n * i;
    double min_dist = -1;
    double min_dist_j = -1;
    for(size_t j=0; j < n; ++j){
      if(i == j)
	continue;
      double d = dist( pt1, points + dim_n * j, dim_n );
      if((min_dist < 0 || min_dist > d) && dens[j] > dens[i]){
	min_dist = d;
	min_dist_j = j;
      }
    }
    ni[i] = 1 + min_dist_j;
    nd[i] = min_dist;
  }
}

// A function to set column names; because it is more complicated than needed
void set_column_names(SEXP matrix, const char **names, int ncol){
  SEXP m_dim = getAttrib(matrix, R_DimSymbol);
  if(length(m_dim) != 2 || INTEGER(m_dim)[1] != ncol){
    warning("points should be a matrix with %d columns, not %d columns", 
	    ncol, INTEGER(m_dim)[1]);
    return;
  }
  SEXP dim_names = PROTECT(allocVector( VECSXP, 2 ));
  SET_VECTOR_ELT(dim_names, 1, allocVector(STRSXP, ncol));
  SEXP colnames = VECTOR_ELT(dim_names, 1);
  for(int i=0; i < ncol; ++i)
    SET_STRING_ELT( colnames, i, mkChar(names[i]));
  setAttrib( matrix, R_DimNamesSymbol, dim_names );
  UNPROTECT(1);
}

// walk the points to define clusters:
void cluster_path(double *nbor_matrix, double max_distance, int *clusters, int nrow){
  double *nbor = nbor_matrix;
  double *dist = nbor_matrix + nrow;
  // A temporary array that holds the path traversed in the current iteration.
  unsigned int *path = malloc(sizeof(unsigned int) * nrow);
  unsigned int path_n = 0;
  // Clear the clusters array
  memset(clusters, 0, sizeof(int) * nrow);
  int current_cluster = 0;
  for(int i=0; i < nrow; ++i){
    // if a cluster identity has already been set, then we skip to the next one
    if(clusters[i] > 0)
      continue;
    // otherwise we set the identity and follow the path until we get to the end
    // Note that the nbor is using 1 based counting so we need to decrement it here
    path_n = 0;
    current_cluster++;
    clusters[i] = current_cluster;
    path[path_n++] = i;
    int j = i;
    while( nbor[j] > 0 && dist[j] <= max_distance ){
      j = nbor[j] - 1;
      if( clusters[j] > 0 ){
	// set identity of full path to clusters[j], decrement current_cluster and break
	for(int k=0; k < path_n; ++k)
	  clusters[ path[k] ] = clusters[j];
	path_n = 0;
	current_cluster--;
	break;
      }
      clusters[j] = current_cluster;
      path[path_n++] = j;
    }
  }
  free(path);
}

SEXP nearest_neighbor(SEXP points_r, SEXP std_r){
  if(TYPEOF(std_r) != REALSXP || length(std_r) != 1)
    error("std_r should give the standard deviation of a normal distribution and be a single real value");
  if(TYPEOF(points_r) != REALSXP)
    error("points should be defined as a double matrix");
  SEXP points_dim = getAttrib(points_r, R_DimSymbol);
  if(length(points_dim) != 2)
    error("points should be a matrix");
  int nrow = INTEGER(points_dim)[0];
  int ncol = INTEGER(points_dim)[1];
  if(nrow < 1 | ncol < 2)
    error("There should be more than one point and at least one dimension");
  if(nrow > ncol)
    warning("There are more rows than columns; points should be in columns");
  double *points = REAL(points_r);
  double std = REAL(std_r)[0];
  // use a third column to give the density of the points
  SEXP ret_value = PROTECT( allocMatrix(REALSXP, ncol, 3) );
  const char *col_names[3] = {"nbor", "dist", "dens"};
  set_column_names( ret_value, col_names, 3 );
  nearest_neighbors( points, ncol, nrow, REAL(ret_value), std );
  UNPROTECT(1);
  return(ret_value);
}

SEXP cluster_by_path(SEXP nbors_r, SEXP max_d_r){
  if(TYPEOF(max_d_r) != REALSXP || length(max_d_r) != 1)
    error("max_d_r should give the maximum distance allowed and should be a single real value");
  if(TYPEOF(nbors_r) != REALSXP)
    error("nbors should be defined as a double matrix");
  SEXP dim = getAttrib(nbors_r, R_DimSymbol);
  if(length(dim) != 2)
    error("nbors should be a matrix");
  int nrow = INTEGER(dim)[0];
  int ncol = INTEGER(dim)[1];
  if(ncol != 3 || nrow < 2)
    error("nbors should be a matrix with 3 columns and at least 2 rows");
  double *nbors = REAL(nbors_r);
  double max_d = REAL(max_d_r)[0];
  // Allocate the return array and call cluster_path
  SEXP clusters_r = PROTECT(allocVector(INTSXP, nrow));
  cluster_path( nbors, max_d, INTEGER(clusters_r), nrow );
  UNPROTECT(1);
  return(clusters_r);
}

static const R_CallMethodDef callMethods[] = {
  {"nearest_neighbor", (DL_FUNC)&nearest_neighbor, 2},
  {"cluster_by_path", (DL_FUNC)&cluster_by_path, 2},
  {NULL, NULL, 0}
};

void R_init_nearest_neighbor(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
