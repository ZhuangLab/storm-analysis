/*
 * A C implementation of the DBSCAN algorithm as described here:
 * http://en.wikipedia.org/wiki/DBSCAN
 * 
 * 11/11
 *
 *
 * Modified to use a kdtree for region queries. This is way way
 * faster than the simple search algorithm I was using before.
 *
 * kdtree code from here:
 *   http://code.google.com/p/kdtree/
 *
 * 11/11
 *
 *
 * Added some additional functions that are somewhat related
 * to clustering.
 *
 * 11/11
 *
 *
 * Hazen
 * 
 * Linux:
 *  gcc -fPIC -g -c -Wall dbscan.c
 *  gcc -shared -Wl,-soname,dbscan.so.1 -o dbscan.so.1.0.1 dbscan.o kdtree.o -lc
 *
 * Windows:
 *  gcc -c dbscan.c
 *  gcc -shared -o dbscan.dll dbscan.o
 */

/* Includes */
#include <stdlib.h>
#include <stdio.h>

#include <assert.h>

#include "kdtree.h"

/* Defines */
#define TESTING 0
#define USEKD 1

#define NOISE -1
#define UNVISITED 0
#define VISITED 1

#define MAXSIZE 10000   // maximum cluster size.


/* Structures */
typedef struct
{
  int end;    // location of the last non-zero element in the array.
  int size;   // total size of the array.
  int *data;  // storage for array data.
} intArr;


/* Function Declarations */

// dbscan specific
void appendN(intArr *, intArr *);
intArr* createIntArr(int);
void dbscan(float *, float *, float *, int *, int *, int, float, int, int);
int expandCluster(intArr *, int, float, int, int);
void mergeN(intArr *, intArr *);
intArr* regionQuery(int *, int, int, int, float);
intArr* regionQueryKD(int *, int, int, int, float);

// additional
void clusterSize(int *, int *, int, int);
void locClSize(int *, int *, int, int);
void recategorize(int *, int *, int, int, int);


/* Global Variables */
int nsize; // size of c,l,x,y,z.
int *c;    // localization category.
int *l;    // localization label.
float *x;  // localization x location (in nm).
float *y;  // localization y location (in nm).
float *z;  // localization z location (in nm).

void *kd;  // kdtree structure pointer.


/* Code */


/*
 * appendN
 *
 * Append the contents of the second intArr to the 
 * the end of the first intArr.
 *
 * N - first intArr.
 * Np - second intArr.
 */
void appendN(intArr *N, intArr *Np)
{
  int i,n_end;
  int *n_data, *np_data;

  // initializations
  n_end = N->end;
  n_data = N->data;
  np_data = Np->data;

  // check size
  if((N->end+Np->end)>N->size){
    printf("N reached maxsize! %d (appendN)\n", N->size);
    return;
  }

  // append data
  for(i=0;i<(Np->end);i++){
    n_data[n_end] = np_data[i];
    n_end++;
  }
  N->end = n_end;
}


/*
 * createIntArr
 *
 * Allocates and initializes a intArr structure.
 *
 * n - size of the array.
 */
intArr* createIntArr(int n)
{
  intArr *ia;

  if(TESTING){
    printf("createIntArr\n");
  }

  ia = (intArr *)malloc(sizeof(intArr));
  ia->end = 0;
  ia->size = n;
  ia->data = (int *)malloc(sizeof(int)*n);

  return ia;
}


/*
 * dbscan
 *
 * Performs the DBSCAN algorithm to identify clusters.
 *
 * db_x - array of x locations.
 * db_y - array of y locations.
 * db_z - array of z locations.
 * db_cat - array localization categories.
 * db_l - pre-allocated storage for the cluster label.
 *        this should be the same size as x,y and z.
 * db_nsize - size of x, y, z and l.
 * eps - cluster size parameter (in nm).
 * min_points - minimum number of points in a cluster.
 * verbose - print cluster information as we go.
 */
void dbscan(float *db_x, float *db_y, float *db_z, int *db_cat, int *db_l, int db_nsize, float eps, int min_points, int verbose)
{
  int i,cluster_size,cn,counts,cur_cat;
  intArr *N;

  if(TESTING){
    printf("dbscan\n");
  }

  // initializations.
  nsize = db_nsize;
  l = db_l;
  x = db_x;
  y = db_y;
  z = db_z;
  c = db_cat;
  cn = 2;
  if(!USEKD){
    eps = eps*eps;
  }

  if(USEKD){
    kd = kd_create(3);
    for(i=0;i<db_nsize;i++){
      assert(kd_insert3(kd, db_x[i], db_y[i], db_z[i], (void *)(long)i) == 0);
    }
  }

  // mark all localizations as unvisited.
  for(i=0;i<nsize;i++){
    l[i] = UNVISITED;
  }

  // cluster localizations.
  for(i=0;i<nsize;i++){
    if(verbose){
      if((i%10000)==0){
	printf("Processing localization %d of %d\n", i, nsize);
      }
    }
    if(l[i]==UNVISITED){
      l[i]=VISITED;
      cur_cat = c[i];
      if(USEKD){
	N = regionQueryKD(&counts,i,min_points,cur_cat,eps);
      }
      else {
	N = regionQuery(&counts,i,min_points,cur_cat,eps);
      }
      if(counts>=min_points){
	l[i] = cn;
	cluster_size = expandCluster(N,cur_cat,eps,min_points,cn);
	if(verbose){
	  printf("  Cluster %d contains %d localizations\n", cn, cluster_size);
	}
	cn++;
      }
      else{
	l[i] = NOISE;
      }
      free(N->data);
      free(N);
    }
  }

  if(USEKD){
    kd_free(kd);
  }
}


/*
 * expandCluster
 *
 * Searches the neighbors of the current cluster to
 * to see if they also have enough neighbors to merit
 * additional analysis.
 * 
 * N - array of not yet visited points in the cluster.
 * cur_cat - localization category.
 * eps - cluster size parameter (in nm).
 * min_points - minimum number of neighbors.
 * cn - current cluster number.
 *
 * returns the cluster size.
 */
int expandCluster(intArr *N, int cur_cat, float eps, int min_points, int cn)
{
  int i,cluster_size,counts;
  int *data;
  intArr *Np;

  if(TESTING){
    printf("expandCluster\n");
  }

  cluster_size = 1;
  data = N->data;
  while((N->end)>0){

    // Mark current position with current cluster number.
    i = data[0];
    l[i] = cn;

    // Find the neighbors of the current position
    // if it has not already been determined to
    // have too few neighbors.
    if(USEKD){
      Np = regionQueryKD(&counts,i,min_points,cur_cat,eps);
    }
    else{
      Np = regionQuery(&counts,i,min_points,cur_cat,eps);
    }


    // If there are enough neighbors, merge N and Np.
    if(counts>=min_points){
      //mergeN(N,Np);
      appendN(N,Np);
    }

    // Free Np.
    free(Np->data);
    free(Np);

    // Remove current position from the N array.
    for(i=0;i<(N->end);i++){
      data[i]=data[i+1];
    }
    N->end--;
    cluster_size++;
  }

  return cluster_size;
}


/*
 * mergeN
 *
 * Merge the contents of a second intArr into the first
 * ignoring duplicates.
 *
 * N - first intArr.
 * Np - second intArr.
 */
void mergeN(intArr *N, intArr *Np)
{
  int i,j,cur,n_end,found;
  int *n_data, *np_data;

  if(TESTING){
    printf("mergeN\n");
  }

  n_end = N->end;
  n_data = N->data;
  np_data = Np->data;

  if(TESTING){
    printf(" N size %d\n", N->end);
    for(i=0;i<(N->end);i++){
      printf("  N %d %d\n", i, n_data[i]);
    }
    printf("\n");
    printf(" Np size %d\n", Np->end);
    for(i=0;i<(Np->end);i++){
      printf("  Np %d %d\n", i, np_data[i]);
    }
    printf("\n");
  }

  for(i=0;i<(Np->end);i++){
    cur = np_data[i];
    found = 0;
    for(j=0;j<n_end;j++){
      if(cur==n_data[j]){
	found = 1;
	break;
      }
    }
    if(!found){
      n_data[N->end] = cur;
      N->end++;
      if((N->end)>(N->size)){
	printf("N reached maxsize! %d (mergeN)\n", N->size);
	break;
      }
    }
  }

  if(TESTING){
    printf(" N size %d\n", N->end);
    for(i=0;i<(N->end);i++){
      printf("  N %d %d\n", i, n_data[i]);
    }
    printf("\n");
  }
}


/*
 * regionQuery
 *
 * Counts how many neighbors surrounding the current point are 
 * within eps of the current point. Returns a intArr structure
 * containing all the unvisited neighbors of the current point.
 *
 * counts - returns the total number of neighbors.
 * index - index of the current point.
 * min_points - minimum cluster size.
 * cur_cat - localization category.
 * eps - cluster size parameter (in nm).
 */
intArr *regionQuery(int *counts, int index, int min_points, int cur_cat, float eps)
{
  int i,cnt;
  float cx,cy,cz,dx,dy,dz;
  intArr *N;

  if(TESTING){
    printf("regionQuery\n");
  }

  // initialization.
  cnt = 0;
  cx = x[index];
  cy = y[index];
  cz = z[index];
  N = createIntArr(MAXSIZE);

  // loop over all the localizations.
  for(i=0;i<nsize;i++){
    if(c[i]==cur_cat){
      dx = cx - x[i];
      dy = cy - y[i];
      dz = cz - z[i];
      if((dx*dx+dy*dy+dz*dz)<eps){
	cnt++;
	if(l[i]==UNVISITED){
	  N->data[N->end] = i;
	  N->end++;
	  if((N->end)==(N->size)){
	    printf("N reached maxsize! %d (regionQuery)\n", N->size);
	    break;
	  }
	  l[i] = VISITED;
	}
      }
    }
  }

  *counts = cnt;
  return N;
}


/*
 * regionQueryKD
 *
 * Counts how many neighbors surrounding the current point are 
 * within eps of the current point. Returns a intArr structure
 * containing all the unvisited neighbors of the current point.
 * Uses a kd-tree to quickly find the points neighbors.
 *
 * counts - returns the total number of neighbors.
 * index - index of the current point.
 * min_points - minimum cluster size.
 * cur_cat - localization category.
 * eps - cluster size parameter (in nm).
 */
intArr *regionQueryKD(int *counts, int index, int min_points, int cur_cat, float eps)
{
  int i,j,cnt;
  void *set;
  intArr *N;

  if(TESTING){
    printf("regionQuery\n");
  }

  cnt = 0;
  N = createIntArr(MAXSIZE);

  // find the set of localizations within eps.
  set = kd_nearest_range3f(kd, x[index], y[index], z[index], eps);
  // printf("rqkd: %d %d %d %d\n", index, kd_res_size(set), c[index], cur_cat);

  // loop over the set to populate N.
  for(i=0;i<kd_res_size(set);i++){
    j = (int)kd_res_item_data(set);
    if(c[j]==cur_cat){
      cnt++;
      if((l[j]==UNVISITED)||(l[j]==NOISE)){
	N->data[N->end] = j;
	N->end++;
	if((N->end)==(N->size)){
	  printf("N reached maxsize! %d (regionQuery)\n", N->size);
	  break;
	}
	//l[j] = VISITED;
      }
    }
    kd_res_next(set);
  }
  kd_res_free(set);

  // if cnts is larger than the min cluster size
  // mark all members of N as visited.
  //
  // FIXME: Marks noise points as visited, which
  //        is not optimal since this means they
  //        will be visited twice.
  if(cnt>=min_points){
    for(i=0;i<(N->end);i++){
      l[N->data[i]] = VISITED;
    }
  }

  *counts = cnt;
  return N;
}


/********************************************************************/

/*
 * clusterSize
 *
 * Returns the size of each cluster.
 *
 * The fiddly adding 2 to max_id & 1 to the cluster id
 * is due to cluster id starting at -1 and going to max_id.
 *
 * counts - preallocated storage for cluster size.
 * cl_id - array of molecule cluster ids.
 * cl_size - size of cl_id.
 * max_id - maximum cluster id number (and size of cl_cnts).
 *
 */
void clusterSize(int *counts, int *cl_id, int cl_size, int max_id)
{
  int i;
  
  for(i=0;i<max_id;i++){
    counts[i] = 0;
  }

  for(i=0;i<cl_size;i++){
    counts[cl_id[i]+1]++;
  }

  // Zero out counts for cluster id numbers that are not valid.
  counts[0] = 0;
  counts[1] = 0; // should be zero anyway..
  counts[2] = 0;
}

/*
 * locClSize
 *
 * For each localization, calculates what size cluster it belongs to.
 *
 * cl_counts - pre-allocated storage for the cluster size.
 * cl_id - array of molecule cluster ids.
 * cl_size - size of cl_counts (and cl_id).
 * max_id - maximum cluster id number.
 */
void locClSize(int *cl_counts, int *cl_id, int cl_size, int max_id)
{
  int i;
  int *counts;

  max_id += 2;
  counts = (int *)malloc(sizeof(int)*max_id);
  clusterSize(counts, cl_id, cl_size, max_id);

  for(i=0;i<cl_size;i++){
    cl_counts[i] = counts[cl_id[i]+1];
  }

  free(counts);
}

/*
 * recategorize
 *
 * Moves all molecules in clusters that are too small to category 0.
 *
 * cl_id - array of molecule cluster ids.
 * category - array of molecule categories.
 * cl_size - size of cluster_id (and category).
 * max_id - maximum cluster id number.
 * min_size - minimum acceptable cluster size.
 */
void recategorize(int *cl_id, int *category, int cl_size, int max_id, int min_size)
{
  int i;
  int *counts;

  max_id += 2;
  counts = (int *)malloc(sizeof(int)*max_id);
  clusterSize(counts, cl_id, cl_size, max_id);

  for(i=0;i<cl_size;i++){
    if(counts[cl_id[i]+1]<min_size){
      category[i]=0;
    }
  }

  free(counts);
}

