#ifndef __UTIL_H
#define __UTIL_H

class MinK_List;

extern timeval g_start_time;		// global parameter: start time
extern timeval g_end_time;			// global parameter: end time

extern float g_runtime;				// global parameter: running time
extern float g_ratio;				// global parameter: overall ratio
extern float g_recall;				// global parameter: recall

// -----------------------------------------------------------------------------
//  struct Result
// -----------------------------------------------------------------------------
struct Result {						// structure for furthest neighbor / hash value
	float key_;							// distance / random projection value
	int   id_;							// object id
};

// -----------------------------------------------------------------------------
int ResultComp(						// compare function for qsort (ascending)
	const void *e1,						// 1st element
	const void *e2);					// 2nd element

// -----------------------------------------------------------------------------
int ResultCompDesc(					// compare function for qsort (descending)
	const void *e1,						// 1st element
	const void *e2);					// 2nd element

// -----------------------------------------------------------------------------
//  uitlity functions
// -----------------------------------------------------------------------------
void create_dir(					// create dir if the path exists
	char *path);						// input path

// -----------------------------------------------------------------------------
int read_data(						// read data/query set from disk
	int   n,							// number of data/query objects
	int   d,			 				// dimensionality
	const char *fname,					// address of data/query set
	float **data);						// data/query objects (return)

// -----------------------------------------------------------------------------
int read_ground_truth(				// read k-NN ground truth 
	int qn,								// number of query objects
	const char *fname,					// address of truth set
	Result **R);						// ground truth results (return)

// -----------------------------------------------------------------------------
float calc_lp_dist(					// calc L_{p} norm
	int   dim,							// dimension
	float p,							// the p value of Lp norm, p in (0, 2]
	const float *vec1,					// 1st point
	const float *vec2);					// 2nd point

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int k,								// top-k value
	const Result *R,					// ground truth results 
	MinK_List *list);					// results returned by algorithms

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int k,								// top-k value
	const Result *R,					// ground truth results 
	const Result *result);				// results returned by algorithms

// -----------------------------------------------------------------------------
int ground_truth(					// find ground truth
	int   n,							// number of data  objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float p,							// the p value of Lp norm, p in (0,2]
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set);			// address of truth set

// -----------------------------------------------------------------------------
void k_nn_search(					// k-NN search
	int   n, 							// cardinality
	int   qn,							// query number
	int   d, 							// dimensionality
	int   k,							// top-k value
	float p,							// the p value of Lp norm, p in (0,2]
	const float **data,					// data objects
	const float **query,				// query objects
	Result **result);					// k-MIP results (return)

#endif // __UTIL_H
