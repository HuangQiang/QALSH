#ifndef __ANN_H
#define __ANN_H

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <unordered_map>

#include "def.h"
#include "util.h"
#include "pri_queue.h"
#include "qalsh.h"
#include "qalsh_plus.h"

struct Result;

// -----------------------------------------------------------------------------
int linear_scan(					// brute-force linear scan (data in disk)
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   B,							// page size
	float p,							// the p value of Lp norm, p in (0,2]
	const float **query,				// query set
	const Result **R,					// truth set
	const char *data_folder,			// data folder
	const char *output_folder);			// output folder

// -----------------------------------------------------------------------------
int indexing_of_qalsh_plus(			// indexing of qalsh+
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	int   leaf,							// leaf size of kd-tree
	int   L,							// number of projection (drusilla)
	int   M,							// number of candidates (drusilla)
	float p,							// the p value of Lp norm, p in (0,2]
	float zeta,							// symmetric factor of p-stable distr.
	float ratio,						// approximation ratio
	const float **data,					// data set
	const char *output_folder);			// output folder

// -----------------------------------------------------------------------------
int knn_of_qalsh_plus(				// k-NN search of qalsh+
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const float **query,				// query set
	const Result **R,					// truth set
	const char *data_folder,			// data folder
	const char *output_folder);			// output folder

// -----------------------------------------------------------------------------
int indexing_of_qalsh(				// indexing of qalsh
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	float p,							// the p value of Lp norm, p in (0,2]
	float zeta,							// symmetric factor of p-stable distr.
	float ratio,						// approximation ratio
	const float **data,					// data set
	const char *output_folder);			// output folder

// -----------------------------------------------------------------------------
int knn_of_qalsh(					// k-NN search of qalsh
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const float **query,				// query set
	const Result **R,					// truth set
	const char *data_folder,			// data folder
	const char *output_folder);			// output folder

#endif // __ANN_H
