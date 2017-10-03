#ifndef __ANN_H
#define __ANN_H

// -----------------------------------------------------------------------------
int ground_truth(					// output ground truth (data in memory)
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	float p,							// the p value of Lp norm, p in (0,2]
	char* data_set,						// address of data set
	char* query_set,					// address of query set
	char* truth_set);					// address of ground truth file

// -----------------------------------------------------------------------------
int indexing(						// build hash tables for the dataset
	int   n,							// number of data points
	int   d,							// dimension of space
	int   B,							// page size
	float p,							// the p value of Lp norm, p in (0,2]
	float zeta,							// symmetric factor of p-stable distr.
	float ratio,						// approximation ratio
	char* data_set,						// address of data set
	char* data_folder,					// folder to store new format of data
	char* output_folder);				// folder to store info of qalsh

// -----------------------------------------------------------------------------
int lshknn(							// k-nn via qalsh (data in disk)
	int   qn,							// number of query points
	int   d,							// dimensionality
	float p,							// the p value of Lp norm, p in (0,2]
	char* query_set,					// path of query set
	char* truth_set,					// groundtrue file
	char* data_folder,					// folder to store new format of data
	char* output_folder);				// folder to store info of qalsh

// -----------------------------------------------------------------------------
int linear_scan(					// brute-force linear scan (data in disk)
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   B,							// page size
	float p,							// the p value of Lp norm, p in (0,2]
	char* query_set,					// address of query set
	char* truth_set,					// address of ground truth file
	char* data_folder,					// folder to store new format of data
	char* output_folder);				// folder to store info of qalsh


#endif
