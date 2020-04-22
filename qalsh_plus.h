#ifndef __QALSH_PLUS_H
#define __QALSH_PLUS_H

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <unordered_map>
#include <vector>

#include "def.h"
#include "util.h"
#include "pri_queue.h"
#include "kd_tree.h"
#include "qalsh.h"

class MinK_List;
class QALSH;

// -----------------------------------------------------------------------------
//  Blocks: an block which stores hash tables for some of data objects
// -----------------------------------------------------------------------------
struct Blocks {
	int   n_pts_;
	int   *index_;
	QALSH *lsh_;

	Blocks() { n_pts_ = -1; index_ = NULL; lsh_ = NULL; }
	~Blocks() { if (lsh_ != NULL) { delete lsh_; lsh_ = NULL; } }
};

// -----------------------------------------------------------------------------
class QALSH_PLUS {
public:
	QALSH_PLUS();					// default constructor
	~QALSH_PLUS();					// destructor

	// -------------------------------------------------------------------------
	int build(						// build index	
		int   n,						// number of data objects
		int   d,						// dimensionality
		int   B,						// page size
		int   leaf,						// leaf size of kd-tree
		int   L,						// number of projection (drusilla)
		int   M,						// number of candidates (drusilla)
		float p,						// l_p distance
		float zeta,						// a parameter of p-stable distr.
		float ratio,					// approximation ratio
		const float **data,				// data objects
		const char *path);				// index path

	// -------------------------------------------------------------------------
	int load(						// load index
		const char *path);				// index path

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	uint64_t knn(					// k-NN search
		int top_k,						// top-k value
		int nb,							// number of blocks for search
		const float *query,				// input query
		const char *data_folder,		// data folder
		MinK_List *list);				// top-k results (return)

	// -------------------------------------------------------------------------
	int   n_pts_;					// number of data objects
	int   dim_;						// dimensionality
	int   B_;						// page size
	int   L_;						// number of projection (drusilla)
	int   M_;						// number of candidates (drusilla)
	float p_;						// l_p distance
	int   num_blocks_;				// number of blocks 
	char  path_[200];				// index path

	int   *sample_id_;				// sample data id
	int   *new_order_id_;			// new order data id after kd-tree partition
	QALSH *lsh_;					// index of sample data objects
	std::unordered_map<int, int> sample_id_to_block_; // sample data id to block
	std::vector<Blocks*> blocks_;	// blocks

protected:
	// -------------------------------------------------------------------------
	void kd_tree_partition(			// kd-tree partition 
		int leaf,						// leaf size of kd-tree
		const float **data,				// data objects
		std::vector<int> &block_size,	// block size (return)
		int *new_order_id);				// new order id (return)

	// -------------------------------------------------------------------------
	void drusilla_select(			// drusilla select
		int   n,						// number of data objects
		const int *new_order_id,		// new order data id
		const float **data,				// data objects
		int   *sample_id);				// sample data id (return)

	// -------------------------------------------------------------------------
	void calc_shift_data(			// calculate shift data objects
		int   n,						// number of data objects
		const int *new_order_id,		// new order data id
		const float **data,				// data objects
		int   &max_id,					// data id with max l2-norm (return)
		float &max_norm,				// max l2-norm (return)
		float *norm,					// l2-norm of shift data (return)
		float **shift_data); 			// shift data (return)

	// -------------------------------------------------------------------------
	int write_params();				// write parameters

	// -------------------------------------------------------------------------
	int read_params();				// read parameters

	// -------------------------------------------------------------------------
	int get_block_order(			// get block order
		int nb,							// number of blocks for search
		MinK_List *list,				// top-t results from sample data
		std::vector<int> &block_order);	// block order (return)
};

#endif // __QALSH_PLUS_H
