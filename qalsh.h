#ifndef __QALSH_H
#define __QALSH_H

#include <vector>
#include <algorithm>
#include <cstring>
#include <vector>
#include <stdio.h>

#include "def.h"
#include "util.h"
#include "random.h"
#include "pri_queue.h"
#include "b_node.h"
#include "b_tree.h"

class BNode;
class BLeafNode;
class BTree;
class MinK_List;

// -----------------------------------------------------------------------------
//  Page: a buffer of one page for c-k-ANN search
// -----------------------------------------------------------------------------
struct Page {
	int size_;						// size for one scan
	int index_pos_;					// cur pos of index key
	int leaf_pos_;					// cur pos of leaf node
	BLeafNode *leaf_node_;			// leaf node (level = 0)
};

// -----------------------------------------------------------------------------
//  QALSH is used to solve the problem of high-dimensional c-Approximate 
//  Nearest Neighbor (c-ANN) search, where the hash tables of QALSH are 
//  indexed by B+ Tree. 
// -----------------------------------------------------------------------------
class QALSH {
public:
	QALSH();						// default constructor
	~QALSH();						// destructor

	// -------------------------------------------------------------------------
	int build(						// build index
		int   n,						// cardinality
		int   d,						// dimensionality
		int   B,						// page size
		float p,						// the p value of L_p norm
		float zeta,						// symmetric factor of p-stable distr.
		float ratio,					// approximation ratio
		const float **data,				// data objects
		const char *path);				// index path

	// -------------------------------------------------------------------------
	float calc_hash_value(			// calc hash value
		int tid,						// hash table id
		const float *point);			// one data object	

	// -------------------------------------------------------------------------
	int load(						// load index
		const char *path);				// index path

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	uint64_t knn(					// k-NN search
		int   top_k,					// top-k value
		const float *query,				// query object
		const char *data_folder,		// data folder
		MinK_List *list);				// k-NN results (return)

	// -------------------------------------------------------------------------
	uint64_t knn(					// k-NN search
		int   top_k,					// top-k value
		float R,						// limited search range
		const float *query,				// query object
		const int   *object_id, 		// object id mapping
		const char  *data_folder,		// data folder
		MinK_List *list);				// k-NN results (return)

	// -------------------------------------------------------------------------
	int   n_pts_;					// cardinality
	int   dim_;						// dimensionality
	int   B_;						// page size
	float p_;						// the p value of L_p norm
	float zeta_;					// symmetric factor of p-stable distr.
	float ratio_;					// approximation ratio
	float w_;						// bucket width
	int   m_;						// number of hashtables
	int   l_;						// collision threshold
	char  path_[300];				// index path

	float **a_;						// hash functions
	BTree **trees_;					// B+ Trees
	uint64_t dist_io_;				// io for computing distance
	uint64_t page_io_;				// io for scanning pages

protected:
	// -------------------------------------------------------------------------
	float calc_l0_prob(				// calc <p1> and <p2> for L_{0.5} distance
		float x);						// x = w / (2.0 * r)

	float calc_l1_prob(				// calc <p1> and <p2> for L_{1.0} distance
		float x);						// x = w / (2.0 * r)

	float calc_l2_prob(				// calc <p1> and <p2> for L_{2.0} distance
		float x);						// x = w / (2.0 * r)
	
	// -------------------------------------------------------------------------
	int bulkload(					// build B+ Trees by bulkloading
		const float **data);			// data set

	// -------------------------------------------------------------------------
	int write_params();				// write parameters to disk

	// -------------------------------------------------------------------------
	int read_params();				// read parameters from disk

	// -------------------------------------------------------------------------
	void get_tree_filename(			// get file name of tree
		int  tree_id,					// tree id
		char *fname);					// file name of tree (return)

	// -------------------------------------------------------------------------
	void init_search_params(		// init parameters for k-NN search
		const float *query,				// query object
		float *q_val,					// hash values of query (return)
		Page  **lptrs,					// left  buffer (return)
		Page  **rptrs);					// right buffer (return)

	// -------------------------------------------------------------------------
	float find_radius(				// find proper radius
		const float *q_val,				// hash value of query
		const Page **lptrs,			// left  buffer
		const Page **rptrs);			// right buffer

	// -------------------------------------------------------------------------
	float update_radius(			// update radius
		float old_radius,				// old radius
		const float *q_val,				// hash value of query
		const Page **lptrs,				// left  buffer
		const Page **rptrs);			// right buffer

	// -------------------------------------------------------------------------
	void update_left_buffer(		// update left buffer
		const Page *rptr,				// right buffer
		Page *lptr);					// left  buffer (return)

	// -------------------------------------------------------------------------
	void update_right_buffer(		// update right buffer
		const Page *lptr,				// left  buffer
		Page* rptr);					// right buffer (return)

	// -------------------------------------------------------------------------
	float calc_dist(				// calc projected distance
		float q_val,					// hash value of query
		const Page *ptr);				// page buffer
	
	// -------------------------------------------------------------------------
	void delete_tree_ptr(			// delete the pointers of B+ Trees
		Page **lptrs,					// left  buffer (return)
		Page **rptrs);					// right buffer (return)
};

#endif // __QALSH_H
