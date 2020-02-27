#ifndef __QALSH_H
#define __QALSH_H

#include <vector>

class BNode;
class BLeafNode;
class BTree;
class MinK_List;

// -----------------------------------------------------------------------------
//  PageBuffer: a buffer of one page for c-k-ANN search
// -----------------------------------------------------------------------------
struct PageBuffer {
	BLeafNode *leaf_node_;			// leaf node (level = 0)
	int index_pos_;					// cur pos of index key
	int leaf_pos_;					// cur pos of leaf node
	int size_;						// size for one scan
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
		const char *index_path);		// index path

	// -------------------------------------------------------------------------
	int load(						// load index
		const char *index_path);		// index path

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	long long knn(					// k-NN search
		int top_k,						// top-k value
		const float *query,				// query object
		const char *data_folder,		// data folder
		MinK_List *list);				// k-NN results (return)

	// -------------------------------------------------------------------------
	long long knn(					// k-NN search
		int top_k,						// top-k value
		float R,						// limited search range
		const float *query,				// query object
		const std::vector<int> &object_id, // object id mapping
		const char *data_folder,		// data folder
		MinK_List *list);				// k-NN results (return)
		
protected:
	int   n_pts_;					// cardinality
	int   dim_;						// dimensionality
	int   B_;						// page size
	float p_;						// the p value of L_p norm
	float zeta_;					// symmetric factor of p-stable distr.
	float appr_ratio_;				// approximation ratio
	char  index_path_[300];			// folder path of index

	float w_;						// bucket width
	float p1_;						// positive probability
	float p2_;						// negative probability
	float alpha_;					// collision threshold percentage
	float beta_;					// false positive percentage
	float delta_;					// error probability
	int   m_;						// number of hashtables
	int   l_;						// collision threshold
	float *a_array_;				// hash functions
	BTree **trees_;					// b-trees

	int   dist_io_;					// io for computing distance
	int   page_io_;					// io for scanning pages
	int   *freq_;					// frequency of data objects
	bool  *checked_;				// whether the data objects are checked
	bool  *bucket_flag_;			// flag of bucket width
	bool  *range_flag_;				// flag of search range
	float *data_;					// one data object
	float *q_val_;					// hash value of query
	
	PageBuffer **lptr_;				// left  pointer of B+ Tree
	PageBuffer **rptr_;				// right pointer of B+ Tree

	// -------------------------------------------------------------------------
	void calc_params();				// calc parameters

	// -------------------------------------------------------------------------
	float calc_l0_prob(				// calc <p1> and <p2> for L_{0.5} distance
		float x);						// x = w / (2.0 * r)

	float calc_l1_prob(				// calc <p1> and <p2> for L_{1.0} distance
		float x);						// x = w / (2.0 * r)

	float calc_l2_prob(				// calc <p1> and <p2> for L_{2.0} distance
		float x);						// x = w / (2.0 * r)

	// -------------------------------------------------------------------------
	void gen_hash_func();			// generate hash functions
	
	// -------------------------------------------------------------------------
	int bulkload(					// build B+ Trees by bulkloading
		const float **data);			// data set

	// -------------------------------------------------------------------------
	float calc_hash_value(			// calc hash value
		int tid,						// hash table id
		const float *point);			// one data object

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
		const float *query);			// query object

	// -------------------------------------------------------------------------
	float find_radius();				// find proper radius

	// -------------------------------------------------------------------------
	float update_radius(			// update radius
		float old_radius);				// old radius

	// -------------------------------------------------------------------------
	void update_left_buffer(		// update left buffer
		const PageBuffer *rptr,			// right buffer
		PageBuffer *lptr);				// left buffer (return)

	// -------------------------------------------------------------------------
	void update_right_buffer(		// update right buffer
		const PageBuffer *lptr,			// left buffer
		PageBuffer* rptr);				// right buffer (return)

	// -------------------------------------------------------------------------
	float calc_dist(				// calc projected distance
		float q_val,					// hash value of query
		const PageBuffer *ptr);			// page buffer
	
	// -------------------------------------------------------------------------
	void delete_tree_ptr();			// delete the pointers of B+ Trees
};

#endif // __QALSH_H
