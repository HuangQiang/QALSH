#ifndef __QALSH_PLUS_H
#define __QALSH_PLUS_H

class QALSH;
class MinK_List;

// -----------------------------------------------------------------------------
//  Blocks: an block which stores hash tables for some of data objects
// -----------------------------------------------------------------------------
struct Blocks {
	int n_pts_;
	vector<int> index_;
	QALSH *lsh_;

	Blocks() { n_pts_ = -1; lsh_ = NULL; }
	~Blocks() { if (lsh_ != NULL) { delete lsh_; lsh_ = NULL; } }
};

// -----------------------------------------------------------------------------
//  QALSH_PLUS: an two-level LSH scheme for high-dimensional c-k-ANN search
// -----------------------------------------------------------------------------
class QALSH_PLUS {
public:
	QALSH_PLUS(						// constructor
		int   n,						// cardinality
		int   d,						// dimensionality
		int   leaf,						// leaf size of kd-tree
		int   L,						// number of projection
		int   M,						// number of candidates for each proj
		float p,						// l_p distance
		float zeta,						// a parameter of p-stable distr.
		float ratio,					// approximation ratio
		const float **data);			// data objects

	// -------------------------------------------------------------------------
	~QALSH_PLUS();					// destructor

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	int knn(						// k-NN seach	
		int   top_k,					// top-k value
		int   nb,						// number of blocks for search
		const float *query,				// input query object
		MinK_List *list);				// k-NN results (return)

protected:
	int   n_pts_;					// cardinality
	int   dim_;						// dimensionality
	int   leaf_;					// leaf size of kd-tree
	int   L_;						// number of projection (drusilla)
	int   M_;						// number of candidates (drusilla)
	float p_;						// l_p distance
	float zeta_;					// a parameter of p-stable distr.
	float appr_ratio_;				// approximation ratio

	int   num_blocks_;				// number of blocks 
	float **new_order_data_;		// new order data objects
	vector<Blocks*> blocks_;		// index of blocks
	
	int   sample_n_pts_;			// number of sample data objects
	int   *sample_id_to_block_;		// sample data id to block
	float **sample_data_;			// sample data objects
	QALSH *lsh_;					// index of sample data objects

	// -------------------------------------------------------------------------
	int bulkload(					// bulkloading
		const float **data);			// data objects

	// -------------------------------------------------------------------------
	int kd_tree_partition(			// kd-tree partition 
		const float **data,				// data objects
		vector<int> &block_size,		// block size (return)
		int *new_order_id);				// new order id (return)

	// -------------------------------------------------------------------------
	int calc_shift_data(			// calculate shift data objects
		int   n,						// number of data objects
		int   d,						// data dimension
		const float **data,				// data objects
		vector<vector<float> > &shift_data); // shift data objects (return)

	// -------------------------------------------------------------------------
	int drusilla_select(			// drusilla select
		int   n,						// number of data objects
		int   d,						// data dimension
		const vector<vector<float> > &shift_data, // shift data objects
		const int *new_order_id,		// new order data id
		int   *sample_id);				// sample data id (return)

	// -------------------------------------------------------------------------
	int get_block_order(			// get block order
		int nb,							// number of blocks for search
		MinK_List *list,				// top-t results from sample data
		vector<int> &block_order);		// block order (return)
};

#endif // QALSH_PLUS
