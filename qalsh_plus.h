#ifndef __QALSH_PLUS_H
#define __QALSH_PLUS_H

class MinK_List;
class QALSH;

// -----------------------------------------------------------------------------
//  Blocks: an block which stores hash tables for some of data objects
// -----------------------------------------------------------------------------
struct Blocks {
	int n_pts_;
	std::vector<int> index_;
	QALSH *lsh_;

	Blocks() { n_pts_ = -1; lsh_ = NULL; }
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
		const char *index_path);		// index path

	// -------------------------------------------------------------------------
	int load(						// load index
		const char *index_path);		// index path

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	long long knn(					// k-NN search
		int top_k,						// top-k value
		int nb,							// number of blocks for search
		const float *query,				// input query
		const char *data_folder,		// data folder
		MinK_List *list);				// top-k results (return)

protected:
	int   n_pts_;					// number of data objects
	int   dim_;						// dimensionality
	int   B_;						// page size
	int   leaf_;					// leaf size of kd-tree
	int   L_;						// number of projection (drusilla)
	int   M_;						// number of candidates (drusilla)
	float p_;						// l_p distance
	float zeta_;					// a parameter of p-stable distr.
	float appr_ratio_;				// approximation ratio
	char  index_path_[200];			// index path

	int   num_blocks_;				// number of blocks 
	std::vector<Blocks*> blocks_;	// blocks

	int   sample_n_pts_;			// number of sample data objects
	QALSH *lsh_;					// index of sample data objects
	std::vector<int> sample_id_;	// sample data id
	std::unordered_map<int, int> sample_id_to_block_; // sample data id to block

	// -------------------------------------------------------------------------
	int bulkload(					// bulkloading for each block
		const float **data);			// original data objects

	// -------------------------------------------------------------------------
	int kd_tree_partition(			// kd-tree partition 
		const float **data,				// data objects
		std::vector<int> &block_size,	// block size (return)
		int *new_order_id);				// new order id (return)

	// -------------------------------------------------------------------------
	int calc_shift_data(			// calculate shift data objects
		int   n,						// number of data objects
		int   d,						// data dimension
		const float **data,				// data objects
		std::vector<std::vector<float> > &shift_data); // shift data objects (return)

	// -------------------------------------------------------------------------
	int drusilla_select(			// drusilla select
		int   n,						// number of data objects
		int   d,						// data dimension
		const std::vector<std::vector<float> > &shift_data, // shift data objects
		const int *new_order_id,		// new order data id
		int   *sample_id);				// sample data id (return)

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
