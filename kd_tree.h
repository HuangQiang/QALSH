#ifndef __KD_TREE_H
#define __KD_TREE_H

class KD_Node;
class MinK_List;

// -----------------------------------------------------------------------------
//	KD_Tree: structure for approximate and exact nearest neighbor search
// -----------------------------------------------------------------------------
class KD_Tree {
public:
	KD_Tree(						// constructor
		int   n,						// number of data objects
		int   d,						// dimensionality
		int   kd_leaf_size,				// leaf size of kd-tree
		const float **data);			// data objects

	~KD_Tree();						// destructor

	// -------------------------------------------------------------------------
	void search(					// k-NN search
		int   top_k,					// top-k value
		float ratio,					// approximation ratio
		const float *query,				// query object
		MinK_List *list);				// k-NN results (return)

	// -------------------------------------------------------------------------
	void traversal(					// traversal kd-tree to get leaf info
		vector<int> &leaf_size,			// leaf size (return)
		int *object_id);				// object id with leaf order (return)

protected:
	int   n_pts_;					// number of data objects
	int   dim_;						// dimensionality
	int   kd_leaf_size_;			// leaf size of kd-tree

	int   *object_id_;				// data objects id
	float *bnd_box_low_;			// bounding box - low  object
	float *bnd_box_high_;			// bounding box - high object
	const float **data_;			// data objects

	KD_Node *root_;					// root of kd-tree

	// -------------------------------------------------------------------------
	void calc_encl_rect(			// calc smallest enclosing rectangle
		KD_Rect &bnds);					// bounding box (return)

	// -------------------------------------------------------------------------
	KD_Node* rkd_tree(				// recursive build kd-tree
		int n,							// number of data objects
		int *object_id,					// object id (return)
		KD_Rect	&bnd_box);				// bounding box for current node (return)

	// -------------------------------------------------------------------------
	void sl_midpt_split(			// sliding mid-object split rule
		int   n,						// number of data objects
		int   *object_id,				// object id (return)
		int   &cut_dim,					// cutting dimension (return)
		float &cut_val,					// cutting value (return)
		int   &num_low);				// num of pts on low side (return)

	// -------------------------------------------------------------------------
	void calc_stat(           		// calc mean, median and variance along dim
		int   n,						// number of data objects
		int   d,						// dimension to check
		const int *object_id,			// object id
		float &mean,                    // mean value (return)
		float &median,					// median value (return)
		float &variance,				// variance (return)
		float &min,						// min value (return)
		float &max);					// max value (return)

	// -------------------------------------------------------------------------
	void plane_split(				// split objects by a plane
		int   n,						// number of data objects
		int   cut_dim,					// cutting dimension
		float cut_val,					// cutting value
		int   *object_id,				// object id (return)
		int   &br1,						// 1st break (values < cv) (return)
		int   &br2);					// 2nd break (values = cv) (return)

	// -------------------------------------------------------------------------
	float calc_box_dist(			// compute distance from object to box
		const float *query);			// query object
};

#endif // __KD_TREE_H
