#include <algorithm>

#include "def.h"
#include "util.h"
#include "pri_queue.h"
#include "kd_node.h"


// -----------------------------------------------------------------------------
//	KD_Leaf: leaf node of kd-tree
// -----------------------------------------------------------------------------
KD_Leaf::KD_Leaf(					// constructor
	int   n,							// number of data objects
	int   dim,							// dimension of data objects
	int   *object_id,					// data objects id
	const float **data)					// data objects
{
	n_pts_     = n;
	dim_       = dim;
	object_id_ = object_id;
	data_      = data;
}

// -----------------------------------------------------------------------------
KD_Leaf::~KD_Leaf()					// destructor
{
}

// -----------------------------------------------------------------------------
void KD_Leaf::search(				// tree search
	float box_dist,						// box distance to query
	float ratio,						// approximation ratio
	const float *query,					// query point
	MinK_List *list)					// k-NN results (return)
{
	for (int i = 0; i < n_pts_; ++i) {
		const float *point = data_[object_id_[i]];
		
		float dist = 0.0F;
		for (int j = 0; j < dim_; ++j) {
			dist += SQR(point[j] - query[j]);
		}
		list->insert(dist, object_id_[i]);
	}
}

// -----------------------------------------------------------------------------
void KD_Leaf::traversal(			// traversal kd-tree
	std::vector<int> &leaf_size)		// leaf size (return)
{
	leaf_size.push_back(n_pts_);
}

// -----------------------------------------------------------------------------
//	KD_Split: split node of kd-tree
// -----------------------------------------------------------------------------
KD_Split::KD_Split(					// constructor
	int     dim,						// dimension of data objectss
	int     cut_dim,					// cutting dimension
	float   cut_val,					// cutting value
	float   low_val,					// low  bound in <dim>
	float   high_val,					// high bound in <dim>
	KD_Node *left,						// left child
	KD_Node *right,						// right child
	const   float **data)				// data objects
{
	dim_        = dim;
	cut_dim_    = cut_dim;
	cut_val_    = cut_val;
	cd_bnds_[0] = low_val;
	cd_bnds_[1] = high_val;
	child_[0]   = left;
	child_[1]   = right;
	data_       = data;
}

// -----------------------------------------------------------------------------
KD_Split::~KD_Split()				// destructor
{
	if (child_[0] != NULL) {
		delete child_[0]; child_[0] = NULL;
	}
	if (child_[1] != NULL) {
		delete child_[1]; child_[1] = NULL;
	}
}

// -----------------------------------------------------------------------------
void KD_Split::search(				// tree search
	float box_dist,						// box distance to query
	float ratio,						// approximation ratio
	const float *query,					// query point
	MinK_List *list)					// k-NN results (return)
{
	float cut_diff = query[cut_dim_] - cut_val_;
	if (cut_diff < 0.0f) {
		// ---------------------------------------------------------------------
		//  query on left of cutting plane
		// ---------------------------------------------------------------------
		child_[0]->search(box_dist, ratio, query, list);

		float box_diff = cd_bnds_[0] - query[cut_dim_];
		if (box_diff < 0.0f) box_diff = 0.0f;

		box_dist += cut_diff * cut_diff - box_diff * box_diff;

		// ---------------------------------------------------------------------
		//  visit right child
		// ---------------------------------------------------------------------
		if (box_dist * ratio < list->max_key()) {
			child_[1]->search(box_dist, ratio, query, list);
		}
	}
	else {
		// ---------------------------------------------------------------------
		//  query on right of cutting plane
		// ---------------------------------------------------------------------
		child_[1]->search(box_dist, ratio, query, list);

		float box_diff = query[cut_dim_] - cd_bnds_[1];
		if (box_diff < 0) box_diff = 0;
		
		box_dist += cut_diff * cut_diff - box_diff * box_diff;

		// ---------------------------------------------------------------------
		//  visit left child
		// ---------------------------------------------------------------------
		if (box_dist * ratio < list->max_key()) {
			child_[0]->search(box_dist, ratio, query, list);
		}
	}
}

// -----------------------------------------------------------------------------
void KD_Split::traversal(			// traversal kd-tree
	std::vector<int> &leaf_size)		// leaf size (return)
{
	child_[0]->traversal(leaf_size);
	child_[1]->traversal(leaf_size);
}
