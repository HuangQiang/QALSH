#ifndef __KD_TREE_H
#define __KD_TREE_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

#include "def.h"
#include "kd_node.h"

class KD_Node;
class MinK_List;

// -----------------------------------------------------------------------------
//	KD_Rect: orthogonal rectangle for bounding rectangle of kd-tree
// -----------------------------------------------------------------------------
class KD_Rect {
public:
	float *low_;					// rectangle lower bound
	float *high_;					// rectangle upper bound

	// -------------------------------------------------------------------------
	KD_Rect(						// constructor
		int   dim,						// dimension
		float l = 0,					// low  boundary, default is zero
		float h = 0);					// high boundary, default is zero

	// -------------------------------------------------------------------------
	KD_Rect(						// copy constructor
		int   dim,						// dimension
		const KD_Rect &rect);			// rectangle to copy

	// -------------------------------------------------------------------------
	KD_Rect(						// construct from points
		int   dim,						// dimension
		const float *low,				// low point
		const float *high);				// high point

	// -------------------------------------------------------------------------
	~KD_Rect();						// destructor

	// -------------------------------------------------------------------------
	bool inside(					// whether a point inside rectangle
		int   dim,						// dimension
		const float *point);			// one point
};

// -----------------------------------------------------------------------------
//	KD_Tree: structure for approximate and exact nearest neighbor search
// -----------------------------------------------------------------------------
class KD_Tree {
public:
	KD_Tree(						// constructor
		int   n,						// number of data objects
		int   d,						// dimensionality
		int   leaf,						// leaf size of kd-tree
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
		std::vector<int> &leaf_size,	// leaf size (return)
		int *object_id);				// object id with leaf order (return)

protected:
	int   n_pts_;					// number of data objects
	int   dim_;						// dimensionality
	int   leaf_;					// leaf size of kd-tree
	const float **data_;			// data objects

	int   *object_id_;				// data objects id
	float *bnd_box_low_;			// bounding box - low  object
	float *bnd_box_high_;			// bounding box - high object
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
