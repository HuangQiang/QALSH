#ifndef __KD_NODE_H
#define __KD_NODE_H

class MinK_List;

// -----------------------------------------------------------------------------
//	KD_Node: node of kd-tree (abstract class)
//  There are two types of KD_Node: KD_Leaf and KD_Split
// -----------------------------------------------------------------------------
class KD_Node {
public:
	virtual ~KD_Node() {}			// virtual destructor

	virtual void search(			// tree search
		float box_dist,					// box distance to query
		float ratio,					// approximation ratio
		const float *query,				// query point
		MinK_List *list) = 0;			// k-NN results (return)

	virtual void traversal(			// traversal kd-tree
		std::vector<int> &leaf_size) = 0; // leaf size (return)

	friend class KD_Tree;			// allow kd-tree to access

protected:
	int   dim_;						// dimension of data objects
	const float **data_;			// data objects
};

// -----------------------------------------------------------------------------
//	KD_Leaf: leaf node of kd-tree
// -----------------------------------------------------------------------------
class KD_Leaf : public KD_Node {
public:
	KD_Leaf(						// constructor
		int   n,						// number of data objects
		int   dim,						// dimension of data objects
		int   *object_id,				// data objects id
		const float **data);			// data objects

	virtual ~KD_Leaf();				// destructor

	virtual void search(			// tree search
		float box_dist,					// box distance to query
		float ratio,					// approximation ratio
		const float *query,				// query point
		MinK_List *list);				// k-NN results (return)

	virtual void traversal(			// traversal kd-tree
		std::vector<int> &leaf_size);	// leaf size (return)

protected:
	int n_pts_;						// number of data objects
	int *object_id_;				// data objects id
};

// -----------------------------------------------------------------------------
//	KD_Split: split node of kd-tree
// -----------------------------------------------------------------------------
class KD_Split : public KD_Node {
public:
	KD_Split(						// constructor
		int     dim,					// dimension of data objectss
		int     cut_dim,				// cutting dimension
		float   cut_val,				// cutting value
		float   low_val,				// low  bound in <dim>
		float   high_val,				// high bound in <dim>
		KD_Node *left,					// left child
		KD_Node *right,					// right child
		const   float **data);			// data objects

	virtual ~KD_Split();			// destructor

	virtual void search(			// tree search
		float box_dist,					// box distance to query
		float ratio,					// approximation ratio
		const float *query,				// query point
		MinK_List *list);				// k-NN results (return)

	virtual void traversal(			// traversal kd-tree
		std::vector<int> &leaf_size);	// leaf size (return)

protected:
	int   cut_dim_;					// cutting dimension
	float cut_val_;					// cutting value
	float cd_bnds_[2];				// cutting bounds
	KD_Node *child_[2];				// children of node
};

#endif // __KD_NODE_H
