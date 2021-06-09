#pragma once

#include <iostream>
#include <algorithm>
#include <vector>
#include <cstring>

#include "def.h"
#include "kd_node.h"

namespace nns {

// -----------------------------------------------------------------------------
//  KD_Rect: orthogonal rectangle for bounding rectangle of kd-tree
// -----------------------------------------------------------------------------
template<class DType>
class KD_Rect {
public:
	DType *low_;					// rectangle lower bound
	DType *high_;					// rectangle upper bound

	// -------------------------------------------------------------------------
	KD_Rect();						// default constructor

	KD_Rect(						// copy constructor
		int   dim,						// dimension
		const KD_Rect<DType> &rect);	// rectangle to copy

	KD_Rect(						// constructor with specific value
		int   dim,						// dimension
		DType l = 0,					// low  boundary, default is zero
		DType h = 0);					// high boundary, default is zero

	KD_Rect(						// constructor with specific points
		int   dim,						// dimension
		const DType *low,				// low point
		const DType *high);				// high point

	// -------------------------------------------------------------------------
	~KD_Rect();						// destructor

	// -------------------------------------------------------------------------
	bool inside(					// whether a point inside rectangle
		int   dim,						// dimension
		const DType *point);			// one point
};

// -----------------------------------------------------------------------------
template<class DType>
KD_Rect<DType>::KD_Rect()			// default constructor
{
	low_ = NULL; high_ = NULL;
}

// -----------------------------------------------------------------------------
template<class DType>
KD_Rect<DType>::KD_Rect(			// copy constructor
	int   dim,							// dimension
	const KD_Rect<DType> &rect)			// copy item
{
	low_  = new DType[dim];
	high_ = new DType[dim];
	for (int i = 0; i < dim; ++i) {
		low_[i]  = rect.low_[i];
		high_[i] = rect.high_[i];
	}
}

// -----------------------------------------------------------------------------
template<class DType>
KD_Rect<DType>::KD_Rect(			// constructor with specific value
	int   dim,							// dimension
	DType l,							// lower bound
	DType h)							// higher bound
{
	low_  = new DType[dim];
	high_ = new DType[dim];
	for (int i = 0; i < dim; ++i) {
		low_[i]  = l;
		high_[i] = h;
	}
}

// -----------------------------------------------------------------------------
template<class DType>
KD_Rect<DType>::KD_Rect(			// constructor with specific points
	int   dim,							// dimension
	const DType *low,					// lower corner point
	const DType *high)					// higher corner point
{
	low_  = new DType[dim];
	high_ = new DType[dim];
	for (int i = 0; i < dim; ++i) {
		low_[i]  = low[i];
		high_[i] = high[i];
	}
}

// -----------------------------------------------------------------------------
template<class DType>
KD_Rect<DType>::~KD_Rect()			// destructor
{
	if (low_  != NULL) { delete[] low_;  low_  = NULL; }
	if (high_ != NULL) { delete[] high_; high_ = NULL; }
}

// -----------------------------------------------------------------------------
template<class DType>
bool KD_Rect<DType>::inside(		// whether a point inside the rectangle
	int   dim,							// dimension
	const DType *point)					// input point
{
	for (int i = 0; i < dim; ++i) {
		if (point[i] < low_[i] || point[i] > high_[i]) return false;
	}
	return true;
}


// -----------------------------------------------------------------------------
//  KD_Tree: structure for exact and approximate k-NN search
// -----------------------------------------------------------------------------
template<class DType>
class KD_Tree {
public:
	KD_Tree(						// constructor
		int   n,						// number of data points
		int   d,						// dimensionality
		int   leaf,						// leaf size of kd-tree
		const DType *data);				// data points

	// -------------------------------------------------------------------------
	~KD_Tree();						// destructor

	// -------------------------------------------------------------------------
	void search(					// k-NN search
		int   top_k,					// top-k value
		float c,						// approximation ratio
		const DType *query,				// query
		MinK_List *list);				// k-NN results (return)

	// -------------------------------------------------------------------------
	void traversal(					// traversal kd-tree to get leaf info
		std::vector<int> &leaf_size,	// leaf size (return)
		int *index);					// data index with leaf order (return)

protected:
	int   n_pts_;					// number of data points
	int   dim_;						// dimensionality
	int   leaf_;					// leaf size of kd-tree
	const DType *data_;				// data points
	
	int   *index_;					// data index
	DType *bnd_box_low_;			// bounding box - low  point
	DType *bnd_box_high_;			// bounding box - high point
	KD_Node<DType> *root_;			// root of kd-tree

	// -------------------------------------------------------------------------
	void calc_encl_rect(			// calc smallest enclosing rectangle
		KD_Rect<DType> &bnds);			// bounding box (return)

	// -------------------------------------------------------------------------
	KD_Node<DType>* rkd_tree(		// recursive build kd-tree
		int n,							// number of data points
		int *index,						// data index (return)
		KD_Rect<DType> &bnd_box);		// bounding box for current node (return)

	// -------------------------------------------------------------------------
	void sl_midpt_split(			// sliding mid-point split rule
		int   n,						// number of data points
		int   *index,					// data index (return)
		int   &cut_dim,					// cutting dimension (return)
		DType &cut_val,					// cutting value (return)
		int   &num_low);				// num of pts on low side (return)

	// -------------------------------------------------------------------------
	void calc_stat(					// calc mean, median and variance along dim
		int   n,						// number of data points
		int   d,						// dimension to check
		const int *index,				// data index
		float &mean,					// mean value (return)
		DType &median,					// median value (return)
		float &variance,				// variance (return)
		DType &min,						// min value (return)
		DType &max);					// max value (return)

	// -------------------------------------------------------------------------
	void plane_split(				// split points by a plane
		int   n,						// number of data points
		int   cut_dim,					// cutting dimension
		DType cut_val,					// cutting value
		int   *index,					// data index (return)
		int   &br1,						// 1st break (values < cv) (return)
		int   &br2);					// 2nd break (values = cv) (return)

	// -------------------------------------------------------------------------
	float calc_box_dist(			// compute distance from point to box
		const DType *query);			// query point
};

// -----------------------------------------------------------------------------
template<class DType>
KD_Tree<DType>::KD_Tree(			// constructor
	int   n,							// number of data points
	int   d,							// dimensionality
	int   leaf,							// leaf size of kd-tree
	const DType *data)					// data points
	: n_pts_(n), dim_(d), leaf_(leaf), data_(data)
{
	index_ = new int[n];
	for (int i = 0; i < n; ++i) index_[i] = i;

	KD_Rect<DType> bnd_box(d);
	calc_encl_rect(bnd_box);

	bnd_box_low_  = new DType[d];
	bnd_box_high_ = new DType[d];
	for (int i = 0; i < d; ++i) {
		bnd_box_low_[i]  = bnd_box.low_[i];
		bnd_box_high_[i] = bnd_box.high_[i];
	}
	root_ = rkd_tree(n, index_, bnd_box);
}

// -----------------------------------------------------------------------------
template<class DType>
KD_Tree<DType>::~KD_Tree()			// destructor
{
	delete root_;

	delete[] index_;
	delete[] bnd_box_low_;
	delete[] bnd_box_high_;
}

// -----------------------------------------------------------------------------
template<class DType>
void KD_Tree<DType>::calc_encl_rect(// calc smallest enclosing rectangle
	KD_Rect<DType> &bnds)				// bounding box (return)
{
	for (int i = 0; i < dim_; ++i) {
		DType low_bnd  = data_[(uint64_t)index_[0]*dim_+i];
		DType high_bnd = data_[(uint64_t)index_[0]*dim_+i];

		// update the low bound and high bound
		for (int j = 1; j < n_pts_; ++j) {
			DType val = data_[(uint64_t)index_[j]*dim_+i];

			if (val < low_bnd) low_bnd = val;
			else if (val > high_bnd) high_bnd = val;
		}
		bnds.low_[i]  = low_bnd;
		bnds.high_[i] = high_bnd;
	}
}

// -----------------------------------------------------------------------------
template<class DType>
KD_Node<DType>* KD_Tree<DType>::rkd_tree(// recursive build kd-tree
	int n,								// number of data points
	int *index,							// data index (return)
	KD_Rect<DType> &bnd_box)			// bounding box for current node (return)
{
	if (n <= leaf_) {
		return new KD_Leaf<DType>(n, dim_, index, data_);
	}
	else {
		int   num_low = -1;			// number of lower side
		int   cut_dim = -1;			// cutting dimension
		DType cut_val = -1;			// cutting value
		
		KD_Node<DType> *left  = NULL; // left child
		KD_Node<DType> *right = NULL; // right child

		sl_midpt_split(n, index, cut_dim, cut_val, num_low);

		// save bounds for cutting dimension
		DType low_val  = bnd_box.low_[cut_dim];
		DType high_val = bnd_box.high_[cut_dim];

		// recursive build left son from index_[0...num_low-1]
		bnd_box.high_[cut_dim] = cut_val;
		left = rkd_tree(num_low, index, bnd_box);
		bnd_box.high_[cut_dim] = high_val;

		// recursive build right son from index_[num_low...n-1]
		bnd_box.low_[cut_dim] = cut_val;
		right = rkd_tree(n-num_low, &index[num_low], bnd_box);
		bnd_box.low_[cut_dim] = low_val;

		return new KD_Split<DType>(dim_, cut_dim, cut_val, low_val, high_val, 
			left, right, data_);
	}
}

// -----------------------------------------------------------------------------
template<class DType>
void KD_Tree<DType>::sl_midpt_split(// sliding mid-point split rule
	int   n,							// number of data points
	int   *index,						// data index (return)
	int   &cut_dim,						// cutting dimension (return)
	DType &cut_val,						// cutting value (return)
	int   &num_low)						// num of pts on low side (return)
{
	cut_val = -1; cut_dim = -1;
	float max_var = MINREAL, mean = -1.0f, var = -1.0f;
	DType median, min, max, cut_min, cut_max;

	for (int i = 0; i < dim_; ++i) {
		calc_stat(n, i, (const int*) index, mean, median, var, min, max);
		if (var > max_var) {
			max_var = var;
			cut_val = median; cut_dim = i; cut_min = min; cut_max = max;
		}
	}
	// permute points according to the <cut_dim> and <cut_val>
	int br1 = -1, br2 = -1;
	plane_split(n, cut_dim, cut_val, index, br1, br2);

	// -------------------------------------------------------------------------
	//	on return:	index[0..br1-1]   < cut_val
	//				index[br1..br2-1] = cut_val
	//				index[br2..n-1]   > cut_val
	//
	//	we can set num_low to any value in the range [br1..br2]
	// -------------------------------------------------------------------------
	if (cut_val < cut_min) num_low = 1;
	else if (cut_val > cut_max) num_low = n - 1;
	else if (br1 > n / 2) num_low = br1;
	else if (br2 < n / 2) num_low = br2;
	else num_low = n / 2;
}

// -----------------------------------------------------------------------------
template<class DType>
void KD_Tree<DType>::calc_stat(		// calc median and variance value
	int   n,							// number of data points
	int   d,							// dimension to check
	const int *index,					// data index
	float &mean,						// mean value (return)
	DType &median,						// median value (return)
	float &var,							// variance (return)
	DType &min,							// min value (return)
	DType &max) 						// max value (return)
{
	// init arr & calc mean, min, and max
	std::vector<DType> arr(n);
	DType val = data_[(uint64_t)index[0]*dim_+d];
	arr[0] = val; mean = (float) val; min = val; max = val;

	for (int i = 1; i < n; ++i) {
		val = data_[(uint64_t)index[i]*dim_+d];

		arr[i] = val; mean += (float) val;
		if (val < min) min = val;
		else if (val > max) max = val;
	}
	mean /= n;

	// calc median
	sort(arr.begin(), arr.end());
	if (n%2 != 0) median = arr[n/2];
	else median = (arr[n/2-1] + arr[n/2]) / 2;

	// calc variance
	var = 0.0f;
	for (int i = 0; i < n; ++i) {
		float tmp = (float) (data_[(uint64_t)index[i]*dim_+d] - mean);
		var += tmp*tmp;
	}
	var /= n;

	arr.clear(); arr.shrink_to_fit();
}

// -----------------------------------------------------------------------------
template<class DType>
void KD_Tree<DType>::plane_split(	// split points by a plane
	int   n,							// number of data points
	int   d,							// cutting dimension
	DType val,							// cutting value
	int   *index,						// data index (return)
	int   &br1,							// 1st break (values < cv) (return)
	int   &br2)							// 2nd break (values = cv) (return)
{
	// partition index[0..n-1] according to cut_val
	int left  = 0;
	int right = n - 1;

	for (;;) {
		while (left < n && data_[(uint64_t)index[left]*dim_+d] < val) {
			left++;
		}
		while (right >= 0 && data_[(uint64_t)index[right]*dim_+d] >= val) {
			right--;
		}
		if (left > right) break;

		SWAP(index[left], index[right]);
		left++; right--;
	}
	// now, we have "index[0..br1-1] < cut_val <= index[br1..n-1]"
	br1 = left;

	// continue partition index[br1..n-1] according to cut_val
	right = n - 1;
	for (;;) {
		while (left < n && data_[(uint64_t)index[left]*dim_+d] <= val) {
			left++;
		}
		while (right >= 0 && data_[(uint64_t)index[right]*dim_+d] > val) {
			right--;
		}
		if (left > right) break;

		SWAP(index[left], index[right]);
		left++; right--;
	}
	// now, we have "index[br1..br2-1] = cut_val < index[br2..n-1]"
	br2 = left;
}

// -----------------------------------------------------------------------------
template<class DType>
void KD_Tree<DType>::search(		// k-NN search
	int   top_k,						// top-k value
	float c,							// approximation ratio
	const DType *query,					// query point
	MinK_List *list)					// k-NN results (return)
{
	// NOTE: the k-NN results with l2-dist sqr NOT l2-dist
	assert(top_k <= n_pts_);
	root_->search(calc_box_dist(query), c*c, query, list);
}

// -----------------------------------------------------------------------------
template<class DType>
float KD_Tree<DType>::calc_box_dist(// compute distance from point to box
	const DType* query)					// query point
{
	float dist = 0.0f;
	for (int i = 0; i < dim_; ++i) {
		if (query[i] < bnd_box_low_[i]) {
			dist += (float) SQR(bnd_box_low_[i] - query[i]);
		}
		else if (query[i] > bnd_box_high_[i]) {
			dist += (float) SQR(query[i] - bnd_box_high_[i]);
		}
	}
	return dist;
}

// -----------------------------------------------------------------------------
template<class DType>
void KD_Tree<DType>::traversal(		// traversal kd-tree to get leaf info
	std::vector<int> &leaf_size,		// leaf size (return)
	int *index)							// data index with leaf order (return)
{
	for (int i = 0; i < n_pts_; ++i) {
		index[i] = index_[i];
	}
	root_->traversal(leaf_size);
}

} // end namespace nns
