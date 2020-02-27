#include <algorithm>
#include <vector>

#include "def.h"
#include "kd_rect.h"
#include "kd_node.h"
#include "kd_tree.h"

// -----------------------------------------------------------------------------
//	KD_Tree: structure for approximate and exact nearest neighbor search
// -----------------------------------------------------------------------------
KD_Tree::KD_Tree(					// constructor
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   kd_leaf_size,					// leaf size of kd-tree
	const float **data)					// data objects
{
	n_pts_        = n;
	dim_          = d;
	kd_leaf_size_ = kd_leaf_size;
	data_         = data;
	object_id_    = new int[n_pts_];
	for (int i = 0; i < n_pts_; ++i) object_id_[i] = i;

	KD_Rect bnd_box(dim_);
	calc_encl_rect(bnd_box);

	bnd_box_low_  = new float[dim_];
	bnd_box_high_ = new float[dim_];
	for (int i = 0; i < dim_; ++i) {
		bnd_box_low_[i] = bnd_box.low_[i];
		bnd_box_high_[i] = bnd_box.high_[i];
	}
	root_ = rkd_tree(n_pts_, object_id_, bnd_box);
}

// -----------------------------------------------------------------------------
KD_Tree::~KD_Tree()					// destructor
{
	if (root_ != NULL) {
		delete root_; root_ = NULL;
	}
	if (object_id_ != NULL) {
		delete[] object_id_; object_id_ = NULL;
	}
	if (bnd_box_low_ != NULL) {
		delete[] bnd_box_low_; bnd_box_low_ = NULL;
	}
	if (bnd_box_high_ != NULL) {
		delete[] bnd_box_high_; bnd_box_high_ = NULL;
	}
}

// -----------------------------------------------------------------------------
void KD_Tree::calc_encl_rect(		// calc smallest enclosing rectangle
	KD_Rect &bnds)						// bounding box (return)
{
	for (int i = 0; i < dim_; ++i) {
		float low_bnd  = data_[object_id_[0]][i];
		float high_bnd = data_[object_id_[0]][i];

		// ---------------------------------------------------------------------
		//  update low and high bound
		// ---------------------------------------------------------------------
		for (int j = 1; j < n_pts_; ++j) {
			float val = data_[object_id_[j]][i];

			if (val < low_bnd) low_bnd = val;
			else if (val > high_bnd) high_bnd = val;
		}
		bnds.low_[i]  = low_bnd;
		bnds.high_[i] = high_bnd;
	}
}

// -----------------------------------------------------------------------------
KD_Node* KD_Tree::rkd_tree(			// recursive build kd-tree
	int n,								// number of data objects
	int *object_id,						// object id (return)
	KD_Rect	&bnd_box)					// bounding box for current node (return)
{
	if (n <= kd_leaf_size_) {
		return new KD_Leaf(n, dim_, object_id, data_);
	}
	else {
		int   num_low = -1;			// number of lower side
		int   cut_dim = -1;			// cutting dimension
		float cut_val = -1.0f;		// cutting value
		
		KD_Node *left  = NULL;		// left child
		KD_Node *right = NULL;		// right child

		sl_midpt_split(n, object_id, cut_dim, cut_val, num_low);

		// ---------------------------------------------------------------------
		//  save bounds for cutting dimension
		// ---------------------------------------------------------------------
		float low_val  = bnd_box.low_[cut_dim];
		float high_val = bnd_box.high_[cut_dim];

		// ---------------------------------------------------------------------
		//  recursive build left son from object_id_[0...num_low-1]
		// ---------------------------------------------------------------------
		bnd_box.high_[cut_dim] = cut_val;
		left = rkd_tree(num_low, object_id, bnd_box);
		bnd_box.high_[cut_dim] = high_val;

		// ---------------------------------------------------------------------
		//  recursive build right son from object_id_[num_low...n-1]
		// ---------------------------------------------------------------------
		bnd_box.low_[cut_dim] = cut_val;
		right = rkd_tree(n - num_low, object_id + num_low, bnd_box);
		bnd_box.low_[cut_dim] = low_val;

		return new KD_Split(dim_, cut_dim, cut_val, low_val, high_val, 
			left, right, data_);
	}
}

// -----------------------------------------------------------------------------
void KD_Tree::sl_midpt_split(		// sliding mid-object split rule
	int   n,							// number of data objects
	int   *object_id,					// object id (return)
	int   &cut_dim,						// cutting dimension (return)
	float &cut_val,						// cutting value (return)
	int   &num_low)						// num of pts on low side (return)
{
	cut_val = -1.0f;
	cut_dim = -1;

	float max_var = MINREAL;
	float cut_min = -1.0f;
	float cut_max = -1.0f;

	float mean    = -1.0f;
	float median  = -1.0f;
	float var     = -1.0f;
	float min     = -1.0f;
	float max     = -1.0f;

	for (int i = 0; i < dim_; ++i) {
		calc_stat(n, i, (const int *) object_id, mean, median, var, min, max);

		if (var > max_var) {
			max_var = var;
			cut_val = median; 
			cut_dim = i;
			cut_min = min;
			cut_max = max;
		}
	}

	// -------------------------------------------------------------------------
	//  permute points according to the <cut_dim> and <cut_val>
	// -------------------------------------------------------------------------
	int br1 = -1;
	int br2 = -1;					
	plane_split(n, cut_dim, cut_val, object_id, br1, br2);

	// -------------------------------------------------------------------------
	//	on return:	object_id[0..br1-1]   < cut_val
	//				object_id[br1..br2-1] = cut_val
	//				object_id[br2..n-1]   > cut_val
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
void KD_Tree::calc_stat(			// calc median and variance value
	int   n,							// number of data objects
	int   d,							// dimension to check
	const int *object_id,				// object id
	float &mean,						// mean value (return)
	float &median,						// median value (return)
	float &variance,					// variance (return)
	float &min,							// min value (return)
	float &max) 						// max value (return)
{	
	// -------------------------------------------------------------------------
	//  calc mean, min, and max
	// -------------------------------------------------------------------------
	std::vector<float> arr(n);
	float val  = data_[object_id[0]][d];

	arr[0] = val;
	min    = val;
	max    = val;
	mean   = val;
	for (int i = 1; i < n; ++i) {
		val = data_[object_id[i]][d];
		arr[i] = val;
		mean += val;

		if (val < min) min = val;
		else if (val > max) max = val;
	}
	mean /= n;

	// -------------------------------------------------------------------------
	//  calc median
	// -------------------------------------------------------------------------
	sort(arr.begin(), arr.end());
	if (n % 2 != 0) median = arr[n / 2];
	else median = (arr[n / 2 - 1] + arr[n / 2]) / 2;

	// -------------------------------------------------------------------------
	//  calc variance
	// -------------------------------------------------------------------------
	variance = 0.0f;
	for (int i = 0; i < n; ++i) {
		variance += SQR(data_[object_id[i]][d] - mean);
	}
	variance /= n;
}

// -----------------------------------------------------------------------------
void KD_Tree::plane_split(			// split points by a plane
	int   n,							// number of data objects
	int   cut_dim,						// cutting dimension
	float cut_val,						// cutting value
	int   *object_id,					// object id (return)
	int   &br1,							// 1st break (values < cv) (return)
	int   &br2)							// 2nd break (values = cv) (return)
{
	// -------------------------------------------------------------------------
	//  partition object_id[0..n-1] according to cut_val
	// -------------------------------------------------------------------------
	int left = 0;
	int right = n - 1;

	for (;;) {
		while (left < n && data_[object_id[left]][cut_dim] < cut_val) {
			left++;
		}
		while (right >= 0 && data_[object_id[right]][cut_dim] >= cut_val) {
			right--;
		}
		if (left > right) break;

		SWAP(object_id[left], object_id[right]);
		left++; right--;
	}
	// -------------------------------------------------------------------------
	//  now, we have "object_id[0..br1-1] < cut_val <= object_id[br1..n-1]"
	// -------------------------------------------------------------------------
	br1 = left;

	// -------------------------------------------------------------------------
	//  continue partition object_id[br1..n-1] according to cut_val
	// -------------------------------------------------------------------------
	right = n - 1;
	for (;;) {
		while (left < n && data_[object_id[left]][cut_dim] <= cut_val) {
			left++;
		}
		while (right >= 0 && data_[object_id[right]][cut_dim] > cut_val) {
			right--;
		}
		if (left > right) break;

		SWAP(object_id[left], object_id[right]);
		left++; right--;
	}
	
	// -------------------------------------------------------------------------
	//  now, we have "object_id[br1..br2-1] = cut_val < object_id[br2..n-1]"
	// -------------------------------------------------------------------------
	br2 = left;
}

// -----------------------------------------------------------------------------
void KD_Tree::search(				// k-NN search
	int   top_k,						// top-k value
	float ratio,						// approximation ratio
	const float *query,					// query object
	MinK_List *list)					// k-NN results (return)
{
	// assert(top_k <= n_pts_);
	root_->search(calc_box_dist(query), SQR(ratio), query, list);
}

// -----------------------------------------------------------------------------
float KD_Tree::calc_box_dist(		// compute distance from point to box
	const float* query)					// query point
{
	float dist = 0.0f;
	for (int i = 0; i < dim_; ++i) {
		if (query[i] < bnd_box_low_[i]) {
			dist += SQR(bnd_box_low_[i] - query[i]);
		}
		else if (query[i] > bnd_box_high_[i]) {
			dist += SQR(query[i] - bnd_box_high_[i]);
		}
	}
	return dist;
}

// -----------------------------------------------------------------------------
void KD_Tree::traversal(			// traversal kd-tree to get leaf info
	std::vector<int> &block_size,			// leaf size (return)
	int *object_id)						// object id with leaf order (return)
{
	for (int i = 0; i < n_pts_; ++i) {
		object_id[i] = object_id_[i];
	}
	root_->traversal(block_size);
}