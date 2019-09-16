#include "headers.h"


// -----------------------------------------------------------------------------
QALSH_PLUS::QALSH_PLUS(				// constructor
	int   n,							// cardinality
	int   d,							// dimensionality
	int   leaf,							// leaf size of kd-tree
	int   L,							// number of projection
	int   M,							// number of candidates for each proj
	float p,							// l_p distance
	float zeta,							// a parameter of p-stable distr.
	float ratio,						// approximation ratio
	const float **data)					// data objects
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_      = n;
	dim_        = d;
	leaf_       = leaf;
	L_          = L;
	M_          = M;
	p_          = p;
	zeta_       = zeta;
	appr_ratio_ = ratio;

	// -------------------------------------------------------------------------
	//  bulkloading
	// -------------------------------------------------------------------------
	bulkload(data);
}

// -----------------------------------------------------------------------------
int QALSH_PLUS::bulkload(			// bulkloading
	const float **data)					// data objects
{
	// -------------------------------------------------------------------------
	//  kd-tree partition
	// -------------------------------------------------------------------------
	int *new_order_id = new int[n_pts_];
	vector<int> block_size;
	kd_tree_partition(data, block_size, new_order_id);

	// -------------------------------------------------------------------------
	//  init <new_order_data_>　for bulkloading data objects 
	// -------------------------------------------------------------------------
	new_order_data_ = new float*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) {
		new_order_data_[i] = new float[dim_];

		int id = new_order_id[i];
		for (int j = 0; j < dim_; ++j) {
			new_order_data_[i][j] = data[id][j];
		}
	}
	
	// -------------------------------------------------------------------------
	//  bulkloading data objects for each block 
	// -------------------------------------------------------------------------
	int sample_size = L_ * M_;
	num_blocks_   = (int) block_size.size();
	sample_n_pts_ = num_blocks_ * sample_size;

	int *sample_id      = new int[sample_n_pts_]; 
	sample_id_to_block_ = new int[sample_n_pts_];
	sample_data_        = new float*[sample_n_pts_];
	for (int i = 0; i < sample_n_pts_; ++i) {
		sample_data_[i] = new float[dim_];
	}

	int start = 0;
	int count = 0;
	for (int i = 0; i < num_blocks_; ++i) {
		// ---------------------------------------------------------------------
		//  calculate shift data
		// ---------------------------------------------------------------------
		int n = block_size[i]; assert(n > sample_size);

		vector<vector<float> > shift_data(n, vector<float>(dim_, 0.0f));
		calc_shift_data(n, dim_, (const float **)new_order_data_+start, shift_data);
		
		// ---------------------------------------------------------------------
		//  select sample data from each blcok 
		// ---------------------------------------------------------------------
		drusilla_select(n, dim_, shift_data, (const int *) new_order_id + start, 
			sample_id + count);

		for (int j = 0; j < sample_size; ++j) {
			sample_id_to_block_[count] = i;

			int id = sample_id[count];
			for (int z = 0; z < dim_; ++z) {
				sample_data_[count][z] = data[id][z];
			}
			++count;
		}

		// ---------------------------------------------------------------------
		//  build index for each blcok 
		// ---------------------------------------------------------------------
		Blocks *block = new Blocks();
		block->n_pts_ = n;
		block->index_.resize(n);
		for (int j = 0; j < n; ++j) {
			block->index_[j] = new_order_id[start + j];
		}

		block->lsh_ = new QALSH(n, dim_, p_, zeta_, appr_ratio_, 
			(const float **) new_order_data_ + start);
		blocks_.push_back(block);

		// ---------------------------------------------------------------------
		//  update parameters
		// ---------------------------------------------------------------------
		start += n;
	}
	assert(count == sample_n_pts_);

	// -------------------------------------------------------------------------
	//  build index for sample data
	// -------------------------------------------------------------------------
	lsh_ = new QALSH(sample_n_pts_, dim_, p_, zeta_, appr_ratio_,
		(const float **) sample_data_);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] new_order_id; new_order_id = NULL;
	delete[] sample_id; sample_id = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int QALSH_PLUS::kd_tree_partition(	// kd-tree partition
	const float **data,					// data objects
	vector<int> &block_size,			// block size
	int *new_order_id)					// new order id
{
	KD_Tree *tree = new KD_Tree(n_pts_, dim_, leaf_, data);
	tree->traversal(block_size, new_order_id);

	delete tree; tree = NULL;
	return 0;
}

// -----------------------------------------------------------------------------
int QALSH_PLUS::calc_shift_data(	// calculate shift data objects
	int   n,							// number of data points
	int   d,							// dimensionality
	const float **data,					// data objects
	vector<vector<float> > &shift_data) // shift data objects (return)
{
	// -------------------------------------------------------------------------
	//  calculate the centroid of data objects
	// -------------------------------------------------------------------------
	vector<float> centroid(d, 0.0f);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < d; ++j) {
			centroid[j] += data[i][j];
		}
	}
	for (int i = 0; i < d; ++i) {
		centroid[i] /= (float) n;
	}

	// -------------------------------------------------------------------------
	//  make a copy of data objects which move to the centroid of data objects
	// -------------------------------------------------------------------------
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < d; ++j) {
			shift_data[i][j] = data[i][j] - centroid[j];
		}
	}
	
	return 0;
}

// -----------------------------------------------------------------------------
int QALSH_PLUS::drusilla_select(	// drusilla select
	int   n,							// number of objects
	int   d,							// data dimension
	const vector<vector<float> > &shift_data, // shift data objects
	const int *new_order_id,			// new order data id
	int   *sample_id)					// sample data id (return)
{
	// -------------------------------------------------------------------------
	//  calc the norm of data objects and find the data object with max norm
	// -------------------------------------------------------------------------
	int   max_id   = -1;
	float max_norm = -1.0f;
	vector<float> norm(n, 0.0f);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < d; ++j) {
			float x = shift_data[i][j];
			norm[i] += x * x;
		}
		norm[i] = sqrt(norm[i]);

		if (norm[i] > max_norm) {
			max_norm = norm[i];
			max_id   = i;
		}
	}

	vector<bool>  close_angle(n);
	vector<float> projection(d);
	Result *score_pair = new Result[n];

	for (int i = 0; i < L_; ++i) {
		// ---------------------------------------------------------------------
		//  determine the projection vector with max norm and normalize it
		// ---------------------------------------------------------------------
		for (int j = 0; j < d; ++j) {
			projection[j] = shift_data[max_id][j] / norm[max_id];
		}

		// ---------------------------------------------------------------------
		//  calculate offsets and distortions
		// ---------------------------------------------------------------------
		for (int j = 0; j < n; ++j) {
			score_pair[j].id_ = j;
			close_angle[j] = false;

			if (norm[j] > 0.0f) {
				float offset = 0.0f;
				for (int k = 0; k < d; ++k) {
					offset += shift_data[j][k] * projection[k];
				}

				float distortion = 0.0f;
				for (int k = 0; k < d; ++k) {
					float x = shift_data[j][k] - offset * projection[k];
					distortion += x * x;
				}

				score_pair[j].key_ = offset * offset - distortion;
				if (atan(distortion / fabs(offset)) < ANGLE) {
					close_angle[j] = true;
				}
			}
			else if (fabs(norm[j]) < FLOATZERO) {
				score_pair[j].key_ = MINREAL + 1.0f;
			}
			else {
				score_pair[j].key_ = MINREAL;
			}
		}

		// ---------------------------------------------------------------------
		//  collect the objects that are well-represented by this projection
		// ---------------------------------------------------------------------
		qsort(score_pair, n, sizeof(Result), ResultCompDesc);

		for (int j = 0; j < M_; ++j) {
			int id = score_pair[j].id_;

			sample_id[i * M_ + j] = new_order_id[id];
			norm[id] = -1.0f;
		}

		// ---------------------------------------------------------------------
		//  find the next largest norm and the corresponding object
		// ---------------------------------------------------------------------
		max_id   = -1;
		max_norm = -1.0f;
		for (int j = 0; j < n; ++j) {
			if (norm[j] > 0.0f && close_angle[j]) {
				norm[j] = 0.0f;
			}

			if (norm[j] > max_norm) {
				max_norm = norm[j];
				max_id   = j;
			}
		}
	}

	// -------------------------------------------------------------------------
	//  initialize parameters for k-FP search
	// -------------------------------------------------------------------------
	delete[] score_pair; score_pair = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
QALSH_PLUS::~QALSH_PLUS()			// destructor
{
	for (int i = 0; i < n_pts_; ++i) {
		delete[] new_order_data_[i]; new_order_data_[i] = NULL;
	}
	delete[] new_order_data_; new_order_data_ = NULL;

	for (int i = 0; i < num_blocks_; ++i) {
		delete blocks_[i]; blocks_[i] = NULL;
	}
	blocks_.clear(); blocks_.shrink_to_fit();

	delete lsh_; lsh_ = NULL;
	delete[] sample_id_to_block_; sample_id_to_block_ = NULL;

	for (int i = 0; i < sample_n_pts_; ++i) {
		delete[] sample_data_[i]; sample_data_[i] = NULL;
	}
	delete[] sample_data_; sample_data_ = NULL;
}

// -----------------------------------------------------------------------------
void QALSH_PLUS::display()			// display parameters
{
	printf("Parameters of QALSH+:\n");
	printf("    n          = %d\n",   n_pts_);
	printf("    d          = %d\n",   dim_);
	printf("    leaf       = %d\n",   leaf_);
	printf("    L          = %d\n",   L_);
	printf("    M          = %d\n",   M_);
	printf("    p          = %.1f\n", p_);
	printf("    zeta       = %.1f\n", zeta_);
	printf("    c          = %.1f\n", appr_ratio_);
	printf("    num_blocks = %d\n",   num_blocks_);
	printf("    sample_n   = %d\n\n", sample_n_pts_);
}

// -----------------------------------------------------------------------------
int QALSH_PLUS::knn(				// c-k-ANN search
	int   top_k,						// top-k value
	int   nb,							// number of blocks for search
	const float *query,					// input query objects
	MinK_List *list)					// k-NN results
{
	// -------------------------------------------------------------------------
	//  use sample data to determine the order of blocks for k-NN search
	// -------------------------------------------------------------------------
	MinK_List *sample_list = new MinK_List(MAXK);
	lsh_->knn(MAXK, query, sample_list);

	vector<int> block_order;
	get_block_order(nb, sample_list, block_order);
	
	// -------------------------------------------------------------------------
	//  use <nb> blocks for k-NN search
	// -------------------------------------------------------------------------
	int radius = MAXREAL;
	int size = (int) block_order.size();
	for (int i = 0; i < size; ++i) {
		int id = block_order[i];
		blocks_[id]->lsh_->knn(top_k, radius, query, blocks_[id]->index_, list);

		radius = list->max_key();
	}
	delete sample_list; sample_list = NULL;
	
	return 0;
}

// -----------------------------------------------------------------------------
int QALSH_PLUS::get_block_order(	// get block order
	int nb,								// number of blocks for search
	MinK_List *list,					// sample results
	vector<int> &block_order)			// block order (return)
{
	assert(nb <= num_blocks_);

	// -------------------------------------------------------------------------
	//  init the counter of each block
	// -------------------------------------------------------------------------
	Result *pair = new Result[num_blocks_];
	for (int i = 0; i < num_blocks_; ++i) {
		pair[i].id_  = i;
		pair[i].key_ = 0.0f;
	}

	// -------------------------------------------------------------------------
	//  select the first <nb> blocks with largest counters
	// -------------------------------------------------------------------------
	int size = (int) list->size();
	for (int i = 0; i < size; ++i) {
		int block_id = sample_id_to_block_[list->ith_id(i)];
		pair[block_id].key_ += 1.0f;
	}

	qsort(pair, num_blocks_, sizeof(Result), ResultCompDesc);
	for (int i = 0; i < nb; ++i) {
		if (fabs(pair[i].key_) < FLOATZERO) break;
		block_order.push_back(pair[i].id_);
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] pair; pair = NULL;

	return 0;
}
