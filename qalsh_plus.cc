#include "headers.h"


// -----------------------------------------------------------------------------
QALSH_Plus::QALSH_Plus()			// constructor
{
	n_pts_              = -1;
	dim_                = -1;
	B_                  = -1;
	kd_leaf_size_       = -1;
	L_                  = -1;
	M_                  = -1;
	p_                  = -1.0f;
	zeta_               = -1.0f;
	appr_ratio_         = -1.0f;
	
	num_blocks_         = -1;
	sample_n_pts_       = -1;
	lsh_                = NULL;
}

// -----------------------------------------------------------------------------
QALSH_Plus::~QALSH_Plus()			// destructor
{
	if (!blocks_.empty()) {
		for (int i = 0; i < num_blocks_; ++i) {
			delete blocks_[i]; blocks_[i] = NULL;
		}
		blocks_.clear(); blocks_.shrink_to_fit();
	}
	
	if (lsh_ != NULL) {
		delete lsh_; lsh_ = NULL;
	}
}

// -----------------------------------------------------------------------------
int QALSH_Plus::build(				// build index		
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	int   kd_leaf_size,					// leaf size of kd-tree
	int   L,							// number of projection (drusilla)
	int   M,							// number of candidates (drusilla)
	float p,							// l_p distance
	float zeta,							// a parameter of p-stable distr.
	float ratio,						// approximation ratio
	const float **data,					// data objects
	const char *index_path)				// index path
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_        = n;
	dim_          = d;
	B_            = B;
	kd_leaf_size_ = kd_leaf_size;
	L_            = L;
	M_            = M;
	p_            = p;
	zeta_         = zeta;
	appr_ratio_   = ratio;

	strcpy(index_path_, index_path);
	create_dir(index_path_);

	// -------------------------------------------------------------------------
	//  build hash tables (bulkloading)
	// -------------------------------------------------------------------------
	bulkload(data);
	display();
	
	return 0;
}

// -----------------------------------------------------------------------------
int QALSH_Plus::bulkload(			// bulkloading
	const float **data)					// data objects
{
	// -------------------------------------------------------------------------
	//  kd-tree partition
	// -------------------------------------------------------------------------
	int *new_order_id = new int[n_pts_];
	vector<int> block_size;
	kd_tree_partition(data, block_size, new_order_id);

	// -------------------------------------------------------------------------
	//  init new_order_data　for bulkloading data objects 
	// -------------------------------------------------------------------------
	float **new_order_data = new float*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) {
		new_order_data[i] = new float[dim_];

		int id = new_order_id[i];
		for (int j = 0; j < dim_; ++j) {
			new_order_data[i][j] = data[id][j];
		}
	}
	
	// -------------------------------------------------------------------------
	//  bulkloading data objects for each block 
	// -------------------------------------------------------------------------
	int sample_size = L_ * M_;
	num_blocks_ = (int) block_size.size();
	sample_n_pts_ = num_blocks_ * sample_size;
	sample_id_.resize(sample_n_pts_); 

	int *sample_id = new int[sample_n_pts_]; 
	float **sample_data = new float*[sample_n_pts_];
	for (int i = 0; i < sample_n_pts_; ++i) {
		sample_data[i] = new float[dim_];
	}

	int start = 0;
	int count = 0;
	for (int i = 0; i < num_blocks_; ++i) {
		// ---------------------------------------------------------------------
		//  calculate shift data
		// ---------------------------------------------------------------------
		int n = block_size[i]; assert(n > sample_size);

		vector<vector<float> > shift_data(n, vector<float>(dim_, 0.0f));
		calc_shift_data(n, dim_, (const float **) new_order_data + start, shift_data);
		
		// ---------------------------------------------------------------------
		//  select sample data from each blcok 
		// ---------------------------------------------------------------------
		drusilla_select(n, dim_, shift_data, (const int *) new_order_id + start, 
			sample_id + count);

		for (int j = 0; j < sample_size; ++j) {
			int id = sample_id[count];

			sample_id_[count] = id;
			sample_id_to_block_[id] = i;
			
			for (int z = 0; z < dim_; ++z) {
				sample_data[count][z] = data[id][z];
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

		char path[200];
		sprintf(path, "%s%d/", index_path_, i);
		create_dir(path);

		block->lsh_ = new QALSH();
		block->lsh_->build(n, dim_, B_, p_, zeta_, appr_ratio_, 
			(const float **) new_order_data + start, path);
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
	char path[200];
	sprintf(path, "%s%s/", index_path_, "sample");
	create_dir(path);

	lsh_ = new QALSH();
	lsh_->build(sample_n_pts_, dim_, B_, p_, zeta_, appr_ratio_, 
		(const float **) sample_data, path);

	// -------------------------------------------------------------------------
	//  write parameters to disk
	// -------------------------------------------------------------------------
	write_params();

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] new_order_id; new_order_id = NULL;
	delete[] sample_id; sample_id = NULL;

	for (int i = 0; i < n_pts_; ++i) {
		delete[] new_order_data[i]; new_order_data[i] = NULL; 
	}
	delete[] new_order_data; new_order_data = NULL;

	for (int i = 0; i < sample_n_pts_; ++i) {
		delete[] sample_data[i]; sample_data[i] = NULL;
	}
	delete[] sample_data; sample_data = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int QALSH_Plus::kd_tree_partition(	// kd-tree partition
	const float **data,					// data objects
	vector<int> &block_size,			// block size
	int *new_order_id)					// new order id
{
	KD_Tree* tree = new KD_Tree(n_pts_, dim_, kd_leaf_size_, data);
	tree->traversal(block_size, new_order_id);

	delete tree; tree = NULL;
	return 0;
}

// -----------------------------------------------------------------------------
int QALSH_Plus::calc_shift_data(	// calculate shift data objects
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
int QALSH_Plus::drusilla_select(	// drusilla select
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

	vector<bool> close_angle(n);
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
int QALSH_Plus::write_params()		// write parameters
{
	char fname[200];
	strcpy(fname, index_path_);
	strcat(fname, "para");

	FILE *fp = fopen(fname, "r");
	if (fp)	{
		printf("Hash tables exist.\n");
		exit(1);
	}

	fp = fopen(fname, "w");
	if (!fp) {
		printf("Could not create %s.\n", fname);
		printf("Perhaps no such folder %s?\n", index_path_);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  write general parameters
	// -------------------------------------------------------------------------
	fprintf(fp, "n = %d\n", n_pts_);
	fprintf(fp, "d = %d\n", dim_);
	fprintf(fp, "B = %d\n", B_);
	fprintf(fp, "kd_leaf_size = %d\n", kd_leaf_size_);
	fprintf(fp, "L = %d\n", L_);
	fprintf(fp, "M = %d\n", M_);
	fprintf(fp, "p = %f\n", p_);
	fprintf(fp, "zeta = %f\n", zeta_);
	fprintf(fp, "c = %f\n", appr_ratio_);
	fprintf(fp, "num_blocks = %d\n", num_blocks_);
	fprintf(fp, "sample_n = %d\n", sample_n_pts_);

	// -------------------------------------------------------------------------
	//  write first level parameters
	// -------------------------------------------------------------------------
	for (int i = 0; i < sample_n_pts_; ++i) {
		int id = sample_id_[i];
		fprintf(fp, "%d,%d\n", id, sample_id_to_block_[id]);
	}

	// -------------------------------------------------------------------------
	//  write second level parameters
	// -------------------------------------------------------------------------
	for (int i = 0; i < num_blocks_; ++i) {
		int n = blocks_[i]->n_pts_;
		
		fprintf(fp, "%d", blocks_[i]->n_pts_);
		for (int j = 0; j < n; ++j) {
			fprintf(fp, " %d", blocks_[i]->index_[j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	
	return 0;
}

// -----------------------------------------------------------------------------
int QALSH_Plus::display()			// display parameters
{
	printf("Parameters of QALSH_Plus:\n");
	printf("    n            = %d\n", n_pts_);
	printf("    d            = %d\n", dim_);
	printf("    B            = %d\n", B_);
	printf("    kd_leaf_size = %d\n", kd_leaf_size_);
	printf("    L            = %d\n", L_);
	printf("    M            = %d\n", M_);
	printf("    p            = %f\n", p_);
	printf("    zeta         = %f\n", zeta_);
	printf("    c            = %f\n", appr_ratio_);
	printf("    num_blocks   = %d\n", num_blocks_);
	printf("    sample_n     = %d\n", sample_n_pts_);
	printf("    index path   = %s\n", index_path_);
	printf("\n");

	return 0;
}

// -----------------------------------------------------------------------------
int QALSH_Plus::load(				// load index
	const char *index_path)				// index path
{
	strcpy(index_path_, index_path);
	return read_params();
}

// -----------------------------------------------------------------------------
int QALSH_Plus::read_params()		// read parameters
{
	char fname[200];
	strcpy(fname, index_path_);
	strcat(fname, "para");

	FILE* fp = fopen(fname, "r");
	if (!fp) {
		printf("Could not open %s.\n", fname);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  read general parameters
	// -------------------------------------------------------------------------
	fscanf(fp, "n = %d\n", &n_pts_);
	fscanf(fp, "d = %d\n", &dim_);
	fscanf(fp, "B = %d\n", &B_);
	fscanf(fp, "kd_leaf_size = %d\n", &kd_leaf_size_);
	fscanf(fp, "L = %d\n", &L_);
	fscanf(fp, "M = %d\n", &M_);
	fscanf(fp, "p = %f\n", &p_);
	fscanf(fp, "zeta = %f\n", &zeta_);
	fscanf(fp, "c = %f\n", &appr_ratio_);
	fscanf(fp, "num_blocks = %d\n", &num_blocks_);
	fscanf(fp, "sample_n = %d\n", &sample_n_pts_);

	// -------------------------------------------------------------------------
	//  laod first level parameters
	// -------------------------------------------------------------------------
	char path[200];
	sprintf(path, "%s%s/", index_path_, "sample");

	lsh_ = new QALSH();
	lsh_->load(path);

	int id = -1, block_id = -1;
	sample_id_.resize(sample_n_pts_);
	for (int i = 0; i < sample_n_pts_; ++i) {
		fscanf(fp, "%d,%d\n", &id, &block_id);

		sample_id_[i] = id;
		sample_id_to_block_[id] = block_id;
	}

	// -------------------------------------------------------------------------
	//  load second level parameters
	// -------------------------------------------------------------------------
	for (int i = 0; i < num_blocks_; ++i) {
		sprintf(path, "%s%d/", index_path_, i);

		Blocks *block = new Blocks();
		block->lsh_ = new QALSH();
		block->lsh_->load(path);

		fscanf(fp, "%d", &block->n_pts_);

		int n = block->n_pts_;
		block->index_.resize(n);
		for (int j = 0; j < n; ++j) {
			fscanf(fp, " %d", &block->index_[j]);
		}
		fscanf(fp, "\n");

		blocks_.push_back(block);
	}
	fclose(fp);
	
	return 0;
}

// -----------------------------------------------------------------------------
int QALSH_Plus::knn(				// k-NN search
	int   top_k,						// top-k value
	int   nb,							// number of blocks for search
	const float *query,					// input query
	const char *data_folder,			// data folder
	MinK_List *list)					// top-k results (return)
{
	long long io_cost = 0;

	// -------------------------------------------------------------------------
	//  use sample data to determine the order of blocks for k-NN search
	// -------------------------------------------------------------------------
	MinK_List *sample_list = new MinK_List(MAXK);	
	io_cost += lsh_->knn(MAXK, MAXREAL, query, sample_id_, data_folder, sample_list);

	vector<int> block_order;
	get_block_order(nb, sample_list, block_order);
	delete sample_list; sample_list = NULL;

	// -------------------------------------------------------------------------
	//  use <nb> blocks for c-k-ANN search
	// -------------------------------------------------------------------------
	float radius = MAXREAL;
	int size = (int) block_order.size();

	for (int i = 0; i < size; ++i) {
		int id = block_order[i];
		io_cost += blocks_[id]->lsh_->knn(top_k, radius, query, 
			blocks_[id]->index_, data_folder, list);

		radius = list->max_key();
	}

	return io_cost;
}

// -----------------------------------------------------------------------------
int QALSH_Plus::get_block_order(	// get block order
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
		assert(block_id >= 0 && block_id < num_blocks_);

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
