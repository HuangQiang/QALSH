#include "qalsh_plus.h"

// -----------------------------------------------------------------------------
QALSH_PLUS::QALSH_PLUS()			// default constructor
{
	n_pts_        = -1;
	dim_          = -1;
	B_            = -1;
	L_            = -1;
	M_            = -1;
	p_            = -1.0f;
	num_blocks_   = -1;
	sample_id_    = NULL;
	new_order_id_ = NULL;
	lsh_          = NULL;
}

// -----------------------------------------------------------------------------
QALSH_PLUS::~QALSH_PLUS()			// destructor
{
	delete lsh_; lsh_ = NULL;

	delete[] sample_id_; sample_id_ = NULL;
	delete[] new_order_id_; new_order_id_ = NULL;
	g_memory -= SIZEINT * n_pts_;
	g_memory -= SIZEINT * (num_blocks_ * L_ * M_);
	
	if (!blocks_.empty()) {
		for (int i = 0; i < num_blocks_; ++i) {
			delete blocks_[i]; blocks_[i] = NULL;
		}
		blocks_.clear(); blocks_.shrink_to_fit();
	}
}

// -----------------------------------------------------------------------------
int QALSH_PLUS::build(				// build index		
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	int   leaf,							// leaf size of kd-tree
	int   L,							// number of projection (drusilla)
	int   M,							// number of candidates (drusilla)
	float p,							// l_p distance
	float zeta,							// a parameter of p-stable distr.
	float ratio,						// approximation ratio
	const float **data,					// data objects
	const char *path)					// index path
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_      = n;
	dim_        = d;
	B_          = B;
	L_          = L;
	M_          = M;
	p_          = p;

	strcpy(path_, path); create_dir(path_);

	// -------------------------------------------------------------------------
	//  kd-tree partition
	// -------------------------------------------------------------------------
	g_memory += SIZEINT * n_pts_;
	new_order_id_ = new int[n_pts_];
	std::vector<int> block_size;
	kd_tree_partition(leaf, data, block_size, new_order_id_);

	num_blocks_ = (int) block_size.size();

	// -------------------------------------------------------------------------
	//  get sample_id and sample_data and build qalsh for each block
	// -------------------------------------------------------------------------
	int sample_size  = L_ * M_;
	int sample_n_pts = num_blocks_ * sample_size;
	sample_id_ = new int[sample_n_pts]; 
	g_memory += SIZEINT * sample_n_pts;

	float **sample_data = new float*[sample_n_pts];
	for (int i = 0; i < sample_n_pts; ++i) sample_data[i] = new float[dim_];

	int start = 0;
	int count = 0;
	for (int i = 0; i < num_blocks_; ++i) {
		// ---------------------------------------------------------------------
		//  get sample id and sample data from each blcok by drusilla select
		// ---------------------------------------------------------------------
		int block_n = block_size[i]; assert(block_n > sample_size);
		drusilla_select(block_n, new_order_id_ + start, data, sample_id_ + count);

		for (int j = 0; j < sample_size; ++j) {
			int did = sample_id_[count];
			sample_id_to_block_[did] = i;
			for (int u = 0; u < dim_; ++u) sample_data[count][u] = data[did][u];
			++count;
		}

		// ---------------------------------------------------------------------
		//  build qalsh for each blcok 
		// ---------------------------------------------------------------------
		Blocks *block = new Blocks();
		block->n_pts_ = block_n;
		block->index_ = new_order_id_ + start;

		char block_path[200];
		sprintf(block_path, "%s%d/", path_, i);
		create_dir(block_path);

		float **block_data = new float*[block_n];
		for (int j = 0; j < block_n; ++j) {
			block_data[j] = new float[dim_];
			int did = block->index_[j];
			for (int u = 0; u < dim_; ++u) block_data[j][u] = data[did][u];
		}

		block->lsh_ = new QALSH();
		block->lsh_->build(block_n, dim_, B_, p_, zeta, ratio, 
			(const float**) block_data, block_path);
		blocks_.push_back(block);

		for (int j = 0; j < block_n; ++j) {
			delete[] block_data[j]; block_data[j] = NULL; 
		}
		delete[] block_data; block_data = NULL;

		// ---------------------------------------------------------------------
		//  update parameters
		// ---------------------------------------------------------------------
		start += block_n;
	}
	assert(start == n_pts_);
	assert(count == sample_n_pts);

	// -------------------------------------------------------------------------
	//  build qalsh for sample data
	// -------------------------------------------------------------------------
	char sample_path[200];
	sprintf(sample_path, "%s%s/", path_, "sample");
	create_dir(sample_path);

	lsh_ = new QALSH();
	lsh_->build(sample_n_pts, dim_, B_, p_, zeta, ratio, 
		(const float**) sample_data, sample_path);

	// -------------------------------------------------------------------------
	//  write parameters to disk
	// -------------------------------------------------------------------------
	if (write_params()) return 1;

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	for (int i = 0; i < sample_n_pts; ++i) {
		delete[] sample_data[i]; sample_data[i] = NULL;
	}
	delete[] sample_data; sample_data = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
void QALSH_PLUS::kd_tree_partition(	// kd-tree partition
	int leaf,							// leaf size of kd-tree
	const float **data,					// data objects
	std::vector<int> &block_size,		// block size (return)
	int *new_order_id)					// new order id (return)
{
	KD_Tree* tree = new KD_Tree(n_pts_, dim_, leaf, data);
	tree->traversal(block_size, new_order_id);

	delete tree; tree = NULL;
}

// -----------------------------------------------------------------------------
void QALSH_PLUS::drusilla_select(	// drusilla select
	int   n,							// number of objects
	const int *new_order_id,			// new order data id
	const float **data,					// data objects
	int   *sample_id)					// sample data id (return)
{
	// -------------------------------------------------------------------------
	//  calc shift data
	// -------------------------------------------------------------------------
	int   max_id   = -1;
	float max_norm = MINREAL;
	float *norm = new float[n];
	float **shift_data = new float*[n];
	for (int i = 0; i < n; ++i) shift_data[i] = new float[dim_];
	
	calc_shift_data(n, new_order_id, data, max_id, max_norm, norm, shift_data);

	// -------------------------------------------------------------------------
	//  drusilla select
	// -------------------------------------------------------------------------
	float  *proj  = new float[dim_];
	Result *score = new Result[n];
	bool   *close_angle = new bool[n];

	for (int i = 0; i < L_; ++i) {
		// ---------------------------------------------------------------------
		//  determine the projection vector with max norm and normalize it
		// ---------------------------------------------------------------------
		for (int j = 0; j < dim_; ++j) {
			proj[j] = shift_data[max_id][j] / norm[max_id];
		}

		// ---------------------------------------------------------------------
		//  calculate offsets and distortions
		// ---------------------------------------------------------------------
		for (int j = 0; j < n; ++j) {
			score[j].id_ = j;
			close_angle[j] = false;

			if (norm[j] > 0.0f) {
				float offset = calc_inner_product(dim_, shift_data[j], proj);

				float distortion = 0.0f;
				for (int u = 0; u < dim_; ++u) {
					distortion += SQR(shift_data[j][u] - offset * proj[u]);
				}

				score[j].key_ = offset * offset - distortion;
				if (atan(sqrt(distortion) / fabs(offset)) < ANGLE) {
					close_angle[j] = true;
				}
			}
			else if (fabs(norm[j]) < FLOATZERO) {
				score[j].key_ = MINREAL + 1.0f;
			}
			else {
				score[j].key_ = MINREAL;
			}
		}
		// ---------------------------------------------------------------------
		//  collect the objects that are well-represented by this projection
		// ---------------------------------------------------------------------
		qsort(score, n, sizeof(Result), ResultCompDesc);

		for (int j = 0; j < M_; ++j) {
			int id = score[j].id_;
			sample_id[i * M_ + j] = new_order_id[id];
			norm[id] = -1.0f;
		}
		// ---------------------------------------------------------------------
		//  find the next largest norm and the corresponding object
		// ---------------------------------------------------------------------
		max_id   = -1;
		max_norm = MINREAL;
		for (int j = 0; j < n; ++j) {
			if (norm[j] > 0.0f && close_angle[j]) { norm[j] = 0.0f; }
			if (norm[j] > max_norm) { max_norm = norm[j]; max_id = j; }
		}
	}
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] norm;        norm        = NULL;
	delete[] close_angle; close_angle = NULL;
	delete[] proj;        proj        = NULL;
	delete[] score;       score       = NULL;

	for (int i = 0; i < n; ++i) {
		delete[] shift_data[i]; shift_data[i] = NULL;
	}
	delete[] shift_data; shift_data = NULL;
}

// -----------------------------------------------------------------------------
void QALSH_PLUS::calc_shift_data(	// calculate shift data objects
	int   n,							// number of data points
	const int *new_order_id,			// new order data id
	const float **data,					// data objects
	int   &max_id,						// data id with max l2-norm (return)
	float &max_norm,					// max l2-norm (return)
	float *norm,						// l2-norm of shift data (return)
	float **shift_data) 				// shift data (return)
{
	// -------------------------------------------------------------------------
	//  calculate the centroid of data objects
	// -------------------------------------------------------------------------
	float *centroid = new float[dim_];
	memset(centroid, 0.0f, dim_ * SIZEFLOAT);
	for (int i = 0; i < n; ++i) {
		int id = new_order_id[i];
		for (int j = 0; j < dim_; ++j) {
			centroid[j] += data[id][j];
		}
	}
	for (int i = 0; i < dim_; ++i) centroid[i] /= n;

	// -------------------------------------------------------------------------
	//  make a copy of data objects which move to the centroid of data objects
	// -------------------------------------------------------------------------
	max_id   = -1;
	max_norm = MINREAL;

	for (int i = 0; i < n; ++i) {
		int id = new_order_id[i];

		norm[i] = 0.0f;
		for (int j = 0; j < dim_; ++j) {
			float tmp = data[id][j] - centroid[j];
			shift_data[i][j] = tmp;
			norm[i] += SQR(tmp);
		}
		norm[i] = sqrt(norm[i]);

		if (norm[i] > max_norm) { max_norm = norm[i]; max_id = i; }
	}
	delete[] centroid; centroid = NULL;
}

// -----------------------------------------------------------------------------
int QALSH_PLUS::write_params()		// write parameters
{
	char fname[200];
	strcpy(fname, path_); strcat(fname, "para");

	FILE *fp = fopen(fname, "r");
	if (fp)	{ printf("Hash tables exist.\n"); exit(1); }
	fp = fopen(fname, "w");
	if (!fp) {
		printf("Could not create %s.\n", fname);
		printf("Perhaps no such folder %s?\n", path_);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  write general parameters
	// -------------------------------------------------------------------------
	fprintf(fp, "n          = %d\n", n_pts_);
	fprintf(fp, "d          = %d\n", dim_);
	fprintf(fp, "B          = %d\n", B_);
	fprintf(fp, "L          = %d\n", L_);
	fprintf(fp, "M          = %d\n", M_);
	fprintf(fp, "p          = %f\n", p_);
	fprintf(fp, "num_blocks = %d\n", num_blocks_);

	// -------------------------------------------------------------------------
	//  write first level parameters
	// -------------------------------------------------------------------------
	int sample_n_pts = num_blocks_ * L_ * M_;
	for (int i = 0; i < sample_n_pts; ++i) {
		int id = sample_id_[i];
		fprintf(fp, "%d,%d\n", id, sample_id_to_block_[id]);
	}

	// -------------------------------------------------------------------------
	//  write second level parameters
	// -------------------------------------------------------------------------
	for (int i = 0; i < num_blocks_; ++i) {		
		fprintf(fp, "%d ", blocks_[i]->n_pts_);
	}
	fprintf(fp, "\n");
	
	for (int i = 0; i < n_pts_; ++i) {
		fprintf(fp, "%d ", new_order_id_[i]);
	}
	fprintf(fp, "\n");
	fclose(fp);
	
	return 0;
}

// -----------------------------------------------------------------------------
void QALSH_PLUS::display()			// display parameters
{
	printf("Parameters of QALSH+:\n");
	printf("    n          = %d\n",   n_pts_);
	printf("    d          = %d\n",   dim_);
	printf("    B          = %d\n",   B_);
	printf("    L          = %d\n",   L_);
	printf("    M          = %d\n",   M_);
	printf("    p          = %.1f\n", p_);
	printf("    num_blocks = %d\n",   num_blocks_);
	printf("    index path = %s\n",   path_);
	printf("\n");
}

// -----------------------------------------------------------------------------
int QALSH_PLUS::load(				// load index
	const char *path)					// index path
{
	strcpy(path_, path);
	return read_params();
}

// -----------------------------------------------------------------------------
int QALSH_PLUS::read_params()		// read parameters
{
	char fname[200];
	strcpy(fname, path_); strcat(fname, "para");

	FILE* fp = fopen(fname, "r");
	if (!fp) {
		printf("Could not open %s\n", fname);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  read general parameters
	// -------------------------------------------------------------------------
	fscanf(fp, "n          = %d\n", &n_pts_);
	fscanf(fp, "d          = %d\n", &dim_);
	fscanf(fp, "B          = %d\n", &B_);
	fscanf(fp, "L          = %d\n", &L_);
	fscanf(fp, "M          = %d\n", &M_);
	fscanf(fp, "p          = %f\n", &p_);
	fscanf(fp, "num_blocks = %d\n", &num_blocks_);

	// -------------------------------------------------------------------------
	//  load first level parameters (lsh, sample_id, and sample_id_to_block)
	// -------------------------------------------------------------------------
	char sample_path[200];
	sprintf(sample_path, "%s%s/", path_, "sample");

	lsh_ = new QALSH();
	lsh_->load(sample_path);

	int sample_n_pts = num_blocks_ * L_ * M_;
	int id = -1, block_id = -1;

	g_memory += SIZEINT * sample_n_pts;
	sample_id_ = new int[sample_n_pts];
	for (int i = 0; i < sample_n_pts; ++i) {
		fscanf(fp, "%d,%d\n", &id, &block_id);

		sample_id_[i] = id;
		sample_id_to_block_[id] = block_id;
	}

	// -------------------------------------------------------------------------
	//  load second level parameters (blocks and new_order_id)
	// -------------------------------------------------------------------------
	char block_path[200];
	for (int i = 0; i < num_blocks_; ++i) {
		sprintf(block_path, "%s%d/", path_, i);

		Blocks *block = new Blocks();
		block->lsh_ = new QALSH();
		block->lsh_->load(block_path);

		fscanf(fp, "%d ", &block->n_pts_);
		blocks_.push_back(block);
	}
	fscanf(fp, "\n");

	g_memory += SIZEINT * n_pts_;
	new_order_id_ = new int[n_pts_];
	for (int i = 0; i < n_pts_; ++i) {
		fscanf(fp, "%d ", &new_order_id_[i]);
	}
	fscanf(fp, "\n");
	fclose(fp);

	int start = 0;
	for (int i = 0; i < num_blocks_; ++i) {
		blocks_[i]->index_ = new_order_id_ + start;
		start += blocks_[i]->n_pts_;
	}
	assert(start == n_pts_);
	
	return 0;
}

// -----------------------------------------------------------------------------
uint64_t QALSH_PLUS::knn(			// k-NN search
	int   top_k,						// top-k value
	int   nb,							// number of blocks for search
	const float *query,					// input query
	const char *data_folder,			// data folder
	MinK_List *list)					// top-k results (return)
{
	assert(nb > 0 && nb <= num_blocks_);
	uint64_t io_cost = 0;

	// -------------------------------------------------------------------------
	//  use sample data to determine the order of blocks for k-NN search
	// -------------------------------------------------------------------------
	MinK_List *sample_list = new MinK_List(MAXK);	
	io_cost += lsh_->knn(MAXK, MAXREAL, query, sample_id_, data_folder, sample_list);

	std::vector<int> block_order;
	get_block_order(nb, sample_list, block_order);

	// -------------------------------------------------------------------------
	//  use <nb> blocks for c-k-ANN search
	// -------------------------------------------------------------------------
	float radius = MAXREAL;
	for (size_t i = 0; i < block_order.size(); ++i) {
		int id = block_order[i];
		Blocks *block = blocks_[id];

		io_cost += block->lsh_->knn(top_k, radius, query, block->index_, 
			data_folder, list);
		radius = list->max_key();
	}
	delete sample_list; sample_list = NULL;

	return io_cost;
}

// -----------------------------------------------------------------------------
int QALSH_PLUS::get_block_order(	// get block order
	int nb,								// number of blocks for search
	MinK_List *list,					// sample results
	std::vector<int> &block_order)		// block order (return)
{
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
		// if (fabs(pair[i].key_) < FLOATZERO) break;
		block_order.push_back(pair[i].id_);
	}
	delete[] pair; pair = NULL;

	return 0;
}
