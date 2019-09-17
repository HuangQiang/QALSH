#include "headers.h"

// -----------------------------------------------------------------------------
QALSH::QALSH()						// constructor
{
	n_pts_       = -1;
	dim_         = -1;
	B_           = -1;
	p_           = -1.0f;
	zeta_        = -1.0f;
	appr_ratio_  = -1.0f;
	
	w_           = -1.0f;
	p1_          = -1.0f;
	p2_          = -1.0f;
	alpha_       = -1.0f;
	beta_        = -1.0f;
	delta_       = -1.0f;
	m_           = -1;
	l_           = -1;
	a_array_     = NULL;
	trees_       = NULL;

	dist_io_     = -1;
	page_io_     = -1;
	freq_        = NULL;
	checked_     = NULL;
	bucket_flag_ = NULL;
	range_flag_  = NULL;
	data_        = NULL;
	q_val_       = NULL;
	lptr_        = NULL;
	rptr_        = NULL;
}

// -----------------------------------------------------------------------------
QALSH::~QALSH()						// destructor
{
	delete[] a_array_;     a_array_     = NULL;
	delete[] freq_;        freq_        = NULL;
	delete[] checked_;     checked_     = NULL;
	delete[] bucket_flag_; bucket_flag_ = NULL;
	delete[] range_flag_;  range_flag_  = NULL;
	delete[] data_;        data_        = NULL;
	delete[] q_val_;       q_val_       = NULL;
	
	for (int i = 0; i < m_; ++i) {
		delete trees_[i]; trees_[i] = NULL;
		if (lptr_[i] != NULL) {
			delete lptr_[i]; lptr_[i] = NULL;
		}
		if (rptr_[i] != NULL) {
			delete rptr_[i]; rptr_[i] = NULL;
		}
	}
	delete[] trees_; trees_ = NULL;
	delete[] lptr_;  lptr_ = NULL;
	delete[] rptr_;  rptr_ = NULL;
}

// -----------------------------------------------------------------------------
int QALSH::build(					// build index
	int   n,							// number of data points
	int   d,							// dimension of space
	int   B,							// page size
	float p,							// the p value of L_p norm
	float zeta,							// symmetric factor of p-stable distr.
	float ratio,						// approximation ratio
	const float **data,					// data objects
	const char *index_path)				// index path
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_      = n;
	dim_        = d;
	B_          = B;
	p_          = p;
	zeta_       = zeta;
	appr_ratio_ = ratio;

	strcpy(index_path_, index_path);
	create_dir(index_path_);

	// -------------------------------------------------------------------------
	//  calc parameters and generate hash functions
	// -------------------------------------------------------------------------
	calc_params();
	gen_hash_func();

	// -------------------------------------------------------------------------
	//  bulkloading
	// -------------------------------------------------------------------------
	bulkload(data);

	return 0;
}

// -----------------------------------------------------------------------------
void QALSH::calc_params()			// calc params of qalsh
{
	delta_ = 1.0f / E;
	beta_  = (float) CANDIDATES / (float) n_pts_;

	// -------------------------------------------------------------------------
	//  init <w_> <p1_> and <p2_> (auto tuning-w)
	//  
	//  w0 ----- best w for L_{0.5} norm to minimize m (auto tuning-w)
	//  w1 ----- best w for L_{1.0} norm to minimize m (auto tuning-w)
	//  w2 ----- best w for L_{2.0} norm to minimize m (auto tuning-w)
	//  other w: use linear combination for interpolation
	// -------------------------------------------------------------------------
	float w0 = (appr_ratio_ - 1.0f) / log(sqrt(appr_ratio_));
	float w1 = 2.0f * sqrt(appr_ratio_);
	float w2 = sqrt((8.0f * appr_ratio_ * appr_ratio_ * log(appr_ratio_))
		/ (appr_ratio_ * appr_ratio_ - 1.0f));

	if (fabs(p_ - 0.5f) < FLOATZERO) {
		w_ = w0;
		p1_ = calc_l0_prob(w_ / 2.0f);
		p2_ = calc_l0_prob(w_ / (2.0f * appr_ratio_));
	}
	else if (fabs(p_ - 1.0f) < FLOATZERO) {
		w_ = w1;
		p1_ = calc_l1_prob(w_ / 2.0f);
		p2_ = calc_l1_prob(w_ / (2.0f * appr_ratio_));
	}
	else if (fabs(p_ - 2.0f) < FLOATZERO) {
		w_ = w2;
		p1_ = calc_l2_prob(w_ / 2.0f);
		p2_ = calc_l2_prob(w_ / (2.0f * appr_ratio_));
	}
	else {
		if (fabs(p_ - 0.8f) < FLOATZERO) {
			w_ = 2.503f;
		}
		else if (fabs(p_ - 1.2f) < FLOATZERO) {
			w_ = 3.151f;
		}
		else if (fabs(p_ - 1.5f) < FLOATZERO) {
			w_ = 3.465f;
		}
		else {
			w_ = (w2 - w1) * p_ + (2.0f * w1 - w2);
		}
		new_stable_prob(p_, zeta_, appr_ratio_, 1.0f, w_, 1000000, p1_, p2_);
	}

	float para1 = sqrt(log(2.0f / beta_));
	float para2 = sqrt(log(1.0f / delta_));
	float para3 = 2.0f * (p1_ - p2_) * (p1_ - p2_);
	float eta   = para1 / para2;

	alpha_ = (eta * p1_ + p2_) / (1.0f + eta);
	m_     = (int) ceil((para1 + para2) * (para1 + para2) / para3);
	l_     = (int) ceil(alpha_ * m_);

	freq_        = new int[n_pts_];
	checked_     = new bool[n_pts_];
	bucket_flag_ = new bool[m_];
	range_flag_  = new bool[m_];
	data_        = new float[dim_];
	q_val_       = new float[m_];
	
	lptr_ = new PageBuffer*[m_];
	rptr_ = new PageBuffer*[m_];
	for (int i = 0; i < m_; ++i) {
		lptr_[i] = new PageBuffer();
		rptr_[i] = new PageBuffer();
	}
}

// -----------------------------------------------------------------------------
inline float QALSH::calc_l0_prob(	// calc prob <p1_> and <p2_> of L1/2 dist
	float x)							// x = w / (2.0 * r)
{
	return new_levy_prob(x);
}

// -----------------------------------------------------------------------------
inline float QALSH::calc_l1_prob(	// calc prob <p1_> and <p2_> of L1 dist
	float x)							// x = w / (2.0 * r)
{
	return new_cauchy_prob(x);
}

// -----------------------------------------------------------------------------
inline float QALSH::calc_l2_prob(	// calc prob <p1_> and <p2_> of L2 dist
	float x)							// x = w / (2.0 * r)
{
	return new_gaussian_prob(x);
}

// -----------------------------------------------------------------------------
void QALSH::display()				// display parameters
{
	printf("Parameters of QALSH:\n");
	printf("    n     = %d\n", n_pts_);
	printf("    d     = %d\n", dim_);
	printf("    B     = %d\n", B_);
	printf("    p     = %.1f\n", p_);
	printf("    zeta  = %.1f\n", zeta_);
	printf("    ratio = %.1f\n", appr_ratio_);
	printf("    w     = %f\n", w_);
	printf("    p1    = %f\n", p1_);
	printf("    p2    = %f\n", p2_);
	printf("    alpha = %f\n", alpha_);
	printf("    beta  = %f\n", beta_);
	printf("    delta = %f\n", delta_);
	printf("    m     = %d\n", m_);
	printf("    l     = %d\n", l_);
	printf("    path  = %s\n", index_path_);
	printf("\n");
}

// -----------------------------------------------------------------------------
void QALSH::gen_hash_func()			// generate hash function <a_array>
{
	int size = m_ * dim_;
	a_array_ = new float[size];
	for (int i = 0; i < size; ++i) {
		if (fabs(p_ - 0.5f) < FLOATZERO) {
			a_array_[i] = levy(1.0f, 0.0f);
		}
		else if (fabs(p_ - 1.0f) < FLOATZERO) {
			a_array_[i] = cauchy(1.0f, 0.0f);
		}
		else if (fabs(p_ - 2.0f) < FLOATZERO) {
			a_array_[i] = gaussian(0.0f, 1.0f);
		}
		else {
			a_array_[i] = p_stable(p_, zeta_, 1.0f, 0.0f);
		}
	}
}

// -----------------------------------------------------------------------------
int QALSH::bulkload(				// build m b-trees by bulkloading
	const float **data)					// data set
{
	// -------------------------------------------------------------------------
	//  write parameters to disk
	// -------------------------------------------------------------------------
	if (write_params()) return 1;

	// -------------------------------------------------------------------------
	//  write hash tables (indexed by B+ Tree) to disk
	// -------------------------------------------------------------------------
	char fname[200];
	trees_ = new BTree*[m_];

	Result *hashtable = new Result[n_pts_];
	for (int i = 0; i < m_; ++i) {
		for (int j = 0; j < n_pts_; ++j) {
			hashtable[j].id_  = j;
			hashtable[j].key_ = calc_hash_value(i, data[j]);
		}
		qsort(hashtable, n_pts_, sizeof(Result), ResultComp);
		
		get_tree_filename(i, fname);
		trees_[i] = new BTree();
		trees_[i]->init(B_, fname);
		if (trees_[i]->bulkload(n_pts_, (const Result *) hashtable)) {
			return 1;
		}
	}
	delete[] hashtable; hashtable = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int QALSH::write_params()			// write parameters to disk
{
	char fname[200];
	strcpy(fname, index_path_);
	strcat(fname, "para");

	FILE *fp = fopen(fname, "r");
	if (fp)	{
		printf("Hash Tables Already Exist\n\n");
		exit(1);
	}

	fp = fopen(fname, "w");
	if (!fp) {
		printf("Could not create %s\n", fname);
		printf("Perhaps no such folder %s?\n", index_path_);
		return 1;
	}

	fprintf(fp, "n     = %d\n", n_pts_);
	fprintf(fp, "d     = %d\n", dim_);
	fprintf(fp, "B     = %d\n", B_);
	fprintf(fp, "p     = %f\n", p_);
	fprintf(fp, "zeta  = %f\n", zeta_);
	fprintf(fp, "ratio = %f\n", appr_ratio_);
	fprintf(fp, "w     = %f\n", w_);
	fprintf(fp, "p1    = %f\n", p1_);
	fprintf(fp, "p2    = %f\n", p2_);
	fprintf(fp, "alpha = %f\n", alpha_);
	fprintf(fp, "beta  = %f\n", beta_);
	fprintf(fp, "delta = %f\n", delta_);
	fprintf(fp, "m     = %d\n", m_);
	fprintf(fp, "l     = %d\n", l_);

	int count = 0;
	for (int i = 0; i < m_; ++i) {
		for (int j = 0; j < dim_; ++j) {
			fprintf(fp, "%f ", a_array_[count++]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);	

	return 0;
}

// -----------------------------------------------------------------------------
inline float QALSH::calc_hash_value( // calc hash value
	int   tid,							// hash table id
	const float *point)					// input data object
{
	float ret  = 0.0f;
	int   base = tid * dim_;
	for (int i = 0; i < dim_; ++i) {
		ret += a_array_[base + i] * point[i];
	}
	return ret;
}

// -----------------------------------------------------------------------------
inline void QALSH::get_tree_filename( // get file name of b-tree
	int  tree_id,						// tree id, from 0 to m-1
	char *fname)						// file name (return)
{
	sprintf(fname, "%s%d.qalsh", index_path_, tree_id);
}

// -----------------------------------------------------------------------------
int QALSH::load(					// load index
	const char *index_path)				// index path
{
	strcpy(index_path_, index_path);
	if (read_params()) return 1;

	char fname[200];
	trees_ = new BTree*[m_];
	for (int i = 0; i < m_; ++i) {
		get_tree_filename(i, fname);

		trees_[i] = new BTree();
		trees_[i]->init_restore(fname);
	}
	return 0;
}

// -----------------------------------------------------------------------------
int QALSH::read_params()			// read parameters from disk
{
	char fname[200];
	strcpy(fname, index_path_);
	strcat(fname, "para");

	FILE *fp = fopen(fname, "r");
	if (!fp) {
		printf("Could not open %s\n", fname);
		return 1;
	}

	fscanf(fp, "n     = %d\n", &n_pts_);
	fscanf(fp, "d     = %d\n", &dim_);
	fscanf(fp, "B     = %d\n", &B_);
	fscanf(fp, "p     = %f\n", &p_);
	fscanf(fp, "zeta  = %f\n", &zeta_);
	fscanf(fp, "ratio = %f\n", &appr_ratio_);
	fscanf(fp, "w     = %f\n", &w_);
	fscanf(fp, "p1    = %f\n", &p1_);
	fscanf(fp, "p2    = %f\n", &p2_);
	fscanf(fp, "alpha = %f\n", &alpha_);
	fscanf(fp, "beta  = %f\n", &beta_);
	fscanf(fp, "delta = %f\n", &delta_);
	fscanf(fp, "m     = %d\n", &m_);
	fscanf(fp, "l     = %d\n", &l_);
	
	a_array_ = new float[m_ * dim_];
	int count = 0;
	for (int i = 0; i < m_; ++i) {
		for (int j = 0; j < dim_; ++j) {
			fscanf(fp, "%f ", &a_array_[count++]);
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	freq_        = new int[n_pts_];
	checked_     = new bool[n_pts_];
	bucket_flag_ = new bool[m_];
	range_flag_  = new bool[m_];
	data_        = new float[dim_];
	q_val_       = new float[m_];
	
	lptr_ = new PageBuffer*[m_];
	rptr_ = new PageBuffer*[m_];
	for (int i = 0; i < m_; ++i) {
		lptr_[i] = new PageBuffer();
		rptr_[i] = new PageBuffer();
	}

	return 0;
}

// -----------------------------------------------------------------------------
long long QALSH::knn(				// k-NN search
	int top_k,							// top-k value
	const float *query,					// query object
	const char *data_folder,			// data folder
	MinK_List *list)					// k-NN results (return)
{
	int candidates = CANDIDATES + top_k - 1; // candidates size
	float kdist = MAXREAL;

	// -------------------------------------------------------------------------
	//  initialize parameters for k-NN search
	// -------------------------------------------------------------------------
	init_search_params(query);

	// -------------------------------------------------------------------------
	//  k-NN search
	// -------------------------------------------------------------------------
	float radius = find_radius();
	float bucket = w_ * radius / 2.0f;

	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialize the stop condition for current round
		// ---------------------------------------------------------------------
		int num_bucket = 0;
		memset(bucket_flag_, true, m_ * SIZEBOOL);

		// ---------------------------------------------------------------------
		//  step 2: find frequent objects
		// ---------------------------------------------------------------------
		while (num_bucket < m_) {
			for (int i = 0; i < m_; ++i) {
				if (!bucket_flag_[i]) continue;

				PageBuffer *lptr = lptr_[i];
				PageBuffer *rptr = rptr_[i];
				// -------------------------------------------------------------
				//  step 2.1: compute <ldist> and <rdist>
				// -------------------------------------------------------------
				float dist  = -1.0f;
				float ldist = MAXREAL;
				float rdist = MAXREAL;
				if (lptr->size_ != -1) ldist = calc_dist(q_val_[i], lptr);
				if (rptr->size_ != -1) rdist = calc_dist(q_val_[i], rptr);

				// -------------------------------------------------------------
				//  step 2.2: determine the closer direction (left or right)
				//  and do collision counting to find frequent objects.
				//
				//  for the frequent object, we calc the L_{p} distance with
				//  query, and update the k-nn result.
				// -------------------------------------------------------------
				if (ldist < bucket && ldist <= rdist) {
					int count = lptr->size_;
					int end   = lptr->leaf_pos_;
					int start = end - count;

					for (int j = end; j > start; --j) {
						int id = lptr->leaf_node_->get_entry_id(j);
						if (++freq_[id] > l_ && !checked_[id]) {
							checked_[id] = true;
							read_data_new_format(id, dim_, B_, data_folder, data_);

							dist  = calc_lp_dist(dim_, p_, kdist, data_, query);
							kdist = list->insert(dist, id);
							if (++dist_io_ >= candidates) break;
						}
					}
					update_left_buffer(rptr, lptr);
				}
				else if (rdist < bucket && ldist > rdist) {
					int count = rptr->size_;
					int start = rptr->leaf_pos_;
					int end   = start + count;

					for (int j = start; j < end; ++j) {
						int id = rptr->leaf_node_->get_entry_id(j);
						if (++freq_[id] > l_ && !checked_[id]) {
							checked_[id] = true;
							read_data_new_format(id, dim_, B_, data_folder, data_);

							dist  = calc_lp_dist(dim_, p_, kdist, data_, query);
							kdist = list->insert(dist, id);
							if (++dist_io_ >= candidates) break;
						}
					}
					update_right_buffer(lptr, rptr);
				}
				else {
					bucket_flag_[i] = false;
					num_bucket++;
				}
				if (num_bucket >= m_ || dist_io_ >= candidates) break;
			}
			if (num_bucket >= m_ || dist_io_ >= candidates) break;
		}
		// ---------------------------------------------------------------------
		//  step 3: stop conditions 1 & 2
		// ---------------------------------------------------------------------
		if (kdist < appr_ratio_ * radius && dist_io_ >= top_k) break;
		if (dist_io_ >= candidates) break;

		// ---------------------------------------------------------------------
		//  step 4: auto-update <radius>
		// ---------------------------------------------------------------------
		radius = update_radius(radius);
		bucket = radius * w_ / 2.0f;
	}
	delete_tree_ptr();	

	return (long long) (page_io_ + dist_io_);
}

// -----------------------------------------------------------------------------
long long QALSH::knn(				// k-NN search
	int top_k,							// top-k value
	float R,							// limited search range
	const float *query,					// query object
	const vector<int> &object_id,		// object id mapping
	const char *data_folder,			// data folder
	MinK_List *list)					// k-NN results (return)
{
	int candidates = CANDIDATES + top_k - 1; // candidates size
	float kdist = R;

	// -------------------------------------------------------------------------
	//  initialize parameters for k-NN search
	// -------------------------------------------------------------------------
	init_search_params(query);

	// -------------------------------------------------------------------------
	//  k-NN search
	// -------------------------------------------------------------------------
	int   num_range = 0;				// used for search range bound
	float radius    = find_radius();
	float bucket    = w_ * radius / 2.0f;

	float range = -1.0f;
	if (R > MAXREAL - 1.0f) range = MAXREAL;
	else range = R * w_ / 2.0f;

	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialize the stop condition for current round
		// ---------------------------------------------------------------------
		int num_bucket = 0;
		memset(bucket_flag_, true, m_ * SIZEBOOL);

		// ---------------------------------------------------------------------
		//  step 2: find frequent objects
		// ---------------------------------------------------------------------
		while (num_bucket < m_ && num_range < m_) {
			for (int i = 0; i < m_; ++i) {
				if (!bucket_flag_[i]) continue;

				PageBuffer *lptr = lptr_[i];
				PageBuffer *rptr = rptr_[i];
				// -------------------------------------------------------------
				//  step 2.1: compute <ldist> and <rdist>
				// -------------------------------------------------------------
				float dist  = -1.0f;
				float ldist = MAXREAL;
				float rdist = MAXREAL;

				if (lptr->size_ != -1) ldist = calc_dist(q_val_[i], lptr);
				if (rptr->size_ != -1) rdist = calc_dist(q_val_[i], rptr);

				// -------------------------------------------------------------
				//  step 2.2: determine the closer direction (left or right)
				//  and do collision counting to find frequent objects.
				//
				//  for the frequent object, we calc the L2 distance with
				//  query, and update the k-nn result.
				// -------------------------------------------------------------
				if (ldist < bucket && ldist < range && ldist <= rdist) {
					int count = lptr->size_;
					int end   = lptr->leaf_pos_;
					int start = end - count;

					for (int j = end; j > start; --j) {
						int id = lptr->leaf_node_->get_entry_id(j);
						if (++freq_[id] > l_ && !checked_[id]) {
							checked_[id] = true;
							int oid = object_id[id];
							read_data_new_format(oid, dim_, B_, data_folder, data_);

							dist  = calc_lp_dist(dim_, p_, kdist, data_, query);
							kdist = list->insert(dist, oid);
							if (++dist_io_ >= candidates) break;
						}
					}
					update_left_buffer(rptr, lptr);
				}
				else if (rdist < bucket && rdist < range && ldist > rdist) {
					int count = rptr->size_;
					int start = rptr->leaf_pos_;
					int end   = start + count;

					for (int j = start; j < end; ++j) {
						int id = rptr->leaf_node_->get_entry_id(j);
						if (++freq_[id] > l_ && !checked_[id]) {
							checked_[id] = true;
							int oid = object_id[id];
							read_data_new_format(oid, dim_, B_, data_folder, data_);

							dist  = calc_lp_dist(dim_, p_, kdist, data_, query);
							kdist = list->insert(dist, oid);
							if (++dist_io_ >= candidates) break;
						}
					}
					update_right_buffer(lptr, rptr);
				}
				else {
					bucket_flag_[i] = false;
					num_bucket++;
					if (ldist >= range && rdist >= range && range_flag_[i]) {
						range_flag_[i] = false;
						num_range++;
					}
				}
				if (num_bucket >= m_ || num_range >= m_) break;
				if (dist_io_ >= candidates) break;
			}
			if (num_bucket >= m_ || num_range >= m_) break;
			if (dist_io_ >= candidates) break;
		}
		// ---------------------------------------------------------------------
		//  step 3: stop conditions 1 & 2
		// ---------------------------------------------------------------------
		if (dist_io_ >= candidates || num_range >= m_) break;

		// ---------------------------------------------------------------------
		//  step 4: auto-update <radius>
		// ---------------------------------------------------------------------
		radius = update_radius(radius);
		bucket = radius * w_ / 2.0f;
	}
	delete_tree_ptr();	

	return (long long) (page_io_ + dist_io_);
}

// -----------------------------------------------------------------------------
void QALSH::init_search_params(		// init parameters for k-NN search
	const float *query)					// query object
{
	page_io_ = 0;
	dist_io_ = 0;

	memset(freq_, 0, n_pts_ * SIZEFLOAT);
	memset(checked_, false, n_pts_ * SIZEBOOL);
	memset(bucket_flag_, true, m_ * SIZEBOOL);
	memset(range_flag_, true, m_ * SIZEBOOL);

	for (int i = 0; i < m_; ++i) {
		lptr_[i]->leaf_node_ = NULL;
		lptr_[i]->index_pos_ = -1;
		lptr_[i]->leaf_pos_  = -1;
		lptr_[i]->size_      = -1;

		rptr_[i]->leaf_node_ = NULL;
		rptr_[i]->index_pos_ = -1;
		rptr_[i]->leaf_pos_  = -1;
		rptr_[i]->size_      = -1;
	}

	BIndexNode *index_node = NULL;
	int  block       = -1;			// variables for index node
	int  follow      = -1;
	bool lescape     = false;
	int  pos         = -1;			// variables for leaf node
	int  increment   = -1;
	int  num_entries = -1;

	for (int i = 0; i < m_; ++i) {
		float q_val = calc_hash_value(i, query);
		BTree *tree = trees_[i];
		PageBuffer *lptr = lptr_[i];
		PageBuffer *rptr = rptr_[i];

		q_val_[i] = q_val;
		block = tree->root_;
		if (block > 1) {
			// ---------------------------------------------------------------------
			//  at least two levels in the B+ Tree: index node and lead node
			// ---------------------------------------------------------------------
			index_node = new BIndexNode();
			index_node->init_restore(tree, block);
			page_io_++;

			// ---------------------------------------------------------------------
			//  find the leaf node whose value is closest and larger than the key
			//  of query q
			// ---------------------------------------------------------------------
			lescape = false;		// locate the position of branch
			while (index_node->get_level() > 1) {
				follow = index_node->find_position_by_key(q_val);

				if (follow == -1) {	// if in the most left branch
					if (lescape) {	// scan the most left branch
						follow = 0;
					}
					else {
						if (block != tree->root_) {
							printf("No branch found\n");
							exit(1);
						}
						else {
							follow = 0;
							lescape = true;
						}
					}
				}
				block = index_node->get_son(follow);
				delete index_node; index_node = NULL;

				index_node = new BIndexNode();
				index_node->init_restore(tree, block);
				page_io_++;			// access a new node (a new page)
			}

			// ---------------------------------------------------------------------
			//  after finding the leaf node whose value is closest to the key of
			//  query, initialize <lptrs[i]> and <rptrs[i]>.
			//
			//  <lescape> = true is that the query has no <lptrs>, the query is 
			//  the smallest value.
			// ---------------------------------------------------------------------
			follow = index_node->find_position_by_key(q_val);
			if (follow < 0) {
				lescape = true;
				follow = 0;
			}

			if (lescape) {			
				// -----------------------------------------------------------------
				//  only init right buffer
				// -----------------------------------------------------------------
				block = index_node->get_son(0);
				rptr->leaf_node_ = new BLeafNode();
				rptr->leaf_node_->init_restore(tree, block);
				rptr->index_pos_ = 0;
				rptr->leaf_pos_ = 0;

				increment = rptr->leaf_node_->get_increment();
				num_entries = rptr->leaf_node_->get_num_entries();
				if (increment > num_entries) {
					rptr->size_ = num_entries;
				}
				else {
					rptr->size_ = increment;
				}
				page_io_++;
			}
			else {					
				// -----------------------------------------------------------------
				//  init left buffer
				// -----------------------------------------------------------------
				block = index_node->get_son(follow);
				lptr->leaf_node_ = new BLeafNode();
				lptr->leaf_node_->init_restore(tree, block);

				pos = lptr->leaf_node_->find_position_by_key(q_val);
				if (pos < 0) pos = 0;
				lptr->index_pos_ = pos;

				increment = lptr->leaf_node_->get_increment();
				if (pos == lptr->leaf_node_->get_num_keys() - 1) {
					num_entries = lptr->leaf_node_->get_num_entries();
					lptr->leaf_pos_ = num_entries - 1;
					lptr->size_ = num_entries - pos * increment;
				}
				else {
					lptr->leaf_pos_ = pos * increment + increment - 1;
					lptr->size_ = increment;
				}
				page_io_++;

				// -----------------------------------------------------------------
				//  init right buffer
				// -----------------------------------------------------------------
				if (pos < lptr->leaf_node_->get_num_keys() - 1) {
					rptr->leaf_node_ = lptr->leaf_node_;
					rptr->index_pos_ = (pos + 1);
					rptr->leaf_pos_  = (pos + 1) * increment;

					if ((pos + 1) == rptr->leaf_node_->get_num_keys() - 1) {
						num_entries = rptr->leaf_node_->get_num_entries();
						rptr->size_ = num_entries - (pos + 1) * increment;
					}
					else {
						rptr->size_ = increment;
					}
				}
				else {
					rptr->leaf_node_ = lptr->leaf_node_->get_right_sibling();
					if (rptr->leaf_node_) {
						rptr->index_pos_ = 0;
						rptr->leaf_pos_ = 0;

						increment = rptr->leaf_node_->get_increment();
						num_entries = rptr->leaf_node_->get_num_entries();
						if (increment > num_entries) {
							rptr->size_ = num_entries;
						}
						else {
							rptr->size_ = increment;
						}
						page_io_++;
					}
				}
			}
		}
		else {
			// ---------------------------------------------------------------------
			//  only one level in the B+ Tree: one lead node
			//  (1) init left buffer
			// ---------------------------------------------------------------------
			lptr->leaf_node_ = new BLeafNode();
			lptr->leaf_node_->init_restore(tree, block);

			pos = lptr->leaf_node_->find_position_by_key(q_val);
			if (pos < 0) pos = 0;
			lptr->index_pos_ = pos;

			increment = lptr->leaf_node_->get_increment();
			if (pos == lptr->leaf_node_->get_num_keys() - 1) {
				num_entries = lptr->leaf_node_->get_num_entries();
				lptr->leaf_pos_ = num_entries - 1;
				lptr->size_ = num_entries - pos * increment;
			}
			else {
				lptr->leaf_pos_ = pos * increment + increment - 1;
				lptr->size_ = increment;
			}
			page_io_++;
			
			// ---------------------------------------------------------------------
			//  (2) init right buffer
			// ---------------------------------------------------------------------
			if (pos < lptr->leaf_node_->get_num_keys() - 1) {
				rptr->leaf_node_ = lptr->leaf_node_;
				rptr->index_pos_ = pos + 1;
				rptr->leaf_pos_  = (pos + 1) * increment;

				if ((pos + 1) == rptr->leaf_node_->get_num_keys() - 1) {
					num_entries = rptr->leaf_node_->get_num_entries();
					rptr->size_ = num_entries - (pos + 1) * increment;
				}
				else {
					rptr->size_ = increment;
				}
			}
			else {
				rptr->leaf_node_ = NULL;
				rptr->index_pos_ = -1;
				rptr->leaf_pos_ = -1;
				rptr->size_ = -1;
			}
		}

		if (index_node != NULL) {
			delete index_node; index_node = NULL;
		}
	}
}

// -----------------------------------------------------------------------------
float QALSH::find_radius()			// find proper radius
{
	float radius = update_radius(1.0f / appr_ratio_);
	if (radius < 1.0f) radius = 1.0f;

	return radius;
}

// -----------------------------------------------------------------------------
float QALSH::update_radius(			// update radius
	float old_radius)					// old radius
{
	// -------------------------------------------------------------------------
	//  find an array of projected distance which is closest to the query in
	//  each of <m> hash tables 
	// -------------------------------------------------------------------------
	vector<float> list;
	for (int i = 0; i < m_; ++i) {
		if (lptr_[i]->size_ != -1) {
			list.push_back(calc_dist(q_val_[i], lptr_[i]));
		}
		if (rptr_[i]->size_ != -1) {
			list.push_back(calc_dist(q_val_[i], rptr_[i]));
		}
	}
	sort(list.begin(), list.end());

	// -------------------------------------------------------------------------
	//  find the median distance and return the new radius
	// -------------------------------------------------------------------------
	int num = (int) list.size();
	if (num == 0) return appr_ratio_ * old_radius;

	float dist = -1.0f;
	if (num % 2 == 0) dist = (list[num / 2 - 1] + list[num / 2]) / 2.0f;
	else dist = list[num / 2];
	
	int kappa = (int) ceil(log(2.0f * dist / w_) / log(appr_ratio_));
	return pow(appr_ratio_, kappa);
}

// -----------------------------------------------------------------------------
void QALSH::update_left_buffer(		// update left buffer
	const PageBuffer *rptr,				// right buffer
	PageBuffer *lptr)					// left buffer (return)
{
	BLeafNode *leaf_node = NULL;
	BLeafNode *old_leaf_node = NULL;

	if (lptr->index_pos_ > 0) {
		lptr->index_pos_--;

		int pos = lptr->index_pos_;
		int increment = lptr->leaf_node_->get_increment();

		lptr->leaf_pos_ = pos * increment + increment - 1;
		lptr->size_ = increment;
	}
	else {
		old_leaf_node = lptr->leaf_node_;
		leaf_node = lptr->leaf_node_->get_left_sibling();

		if (leaf_node) {
			lptr->leaf_node_ = leaf_node;
			lptr->index_pos_ = lptr->leaf_node_->get_num_keys() - 1;

			int pos = lptr->index_pos_;
			int increment = lptr->leaf_node_->get_increment();
			int num_entries = lptr->leaf_node_->get_num_entries();
			
			lptr->leaf_pos_ = num_entries - 1;
			lptr->size_ = num_entries - pos * increment;
			page_io_++;
		}
		else {
			lptr->leaf_node_ = NULL;
			lptr->index_pos_ = -1;
			lptr->leaf_pos_ = -1;
			lptr->size_ = -1;
		}

		if (rptr->leaf_node_ != old_leaf_node) {
			delete old_leaf_node; old_leaf_node = NULL;
		}
	}
}

// -----------------------------------------------------------------------------
void QALSH::update_right_buffer(	// update right buffer
	const PageBuffer *lptr,				// left buffer
	PageBuffer *rptr)					// right buffer (return)
{
	BLeafNode *leaf_node = NULL;
	BLeafNode *old_leaf_node = NULL;

	if (rptr->index_pos_ < rptr->leaf_node_->get_num_keys() - 1) {
		rptr->index_pos_++;

		int pos = rptr->index_pos_;
		int increment = rptr->leaf_node_->get_increment();
		
		rptr->leaf_pos_ = pos * increment;
		if (pos == rptr->leaf_node_->get_num_keys() - 1) {
			int num_entries = rptr->leaf_node_->get_num_entries();
			rptr->size_ = num_entries - pos * increment;
		}
		else {
			rptr->size_ = increment;
		}
	}
	else {
		old_leaf_node = rptr->leaf_node_;
		leaf_node = rptr->leaf_node_->get_right_sibling();

		if (leaf_node) {
			rptr->leaf_node_ = leaf_node;
			rptr->index_pos_ = 0;
			rptr->leaf_pos_ = 0;

			int increment = rptr->leaf_node_->get_increment();
			int num_entries = rptr->leaf_node_->get_num_entries();
			if (increment > num_entries) {
				rptr->size_ = num_entries;
			}
			else {
				rptr->size_ = increment;
			}
			page_io_++;
		}
		else {
			rptr->leaf_node_ = NULL;
			rptr->index_pos_ = -1;
			rptr->leaf_pos_ = -1;
			rptr->size_ = -1;
		}

		if (lptr->leaf_node_ != old_leaf_node) {
			delete old_leaf_node; old_leaf_node = NULL;
		}
	}
}

// -----------------------------------------------------------------------------
inline float QALSH::calc_dist(		// calc projected distance
	float q_val,						// hash value of query
	const PageBuffer *ptr)				// page buffer
{
	int   pos  = ptr->index_pos_;
	float key  = ptr->leaf_node_->get_key(pos);
	float dist = fabs(key - q_val);

	return dist;
}

// -----------------------------------------------------------------------------
void QALSH::delete_tree_ptr()		// delete the pointers of B+ Trees
{
	for (int i = 0; i < m_; ++i) {
		// ---------------------------------------------------------------------
		//  CANNOT remove the condition
		//              lptrs[i]->leaf_node != rptrs[i]->leaf_node
		//  because "lptrs[i]->leaf_node" and "rptrs[i]->leaf_node" may point 
		//  to the same address, then we would delete it twice and receive the 
		//  runtime error or segmentation fault.
		// ---------------------------------------------------------------------
		if (lptr_[i]->leaf_node_ && lptr_[i]->leaf_node_ != rptr_[i]->leaf_node_) {
			delete lptr_[i]->leaf_node_; lptr_[i]->leaf_node_ = NULL;
		}
		if (rptr_[i]->leaf_node_) {
			delete rptr_[i]->leaf_node_; rptr_[i]->leaf_node_ = NULL;
		}
	}
}
