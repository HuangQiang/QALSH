#include "qalsh.h"

// -----------------------------------------------------------------------------
QALSH::QALSH()						// constructor
{
	n_pts_ = -1;
	dim_   = -1;
	B_     = -1;
	p_     = -1.0f;
	zeta_  = -1.0f;
	ratio_ = -1.0f;
	w_     = -1.0f;
	m_     = -1;
	l_     = -1;
	a_     = NULL;
	trees_ = NULL;

	dist_io_ = -1;
	page_io_ = -1;
}

// -----------------------------------------------------------------------------
QALSH::~QALSH()						// destructor
{
	for (int i = 0; i < m_; ++i) {
		delete[] a_[i]; a_[i] = NULL;
		delete trees_[i]; trees_[i] = NULL;
	}
	delete[] a_; a_ = NULL;
	delete[] trees_; trees_ = NULL;

	g_memory -= SIZEFLOAT * m_ * dim_;
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
	const char *path)					// index path
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_ = n;
	dim_   = d;
	B_     = B;
	p_     = p;
	zeta_  = zeta;
	ratio_ = ratio;

	strcpy(path_, path); 
	create_dir(path_);

	// -------------------------------------------------------------------------
	//  init <w_> <m_> and <l_> (auto tuning-w)
	//  
	//  w0 ----- best w for L_{0.5} norm to minimize m (auto tuning-w)
	//  w1 ----- best w for L_{1.0} norm to minimize m (auto tuning-w)
	//  w2 ----- best w for L_{2.0} norm to minimize m (auto tuning-w)
	//  other w: use linear combination for interpolation
	// -------------------------------------------------------------------------
	float delta = 1.0f / E;
	float beta  = (float) CANDIDATES / (float) n_pts_;

	float w0 = (ratio_ - 1.0f) / log(sqrt(ratio_));
	float w1 = 2.0f * sqrt(ratio_);
	float w2 = sqrt((8.0f * SQR(ratio_) * log(ratio_)) / (SQR(ratio_) - 1.0f));
	float p1 = -1.0f, p2 = -1.0f;

	if (fabs(p_ - 0.5f) < FLOATZERO) {
		w_ = w0;
		p1 = calc_l0_prob(w_ / 2.0f);
		p2 = calc_l0_prob(w_ / (2.0f * ratio_));
	}
	else if (fabs(p_ - 1.0f) < FLOATZERO) {
		w_ = w1;
		p1 = calc_l1_prob(w_ / 2.0f);
		p2 = calc_l1_prob(w_ / (2.0f * ratio_));
	}
	else if (fabs(p_ - 2.0f) < FLOATZERO) {
		w_ = w2;
		p1 = calc_l2_prob(w_ / 2.0f);
		p2 = calc_l2_prob(w_ / (2.0f * ratio_));
	}
	else {
		if (fabs(p_-0.8f) < FLOATZERO) w_ = 2.503f;
		else if (fabs(p_-1.2f) < FLOATZERO) w_ = 3.151f;
		else if (fabs(p_-1.5f) < FLOATZERO) w_ = 3.465f;
		else w_ = (w2 - w1) * p_ + (2.0f * w1 - w2);

		new_stable_prob(p_, zeta_, ratio_, 1.0f, w_, 1000000, p1, p2);
	}

	float para1 = sqrt(log(2.0f / beta));
	float para2 = sqrt(log(1.0f / delta));
	float para3 = 2.0f * (p1 - p2) * (p1 - p2);
	float eta   = para1 / para2;
	float alpha = (eta * p1 + p2) / (1.0f + eta);

	m_ = (int) ceil((para1 + para2) * (para1 + para2) / para3);
	l_ = (int) ceil(alpha * m_);

	// -------------------------------------------------------------------------
	//  generate hash functions
	// -------------------------------------------------------------------------
	g_memory += SIZEFLOAT * m_ * dim_;
	a_ = new float*[m_];
	for (int i = 0; i < m_; ++i) {
		a_[i] = new float[dim_];
		for (int j = 0; j < dim_; ++j) {
			if (fabs(p_-0.5f) < FLOATZERO) a_[i][j] = levy(1.0f, 0.0f);
			else if (fabs(p_-1.0f) < FLOATZERO) a_[i][j] = cauchy(1.0f, 0.0f);
			else if (fabs(p_-2.0f) < FLOATZERO) a_[i][j] = gaussian(0.0f, 1.0f);
			else a_[i][j] = p_stable(p_, zeta_, 1.0f, 0.0f);
		}
	}

	// -------------------------------------------------------------------------
	//  write parameters to disk
	// -------------------------------------------------------------------------
	if (write_params()) return 1;

	// -------------------------------------------------------------------------
	//  bulkloading
	// -------------------------------------------------------------------------
	if (bulkload(data)) return 1;

	return 0;
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
int QALSH::write_params()			// write parameters to disk
{
	char fname[200];
	strcpy(fname, path_); strcat(fname, "para");

	FILE *fp = fopen(fname, "r");
	if (fp)	{ printf("Hash Tables Already Exist\n\n"); exit(1); } // todo modify
	fp = fopen(fname, "w");
	if (!fp) {
		printf("Could not create %s\n", fname);
		printf("Perhaps no such folder %s?\n", path_);
		return 1;
	}

	fprintf(fp, "n     = %d\n", n_pts_);
	fprintf(fp, "d     = %d\n", dim_);
	fprintf(fp, "B     = %d\n", B_);
	fprintf(fp, "p     = %f\n", p_);
	fprintf(fp, "zeta  = %f\n", zeta_);
	fprintf(fp, "ratio = %f\n", ratio_);
	fprintf(fp, "w     = %f\n", w_);
	fprintf(fp, "m     = %d\n", m_);
	fprintf(fp, "l     = %d\n", l_);

	for (int i = 0; i < m_; ++i) {
		for (int j = 0; j < dim_; ++j) fprintf(fp, "%f ", a_[i][j]);
		fprintf(fp, "\n");
	}
	fclose(fp);

	return 0;
}

// -----------------------------------------------------------------------------
int QALSH::bulkload(				// build m b-trees by bulkloading
	const float **data)					// data set
{
	char fname[200];
	Result *table = new Result[n_pts_];

	trees_ = new BTree*[m_];
	for (int i = 0; i < m_; ++i) {
		for (int j = 0; j < n_pts_; ++j) {
			table[j].id_  = j;
			table[j].key_ = calc_hash_value(i, data[j]);
		}
		qsort(table, n_pts_, sizeof(Result), ResultComp);
		
		get_tree_filename(i, fname);
		trees_[i] = new BTree();
		trees_[i]->init(B_, fname);
		if (trees_[i]->bulkload(n_pts_, table)) return 1;
	}
	delete[] table; table = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
float QALSH::calc_hash_value( 		// calc hash value
	int   tid,							// hash table id
	const float *data)					// input data object
{
	return calc_inner_product(dim_, a_[tid], data);
}

// -----------------------------------------------------------------------------
inline void QALSH::get_tree_filename( // get file name of b-tree
	int  tid,							// tree id, from 0 to m-1
	char *fname)						// file name (return)
{
	sprintf(fname, "%s%d.qalsh", path_, tid);
}

// -----------------------------------------------------------------------------
void QALSH::display()				// display parameters
{
	printf("Parameters of QALSH:\n");
	printf("    n     = %d\n",   n_pts_);
	printf("    d     = %d\n",   dim_);
	printf("    B     = %d\n",   B_);
	printf("    p     = %.1f\n", p_);
	printf("    zeta  = %.1f\n", zeta_);
	printf("    ratio = %.1f\n", ratio_);
	printf("    w     = %f\n",   w_);
	printf("    m     = %d\n",   m_);
	printf("    l     = %d\n",   l_);
	printf("    path  = %s\n",   path_);
	printf("\n");
}

// -----------------------------------------------------------------------------
int QALSH::load(					// load index
	const char *path)					// index path
{
	// -------------------------------------------------------------------------
	//  read parameters from disk
	// -------------------------------------------------------------------------
	strcpy(path_, path);
	if (read_params()) return 1;

	// -------------------------------------------------------------------------
	//  init b+ trees for k-NN search
	// -------------------------------------------------------------------------
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
	strcpy(fname, path_); strcat(fname, "para");

	FILE *fp = fopen(fname, "r");
	if (!fp) { printf("Could not open %s\n", fname); return 1; }

	fscanf(fp, "n     = %d\n", &n_pts_);
	fscanf(fp, "d     = %d\n", &dim_);
	fscanf(fp, "B     = %d\n", &B_);
	fscanf(fp, "p     = %f\n", &p_);
	fscanf(fp, "zeta  = %f\n", &zeta_);
	fscanf(fp, "ratio = %f\n", &ratio_);
	fscanf(fp, "w     = %f\n", &w_);
	fscanf(fp, "m     = %d\n", &m_);
	fscanf(fp, "l     = %d\n", &l_);
	
	g_memory += SIZEFLOAT * m_ * dim_;
	a_ = new float*[m_];
	for (int i = 0; i < m_; ++i) {
		a_[i] = new float[dim_];
		for (int j = 0; j < dim_; ++j) fscanf(fp, "%f ", &a_[i][j]);
		fscanf(fp, "\n");
	}
	fclose(fp);

	return 0;
}

// -----------------------------------------------------------------------------
uint64_t QALSH::knn(				// k-NN search
	int   top_k,						// top-k value
	const float *query,					// query object
	const char *data_folder,			// data folder
	MinK_List *list)					// k-NN results (return)
{
	int   *freq    = new int[n_pts_];
	bool  *checked = new bool[n_pts_];
	bool  *flag    = new bool[m_];
	float *q_val   = new float[m_];
	float *data    = new float[dim_];

	Page **lptrs = new Page*[m_];
	Page **rptrs = new Page*[m_];
	for (int i = 0; i < m_; ++i) {
		lptrs[i] = new Page();
		rptrs[i] = new Page();
	}

	// -------------------------------------------------------------------------
	//  initialize parameters for k-NN search
	// -------------------------------------------------------------------------
	memset(freq, 0, n_pts_ * SIZEFLOAT);
	memset(checked, false, n_pts_ * SIZEBOOL);

	init_search_params(query, q_val, lptrs, rptrs);

	// -------------------------------------------------------------------------
	//  k-NN search
	// -------------------------------------------------------------------------
	int   candidates = CANDIDATES + top_k - 1; // candidates size
	float kdist = MAXREAL;
	float radius = find_radius(q_val, (const Page**) lptrs, (const Page**) rptrs);
	float bucket = w_ * radius / 2.0f;

	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialize the stop condition for current round
		// ---------------------------------------------------------------------
		int num_flag = 0;
		memset(flag, true, m_ * SIZEBOOL);

		// ---------------------------------------------------------------------
		//  step 2: find frequent objects
		// ---------------------------------------------------------------------
		while (num_flag < m_) {
			for (int i = 0; i < m_; ++i) {
				if (!flag[i]) continue;

				// -------------------------------------------------------------
				//  step 2.1: compute <ldist> and <rdist>
				// -------------------------------------------------------------
				Page *lptr = lptrs[i];
				Page *rptr = rptrs[i];

				float ldist = MAXREAL;
				float rdist = MAXREAL;
				if (lptr->size_ != -1) ldist = calc_dist(q_val[i], lptr);
				if (rptr->size_ != -1) rdist = calc_dist(q_val[i], rptr);

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
						if (++freq[id] > l_ && !checked[id]) {
							checked[id] = true;
							read_data_new_format(id, dim_, B_, data_folder, data);

							float dist = calc_lp_dist(dim_, p_, kdist, data, query);
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
						if (++freq[id] > l_ && !checked[id]) {
							checked[id] = true;
							read_data_new_format(id, dim_, B_, data_folder, data);

							float dist = calc_lp_dist(dim_, p_, kdist, data, query);
							kdist = list->insert(dist, id);
							if (++dist_io_ >= candidates) break;
						}
					}
					update_right_buffer(lptr, rptr);
				}
				else {
					flag[i] = false;
					++num_flag;
				}
				if (num_flag >= m_ || dist_io_ >= candidates) break;
			}
			if (num_flag >= m_ || dist_io_ >= candidates) break;
		}
		// ---------------------------------------------------------------------
		//  step 3: stop conditions 1 & 2
		// ---------------------------------------------------------------------
		if (kdist < ratio_ * radius && dist_io_ >= top_k) break;
		if (dist_io_ >= candidates) break;

		// ---------------------------------------------------------------------
		//  step 4: auto-update <radius>
		// ---------------------------------------------------------------------
		radius = update_radius(radius, q_val, (const Page**) lptrs, 
			(const Page**) rptrs);
		bucket = radius * w_ / 2.0f;
	}
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete_tree_ptr(lptrs, rptrs);

	delete[] freq;    freq    = NULL;
	delete[] checked; checked = NULL;
	delete[] flag;    flag    = NULL;
	delete[] q_val;   q_val   = NULL;
	delete[] data;    data    = NULL;

	return page_io_ + dist_io_;
}

// -----------------------------------------------------------------------------
uint64_t QALSH::knn(				// k-NN search
	int   top_k,						// top-k value
	float R,							// limited search range
	const float *query,					// query object
	const int   *object_id, 			// object id mapping
	const char *data_folder,			// data folder
	MinK_List *list)					// k-NN results (return)
{
	int   *freq        = new int[n_pts_];
	bool  *checked     = new bool[n_pts_];
	bool  *bucket_flag = new bool[m_];
	bool  *range_flag  = new bool[m_];
	float *q_val       = new float[m_];
	float *data        = new float[dim_];

	Page **lptrs = new Page*[m_];
	Page **rptrs = new Page*[m_];
	for (int i = 0; i < m_; ++i) {
		lptrs[i] = new Page();
		rptrs[i] = new Page();
	}

	// -------------------------------------------------------------------------
	//  initialize parameters for k-NN search
	// -------------------------------------------------------------------------
	memset(freq, 0, n_pts_ * SIZEFLOAT);
	memset(checked, false, n_pts_ * SIZEBOOL);
	memset(range_flag, true, m_ * SIZEBOOL);

	init_search_params(query, q_val, lptrs, rptrs);

	// -------------------------------------------------------------------------
	//  k-NN search
	// -------------------------------------------------------------------------
	int candidates = CANDIDATES + top_k - 1; // candidates size
	int num_range = 0;				// used for search range bound
	
	float kdist  = R;
	float radius = find_radius(q_val, (const Page**) lptrs, (const Page**) rptrs);
	float bucket = w_ * radius / 2.0f;
	float range  = R > MAXREAL - 1.0f ? MAXREAL : R * w_ / 2.0f;

	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialize the stop condition for current round
		// ---------------------------------------------------------------------
		int num_bucket = 0;
		memset(bucket_flag, true, m_ * SIZEBOOL);

		// ---------------------------------------------------------------------
		//  step 2: find frequent objects
		// ---------------------------------------------------------------------
		while (num_bucket < m_ && num_range < m_) {
			for (int i = 0; i < m_; ++i) {
				if (!bucket_flag[i]) continue;

				// -------------------------------------------------------------
				//  step 2.1: compute <ldist> and <rdist>
				// -------------------------------------------------------------
				Page *lptr = lptrs[i];
				Page *rptr = rptrs[i];

				float ldist = MAXREAL;
				float rdist = MAXREAL;
				if (lptr->size_ != -1) ldist = calc_dist(q_val[i], lptr);
				if (rptr->size_ != -1) rdist = calc_dist(q_val[i], rptr);

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
						if (++freq[id] > l_ && !checked[id]) {
							checked[id] = true;
							int oid = object_id[id];
							read_data_new_format(oid, dim_, B_, data_folder, data);

							float dist  = calc_lp_dist(dim_, p_, kdist, data, query);
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
						if (++freq[id] > l_ && !checked[id]) {
							checked[id] = true;
							int oid = object_id[id];
							read_data_new_format(oid, dim_, B_, data_folder, data);

							float dist  = calc_lp_dist(dim_, p_, kdist, data, query);
							kdist = list->insert(dist, oid);
							if (++dist_io_ >= candidates) break;
						}
					}
					update_right_buffer(lptr, rptr);
				}
				else {
					bucket_flag[i] = false;
					num_bucket++;
					if (ldist >= range && rdist >= range && range_flag[i]) {
						range_flag[i] = false;
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
		radius = update_radius(radius, q_val, (const Page**) lptrs, 
			(const Page**) rptrs);
		bucket = radius * w_ / 2.0f;
	}
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete_tree_ptr(lptrs, rptrs);

	delete[] freq; freq = NULL;
	delete[] checked; checked = NULL;
	delete[] bucket_flag; bucket_flag = NULL;
	delete[] range_flag; range_flag = NULL;
	delete[] q_val; q_val = NULL;
	delete[] data; data = NULL;

	return page_io_ + dist_io_;
}

// -----------------------------------------------------------------------------
void QALSH::init_search_params(		// init parameters for k-NN search
	const float *query,					// query object
	float *q_val,						// hash values of query (return)
	Page  **lptrs,						// left buffer (return)
	Page  **rptrs)						// right buffer (return)
{
	page_io_ = 0;
	dist_io_ = 0;

	for (int i = 0; i < m_; ++i) {
		lptrs[i]->leaf_node_ = NULL;
		lptrs[i]->index_pos_ = -1;
		lptrs[i]->leaf_pos_  = -1;
		lptrs[i]->size_      = -1;

		rptrs[i]->leaf_node_ = NULL;
		rptrs[i]->index_pos_ = -1;
		rptrs[i]->leaf_pos_  = -1;
		rptrs[i]->size_      = -1;
	}

	BIndexNode *index_node = NULL;
	int  block       = -1;			// variables for index node
	int  follow      = -1;
	bool lescape     = false;
	int  pos         = -1;			// variables for leaf node
	int  increment   = -1;
	int  num_entries = -1;

	for (int i = 0; i < m_; ++i) {
		float q_v = calc_hash_value(i, query);
		BTree *tree = trees_[i];
		Page  *lptr = lptrs[i];
		Page  *rptr = rptrs[i];

		q_val[i] = q_v;
		block = tree->root_;
		if (block > 1) {
			// ---------------------------------------------------------------------
			//  at least two levels in the B+ Tree: index node and lead node
			// ---------------------------------------------------------------------
			index_node = new BIndexNode();
			index_node->init_restore(tree, block);
			++page_io_;

			// ---------------------------------------------------------------------
			//  find the leaf node whose value is closest and larger than the key
			//  of query q
			// ---------------------------------------------------------------------
			lescape = false;		// locate the position of branch
			while (index_node->get_level() > 1) {
				follow = index_node->find_position_by_key(q_v);

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
				++page_io_;			// access a new node (a new page)
			}

			// ---------------------------------------------------------------------
			//  after finding the leaf node whose value is closest to the key of
			//  query, initialize <lptrs[i]> and <rptrs[i]>.
			//
			//  <lescape> = true is that the query has no <lptrs>, the query is 
			//  the smallest value.
			// ---------------------------------------------------------------------
			follow = index_node->find_position_by_key(q_v);
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
				++page_io_;
			}
			else {					
				// -----------------------------------------------------------------
				//  init left buffer
				// -----------------------------------------------------------------
				block = index_node->get_son(follow);
				lptr->leaf_node_ = new BLeafNode();
				lptr->leaf_node_->init_restore(tree, block);

				pos = lptr->leaf_node_->find_position_by_key(q_v);
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
				++page_io_;

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
						++page_io_;
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

			pos = lptr->leaf_node_->find_position_by_key(q_v);
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
			++page_io_;
			
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
float QALSH::find_radius(			// find proper radius
	const float *q_val,					// hash value of query
	const Page **lptrs,					// left buffer
	const Page **rptrs)					// right buffer
{
	float radius = update_radius(1.0f / ratio_, q_val, lptrs, rptrs);
	if (radius < 1.0f) radius = 1.0f;

	return radius;
}

// -----------------------------------------------------------------------------
float QALSH::update_radius(			// update radius
	float old_radius,					// old radius
	const float *q_val,					// hash value of query
	const Page **lptrs,					// left buffer
	const Page **rptrs)					// right buffer
{
	// -------------------------------------------------------------------------
	//  find an array of projected distance which is closest to the query in
	//  each of <m> hash tables 
	// -------------------------------------------------------------------------
	std::vector<float> list;
	for (int i = 0; i < m_; ++i) {
		if (lptrs[i]->size_ != -1) {
			list.push_back(calc_dist(q_val[i], lptrs[i]));
		}
		if (rptrs[i]->size_ != -1) {
			list.push_back(calc_dist(q_val[i], rptrs[i]));
		}
	}
	std::sort(list.begin(), list.end());

	// -------------------------------------------------------------------------
	//  find the median distance and return the new radius
	// -------------------------------------------------------------------------
	int num = (int) list.size();
	if (num == 0) return ratio_ * old_radius;

	float dist = -1.0f;
	if (num % 2 == 0) dist = (list[num / 2 - 1] + list[num / 2]) / 2.0f;
	else dist = list[num / 2];
	
	int kappa = (int) ceil(log(2.0f * dist / w_) / log(ratio_));
	return pow(ratio_, kappa);
}

// -----------------------------------------------------------------------------
void QALSH::update_left_buffer(		// update left buffer
	const Page *rptr,					// right buffer
	Page *lptr)							// left buffer (return)
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
			++page_io_;
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
	const Page *lptr,					// left buffer
	Page *rptr)							// right buffer (return)
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
			++page_io_;
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
	const Page *ptr)					// page buffer
{
	int   pos  = ptr->index_pos_;
	float key  = ptr->leaf_node_->get_key(pos);

	return fabs(key - q_val);
}

// -----------------------------------------------------------------------------
void QALSH::delete_tree_ptr(		// delete the pointers of B+ Trees
	Page **lptrs,						// left buffer (return)
	Page **rptrs)						// right buffer (return)
{
	for (int i = 0; i < m_; ++i) {
		// ---------------------------------------------------------------------
		//  CANNOT remove the condition
		//              lptrs[i]->leaf_node != rptrs[i]->leaf_node
		//  because "lptrs[i]->leaf_node" and "rptrs[i]->leaf_node" may point 
		//  to the same address, then we would delete it twice and receive the 
		//  runtime error or segmentation fault.
		// ---------------------------------------------------------------------
		if (lptrs[i]->leaf_node_ && lptrs[i]->leaf_node_ != rptrs[i]->leaf_node_) {
			delete lptrs[i]->leaf_node_; lptrs[i]->leaf_node_ = NULL;
		}
		if (rptrs[i]->leaf_node_) {
			delete rptrs[i]->leaf_node_; rptrs[i]->leaf_node_ = NULL;
		}

		delete lptrs[i]; lptrs[i] = NULL;
		delete rptrs[i]; rptrs[i] = NULL;
	}
	delete[] lptrs; lptrs = NULL;
	delete[] rptrs; rptrs = NULL;
}
