#include "headers.h"


// -----------------------------------------------------------------------------
QALSH::QALSH(						// constructor
	int   n,							// cardinality
	int   d,							// dimensionality
	float p,							// l_p distance
	float zeta,							// a parameter of p-stable distr.
	float ratio,						// approximation ratio
	const float **data)					// data objects
{
	// -------------------------------------------------------------------------
	//  init basic parameters
	// -------------------------------------------------------------------------
	n_pts_      = n;
	dim_        = d;
	p_          = p;
	zeta_       = zeta;
	appr_ratio_ = ratio;
	data_       = data;
	delta_      = 1.0f / E;
	beta_       = (float) CANDIDATES / (float) n;

	// -------------------------------------------------------------------------
	//  init <w_> <p1_> and <p2_> (auto tuning-w)
	//  
	//  w0 ----- best w for L_{0.5} norm to minimize m (auto tuning-w)
	//  w1 ----- best w for L_{1.0} norm to minimize m (auto tuning-w)
	//  w2 ----- best w for L_{2.0} norm to minimize m (auto tuning-w)
	//  other w: use linear combination for interpolation
	// -------------------------------------------------------------------------
	float w0 = (ratio - 1.0f) / log(sqrt(ratio));
	float w1 = 2.0f * sqrt(ratio);
	float w2 = sqrt((8.0f*ratio*ratio*log(ratio)) / (ratio*ratio- 1.0f));

	if (fabs(p_ - 0.5f) < FLOATZERO) {
		w_  = w0;
		p1_ = calc_l0_prob(w_ / 2.0f);
		p2_ = calc_l0_prob(w_ / (2.0f * appr_ratio_));
	}
	else if (fabs(p_ - 1.0f) < FLOATZERO) {
		w_  = w1;
		p1_ = calc_l1_prob(w_ / 2.0f);
		p2_ = calc_l1_prob(w_ / (2.0f * appr_ratio_));
	}
	else if (fabs(p_ - 2.0f) < FLOATZERO) {
		w_  = w2;
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

	// -------------------------------------------------------------------------
	//  generate hash functions
	// -------------------------------------------------------------------------
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

	// -------------------------------------------------------------------------
	//  bulkloading
	// -------------------------------------------------------------------------
	freq_        = new int[n_pts_];
	lpos_        = new int[m_];
	rpos_        = new int[m_];
	checked_     = new bool[n_pts_];
	bucket_flag_ = new bool[m_];
	range_flag_  = new bool[m_];
	q_val_       = new float[m_];
	tables_      = new Result*[m_];

	for (int i = 0; i < m_; ++i) {
		tables_[i] = new Result[n_pts_];
		for (int j = 0; j < n_pts_; ++j) {
			tables_[i][j].id_  = j;
			tables_[i][j].key_ = calc_hash_value(i, data_[j]);
		}
		qsort(tables_[i], n_pts_, sizeof(Result), ResultComp);
	}
}

// -----------------------------------------------------------------------------
inline float QALSH::calc_l0_prob(	// calc prob of L1/2 dist
	float x)							// x = w / (2.0 * r)
{
	return new_levy_prob(x);
}

// -----------------------------------------------------------------------------
inline float QALSH::calc_l1_prob(	// calc prob of L1 dist
	float x)							// x = w / (2.0 * r)
{
	return new_cauchy_prob(x);
}

// -----------------------------------------------------------------------------
inline float QALSH::calc_l2_prob(	// calc prob of L2 dist
	float x)							// x = w / (2.0 * r)
{
	return new_gaussian_prob(x);
}

// -----------------------------------------------------------------------------
inline float QALSH::calc_hash_value( // calc hash value
	int   table_id,						// hash table id
	const float *data)					// one data/query object
{
	float ret  = 0.0f;
	int   base = table_id * dim_;
	for (int i = 0; i < dim_; ++i) {
		ret += (a_array_[base + i] * data[i]);
	}
	return ret;
}

// -----------------------------------------------------------------------------
QALSH::~QALSH()						// destructor
{
	delete[] a_array_;     a_array_     = NULL;
	delete[] freq_;        freq_        = NULL;
	delete[] lpos_;        lpos_        = NULL;
	delete[] rpos_;        rpos_        = NULL;
	delete[] checked_;     checked_     = NULL;
	delete[] bucket_flag_; bucket_flag_ = NULL;
	delete[] range_flag_;  range_flag_  = NULL;
	delete[] q_val_;       q_val_       = NULL;

	for (int i = 0; i < m_; ++i) {
		delete[] tables_[i]; tables_[i] = NULL;
	}
	delete[] tables_; tables_ = NULL;
}

// -----------------------------------------------------------------------------
void QALSH::display()				// display parameters
{
	printf("Parameters of QALSH:\n");
	printf("    n     = %d\n",   n_pts_);
	printf("    d     = %d\n",   dim_);
	printf("    p     = %.1f\n", p_);
	printf("    zeta  = %f\n",   zeta_);
	printf("    ratio = %.1f\n", appr_ratio_);
	printf("    w     = %f\n",   w_);
	printf("    p1    = %f\n",   p1_);
	printf("    p2    = %f\n",   p2_);
	printf("    alpha = %f\n",   alpha_);
	printf("    beta  = %f\n",   beta_);
	printf("    delta = %f\n",   delta_);
	printf("    m     = %d\n",   m_);
	printf("    l     = %d\n\n", l_);
}

// -----------------------------------------------------------------------------
int QALSH::knn(						// k-nn search
	int   top_k,						// top-k value
	const float *query,					// input query object
	MinK_List *list)					// k-NN results (return)
{
	int candidates = CANDIDATES + top_k - 1; // candidate size
	float kdist = MAXREAL;

	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	memset(freq_, 0, n_pts_ * SIZEFLOAT);
	memset(checked_, false, n_pts_ * SIZEBOOL);
	memset(bucket_flag_, true, m_ * SIZEBOOL);

	for (int i = 0; i < m_; ++i) {
		q_val_[i] = calc_hash_value(i, query);
		int pos = binary_search_pos(i, q_val_[i]);
		if (pos == 0) {
			lpos_[i] = -1;  rpos_[i] = pos;
		} 
		else {
			lpos_[i] = pos; rpos_[i] = pos + 1;
		}
	}

	// -------------------------------------------------------------------------
	//  c-k-ANN search
	// -------------------------------------------------------------------------
	int   dist_cnt     = 0;			// number of candidates computation
	float radius       = 1.0f;		// search radius
	float bucket_width = radius * w_ / 2.0f; // bucket width

	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialize the stop condition for current round
		// ---------------------------------------------------------------------
		int num_bucket = 0;
		memset(bucket_flag_, true, m_ * SIZEBOOL);

		// ---------------------------------------------------------------------
		//  step 2: (R,c)-NN search
		// ---------------------------------------------------------------------
		int count = 0;
		while (num_bucket < m_) {
			float  ldist  = -1.0f;	// left  proj dist with query
			float  rdist  = -1.0f;	// right proj dist with query
			Result *table = NULL;

			for (int j = 0; j < m_; ++j) {
				if (!bucket_flag_[j]) continue;

				table = tables_[j];
				// -------------------------------------------------------------
				//  step 2.1: scan the left part of hash table
				// -------------------------------------------------------------
				count = 0;
				while (count < SCAN_SIZE) {
					ldist = MAXREAL;
					if (lpos_[j] >= 0) {
						ldist = fabs(q_val_[j] - table[lpos_[j]].key_);
					}
					if (ldist > bucket_width) break;

					int id = table[lpos_[j]].id_;
					if (++freq_[id] >= l_ && !checked_[id]) {
						checked_[id] = true;
						float dist = calc_lp_dist(dim_, p_, data_[id], query);
						kdist = list->insert(dist, id);

						if (++dist_cnt >= candidates) break;
					}
					--lpos_[j];
					++count;
				}
				if (dist_cnt >= candidates) break;

				// -------------------------------------------------------------
				//  step 2.2: scan the right part of hash table
				// -------------------------------------------------------------
				count = 0;
				while (count < SCAN_SIZE) {
					rdist = MAXREAL;
					if (rpos_[j] < n_pts_) {
						rdist = fabs(q_val_[j] - table[rpos_[j]].key_);
					}
					if (rdist > bucket_width) break;

					int id = table[rpos_[j]].id_;
					if (++freq_[id] >= l_ && !checked_[id]) {
						checked_[id] = true;
						float dist = calc_lp_dist(dim_, p_, data_[id], query);
						kdist = list->insert(dist, id);

						if (++dist_cnt >= candidates) break;
					}
					++rpos_[j];
					++count;
				}
				if (dist_cnt >= candidates) break;

				// -------------------------------------------------------------
				//  step 2.3: check whether this bucket_width is finished scanned
				// -------------------------------------------------------------
				if (ldist > bucket_width && rdist > bucket_width) {
					bucket_flag_[j] = false;
					++num_bucket;
					if (num_bucket >= m_) break;
				}
			}
			if (num_bucket >= m_ || dist_cnt >= candidates) break;
		}
		// ---------------------------------------------------------------------
		//  step 3: stop conditions t1 and t2
		// ---------------------------------------------------------------------
		if (kdist < appr_ratio_ * radius && dist_cnt >= top_k) break;
		if (dist_cnt >= candidates) break;

		// ---------------------------------------------------------------------
		//  step 4: auto-update radius
		// ---------------------------------------------------------------------
		radius = appr_ratio_ * radius;
		bucket_width = radius * w_ / 2.0f;
	}
	return 0;
}

// -----------------------------------------------------------------------------
int QALSH::knn(						// k-NN search
	int   top_k,						// top-k value
	float R,							// limited search range
	const float *query,					// input query object
	const vector<int> &object_id,		// object id mapping
	MinK_List *list)					// k-NN results (return)
{
	int candidates = CANDIDATES + top_k - 1; // candidate size
	float kdist = MAXREAL;			// k-th ANN distance

	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	memset(freq_, 0, n_pts_ * SIZEFLOAT);
	memset(checked_, false, n_pts_ * SIZEBOOL);
	memset(bucket_flag_, true, m_ * SIZEBOOL);
	memset(range_flag_, true, m_ * SIZEBOOL);

	for (int i = 0; i < m_; ++i) {
		q_val_[i] = calc_hash_value(i, query);
		int pos = binary_search_pos(i, q_val_[i]);
		if (pos == 0) {
			lpos_[i] = -1;  rpos_[i] = pos;
		} else {
			lpos_[i] = pos; rpos_[i] = pos + 1;
		}
	}

	// -------------------------------------------------------------------------
	//  k-nn search via dynamic collision counting
	// -------------------------------------------------------------------------
	int   dist_cnt     = 0;			// number of candidates computation
	int   num_range    = 0;			// number of search range flag
	float radius       = 1.0f;		// search radius
	float bucket_width = radius * w_ / 2.0f;  // bucket width
	
	float range = -1.0f;			// limited search range
	if (R > MAXREAL - 1.0f) range = MAXREAL;
	else range = R * w_ / 2.0f;

	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialize the stop condition for current round
		// ---------------------------------------------------------------------
		int num_bucket = 0;
		memset(bucket_flag_, true, m_ * SIZEBOOL);

		// ---------------------------------------------------------------------
		//  step 2: (R,c)-NN search
		// ---------------------------------------------------------------------
		int count = 0;
		while (num_bucket < m_ && num_range < m_) {
			float  ldist  = -1.0f;	// left  proj dist to query
			float  rdist  = -1.0f;	// right proj dist to query
			Result *table = NULL;

			for (int j = 0; j < m_; ++j) {
				if (!bucket_flag_[j]) continue;

				table = tables_[j];
				// -------------------------------------------------------------
				//  step 2.1: scan the left part of hash table
				// -------------------------------------------------------------
				count = 0;
				while (count < SCAN_SIZE) {
					ldist = MAXREAL;
					if (lpos_[j] >= 0) {
						ldist = fabs(q_val_[j] - table[lpos_[j]].key_);
					}
					if (ldist > bucket_width || ldist > range) break;

					int id = table[lpos_[j]].id_;
					if (++freq_[id] >= l_ && !checked_[id]) {
						checked_[id] = true;
						float dist = calc_lp_dist(dim_, p_, data_[id], query);
						kdist = list->insert(dist, object_id[id]);

						if (++dist_cnt >= candidates) break;
					}
					--lpos_[j];
					++count;
				}
				if (dist_cnt >= candidates) break;

				// -------------------------------------------------------------
				//  step 2.2: scan the right part of hash table
				// -------------------------------------------------------------
				count = 0;
				while (count < SCAN_SIZE) {
					rdist = MAXREAL;
					if (rpos_[j] < n_pts_) {
						rdist = fabs(q_val_[j] - table[rpos_[j]].key_);
					}
					if (rdist > bucket_width || rdist > range) break;

					int id = tables_[j][rpos_[j]].id_;
					if (++freq_[id] >= l_ && !checked_[id]) {
						checked_[id] = true;
						float dist = calc_lp_dist(dim_, p_, data_[id], query);
						list->insert(dist, object_id[id]);

						if (++dist_cnt >= candidates) break;
					}
					++rpos_[j];
					++count;
				}
				if (dist_cnt >= candidates) break;

				// -------------------------------------------------------------
				//  step 2.3: check whether this bucket_width is finished scanned
				// -------------------------------------------------------------
				if (ldist > bucket_width && rdist > bucket_width) {
					bucket_flag_[j] = false;
					++num_bucket;
				}
				if (ldist > range && rdist > range) {
					bucket_flag_[j] = false;
					++num_bucket;
					if (range_flag_[j]) {
						range_flag_[j] = false;
						++num_range;
					}
				}
				if (num_bucket >= m_ || num_range >= m_) break;
			}
			if (num_bucket >= m_ || num_range >= m_) break;
			if (dist_cnt >= candidates) break;
		}

		// ---------------------------------------------------------------------
		//  step 3: terminating conditions
		// ---------------------------------------------------------------------
		if (num_range >= m_ || dist_cnt >= candidates) break;

		// ---------------------------------------------------------------------
		//  step 4: auto-update <radius>
		// ---------------------------------------------------------------------
		radius = appr_ratio_ * radius;
		bucket_width = radius * w_ / 2.0f;
	}
	return 0;
}

// -----------------------------------------------------------------------------
int QALSH::binary_search_pos(		// binary search position
	int table_id,						// hash table is
	float value)						// hash value
{
	int left = 0;
	int right = n_pts_ - 1;
	int mid = 0;

	Result *table = tables_[table_id];
	while (left < right) {
		mid = (left + right + 1) / 2;
		if (fabs(table[mid].key_ - value) < FLOATZERO) {
			return mid;
		}
		if (table[mid].key_ < value) left = mid;
		else right = mid - 1;
	}
	// assert(left >= 0 && left < n_pts_);

	return left;
}











