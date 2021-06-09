#pragma once

#include <algorithm>
#include <cstring>
#include <vector>

#include "def.h"
#include "util.h"
#include "random.h"
#include "pri_queue.h"
#include "b_node.h"
#include "b_tree.h"

namespace nns {

// -----------------------------------------------------------------------------
//  Page: a buffer of one page for c-ANNS
// -----------------------------------------------------------------------------
struct Page {
	int size_;						// size for one scan
	int key_pos_;					// current pos of key_ in this leaf node
	int idx_pos_;					// current pos of id_  in this leaf node
	BLeafNode *node_;				// leaf node (level = 0)
};

// -----------------------------------------------------------------------------
//  Query-Aware Locality-Sensitive Hashing (QALSH) is designed to deal with the 
//  problem of c-Approximate Nearest Neighbor Search (c-ANNS). This is an 
//  external memory implementation. We design a new variant of B+ Tree to index 
//  the hash tables. The work was published in PVLDB 2015 as follows.
//
//  Qiang Huang, Jianlin Feng, Yikai Zhang, Qiong Fang, and Wilfred Ng. 
//  Query-aware locality-sensitive hashing for approximate nearest neighbor 
//  search, Proceedings of the VLDB Endowment (PVLDB), 9(1), pages 1–12, 2015.
// -----------------------------------------------------------------------------
template<class DType>
class QALSH {
public:
	int   n_pts_;					// number of data points
	int   dim_;						// data dimension
	int   B_;						// page size
	float p_;						// l_p distance, p \in (0,2]
	float zeta_;					// symmetric factor of p-stable distr.
	float c_;						// approximation ratio
	char  path_[300];				// index path
	const int *index_;				// data index

	float w_;						// bucket width
	int   m_;						// number of hash tables
	int   l_;						// collision threshold
	float *a_;						// query-aware lsh hash functions
	BTree **trees_;					// B+ Trees
	uint64_t dist_io_;				// io for computing distance
	uint64_t page_io_;				// io for scanning pages

	// -------------------------------------------------------------------------
	QALSH(							// constructor (build lsh index)
		int   n,						// number of data points
		int   d,						// data dimension
		int   B,						// page size
		float p,						// l_p distance, p \in (0,2]
		float zeta,						// symmetric factor of p-stable distr.
		float c,						// approximation ratio
		const DType *data,				// data points
		const char *path,				// index path
		const int *index = NULL);		// data index

	// -------------------------------------------------------------------------
	QALSH(							// constructor (load lsh index)
		const char *path,				// index path
		const int  *index = NULL);		// data index

	// -------------------------------------------------------------------------
	~QALSH();						// destructor

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	uint64_t get_memory_usage() {	// get estimated memory usage
		uint64_t ret = 0ULL;
		ret += sizeof(*this);
		ret += SIZEFLOAT * m_ * dim_; // a_
		for (int i = 0; i < m_; ++i) { // trees_
			ret += B_; // each tree only allocates B_ bytes
		}
		return ret;
	}

	// -------------------------------------------------------------------------
	uint64_t knn(					// k-NN search
		int   top_k,					// top-k value
		const DType *query,				// query point
		const char *dfolder,			// data folder
		MinK_List *list);				// k-NN results (return)

	// -------------------------------------------------------------------------
	uint64_t knn2(					// k-NN search (assis func for QALSH_PLUS)
		int   top_k,					// top-k value
		const DType *query,				// query point
		const char *dfolder,			// data folder
		MinK_List *list);				// k-NN results (return)

protected:
	// -------------------------------------------------------------------------
	inline float calc_l0_prob(float x) { return new_levy_prob(x); }

	inline float calc_l1_prob(float x) { return new_cauchy_prob(x); }
	
	inline float calc_l2_prob(float x) { return new_gaussian_prob(x); }

	// -------------------------------------------------------------------------
	int write_params();				// write parameters to disk

	// -------------------------------------------------------------------------
	int bulkload(					// build B+ Trees by bulkloading
		const DType *data);				// data points
	
	// -------------------------------------------------------------------------
	inline float calc_hash_value(int tid, const DType *data) { 
		return calc_inner_product<DType>(dim_, &a_[tid*dim_], data);
	}
	
	// -------------------------------------------------------------------------
	inline void get_tree_filename(int tid, char *fname) { // get fname of b+tree
		sprintf(fname, "%s%d.qalsh", path_, tid);
	}

	// -------------------------------------------------------------------------
	int read_params();				// read parameters from disk

	// -------------------------------------------------------------------------
	void init_search_params(		// init parameters for k-NN search
		const DType *query,				// query point
		float *q_val,					// hash values of query (return)
		Page  **lptrs,					// left  buffer (return)
		Page  **rptrs);					// right buffer (return)

	// -------------------------------------------------------------------------
	float find_radius(				// find proper radius
		const float *q_val,				// hash value of query
		const Page **lptrs,			// left  buffer
		const Page **rptrs);			// right buffer

	float update_radius(			// update radius
		float old_radius,				// old radius
		const float *q_val,				// hash value of query
		const Page **lptrs,				// left  buffer
		const Page **rptrs);			// right buffer

	// -------------------------------------------------------------------------
	void update_left_buffer(		// update left buffer
		const Page *rptr,				// right buffer
		Page *lptr);					// left  buffer (return)

	void update_right_buffer(		// update right buffer
		const Page *lptr,				// left  buffer
		Page* rptr);					// right buffer (return)

	// -------------------------------------------------------------------------
	float calc_dist(				// calc projected distance
		float q_val,					// hash value of query
		const Page *ptr);				// page buffer
	
	// -------------------------------------------------------------------------
	void delete_tree_ptr(			// delete the pointers of B+ Trees
		Page **lptrs,					// left  buffer (return)
		Page **rptrs);					// right buffer (return)
};

// -----------------------------------------------------------------------------
template<class DType>
QALSH<DType>::QALSH(				// constructor (build lsh index)
	int   n,							// number of data points
	int   d,							// dimension of space
	int   B,							// page size
	float p,							// l_p distance, p \in (0,2]
	float zeta,							// symmetric factor of p-stable distr.
	float c,							// approximation ratio
	const DType *data,					// data points
	const char *path,					// index path
	const int *index)					// data index
	: n_pts_(n), dim_(d), B_(B), p_(p), zeta_(zeta), c_(c), index_(index)
{
	dist_io_ = 0;
	page_io_ = 0;
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

	float w0 = (c_ - 1.0f) / log(sqrt(c_));
	float w1 = 2.0f * sqrt(c_);
	float w2 = sqrt((8.0f * SQR(c_) * log(c_)) / (SQR(c_) - 1.0f));
	float p1 = -1.0f, p2 = -1.0f;

	if (fabs(p_ - 0.5f) < FLOATZERO) {
		w_ = w0;
		p1 = calc_l0_prob(w_ / 2.0f);
		p2 = calc_l0_prob(w_ / (2.0f * c_));
	}
	else if (fabs(p_ - 1.0f) < FLOATZERO) {
		w_ = w1;
		p1 = calc_l1_prob(w_ / 2.0f);
		p2 = calc_l1_prob(w_ / (2.0f * c_));
	}
	else if (fabs(p_ - 2.0f) < FLOATZERO) {
		w_ = w2;
		p1 = calc_l2_prob(w_ / 2.0f);
		p2 = calc_l2_prob(w_ / (2.0f * c_));
	}
	else {
		if (fabs(p_-0.8f) < FLOATZERO) w_ = 2.503f;
		else if (fabs(p_-1.2f) < FLOATZERO) w_ = 3.151f;
		else if (fabs(p_-1.5f) < FLOATZERO) w_ = 3.465f;
		else w_ = (w2 - w1) * p_ + (2.0f * w1 - w2);

		new_stable_prob(p_, zeta_, c_, 1.0f, w_, 1000000, p1, p2);
	}

	float para1 = sqrt(log(2.0f / beta));
	float para2 = sqrt(log(1.0f / delta));
	float para3 = 2.0f * (p1 - p2) * (p1 - p2);
	float eta   = para1 / para2;
	float alpha = (eta * p1 + p2) / (1.0f + eta);

	m_ = (int) ceil((para1 + para2) * (para1 + para2) / para3);
	l_ = (int) ceil(alpha * m_);

	// generate hash functions
	a_ = new float[m_*dim_];
	for (int i = 0; i < m_*dim_; ++i) {
		if (fabs(p_-0.5f) < FLOATZERO) a_[i] = levy(1.0f, 0.0f);
		else if (fabs(p_-1.0f) < FLOATZERO) a_[i] = cauchy(1.0f, 0.0f);
		else if (fabs(p_-2.0f) < FLOATZERO) a_[i] = gaussian(0.0f, 1.0f);
		else a_[i] = p_stable(p_, zeta_, 1.0f, 0.0f);
	}

	// write parameters to disk
	if (write_params()) exit(1);

	//  bulkloading
	if (bulkload(data)) exit(1);
}

// -----------------------------------------------------------------------------
template<class DType>
int QALSH<DType>::write_params()	// write parameters to disk
{
	char fname[200]; sprintf(fname, "%spara", path_);
	FILE *fp = fopen(fname, "rb");
	if (fp) { printf("Hash Tables Already Exist\n\n"); exit(1); }

	fp = fopen(fname, "wb");
	if (!fp) {
		printf("Could not create %s\n", fname);
		printf("Perhaps no such folder %s?\n", path_);
		return 1;
	}

	fwrite(&n_pts_, SIZEINT,   1, fp);
	fwrite(&dim_,   SIZEINT,   1, fp);
	fwrite(&B_,     SIZEINT,   1, fp);
	fwrite(&m_,     SIZEINT,   1, fp);
	fwrite(&l_,     SIZEINT,   1, fp);
	fwrite(&p_,     SIZEFLOAT, 1, fp);
	fwrite(&zeta_,  SIZEFLOAT, 1, fp);
	fwrite(&c_,     SIZEFLOAT, 1, fp);
	fwrite(&w_,     SIZEFLOAT, 1, fp);
	
	fwrite(a_, SIZEFLOAT, m_*dim_, fp);
	fclose(fp);
	return 0;
}

// -----------------------------------------------------------------------------
template<class DType>
int QALSH<DType>::bulkload(			// build B+ Trees by bulkloading
	const DType *data)					// data set
{
	Result *table = new Result[n_pts_];
	trees_ = new BTree*[m_];
	for (int i = 0; i < m_; ++i) {
		// calc hash value and init the hash table
		for (int j = 0; j < n_pts_; ++j) {
			table[j].id_  = j;
			table[j].key_ = calc_hash_value(i, &data[(uint64_t)j*dim_]);
		}
		qsort(table, n_pts_, sizeof(Result), ResultComp);
		
		// use B+ tree to index the hash table
		char fname[200]; get_tree_filename(i, fname);
		trees_[i] = new BTree();
		trees_[i]->init(B_, fname);

		if (trees_[i]->bulkload(n_pts_, table)) return 1;
	}
	delete[] table;
	return 0;
}

// -----------------------------------------------------------------------------
template<class DType>
QALSH<DType>::~QALSH()				// destructor
{
	for (int i = 0; i < m_; ++i) {
		delete trees_[i]; trees_[i] = NULL;
	}
	delete[] trees_;
	delete[] a_;
}

// -----------------------------------------------------------------------------
template<class DType>
QALSH<DType>::QALSH(				// constructor (load lsh index)
	const char *path,					// index path
	const int  *index)					// data index
	: index_(index)
{
	dist_io_ = 0;
	page_io_ = 0;
	strcpy(path_, path);

	// read parameters from disk
	if (read_params()) exit(1);

	// init b+ trees for k-NN search
	trees_ = new BTree*[m_];
	for (int i = 0; i < m_; ++i) {
		char fname[200]; get_tree_filename(i, fname);
		trees_[i] = new BTree();
		trees_[i]->init_restore(fname);
	}
}

// -----------------------------------------------------------------------------
template<class DType>
int QALSH<DType>::read_params()		// read parameters from disk
{
	char fname[200]; sprintf(fname, "%spara", path_);
	FILE *fp = fopen(fname, "rb");
	if (!fp) { printf("Could not open %s\n", fname); return 1; }

	fread(&n_pts_, SIZEINT,   1, fp);
	fread(&dim_,   SIZEINT,   1, fp);
	fread(&B_,     SIZEINT,   1, fp);
	fread(&m_,     SIZEINT,   1, fp);
	fread(&l_,     SIZEINT,   1, fp);
	fread(&p_,     SIZEFLOAT, 1, fp);
	fread(&zeta_,  SIZEFLOAT, 1, fp);
	fread(&c_,     SIZEFLOAT, 1, fp);
	fread(&w_,     SIZEFLOAT, 1, fp);
	
	a_ = new float[m_*dim_];
	fread(a_, SIZEFLOAT, m_*dim_, fp);
	fclose(fp);
	return 0;
}

// -----------------------------------------------------------------------------
template<class DType>
void QALSH<DType>::display()		// display parameters
{
	printf("Parameters of QALSH:\n");
	printf("n    = %d\n",   n_pts_);
	printf("d    = %d\n",   dim_);
	printf("B    = %d\n",   B_);
	printf("p    = %.1f\n", p_);
	printf("zeta = %.1f\n", zeta_);
	printf("c    = %.1f\n", c_);
	printf("w    = %f\n",   w_);
	printf("m    = %d\n",   m_);
	printf("l    = %d\n",   l_);
	printf("path = %s\n\n", path_);
}

// -----------------------------------------------------------------------------
template<class DType>
uint64_t QALSH<DType>::knn(			// k-NN search
	int   top_k,						// top-k value
	const DType *query,					// query point
	const char *dfolder,				// data folder
	MinK_List *list)					// k-NN results (return)
{
	list->reset();

	// initialize parameters for c-k-ANNS
	int  *freq    = new int[n_pts_]; memset(freq, 0, n_pts_*SIZEFLOAT);
	bool *checked = new bool[n_pts_]; memset(checked, false, n_pts_*SIZEBOOL);
	bool *flag    = new bool[m_]; memset(flag, true, m_ * SIZEBOOL);
	
	DType *data  = new DType[dim_];
	float *q_val = new float[m_];
	Page **lptrs = new Page*[m_];
	Page **rptrs = new Page*[m_];
	for (int i = 0; i < m_; ++i) {
		lptrs[i] = new Page();
		rptrs[i] = new Page();
	}
	init_search_params(query, q_val, lptrs, rptrs);

	// c-k-ANNS via dynamic collision counting framework
	int   candidates = CANDIDATES + top_k - 1; // candidates size
	float kdist  = MAXREAL;
	float radius = find_radius(q_val, (const Page**)lptrs, (const Page**)rptrs);
	float bucket = w_ * radius / 2.0f;

	while (true) {
		// step 1: initialize the stop condition for current round
		int num_flag = 0;
		memset(flag, true, m_*SIZEBOOL);

		// step 2: (R,c)-NN search (find frequent data points)
		while (num_flag < m_) {
			for (int i = 0; i < m_; ++i) {
				if (!flag[i]) continue;

				// step 2.1: compute <ldist> and <rdist>
				Page *lptr = lptrs[i];
				Page *rptr = rptrs[i];

				float dist, ldist = MAXREAL, rdist = MAXREAL;
				if (lptr->size_ != -1) ldist = calc_dist(q_val[i], lptr);
				if (rptr->size_ != -1) rdist = calc_dist(q_val[i], rptr);

				// step 2.2: determine the closer direction (left or right)
				// and do collision counting to find frequent points
				if (ldist < bucket && ldist <= rdist) {
					int count = lptr->size_;
					int end   = lptr->idx_pos_;
					int start = end - count;

					for (int j = end; j > start; --j) {
						int id = lptr->node_->get_entry_id(j);
						if (++freq[id] > l_ && !checked[id]) {
							checked[id] = true;
							read_data_new_format<DType>(id, dim_, B_, dfolder, data);
							dist  = calc_lp_dist<DType>(dim_, p_, kdist, data, query);
							kdist = list->insert(dist, id);
							if (++dist_io_ >= candidates) break;
						}
					}
					update_left_buffer(rptr, lptr);
				}
				else if (rdist < bucket && ldist > rdist) {
					int count = rptr->size_;
					int start = rptr->idx_pos_;
					int end   = start + count;

					for (int j = start; j < end; ++j) {
						int id = rptr->node_->get_entry_id(j);
						if (++freq[id] > l_ && !checked[id]) {
							checked[id] = true;
							read_data_new_format<DType>(id, dim_, B_, dfolder, data);
							dist  = calc_lp_dist<DType>(dim_, p_, kdist, data, query);
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
		// step 3: stop conditions 1 & 2
		if (kdist < c_*radius && dist_io_ >= top_k) break;
		if (dist_io_ >= candidates) break;

		// step 4: auto-update <radius>
		radius = update_radius(radius, q_val, (const Page**) lptrs, 
			(const Page**) rptrs);
		bucket = radius * w_ / 2.0f;
	}
	// release space
	delete_tree_ptr(lptrs, rptrs);
	delete[] freq;
	delete[] checked;
	delete[] flag;
	delete[] q_val;
	delete[] data;

	return page_io_ + dist_io_;
}

// -----------------------------------------------------------------------------
template<class DType>
uint64_t QALSH<DType>::knn2(		// k-NN search
	int   top_k,						// top-k value
	const DType *query,					// query point
	const char *dfolder,				// data folder
	MinK_List *list)					// k-NN results (return)
{
	// initialize parameters for c-k-ANNS
	int  *freq = new int[n_pts_]; memset(freq, 0, n_pts_*SIZEFLOAT);
	bool *checked = new bool[n_pts_]; memset(checked, false, n_pts_*SIZEBOOL);
	bool *bucket_flag = new bool[m_]; memset(bucket_flag, true, m_*SIZEBOOL);
	bool *range_flag = new bool[m_]; memset(range_flag, true, m_*SIZEBOOL);
	
	DType *data  = new DType[dim_];
	float *q_val = new float[m_];
	Page **lptrs = new Page*[m_];
	Page **rptrs = new Page*[m_];
	for (int i = 0; i < m_; ++i) {
		lptrs[i] = new Page();
		rptrs[i] = new Page();
	}
	init_search_params(query, q_val, lptrs, rptrs);

	// c-k-ANNS via dynamic collision counting framework
	int candidates = CANDIDATES+top_k-1; // candidates size
	int num_range  = 0;				// used for search range bound
	
	float kdist  = list->max_key();
	float radius = find_radius(q_val, (const Page**)lptrs, (const Page**)rptrs);
	float bucket = w_ * radius / 2.0f;
	float range  = kdist > MAXREAL-1.0f ? MAXREAL : kdist*w_/2.0f;

	while (true) {
		// step 1: initialize the stop condition for current round
		int num_bucket = 0;
		memset(bucket_flag, true, m_*SIZEBOOL);

		// step 2: (R,c)-NN search (find frequent data points)
		while (num_bucket < m_ && num_range < m_) {
			for (int i = 0; i < m_; ++i) {
				if (!bucket_flag[i]) continue;

				// step 2.1: compute <ldist> and <rdist>
				Page *lptr = lptrs[i];
				Page *rptr = rptrs[i];

				float dist, ldist = MAXREAL, rdist = MAXREAL;
				if (lptr->size_ != -1) ldist = calc_dist(q_val[i], lptr);
				if (rptr->size_ != -1) rdist = calc_dist(q_val[i], rptr);

				// step 2.2: determine the closer direction (left or right)
				// and do collision counting to find frequent points
				if (ldist < bucket && ldist < range && ldist <= rdist) {
					int count = lptr->size_;
					int end   = lptr->idx_pos_;
					int start = end - count;

					for (int j = end; j > start; --j) {
						int id = lptr->node_->get_entry_id(j);
						if (++freq[id] > l_ && !checked[id]) {
							checked[id] = true;
							int oid = index_[id];
							read_data_new_format<DType>(oid, dim_, B_, dfolder, data);
							dist  = calc_lp_dist<DType>(dim_, p_, kdist, data, query);
							kdist = list->insert(dist, oid);
							if (++dist_io_ >= candidates) break;
						}
					}
					update_left_buffer(rptr, lptr);
				}
				else if (rdist < bucket && rdist < range && ldist > rdist) {
					int count = rptr->size_;
					int start = rptr->idx_pos_;
					int end   = start + count;

					for (int j = start; j < end; ++j) {
						int id = rptr->node_->get_entry_id(j);
						if (++freq[id] > l_ && !checked[id]) {
							checked[id] = true;
							int oid = index_[id];
							read_data_new_format<DType>(oid, dim_, B_, dfolder, data);
							dist  = calc_lp_dist<DType>(dim_, p_, kdist, data, query);
							kdist = list->insert(dist, oid);
							if (++dist_io_ >= candidates) break;
						}
					}
					update_right_buffer(lptr, rptr);
				}
				else {
					bucket_flag[i] = false;
					++num_bucket;
					if (ldist >= range && rdist >= range && range_flag[i]) {
						range_flag[i] = false;
						++num_range;
					}
				}
				if (num_bucket >= m_ || num_range >= m_) break;
				if (dist_io_ >= candidates) break;
			}
			if (num_bucket >= m_ || num_range >= m_) break;
			if (dist_io_ >= candidates) break;
		}
		// step 3: stop conditions 1 & 2
		if (dist_io_ >= candidates || num_range >= m_) break;

		// step 4: auto-update <radius>
		radius = update_radius(radius, q_val, (const Page**) lptrs, 
			(const Page**) rptrs);
		bucket = radius * w_ / 2.0f;
	}
	// release space
	delete_tree_ptr(lptrs, rptrs);
	delete[] freq;
	delete[] checked;
	delete[] bucket_flag;
	delete[] range_flag;
	delete[] q_val;
	delete[] data;

	return page_io_ + dist_io_;
}

// -----------------------------------------------------------------------------
template<class DType>
void QALSH<DType>::init_search_params(// init parameters for k-NN search
	const DType *query,					// query point
	float *q_val,						// hash values of query (return)
	Page  **lptrs,						// left buffer (return)
	Page  **rptrs)						// right buffer (return)
{
	page_io_ = 0;
	dist_io_ = 0;

	for (int i = 0; i < m_; ++i) {
		lptrs[i]->node_    = NULL;
		lptrs[i]->key_pos_ = -1;
		lptrs[i]->idx_pos_ = -1;
		lptrs[i]->size_    = -1;

		rptrs[i]->node_    = NULL;
		rptrs[i]->key_pos_ = -1;
		rptrs[i]->idx_pos_ = -1;
		rptrs[i]->size_    = -1;
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
			// -----------------------------------------------------------------
			//  at least two levels in the B+ Tree: index node and lead node
			// -----------------------------------------------------------------
			index_node = new BIndexNode();
			index_node->init_restore(tree, block);
			++page_io_;

			// -----------------------------------------------------------------
			//  find the leaf node whose value is closest and larger than the 
			//  key of query q
			// -----------------------------------------------------------------
			lescape = false;		// locate the position of branch
			while (index_node->get_level() > 1) {
				follow = index_node->find_position_by_key(q_v);

				if (follow == -1) { // scan the most left branch
					if (lescape) { 
						follow = 0;
					} else {
						if (block != tree->root_) { 
							printf("No branch found\n"); exit(1);
						} else {
							follow = 0; lescape = true;
						}
					}
				}
				block = index_node->get_son(follow);
				delete index_node; index_node = NULL;

				index_node = new BIndexNode();
				index_node->init_restore(tree, block);
				++page_io_;			// access a new node (a new page)
			}

			// -----------------------------------------------------------------
			//  after finding the leaf node whose value is closest to the key of
			//  query, initialize <lptrs[i]> and <rptrs[i]>.
			//
			//  <lescape> = true is that the query has no <lptrs>, the query is 
			//  the smallest value.
			// -----------------------------------------------------------------
			follow = index_node->find_position_by_key(q_v);
			if (follow < 0) {
				lescape = true; follow = 0;
			}

			if (lescape) {
				// -------------------------------------------------------------
				//  only init right buffer
				// -------------------------------------------------------------
				block = index_node->get_son(0);
				rptr->node_ = new BLeafNode();
				rptr->node_->init_restore(tree, block);
				rptr->key_pos_ = 0;
				rptr->idx_pos_ = 0;

				increment = rptr->node_->get_increment();
				num_entries = rptr->node_->get_num_entries();
				if (increment > num_entries) rptr->size_ = num_entries;
				else rptr->size_ = increment;

				++page_io_;
			}
			else {					
				// -------------------------------------------------------------
				//  init left buffer
				// -------------------------------------------------------------
				block = index_node->get_son(follow);
				lptr->node_ = new BLeafNode();
				lptr->node_->init_restore(tree, block);

				pos = lptr->node_->find_position_by_key(q_v);
				if (pos < 0) pos = 0;
				lptr->key_pos_ = pos;

				increment = lptr->node_->get_increment();
				if (pos == lptr->node_->get_num_keys()-1) {
					num_entries = lptr->node_->get_num_entries();

					lptr->idx_pos_ = num_entries - 1;
					lptr->size_ = num_entries - pos*increment;
				}
				else {
					lptr->idx_pos_ = pos*increment + increment - 1;
					lptr->size_ = increment;
				}
				++page_io_;

				// -------------------------------------------------------------
				//  init right buffer
				// -------------------------------------------------------------
				if (pos < lptr->node_->get_num_keys() - 1) {
					rptr->node_    = lptr->node_;
					rptr->key_pos_ = pos + 1;
					rptr->idx_pos_ = (pos+1) * increment;

					if ((pos+1) == rptr->node_->get_num_keys()-1) {
						num_entries = rptr->node_->get_num_entries();
						rptr->size_ = num_entries - (pos+1)*increment;
					} else {
						rptr->size_ = increment;
					}
				}
				else {
					rptr->node_ = lptr->node_->get_right_sibling();
					if (rptr->node_) {
						rptr->key_pos_ = 0;
						rptr->idx_pos_ = 0;

						increment = rptr->node_->get_increment();
						num_entries = rptr->node_->get_num_entries();
						if (increment > num_entries) rptr->size_ = num_entries;
						else rptr->size_ = increment;
						
						++page_io_;
					}
				}
			}
		}
		else {
			// -----------------------------------------------------------------
			//  only one level in the B+ Tree: one lead node
			//  (1) init left buffer
			// -----------------------------------------------------------------
			lptr->node_ = new BLeafNode();
			lptr->node_->init_restore(tree, block);

			pos = lptr->node_->find_position_by_key(q_v);
			if (pos < 0) pos = 0;
			lptr->key_pos_ = pos;

			increment = lptr->node_->get_increment();
			if (pos == lptr->node_->get_num_keys() - 1) {
				num_entries = lptr->node_->get_num_entries();

				lptr->idx_pos_ = num_entries - 1;
				lptr->size_ = num_entries - pos*increment;
			}
			else {
				lptr->idx_pos_ = pos*increment + increment - 1;
				lptr->size_ = increment;
			}
			++page_io_;
			
			// -----------------------------------------------------------------
			//  (2) init right buffer
			// -----------------------------------------------------------------
			if (pos < lptr->node_->get_num_keys() - 1) {
				rptr->node_    = lptr->node_;
				rptr->key_pos_ = pos + 1;
				rptr->idx_pos_ = (pos+1) * increment;

				if ((pos+1) == rptr->node_->get_num_keys()-1) {
					num_entries = rptr->node_->get_num_entries();
					rptr->size_ = num_entries - (pos+1)*increment;
				} else {
					rptr->size_ = increment;
				}
			}
			else {
				rptr->node_    = NULL;
				rptr->key_pos_ = -1;
				rptr->idx_pos_ = -1;
				rptr->size_    = -1;
			}
		}
		if (index_node != NULL) { delete index_node; index_node = NULL; }
	}
}

// -----------------------------------------------------------------------------
template<class DType>
float QALSH<DType>::find_radius(	// find proper radius
	const float *q_val,					// hash value of query
	const Page **lptrs,					// left buffer
	const Page **rptrs)					// right buffer
{
	float radius = update_radius(1.0f / c_, q_val, lptrs, rptrs);
	if (radius < 1.0f) radius = 1.0f;

	return radius;
}

// -----------------------------------------------------------------------------
template<class DType>
float QALSH<DType>::update_radius(	// update radius
	float old_radius,					// old radius
	const float *q_val,					// hash value of query
	const Page **lptrs,					// left buffer
	const Page **rptrs)					// right buffer
{
	// calc the projected distance which is closest to q among m hash tables 
	std::vector<float> list;
	for (int i = 0; i < m_; ++i) {
		if (lptrs[i]->size_ != -1) list.push_back(calc_dist(q_val[i], lptrs[i]));
		if (rptrs[i]->size_ != -1) list.push_back(calc_dist(q_val[i], rptrs[i]));
	}
	std::sort(list.begin(), list.end());

	// find the median distance and return the new radius
	int num = (int) list.size();
	if (num == 0) return c_ * old_radius;

	float dist = -1.0f;
	if (num % 2 == 0) dist = (list[num/2-1] + list[num/2]) / 2.0f;
	else dist = list[num / 2];
	
	int kappa = (int) ceil(log(2.0f*dist/w_) / log(c_));
	list.clear(); list.shrink_to_fit();

	return pow(c_, kappa);
}

// -----------------------------------------------------------------------------
template<class DType>
void QALSH<DType>::update_left_buffer(// update left buffer
	const Page *rptr,					// right buffer
	Page *lptr)							// left buffer (return)
{
	BLeafNode *leaf_node = NULL;
	BLeafNode *old_leaf_node = NULL;

	if (lptr->key_pos_ > 0) {
		lptr->key_pos_--;

		int pos        = lptr->key_pos_;
		int increment  = lptr->node_->get_increment();
		lptr->idx_pos_ = pos*increment + increment - 1;
		lptr->size_    = increment;
	}
	else {
		old_leaf_node = lptr->node_;
		leaf_node = lptr->node_->get_left_sibling();

		if (leaf_node) {
			lptr->node_     = leaf_node;
			lptr->key_pos_  = lptr->node_->get_num_keys() - 1;

			int pos         = lptr->key_pos_;
			int increment   = lptr->node_->get_increment();
			int num_entries = lptr->node_->get_num_entries();
			lptr->idx_pos_  = num_entries - 1;
			lptr->size_     = num_entries - pos*increment;
			++page_io_;
		}
		else {
			lptr->node_    = NULL;
			lptr->key_pos_ = -1;
			lptr->idx_pos_ = -1;
			lptr->size_    = -1;
		}

		if (rptr->node_ != old_leaf_node) {
			delete old_leaf_node; old_leaf_node = NULL;
		}
	}
}

// -----------------------------------------------------------------------------
template<class DType>
void QALSH<DType>::update_right_buffer(// update right buffer
	const Page *lptr,					// left buffer
	Page *rptr)							// right buffer (return)
{
	BLeafNode *leaf_node = NULL;
	BLeafNode *old_leaf_node = NULL;

	if (rptr->key_pos_ < rptr->node_->get_num_keys()-1) {
		rptr->key_pos_++;

		int pos       = rptr->key_pos_;
		int increment = rptr->node_->get_increment();
		
		rptr->idx_pos_ = pos * increment;
		if (pos == rptr->node_->get_num_keys()-1) {
			int num_entries = rptr->node_->get_num_entries();
			rptr->size_ = num_entries - pos*increment;
		} else {
			rptr->size_ = increment;
		}
	}
	else {
		old_leaf_node = rptr->node_;
		leaf_node = rptr->node_->get_right_sibling();

		if (leaf_node) {
			rptr->node_    = leaf_node;
			rptr->key_pos_ = 0;
			rptr->idx_pos_ = 0;

			int increment   = rptr->node_->get_increment();
			int num_entries = rptr->node_->get_num_entries();
			if (increment > num_entries) rptr->size_ = num_entries;
			else rptr->size_ = increment;

			++page_io_;
		}
		else {
			rptr->node_    = NULL;
			rptr->key_pos_ = -1;
			rptr->idx_pos_ = -1;
			rptr->size_    = -1;
		}

		if (lptr->node_ != old_leaf_node) {
			delete old_leaf_node; old_leaf_node = NULL;
		}
	}
}

// -----------------------------------------------------------------------------
template<class DType>
inline float QALSH<DType>::calc_dist(// calc projected distance
	float q_val,						// hash value of query
	const Page *ptr)					// page buffer
{
	int   pos = ptr->key_pos_;
	float key = ptr->node_->get_key(pos);

	return fabs(key - q_val);
}

// -----------------------------------------------------------------------------
template<class DType>
void QALSH<DType>::delete_tree_ptr(	// delete the pointers of B+ Trees
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
		if (lptrs[i]->node_ && lptrs[i]->node_ != rptrs[i]->node_) {
			delete lptrs[i]->node_; lptrs[i]->node_ = NULL;
		}
		if (rptrs[i]->node_) {
			delete rptrs[i]->node_; rptrs[i]->node_ = NULL;
		}

		delete lptrs[i]; lptrs[i] = NULL;
		delete rptrs[i]; rptrs[i] = NULL;
	}
	delete[] lptrs; lptrs = NULL;
	delete[] rptrs; rptrs = NULL;
}

} // end namespace nns
