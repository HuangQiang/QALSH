#pragma once

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <unordered_map>
#include <vector>

#include "def.h"
#include "util.h"
#include "pri_queue.h"
#include "kd_tree.h"
#include "qalsh.h"

namespace nns {

// -----------------------------------------------------------------------------
template<class DType>
class QALSH_PLUS {
public:
	QALSH_PLUS(						// constructor (build index)
		int   n,						// number of data points
		int   d,						// dimensionality
		int   B,						// page size
		int   leaf,						// leaf size of kd-tree
		int   L,						// number of projection (drusilla)
		int   M,						// number of candidates (drusilla)
		float p,						// l_p distance
		float zeta,						// a parameter of p-stable distr.
		float c,						// approximation ratio
		const DType *data,				// data points
		const char *path);				// index path

	// -------------------------------------------------------------------------
	QALSH_PLUS(						// constructor (load index)
		const char *path);				// index path

	// -------------------------------------------------------------------------
	~QALSH_PLUS();					// destructor

	// -------------------------------------------------------------------------
	inline int get_num_blocks() { return n_blocks_; }

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	uint64_t get_memory_usage() {	// get estimated memory usage
		uint64_t ret = 0ULL;
		ret += sizeof(*this);
		ret += SIZEINT * n_blocks_; // block_size_
		ret += SIZEINT * n_pts_;	// index_
		ret += SIZEINT * n_blocks_*n_samples_; // sample_index_
		ret += SIZEINT * n_pts_;	// sample_index_to_block_
		ret += lsh_->get_memory_usage(); // first level lsh
		for (int i = 0; i < n_blocks_; ++i) { // second level lsh
			ret += blocks_[i]->get_memory_usage();
		}
		return ret;
	}

	// -------------------------------------------------------------------------
	uint64_t knn(					// k-NN search
		int   top_k,					// top-k value
		int   nb,						// number of blocks for search
		const DType *query,				// query point
		const char *dfolder,			// data folder
		MinK_List *list);				// top-k results (return)

protected:
	int  n_pts_;					// number of data points
	int  dim_;						// data dimension
	int  n_samples_;				// number of samples for drusilla-select
	char path_[200];				// index path

	int  n_blocks_;					// number of blocks 
	int  *block_size_;				// an array the block size of n_blocks_
	int  *index_;					// data index after kd-tree partition
	int  *sample_index_;			// sample data index
	int  *sample_index_to_block_;	// sample data id to block
	QALSH<DType> *lsh_;				// first level lsh index for sample data
	std::vector<QALSH<DType>*> blocks_; // second level lsh index for blocks

	// -------------------------------------------------------------------------
	void kd_tree_partition(			// kd-tree partition
		int   leaf,						// leaf size of kd-tree
		const DType *data);				// data points

	// -------------------------------------------------------------------------
	void copy(						// copy original data to destination data
		const DType *orig_data, 		// original data
		DType *dest_data);				// destination data (return)

	// -------------------------------------------------------------------------
	void drusilla_select(			// drusilla select for a block
		int   n,						// number of data points in this block
		int   L,						// #projection for drusilla-select
		int   M,						// #candidates for drusilla-select
		const int *index,				// data index for this block
		const DType *data,				// data points in this block
		int   *sample_index,			// sample data index (return)
		DType *sample_data);			// sample data (return)

	// -------------------------------------------------------------------------
	void calc_shift_data(			// calculate shift data points
		int   n,						// number of data points in this block
		const DType *data,				// data points in this block
		int   &max_id,					// locl id with max l2-norm (return)
		float &max_norm,				// max l2-norm (return)
		float *norm,					// l2-norm of shift data (return)
		float *shift_data); 			// shift data (return)

	// -------------------------------------------------------------------------
	void shift(						// shift the original data by centroid
		const DType *data,				// original data point
		const float *centroid,			// centroid
		float &norm,					// l2-norm of shifted data (return)
		float *shift_data);				// shifted data (return)

	// -------------------------------------------------------------------------
	void select_proj(				// select project vector
		float norm,						// max l2-norm
		const float *shift_data,		// shift data with max l2-norm
		float *proj);					// projection vector

	// -------------------------------------------------------------------------
	float calc_distortion(			// calc distortion
		float offset,					// offset
		const float *proj,				// projection vector
		const float *shift_data);		// input shift data

	// -------------------------------------------------------------------------
	int write_params();				// write parameters

	// -------------------------------------------------------------------------
	int read_params();				// read parameters

	// -------------------------------------------------------------------------
	uint64_t get_block_order(		// get block order
		int nb,							// number of blocks for search
		const DType *query,				// query point
		const char *dfolder,			// data folder
		std::vector<int> &block_order);	// block order (return)
};

// -----------------------------------------------------------------------------
template<class DType>
QALSH_PLUS<DType>::QALSH_PLUS(		// constructor (build index)
	int   n,							// number of data points
	int   d,							// dimensionality
	int   B,							// page size
	int   leaf,							// leaf size of kd-tree
	int   L,							// number of projection (drusilla)
	int   M,							// number of candidates (drusilla)
	float p,							// l_p distance
	float zeta,							// a parameter of p-stable distr.
	float c,							// approximation ratio
	const DType *data,					// data points
	const char *path)					// index path
	: n_pts_(n), dim_(d), n_samples_(L*M)
{
	strcpy(path_, path);
	create_dir(path_);

	index_ = new int[n_pts_];
	sample_index_to_block_ = new int[n_pts_];
	memset(sample_index_to_block_, -1, n_pts_);

	// kd-tree partition (get index_, n_blocks_, and block_size_)
	kd_tree_partition(leaf, data);

	// init sample_index_ and build qalsh for each block
	int n_sample_pts = n_blocks_*n_samples_;
	DType *sample_data = new DType[(uint64_t) n_sample_pts*dim_];
	sample_index_ = new int[n_sample_pts];

	int start = 0;
	int count = 0;
	for (int i = 0; i < n_blocks_; ++i) {
		int n_blk = block_size_[i];
		const int *index  = (const int*) &index_[start];
		int *sample_index = &sample_index_[count];

		// get the block data from index
		DType *blk_data = new DType[(uint64_t)n_blk*dim_];
		for (int j = 0; j < n_blk; ++j) {
			copy(&data[(uint64_t)index[j]*dim_], &blk_data[(uint64_t)j*dim_]);
		}

		// get sample data index (representative data) by drusilla select
		assert(n_blk > n_samples_);
		drusilla_select(n_blk, L, M, index, (const DType*) blk_data, 
			sample_index, &sample_data[(uint64_t)count*dim_]);

		for (int j = 0; j < n_samples_; ++j) {
			sample_index_to_block_[sample_index[j]] = i;
		}
		
		// build qalsh for each blcok 
		char block_path[200]; sprintf(block_path, "%s%d/", path_, i);
		create_dir(block_path);

		QALSH<DType> *lsh = new QALSH<DType>(n_blk, dim_, B, p, zeta, c, 
			(const DType*) blk_data, block_path, index);
		blocks_.push_back(lsh);
		delete[] blk_data;

		// update parameters
		start += n_blk;
		count += n_samples_;
	}
	assert(start == n_pts_ && count == n_sample_pts);

	// build qalsh for sample data
	char sample_path[200]; sprintf(sample_path, "%ssample/", path_);
	create_dir(sample_path);

	lsh_ = new QALSH<DType>(n_sample_pts, dim_, B, p, zeta, c,
		(const DType*) sample_data, sample_path, (const int*) sample_index_);

	// write parameters to disk
	if (write_params()) exit(1);
	delete[] sample_data;
}

// -----------------------------------------------------------------------------
template<class DType>
void QALSH_PLUS<DType>::kd_tree_partition(// kd-tree partition
	int   leaf,							// leaf size of kd-tree
	const DType *data)					// data points
{
	// build a kd-tree for input data with specific leaf size
	KD_Tree<DType>* tree = new KD_Tree<DType>(n_pts_, dim_, leaf, data);

	// init index_
	std::vector<int> block_size;
	tree->traversal(block_size, index_);
	
	// init n_blocks_, and block_size_
	n_blocks_ = (int) block_size.size(); assert(n_blocks_ > 0);
	block_size_ = new int[n_blocks_];
	for (int i = 0; i < n_blocks_; ++i) {
		block_size_[i] = block_size[i];
	}
	// release space
	delete tree;
	block_size.clear(); block_size.shrink_to_fit();
}

// -----------------------------------------------------------------------------
template<class DType>
void QALSH_PLUS<DType>::copy(		// copy original data to destination data
	const DType *orig_data, 			// original data
	DType *dest_data)					// destination data (return)
{
	for (int i = 0; i < dim_; ++i) {
		dest_data[i] = orig_data[i];
	}
}

// -----------------------------------------------------------------------------
template<class DType>
void QALSH_PLUS<DType>::drusilla_select(// drusilla select
	int   n,							// number of data points in this block
	int   L,							// #projection for drusilla-select
	int   M,							// #candidates for drusilla-select
	const int *index,					// data index for this block
	const DType *data,					// data points in this block
	int   *sample_index,				// sample data index (return)
	DType *sample_data)					// sample data (return)
{
	// calc shift data
	int   max_id      = -1;
	float max_norm    = MINREAL;
	float *norm       = new float[n];
	float *shift_data = new float[(uint64_t)n*dim_];
	
	calc_shift_data(n, data, max_id, max_norm, norm, shift_data);

	// drusilla select
	float  *proj        = new float[dim_];
	Result *score       = new Result[n];
	bool   *close_angle = new bool[n];
	float  offset       = -1.0f;
	float  distortion   = -1.0f;

	for (int i = 0; i < L; ++i) {
		// select the projection vector with largest norm and normalize it
		select_proj(norm[max_id], &shift_data[(uint64_t)max_id*dim_], proj);

		// calculate offsets and distortions
		for (int j = 0; j < n; ++j) {
			close_angle[j] = false;
			score[j].id_   = j;

			if (norm[j] > 0.0f) {
				const float *tmp = &shift_data[(uint64_t)j*dim_];
				offset = calc_inner_product<float>(dim_, (const float*)proj, tmp);
				distortion = calc_distortion(offset, (const float*)proj, tmp);
				score[j].key_ = offset*offset - distortion;

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
		// collect the points that are well-represented by this projection
		qsort(score, n, sizeof(Result), ResultCompDesc);
		for (int j = 0; j < M; ++j) {
			int id  = score[j].id_;
			int loc = i * M + j;

			sample_index[loc] = index[id];
			copy(&data[(uint64_t)id*dim_], &sample_data[(uint64_t)loc*dim_]);
			norm[id] = -1.0f;
		}
		//  find the next largest norm and the corresponding point
		max_id = -1; max_norm = MINREAL;
		for (int j = 0; j < n; ++j) {
			if (norm[j] > 0.0f && close_angle[j]) { norm[j] = 0.0f; }
			if (norm[j] > max_norm) { max_norm = norm[j]; max_id = j; }
		}
	}
	// release space
	delete[] close_angle;
	delete[] score;
	delete[] proj;
	delete[] norm;
	delete[] shift_data;
}

// -----------------------------------------------------------------------------
template<class DType>
void QALSH_PLUS<DType>::calc_shift_data(// calculate shift data points
	int   n,							// number of data points in this block
	const DType *data,					// data points in this block
	int   &max_id,						// locl id with max l2-norm (return)
	float &max_norm,					// max l2-norm (return)
	float *norm,						// l2-norm of shift data (return)
	float *shift_data)					// shift data (return)
{
	// calculate the centroid of data points
	float *centroid = new float[dim_]; memset(centroid, 0.0f, dim_*SIZEFLOAT);
	for (int i = 0; i < n; ++i) {
		const DType *tmp = &data[(uint64_t) i*dim_];
		for (int j = 0; j < dim_; ++j) {
			centroid[j] += (float) tmp[j];
		}
	}
	for (int i = 0; i < dim_; ++i) centroid[i] /= n;

	// make a copy of data points which move to the centroid of data points
	max_id = -1; max_norm = MINREAL;
	for (int i = 0; i < n; ++i) {
		uint64_t loc = (uint64_t) i * dim_;
		shift(&data[loc], centroid, norm[i], &shift_data[loc]);
		if (norm[i] > max_norm) { max_norm = norm[i]; max_id = i; }
	}
	delete[] centroid;
}

// -----------------------------------------------------------------------------
template<class DType>
void QALSH_PLUS<DType>::shift(		// shift the original data by centroid
	const DType *data,					// original data point
	const float *centroid,				// centroid
	float &norm,						// l2-norm of shifted data (return)
	float *shift_data)					// shifted data (return)
{
	norm = 0.0f;
	for (int j = 0; j < dim_; ++j) {
		float tmp = (float) data[j] - centroid[j];
		shift_data[j] = tmp; norm += tmp*tmp;
	}
	norm = sqrt(norm);
}

// -----------------------------------------------------------------------------
template<class DType>
void QALSH_PLUS<DType>::select_proj(// select project vector
	float norm,							// max l2-norm
	const float *shift_data,			// shift data with max l2-norm
	float *proj)						// projection vector
{
	for (int j = 0; j < dim_; ++j) {
		proj[j] = shift_data[j] / norm;
	}
}

// -----------------------------------------------------------------------------
template<class DType>
float QALSH_PLUS<DType>::calc_distortion(// calc distortion
	float offset,							// offset
	const float *proj,						// projection vector
	const float *shift_data)				// input shift data
{
	float distortion = 0.0f;
	for (int j = 0; j < dim_; ++j) {
		float tmp = shift_data[j] - offset * proj[j];
		distortion += tmp * tmp;
	}
	return distortion;
}

// -----------------------------------------------------------------------------
template<class DType>
int QALSH_PLUS<DType>::write_params() // write parameters
{
	char fname[200]; sprintf(fname, "%spara", path_);
	FILE *fp = fopen(fname, "rb");
	if (fp)	{ printf("Hash Tables Already Exist\n"); exit(1); }

	fp = fopen(fname, "wb");
	if (!fp) {
		printf("Could not create %s\n", fname);
		printf("Perhaps no such folder %s?\n", path_);
		return 1;
	}

	// write general parameters
	fwrite(&n_pts_,     SIZEINT, 1, fp);
	fwrite(&dim_,       SIZEINT, 1, fp);
	fwrite(&n_samples_, SIZEINT, 1, fp);
	fwrite(&n_blocks_,  SIZEINT, 1, fp);

	// write first level parameters
	fwrite(sample_index_, SIZEINT, n_blocks_*n_samples_, fp);

	// write second level parameters
	fwrite(block_size_, SIZEINT, n_blocks_, fp);
	fwrite(index_,      SIZEINT, n_pts_,    fp);
	fclose(fp);
	return 0;
}

// -----------------------------------------------------------------------------
template<class DType>
QALSH_PLUS<DType>::QALSH_PLUS(		// load index
	const char *path)					// index path
{
	strcpy(path_, path);

	// read parameters from disk
	if (read_params()) exit(1);

	// load first level lsh index (lsh_)
	char sample_path[200]; sprintf(sample_path, "%ssample/", path_);
	lsh_ = new QALSH<DType>(sample_path, sample_index_);

	// load second level lsh index (blocks_)
	int start = 0;
	for (int i = 0; i < n_blocks_; ++i) {
		char block_path[200]; sprintf(block_path, "%s%d/", path_, i);
		QALSH<DType> *lsh = new QALSH<DType>(block_path, 
			(const int*) &index_[start]);
		
		blocks_.push_back(lsh);
		start += block_size_[i];
	}
}

// -----------------------------------------------------------------------------
template<class DType>
int QALSH_PLUS<DType>::read_params()// read parameters
{
	char fname[200]; sprintf(fname, "%spara", path_);
	FILE* fp = fopen(fname, "rb");
	if (!fp) { printf("Could not open %s\n", fname); return 1; }

	// read general parameters
	fread(&n_pts_,     SIZEINT, 1, fp);
	fread(&dim_,       SIZEINT, 1, fp);
	fread(&n_samples_, SIZEINT, 1, fp);
	fread(&n_blocks_,  SIZEINT, 1, fp);

	// load first level parameters, sample_index_ and sample_index_to_block_
	int n_sample_pts = n_blocks_*n_samples_;
	sample_index_ = new int[n_sample_pts];
	fread(sample_index_, SIZEINT, n_sample_pts, fp);
	
	sample_index_to_block_ = new int[n_pts_];
	memset(sample_index_to_block_, -1, n_pts_);
	int bid = 0; // block id
	for (int i = 0; i < n_sample_pts; ++i) {
		int id = sample_index_[i];
		sample_index_to_block_[id] = bid;

		if ((i+1)%n_samples_ == 0) ++bid;
	}
	assert(bid == n_blocks_);

	// load second level parameters (block_size_, index_)
	block_size_ = new int[n_blocks_];
	index_ = new int[n_pts_];

	fread(block_size_, SIZEINT, n_blocks_, fp);
	fread(index_,      SIZEINT, n_pts_,    fp);
	fclose(fp);
	return 0;
}

// -----------------------------------------------------------------------------
template<class DType>
QALSH_PLUS<DType>::~QALSH_PLUS()	// destructor
{
	for (int i = 0; i < n_blocks_; ++i) {
		delete blocks_[i]; blocks_[i] = NULL;
	}
	blocks_.clear(); blocks_.shrink_to_fit();
	delete lsh_;

	delete[] block_size_;
	delete[] sample_index_to_block_;
	delete[] sample_index_;
	delete[] index_;
}

// -----------------------------------------------------------------------------
template<class DType>
void QALSH_PLUS<DType>::display()	// display parameters
{
	printf("Parameters of QALSH+:\n");
	printf("n         = %d\n", n_pts_);
	printf("d         = %d\n", dim_);
	printf("n_samples = %d\n", n_samples_);
	printf("n_blocks  = %d\n", n_blocks_);
	printf("path      = %s\n", path_);
	printf("\n");
}

// -----------------------------------------------------------------------------
template<class DType>
uint64_t QALSH_PLUS<DType>::knn(	// k-NN search
	int   top_k,						// top-k value
	int   nb,							// number of blocks for search
	const DType *query,					// input query
	const char *dfolder,				// data folder
	MinK_List *list)					// top-k results (return)
{
	assert(nb > 0 && nb <= n_blocks_);
	list->reset();

	// use sample data to determine the order of blocks for c-k-ANNS
	uint64_t page_io = 0;
	std::vector<int> block_order;
	page_io += get_block_order(nb, query, dfolder, block_order);

	// use <nb> blocks for c-k-ANNS
	for (int bid : block_order) {
		page_io += blocks_[bid]->knn2(top_k, query, dfolder, list);
	}
	block_order.clear(); block_order.shrink_to_fit();

	return page_io;
}

// -----------------------------------------------------------------------------
template<class DType>
uint64_t QALSH_PLUS<DType>::get_block_order(// get block order
	int   nb,							// number of blocks for search
	const DType *query,					// query point
	const char *dfolder,				// data folder
	std::vector<int> &block_order)		// block order (return)
{
	MinK_List *list = new MinK_List(MAXK);
	uint64_t page_io = lsh_->knn2(MAXK, query, dfolder, list);

	// init the counter of each block
	Result *pair = new Result[n_blocks_];
	for (int i = 0; i < n_blocks_; ++i) {
		pair[i].id_  = i;
		pair[i].key_ = 0.0f;
	}
	// select the first <nb> blocks with largest counters
	for (int i = 0; i < list->size(); ++i) {
		int bid = sample_index_to_block_[list->ith_id(i)];
		pair[bid].key_ += 1.0f;
	}
	qsort(pair, n_blocks_, sizeof(Result), ResultCompDesc);
	
	for (int i = 0; i < nb; ++i) {
		// if (fabs(pair[i].key_) < FLOATZERO) break;
		block_order.push_back(pair[i].id_);
	}
	delete[] pair;
	delete list;

	return page_io;
}

} // end namespace nns
