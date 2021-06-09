#pragma once

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>

#include <unistd.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>

#include "def.h"
#include "pri_queue.h"

namespace nns {

extern timeval  g_start_time;		// global parameter: start time
extern timeval  g_end_time;			// global parameter: end time

extern float    g_indexing_time;	// global parameter: indexing time
extern float    g_estimated_mem;	// global parameter: estimated memory

extern float    g_runtime;			// global parameter: running time
extern float    g_ratio;			// global parameter: overall ratio
extern float    g_recall;			// global parameter: recall
extern uint64_t g_page_io;			// global parameter: page i/o

// -------------------------------------------------------------------------
void create_dir(					// create directory
	char *path);						// input path

// -----------------------------------------------------------------------------
int write_buffer_to_page(			// write buffer to one page
	int   B,							// page size
	const char *fname,					// file name of data
	const char *buffer);				// buffer to store data

// -----------------------------------------------------------------------------
int read_buffer_from_page(			// read buffer from page
	int   B,							// page size
	const char *fname,					// file name of data
	char  *buffer);						// buffer to store data (return)

// -----------------------------------------------------------------------------
int write_ground_truth(				// write ground truth to disk
	int   n,							// number of ground truth results
	int   d,			 				// dimension of ground truth results
	float p,							// l_p distance
	const char *prefix,					// prefix of truth set
	const Result *truth);				// ground truth

// -----------------------------------------------------------------------------
float calc_ratio(					// calc overall ratio [1,\infinity)
	int   k,							// top-k value
	const Result *truth,				// ground truth results 
	MinK_List *list);					// top-k approximate results

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int   k,							// top-k value
	const Result *truth,				// ground truth results 
	MinK_List *list);					// results returned by algorithms

// -----------------------------------------------------------------------------
template<class DType>
int read_data(						// read data (binary) from disk
	int   n,							// cardinality
	int   d,			 				// dimensionality
	int   sign,							// 0-data; 1-query; 2-truth
	float p,							// L_p distance
	const char *prefix,					// prefix of data set
	DType *data)						// data (return)
{
	char fname[200]; 
	switch (sign) {
		case 0: sprintf(fname, "%s.ds", prefix); break;
		case 1: sprintf(fname, "%s.q", prefix); break;
		case 2: sprintf(fname, "%s.gt%3.1f", prefix, p); break;
		default: printf("Parameters error!\n"); return 1;
	}
	
	FILE *fp = fopen(fname, "rb");
	if (!fp) { printf("Could not open %s\n", fname); return 1; }

	uint64_t size = (uint64_t) n*d;
	uint32_t max_uint = 4294967295U;
	if (size < max_uint) {
		// if no longer than size_t, read the whole array directly
		fread(data, sizeof(DType), (uint32_t) size, fp);
	}
	else {
		// we store data points in linear order, n*d == d*n (save multi times)
		for (int i = 0; i < d; ++i) {
			fread(&data[(uint64_t)i*n], sizeof(DType), n, fp);
		}
	}
	fclose(fp);
	return 0;
}

// -----------------------------------------------------------------------------
template<class DType>
inline void write_data_to_buffer(	// write data to buffer
	int   n,							// number of data points
	int   d,							// dimensionality
	const DType *data,					// data set
	char  *buffer)						// buffer to store data (return)
{
	size_t start = 0;
	size_t size  = sizeof(DType);
	for (int i = 0; i < n*d; ++i) {
		memcpy(&buffer[start], &data[i], size);
		start += size;
	}
}

// -----------------------------------------------------------------------------
template<class DType>
int write_data_new_form(			// write dataset with new format
	int   n,							// number of data points
	int   d,							// data dimension
	int   B,							// page size
	const DType *data,					// data points
	const char *dfolder) 				// data folder
{
	char dpath[200]; sprintf(dpath, "%sdata/", dfolder);

	// check whether the directory exists
	if (access(dpath, F_OK) == 0) { 
		printf("New format data exist. No need writing data again.\n");
		return 0; 
	}
	create_dir(dpath);
	
	// compute num of data in one page and total number of data file
	int num = (int) floor((float) B / (d*sizeof(DType)));
	int total_file = (int) ceil((float) n / num); assert(total_file > 0);
	
	// write new format data for qalsh
	char *buffer = new char[B]; memset(buffer, 0, B*SIZECHAR); // one page
	int  start = 0;
	for (int i = 0; i < total_file; ++i) {
		// write data to buffer
		if (start + num > n) num = n - start;
		write_data_to_buffer<DType>(num, d, &data[(uint64_t)start*d], buffer);

		// write one page of data to disk
		char fname[200]; sprintf(fname, "%s%d.data", dpath, i);
		write_buffer_to_page(B, fname, (const char*) buffer);
		start += num;
	}
	assert(start == n);
	delete[] buffer;

	return 0;
}

// -----------------------------------------------------------------------------
template<class DType>
inline void read_data_from_buffer(	// read data from buffer
	int   id,							// index of data in buffer
	int   d,							// dimensionality
	const char *buffer,					// buffer to store data
	DType *data)						// data point (return)
{
	size_t size  = sizeof(DType);
	size_t start = id * d * size;
	for (int i = 0; i < d; ++i) {
		memcpy(&data[i], &buffer[start], size);
		start += size;
	}
}

// -----------------------------------------------------------------------------
template<class DType>
int read_data_new_format(			// read data with new format from disk
	int   id,							// index of data
	int   d,							// dimensionality
	int   B,							// page size
	const char *dfolder,				// data folder
	DType *data)						// data point (return)
{
	// compute num of data in one page and the data file id
	int num = (int) floor((float) B / (d*sizeof(DType)));
	int fid = (int) floor((float) id / num);

	// read buffer (one page of data) in new format from disk
	char *buffer = new char[B]; memset(buffer, 0, B*SIZECHAR);
	char fname[200]; sprintf(fname, "%sdata/%d.data", dfolder, fid);
	
	if (read_buffer_from_page(B, fname, buffer)) return 1;

	// read the specific data point from buffer
	read_data_from_buffer<DType>(id%num, d, (const char*) buffer, data);

	delete[] buffer;
	return 0;
}

// -----------------------------------------------------------------------------
template<class DType>
float calc_inner_product(			// calc inner product
	int   dim,							// dimension
	const float *p1,					// 1st point
	const DType *p2)					// 2nd point
{
	float r = 0.0f;
	for (int i = 0; i < dim; ++i) {
		r += (float) (p1[i] * p2[i]);
	}
	return r;
}

// -----------------------------------------------------------------------------
template<class DType>
float calc_l2_sqr(					// calc l2 square distance
	int   dim,							// dimension
	float threshold,					// threshold
	const DType *p1,					// 1st point
	const DType *p2)					// 2nd point
{
	unsigned d = dim & ~unsigned(7);
	const DType *aa = p1, *end_a = aa + d;
	const DType *bb = p2, *end_b = bb + d;

	// -------------------------------------------------------------------------
	// __builtin_prefetch (const void *addr[, rw[, locality]])
	// addr (required): Represents the address of the memory.
	//
	// rw (optional): A compile-time constant which can take the values:
	//   0 (default): prepare the prefetch for a read
	//   1 : prepare the prefetch for a write to the memory
	// 
	// locality (optional): A compile-time constant integer which can take the
	// following temporal locality (L) values:
	//   0: None, the data can be removed from the cache after the access
	//   1: Low, L3 cache, leave the data in L3 cache after access
	//   2: Moderate, L2 cache, leave the data in L2, L3 cache after access
	//   3 (default): High, L1 cache, leave the data in L1, L2, and L3 cache
	// -------------------------------------------------------------------------
	const int SHIFT = 8 * sizeof(DType);
	__builtin_prefetch(aa, 0, 3);
	__builtin_prefetch(bb, 0, 0);

	float r = 0.0f;
	float r0, r1, r2, r3, r4, r5, r6, r7;

	const DType *a = end_a, *b = end_b;

	r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = 0.0f;
	switch (dim & 7) {
		case 7: r6 = (float) SQR(a[6] - b[6]);
		case 6: r5 = (float) SQR(a[5] - b[5]);
		case 5: r4 = (float) SQR(a[4] - b[4]);
		case 4: r3 = (float) SQR(a[3] - b[3]);
		case 3: r2 = (float) SQR(a[2] - b[2]);
		case 2: r1 = (float) SQR(a[1] - b[1]);
		case 1: r0 = (float) SQR(a[0] - b[0]);
	}

	a = aa; b = bb;
	for (; a < end_a; a += 8, b += 8) {
		__builtin_prefetch(a+SHIFT, 0, 3);
		__builtin_prefetch(b+SHIFT, 0, 0);

		r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;
		if (r > threshold) return r;

		r0 = (float) SQR(a[0] - b[0]);
		r1 = (float) SQR(a[1] - b[1]);
		r2 = (float) SQR(a[2] - b[2]);
		r3 = (float) SQR(a[3] - b[3]);
		r4 = (float) SQR(a[4] - b[4]);
		r5 = (float) SQR(a[5] - b[5]);
		r6 = (float) SQR(a[6] - b[6]);
		r7 = (float) SQR(a[7] - b[7]);
	}
	r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;
	
	return r;
}

// -----------------------------------------------------------------------------
template<class DType>
float calc_l1_dist(					// calc Manhattan distance (l_1)
	int   dim,							// dimension
	float threshold,					// threshold
	const DType *p1,					// 1st point
	const DType *p2)					// 2nd point
{
	unsigned d = dim & ~unsigned(7);
	const DType *aa = p1, *end_a = aa + d;
	const DType *bb = p2, *end_b = bb + d;

	const int SHIFT = 8 * sizeof(DType);
	__builtin_prefetch(aa, 0, 3);
	__builtin_prefetch(bb, 0, 0);

	float r = 0.0f;
	float r0, r1, r2, r3, r4, r5, r6, r7;

	const DType *a = end_a, *b = end_b;

	r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = 0.0f;
	switch (dim & 7) {
		case 7: r6 = (float) fabs(a[6] - b[6]);
		case 6: r5 = (float) fabs(a[5] - b[5]);
		case 5: r4 = (float) fabs(a[4] - b[4]);
		case 4: r3 = (float) fabs(a[3] - b[3]);
		case 3: r2 = (float) fabs(a[2] - b[2]);
		case 2: r1 = (float) fabs(a[1] - b[1]);
		case 1: r0 = (float) fabs(a[0] - b[0]);
	}

	a = aa; b = bb;
	for (; a < end_a; a += 8, b += 8) {
		__builtin_prefetch(a+SHIFT, 0, 3);
		__builtin_prefetch(b+SHIFT, 0, 0);

		r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;
		if (r > threshold) return r;

		r0 = (float) fabs(a[0] - b[0]);
		r1 = (float) fabs(a[1] - b[1]);
		r2 = (float) fabs(a[2] - b[2]);
		r3 = (float) fabs(a[3] - b[3]);
		r4 = (float) fabs(a[4] - b[4]);
		r5 = (float) fabs(a[5] - b[5]);
		r6 = (float) fabs(a[6] - b[6]);
		r7 = (float) fabs(a[7] - b[7]);
	}
	r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;
	
	return r;
}

// -----------------------------------------------------------------------------
template<class DType>
float calc_l0_sqrt(					// calc l_{0.5} sqrt distance
	int   dim,							// dimension
	float threshold,					// threshold
	const DType *p1,					// 1st point
	const DType *p2)					// 2nd point
{
	unsigned d = dim & ~unsigned(7);
	const DType *aa = p1, *end_a = aa + d;
	const DType *bb = p2, *end_b = bb + d;

	const int SHIFT = 8 * sizeof(DType);
	__builtin_prefetch(aa, 0, 3);
	__builtin_prefetch(bb, 0, 0);

	float r = 0.0f;
	float r0, r1, r2, r3, r4, r5, r6, r7;

	const DType *a = end_a, *b = end_b;

	r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = 0.0f;
	switch (dim & 7) {
		case 7: r6 = sqrt((float) fabs(a[6] - b[6]));
		case 6: r5 = sqrt((float) fabs(a[5] - b[5]));
		case 5: r4 = sqrt((float) fabs(a[4] - b[4]));
		case 4: r3 = sqrt((float) fabs(a[3] - b[3]));
		case 3: r2 = sqrt((float) fabs(a[2] - b[2]));
		case 2: r1 = sqrt((float) fabs(a[1] - b[1]));
		case 1: r0 = sqrt((float) fabs(a[0] - b[0]));
	}

	a = aa; b = bb;
	for (; a < end_a; a += 8, b += 8) {
		__builtin_prefetch(a+SHIFT, 0, 3);
		__builtin_prefetch(b+SHIFT, 0, 0);

		r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;
		if (r > threshold) return r;

		r0 = sqrt((float) fabs(a[0] - b[0]));
		r1 = sqrt((float) fabs(a[1] - b[1]));
		r2 = sqrt((float) fabs(a[2] - b[2]));
		r3 = sqrt((float) fabs(a[3] - b[3]));
		r4 = sqrt((float) fabs(a[4] - b[4]));
		r5 = sqrt((float) fabs(a[5] - b[5]));
		r6 = sqrt((float) fabs(a[6] - b[6]));
		r7 = sqrt((float) fabs(a[7] - b[7]));
	}
	r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;
	
	return r;
}

// -----------------------------------------------------------------------------
template<class DType>
float calc_lp_pow(					// calc l_p pow_p distance
	int   dim,							// dimension
	float p,							// l_p distance, p \in (0,2]
	float threshold,					// threshold
	const DType *p1,					// 1st point
	const DType *p2)					// 2nd point
{
	unsigned d = dim & ~unsigned(7);
	const DType *aa = p1, *end_a = aa + d;
	const DType *bb = p2, *end_b = bb + d;

	const int SHIFT = 8 * sizeof(DType);
	__builtin_prefetch(aa, 0, 3);
	__builtin_prefetch(bb, 0, 0);

	float r = 0.0f;
	float r0, r1, r2, r3, r4, r5, r6, r7;

	const DType *a = end_a, *b = end_b;

	r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = 0.0f;
	switch (dim & 7) {
		case 7: r6 = pow((float) fabs(a[6] - b[6]), p);
		case 6: r5 = pow((float) fabs(a[5] - b[5]), p);
		case 5: r4 = pow((float) fabs(a[4] - b[4]), p);
		case 4: r3 = pow((float) fabs(a[3] - b[3]), p);
		case 3: r2 = pow((float) fabs(a[2] - b[2]), p);
		case 2: r1 = pow((float) fabs(a[1] - b[1]), p);
		case 1: r0 = pow((float) fabs(a[0] - b[0]), p);
	}

	a = aa; b = bb;
	for (; a < end_a; a += 8, b += 8) {
		__builtin_prefetch(a+SHIFT, 0, 3);
		__builtin_prefetch(b+SHIFT, 0, 0);

		r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;
		if (r > threshold) return r;

		r0 = pow((float) fabs(a[0] - b[0]), p);
		r1 = pow((float) fabs(a[1] - b[1]), p);
		r2 = pow((float) fabs(a[2] - b[2]), p);
		r3 = pow((float) fabs(a[3] - b[3]), p);
		r4 = pow((float) fabs(a[4] - b[4]), p);
		r5 = pow((float) fabs(a[5] - b[5]), p);
		r6 = pow((float) fabs(a[6] - b[6]), p);
		r7 = pow((float) fabs(a[7] - b[7]), p);
	}
	r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;
	
	return r;
}

// -----------------------------------------------------------------------------
template<class DType>
float calc_lp_dist(					// calc l_p distance
	int   dim,							// dimension
	float p,							// l_p distance, p \in (0,2]
	float threshold,					// threshold
	const DType *p1,					// 1st point
	const DType *p2)					// 2nd point
{
	if (fabs(p - 2.0f) < FLOATZERO) {
		return sqrt(calc_l2_sqr<DType>(dim, SQR(threshold), p1, p2));
	}
	else if (fabs(p - 1.0f) < FLOATZERO) {
		return calc_l1_dist<DType>(dim, threshold, p1, p2);
	}
	else if (fabs(p - 0.5f) < FLOATZERO) {
		float ret = calc_l0_sqrt<DType>(dim, sqrt(threshold), p1, p2);
		return SQR(ret);
	}
	else {
		float ret = calc_lp_pow<DType>(dim, p, pow(threshold, p), p1, p2);
		return pow(ret, 1.0f / p);
	}
}

// -----------------------------------------------------------------------------
template<class DType>
void kNN_search(					// k-NN search
	int   n, 							// cardinality
	int   d, 							// dimensionality
	int   k,							// top-k value
	float p,							// l_p distance, p \in (0,2]
	const DType *data,					// data points
	const DType *query,					// query point
	MinK_List *list)					// top-k results (return)
{
	float dist = -1.0f, kdist = MAXREAL;
	list->reset();
	for (int j = 0; j < n; ++j) {
		// data ID starts from 0
		dist = calc_lp_dist<DType>(d, p, kdist, &data[(uint64_t)j*d], query);
		kdist = list->insert(dist, j);
	}
}

// -----------------------------------------------------------------------------
template<class DType>
uint64_t linear(					// linear scan search
	int   n,							// number of data points
	int   d,							// dimensionality
	int   B,							// page size
	int   p,							// l_p distance, p \in (0,2]
	int   top_k,						// top-k value
	const DType *query,					// query point
	const char *dfolder,				// data folder
	MinK_List *list)					// k-NN results (return)
{
	list->reset();

	// assume data in disk, every time read a page of data ONLY
	char dpath[200]; sprintf(dpath, "%sdata/", dfolder);
	char *buffer = new char[B];		// one page buffer
	DType *data  = new DType[d];	// one data point

	// compute num of data in one page and total number of data file
	int num = (int) floor((float) B / (d*sizeof(DType)));
	int total_file = (int) ceil((float) n / num); assert(total_file > 0);
	
	// linear scan to find the k-NN of query
	int   id = 0, start = 0;
	float dist, kdist = MAXREAL;

	for (int i = 0; i < total_file; ++i) {
		// read one page of data into buffer
		char fname[200]; sprintf(fname, "%s%d.data", dpath, i);
		read_buffer_from_page(B, fname, buffer);

		// linear scan data points in one page buffer
		if (start + num > n) num = n - start;
		for (int j = 0; j < num; ++j) {
			read_data_from_buffer<DType>(j, d, (const char*) buffer, data);
			dist = calc_lp_dist<DType>(d, p, kdist, (const DType*) data, query);
			
			// data ID starts from 0
			kdist = list->insert(dist, id++);
		}
		start += num;
	}
	assert(start == n && id == n);
	delete[] buffer;
	delete[] data;
	
	return (uint64_t) total_file;
}

} // end namespace nns
