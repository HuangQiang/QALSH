#pragma once

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstring>

#include "def.h"
#include "util.h"
#include "qalsh.h"
#include "qalsh_plus.h"

namespace nns {

// -----------------------------------------------------------------------------
template<class DType>
int ground_truth(					// find ground truth
	int   n,							// number of data  points
	int   qn,							// number of query points
	int   d,							// dimensionality
	float p,							// l_p distance, p \in (0,2]
	const char *prefix,					// prefix of ground truth results
	const DType *data,					// data points
	const DType *query)					// query points
{
	gettimeofday(&g_start_time, NULL);
	Result *truth = new Result[qn*MAXK];
	MinK_List *list = new MinK_List(MAXK);

	for (int i = 0; i < qn; ++i) {
		kNN_search(n, d, MAXK, p, data, &query[(uint64_t)i*d], list);

		for (int j = 0; j < MAXK; ++j) {
			truth[i*MAXK+j].id_  = list->ith_id(j);
			truth[i*MAXK+j].key_ = list->ith_key(j);
		}
	}
	write_ground_truth(qn, MAXK, p, prefix, (const Result*) truth);
	delete   list;
	delete[] truth;

	gettimeofday(&g_end_time, NULL);
	float truth_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Ground Truth: %f Seconds\n\n", truth_time);

	return 0;
}

// -----------------------------------------------------------------------------
template<class DType>
int linear_scan(					// brute-force linear scan (data on disk)
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimensionality
	int   B,							// page size
	float p,							// l_p distance, p \in (0,2]
	const DType *query,					// query points
	const Result *truth,				// ground truth results
	const char *dfolder,				// data folder
	const char *ofolder)				// output folder
{
	char fname[200]; sprintf(fname, "%slinear.out", ofolder);
	FILE *fp = fopen(fname, "a+");
	if (!fp) { printf("Could not create %s\n", fname); return 1; }
	
	//  k-NN search by Linear Scan (assume data on disk)
	printf("k-NN Search by Linear Scan:\n");
	printf("Top-k\t\tRatio\t\tI/O\t\tTime (ms)\tRecall\n");
	for (int top_k : TOPKs) {
		gettimeofday(&g_start_time, NULL);
		MinK_List *list = new MinK_List(top_k);
		g_ratio   = 0.0f;
		g_recall  = 0.0f;
		g_page_io = 0;

		for (int i = 0; i < qn; ++i) {
			g_page_io += linear<DType>(n, d, B, p, top_k, 
				&query[(uint64_t)i*d], dfolder, list);
			g_ratio   += calc_ratio(top_k,  &truth[(uint64_t)i*MAXK], list);
			g_recall  += calc_recall(top_k, &truth[(uint64_t)i*MAXK], list);
		}
		delete list;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime*1000.0f) / qn;
		g_page_io = (uint64_t) ceil((double) g_page_io/qn);

		printf("%d\t\t%.4f\t\t%llu\t\t%.2f\t\t%.2f\n", top_k, g_ratio, 
			g_page_io, g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%llu\t%f\t%f\n", top_k, g_ratio, g_page_io, 
			g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");

	fclose(fp);
	return 0;
}

// -----------------------------------------------------------------------------
template<class DType>
int indexing_of_qalsh_plus(			// indexing of qalsh+
	int   n,							// number of data points
	int   d,							// dimensionality
	int   B,							// page size
	int   leaf,							// leaf size of kd-tree
	int   L,							// number of projection (drusilla)
	int   M,							// number of candidates (drusilla)
	float p,							// l_p distance, p \in (0,2]
	float zeta,							// symmetric factor of p-stable distr.
	float c,							// approximation ratio
	const DType *data,					// data points
	const char *ofolder)				// output folder
{
	char fname[200]; sprintf(fname, "%sqalsh_plus.out", ofolder);
	FILE *fp = fopen(fname, "a+");
	if (!fp) { printf("Could not create %s\n", fname); return 1; }

	//  indexing of QALSH+
	gettimeofday(&g_start_time, NULL);
	char path[200]; sprintf(path, "%sqalsh_plus/", ofolder);
	QALSH_PLUS<DType> *lsh = new QALSH_PLUS<DType>(n, d, B, leaf, L, M, p, 
		zeta, c, data, path);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	g_indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	g_estimated_mem = lsh->get_memory_usage() / 1048576.0f;

	printf("Indexing Time = %f Seconds\n", g_indexing_time);
	printf("Estimated Mem = %f MB\n\n", g_estimated_mem);
	fprintf(fp, "Indexing Time = %f Seconds\n", g_indexing_time);
	fprintf(fp, "Estimated Mem = %f MB\n\n", g_estimated_mem);

	fclose(fp);
	delete lsh;
	return 0;
}

// -----------------------------------------------------------------------------
template<class DType>
int knn_of_qalsh_plus(				// k-NN search of qalsh+
	int   qn,							// number of query points
	int   d,							// dimensionality
	const DType *query,					// query points
	const Result *truth,				// ground truth
	const char *dfolder,				// data folder
	const char *ofolder)				// output folder
{
	char fname[200]; sprintf(fname, "%sqalsh_plus.out", ofolder);
	FILE *fp = fopen(fname, "a+");
	if (!fp) { printf("Could not create %s\n", fname); return 1; }

	// load QALSH+
	gettimeofday(&g_start_time, NULL);
	char path[200]; sprintf(path, "%sqalsh_plus/", ofolder);
	QALSH_PLUS<DType> *lsh = new QALSH_PLUS<DType>(path);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	g_indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Load QALSH+ Index = %f Seconds\n\n", g_indexing_time);

	// c-k-ANNS by QALSH+
	printf("k-NN Search by QALSH+: \n");
	for (int nb = 1; nb <= lsh->get_num_blocks(); ++nb) {
		printf("nb = %d\n", nb);
		fprintf(fp, "nb = %d\n", nb);

		printf("Top-k\t\tRatio\t\tI/O\t\tTime (ms)\tRecall\n");
		for (int top_k : TOPKs) {
			gettimeofday(&g_start_time, NULL);
			MinK_List *list = new MinK_List(top_k);
			g_ratio   = 0.0f;
			g_recall  = 0.0f;
			g_page_io = 0;

			for (int i = 0; i < qn; ++i) {
				g_page_io += lsh->knn(top_k, nb, &query[(uint64_t)i*d], dfolder, list);
				g_ratio   += calc_ratio(top_k, &truth[(uint64_t)i*MAXK], list);
				g_recall  += calc_recall(top_k, &truth[(uint64_t)i*MAXK], list);
			}
			delete list;
			gettimeofday(&g_end_time, NULL);
			g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
				(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

			g_ratio   = g_ratio / qn;
			g_recall  = g_recall / qn;
			g_runtime = (g_runtime*1000.0f) / qn;
			g_page_io = (uint64_t) ceil((double) g_page_io/qn);

			printf("%d\t\t%.4f\t\t%llu\t\t%.2f\t\t%.2f\n", top_k, g_ratio, 
				g_page_io, g_runtime, g_recall);
			fprintf(fp, "%d\t%f\t%llu\t%f\t%f\n", top_k, g_ratio, g_page_io, 
				g_runtime, g_recall);
		}
		printf("\n");
		fprintf(fp, "\n");
	}
	fclose(fp);
	delete lsh;
	return 0;
}

// -----------------------------------------------------------------------------
template<class DType>
int indexing_of_qalsh(				// indexing of qalsh
	int   n,							// number of data points
	int   d,							// dimensionality
	int   B,							// page size
	float p,							// l_p distance, p \in (0,2]
	float zeta,							// symmetric factor of p-stable distr.
	float c,							// approximation ratio
	const DType *data,					// data points
	const char *ofolder)				// output folder
{
	char fname[200]; sprintf(fname, "%sqalsh.out", ofolder);
	FILE *fp = fopen(fname, "a+");
	if (!fp) { printf("Could not create %s\n", fname); return 1; }

	// indexing of QALSH
	gettimeofday(&g_start_time, NULL);
	char path[200]; sprintf(path, "%sqalsh/", ofolder);
	QALSH<DType> *lsh = new QALSH<DType>(n, d, B, p, zeta, c, data, path);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	g_indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	g_estimated_mem = lsh->get_memory_usage() / 1048576.0f;

	printf("Indexing Time = %f Seconds\n", g_indexing_time);
	printf("Estimated Mem = %f MB\n\n", g_estimated_mem);
	fprintf(fp, "Indexing Time = %f Seconds\n", g_indexing_time);
	fprintf(fp, "Estimated Mem = %f MB\n\n", g_estimated_mem);

	fclose(fp);
	delete lsh;
	return 0;
}

// -----------------------------------------------------------------------------
template<class DType>
int knn_of_qalsh(					// k-NN search of qalsh
	int   qn,							// number of query points
	int   d,							// dimensionality
	const DType *query,					// query points
	const Result *truth,				// ground truth
	const char *dfolder,				// data folder
	const char *ofolder)				// output folder
{
	char fname[200]; sprintf(fname, "%sqalsh.out", ofolder);
	FILE *fp = fopen(fname, "a+");
	if (!fp) { printf("Could not create %s\n", fname); return 1; }

	// load QALSH
	gettimeofday(&g_start_time, NULL);
	char path[200]; sprintf(path, "%sqalsh/", ofolder);
	QALSH<DType> *lsh = new QALSH<DType>(path);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	g_indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Load QALSH Index = %f Seconds\n\n", g_indexing_time);

	// c-k-ANNS by QALSH
	printf("k-NN Search by QALSH: \n");
	printf("Top-k\t\tRatio\t\tI/O\t\tTime (ms)\tRecall\n");
	for (int top_k : TOPKs) {
		gettimeofday(&g_start_time, NULL);
		MinK_List *list = new MinK_List(top_k);
		g_ratio   = 0.0f;
		g_recall  = 0.0f;
		g_page_io = 0;

		for (int i = 0; i < qn; ++i) {
			g_page_io += lsh->knn(top_k, &query[(uint64_t)i*d], dfolder, list);
			g_ratio   += calc_ratio(top_k,  &truth[(uint64_t)i*MAXK], list);
			g_recall  += calc_recall(top_k, &truth[(uint64_t)i*MAXK], list);
		}
		delete list;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime*1000.0f) / qn;
		g_page_io = (uint64_t) ceil((double) g_page_io/qn);

		printf("%d\t\t%.4f\t\t%llu\t\t%.2f\t\t%.2f\n", top_k, g_ratio, 
			g_page_io, g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%llu\t%f\t%f\n", top_k, g_ratio, g_page_io, 
			g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	
	fclose(fp);
	delete lsh;
	return 0;
}

} // end namespace nns
