#include "headers.h"

// -----------------------------------------------------------------------------
int linear_scan(					// k-NN search by linear scan
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimensionality
	float p,							// the p value of Lp norm, p in (0,2]
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *out_path)				// output path
{
	char output_set[200];
	sprintf(output_set, "%slinear_mem.out", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s.\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  k-NN search by linear scan
	// -------------------------------------------------------------------------
	printf("Top-k NN by Linear Scan:\n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		
		int top_k = TOPK[num];
		Result **result = new Result*[qn];
		for (int i = 0; i < qn; ++i) {
			result[i] = new Result[top_k];
		}
		k_nn_search(n, qn, d, top_k, p, data, query, result);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			g_recall += calc_recall(top_k, (const Result*) R[i], 
				(const Result *) result[i]);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += result[i][j].key_ / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		for (int i = 0; i < qn; ++i) {
			delete[] result[i]; result[i] = NULL;
		}
		delete[] result; result = NULL;

		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.2f\t\t%.2f\n", top_k, g_ratio, g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);
	
	return 0;
}

// -----------------------------------------------------------------------------
int qalsh_plus(						// k-NN search by qalsh+
	int   n,							// number of data  objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   leaf,							// leaf size of kd-tree
	int   L,							// number of projection (drusilla)
	int   M,							// number of candidates (drusilla)
	float p,							// the p value of Lp norm, p in (0,2]
	float zeta,							// symmetric factor of p-stable distr.
	float ratio,						// approximation ratio
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *out_path)				// output path
{
	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%sqalsh_plus_mem.out", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s.\n", output_set);
		return 1;
	}

	gettimeofday(&g_start_time, NULL);
	QALSH_PLUS *lsh = new QALSH_PLUS(n, d, leaf, L, M, p, zeta, ratio, data);
	lsh->display();
	
	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time = %f Seconds\n\n", indexing_time);
	fprintf(fp, "Indexing Time = %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  k-NN search by QALSH+
	// -------------------------------------------------------------------------
	int start = 2, end = 10;
	if (leaf >= 40000) {
		start = 5; end = 20;
	}
	printf("Top-k NN Search by QALSH+:\n");
	for (int nb = start; nb <= end; ++nb) {
		printf("  nb = %d\n", nb);
		fprintf(fp, "nb = %d\n", nb);

		printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
		for (int num = 0; num < MAX_ROUND; ++num) {
			gettimeofday(&g_start_time, NULL);
			int top_k = TOPK[num];
			MinK_List *list = new MinK_List(top_k);

			g_ratio  = 0.0f;
			g_recall = 0.0f;
			for (int i = 0; i < qn; ++i) {
				list->reset();
				lsh->knn(top_k, nb, query[i], list);
				g_recall += calc_recall(top_k, (const Result*) R[i], list);

				float ratio = 0.0f;
				for (int j = 0; j < top_k; ++j) {
					ratio += list->ith_key(j) / R[i][j].key_;
				}
				g_ratio += ratio / top_k;
			}
			delete list; list = NULL;
			gettimeofday(&g_end_time, NULL);
			g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
				(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

			g_ratio   = g_ratio / qn;
			g_recall  = g_recall / qn;
			g_runtime = (g_runtime * 1000.0f) / qn;

			printf("  %3d\t\t%.4f\t\t%.2f\t\t%.2f\n", top_k, g_ratio, g_runtime, g_recall);
			fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, g_recall);
		}
		printf("\n");
		fprintf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	
	return 0;
}

// -----------------------------------------------------------------------------
int qalsh(							// k-NN search by qalsh
	int   n,							// number of data  objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float p,							// the p value of Lp norm, p in (0,2]
	float zeta,							// symmetric factor of p-stable distr.
	float ratio,						// approximation ratio
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *out_path)				// output path
{	
	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	char output_set[200];
	sprintf(output_set, "%sqalsh_mem.out", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s.\n", output_set);
		return 1;
	}

	gettimeofday(&g_start_time, NULL);
	QALSH *lsh = new QALSH(n, d, p, zeta, ratio, data);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time = %f Seconds\n\n", indexing_time);
	fprintf(fp, "Indexing Time = %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  c-k-ANN search
	// -------------------------------------------------------------------------
	printf("Top-k NN Search by QALSH:\n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MinK_List *list = new MinK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			lsh->knn(top_k, query[i], list);
			g_recall += calc_recall(top_k, (const Result*) R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&g_end_time, NULL);
		g_runtime = g_end_time.tv_sec - g_start_time.tv_sec + 
			(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime = (g_runtime * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.2f\t\t%.2f\n", top_k, g_ratio, g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, g_ratio, g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;

	return 0;
}

