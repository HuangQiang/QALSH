#include "ann.h"

// -----------------------------------------------------------------------------
int linear_scan(					// brute-force linear scan (data in disk)
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   B,							// page size
	float p,							// the p value of Lp norm, p in (0,2]
	const float **query,				// query set
	const Result **R,					// truth set
	const char *data_folder,			// data folder
	const char *output_folder)			// output folder
{
	char output_set[200];
	strcpy(output_set, output_folder);
	strcat(output_set, "linear.out");

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s.\n", output_set);
		return 1;
	}
	
	// -------------------------------------------------------------------------
	//  k-NN search by Linear Scan
	// -------------------------------------------------------------------------
	printf("Top-k NN Search by Linear Scan:\n");
	printf("  Top-k\t\tRatio\t\tI/O\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MinK_List *list = new MinK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		g_io     = 0;		
		for (int i = 0; i < qn; ++i) {
			list->reset();
			g_io += linear(n, d, B, p, top_k, query[i], data_folder, list);
			g_recall += calc_recall(top_k, R[i], list);

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
		g_io      = (long long) ceil((float) g_io / (float) qn);

		printf("  %3d\t\t%.4f\t\t%lld\t\t%.2f\t\t%.2f\n", top_k, g_ratio, 
			g_io, g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%lld\t%f\t%f\n", top_k, g_ratio, g_io, 
			g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);
	
	return 0;
}

// -----------------------------------------------------------------------------
int indexing_of_qalsh_plus(			// indexing of qalsh+
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	int   leaf,							// leaf size of kd-tree
	int   L,							// number of projection (drusilla)
	int   M,							// number of candidates (drusilla)
	float p,							// the p value of Lp norm, p in (0,2]
	float zeta,							// symmetric factor of p-stable distr.
	float ratio,						// approximation ratio
	const float **data,					// data set
	const char *output_folder)			// output folder
{
	char output_set[200];
	strcpy(output_set, output_folder);
	strcat(output_set, "qalsh_plus.out");

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  indexing of qalsh+
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	char index_path[200];
	strcpy(index_path, output_folder);
	strcat(index_path, "qalsh_plus/");
	
	QALSH_PLUS *lsh = new QALSH_PLUS();
	lsh->build(n, d, B, leaf, L, M, p, zeta, ratio, data, index_path);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time = %f Seconds\n", indexing_time);
	printf("Memory = %f MB\n\n", g_memory / 1048576.0f);
	
	fprintf(fp, "index_time = %f Seconds\n", indexing_time);
	fprintf(fp, "memory     = %f MB\n\n", g_memory / 1048576.0f);
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	assert(g_memory == 0);

	return 0;
}

// -----------------------------------------------------------------------------
int knn_of_qalsh_plus(				// k-NN search of qalsh+
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const float **query,				// query set
	const Result **R,					// truth set
	const char *data_folder,			// data folder
	const char *output_folder)			// output folder
{
	char output_set[200];
	strcpy(output_set, output_folder);
	strcat(output_set, "qalsh_plus.out");

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  load qalsh+ 
	// -------------------------------------------------------------------------
	char index_path[200];
	strcpy(index_path, output_folder);
	strcat(index_path, "qalsh_plus/");

	QALSH_PLUS *lsh = new QALSH_PLUS();
	if (lsh->load(index_path)) {
		printf("Could not load QALSH_Plus\n");
		return 1;
	}
	lsh->display();

	// -------------------------------------------------------------------------
	//  c-k-ANN search by QALSH+
	// -------------------------------------------------------------------------
	int start = 1;
	int end = lsh->num_blocks_;
	assert(end >= start);

	printf("Top-k NN Search by QALSH+: \n");
	for (int nb = start; nb <= end; ++nb) {
		printf("  nb = %d\n", nb);
		fprintf(fp, "nb = %d\n", nb);

		printf("  Top-k\t\tRatio\t\tI/O\t\tTime (ms)\tRecall\n");
		for (int num = 0; num < MAX_ROUND; ++num) {
			gettimeofday(&g_start_time, NULL);
			int top_k = TOPK[num];
			MinK_List *list = new MinK_List(top_k);

			g_ratio  = 0.0f;
			g_recall = 0.0f;
			g_io     = 0;
			for (int i = 0; i < qn; ++i) {
				list->reset();
				g_io += lsh->knn(top_k, nb, query[i], data_folder, list);
				g_recall += calc_recall(top_k, R[i], list);

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
			g_io      = (long long) ceil((float) g_io / (float) qn);

			printf("  %3d\t\t%.4f\t\t%lld\t\t%.2f\t\t%.2f\n", top_k, g_ratio, 
				g_io, g_runtime, g_recall);
			fprintf(fp, "%d\t%f\t%lld\t%f\t%f\n", top_k, g_ratio, g_io, 
				g_runtime, g_recall);
		}
		printf("\n");
		fprintf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	assert(g_memory == 0);

	return 0;
}

// -----------------------------------------------------------------------------
int indexing_of_qalsh(				// indexing of qalsh
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	float p,							// the p value of Lp norm, p in (0,2]
	float zeta,							// symmetric factor of p-stable distr.
	float ratio,						// approximation ratio
	const float **data,					// data set
	const char *output_folder)			// output folder
{
	char output_set[200];
	strcpy(output_set, output_folder);
	strcat(output_set, "qalsh.out");
	
	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  indexing of QALSH
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	char index_path[200];
	strcpy(index_path, output_folder);
	strcat(index_path, "qalsh/");

	QALSH *lsh = new QALSH();
	lsh->build(n, d, B, p, zeta, ratio, data, index_path);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time = %f Seconds\n", indexing_time);
	printf("Memory = %f MB\n\n", g_memory / 1048576.0f);
	
	fprintf(fp, "index_time = %f Seconds\n", indexing_time);
	fprintf(fp, "memory     = %f MB\n\n", g_memory / 1048576.0f);
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	assert(g_memory == 0);

	return 0;
}

// -----------------------------------------------------------------------------
int knn_of_qalsh(					// k-NN search of qalsh
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const float **query,				// query set
	const Result **R,					// truth set
	const char *data_folder,			// data folder
	const char *output_folder)			// output folder
{
	// -------------------------------------------------------------------------
	//  load QALSH
	// -------------------------------------------------------------------------
	char index_path[200];
	strcpy(index_path, output_folder);
	strcat(index_path, "qalsh/");

	QALSH *lsh = new QALSH();
	if (lsh->load(index_path)) {
		printf("Could not load QALSH\n");
		return 1;
	}

	// -------------------------------------------------------------------------
	//  c-k-ANN search by QALSH
	// -------------------------------------------------------------------------
	char output_set[200];
	strcpy(output_set, output_folder);
	strcat(output_set, "qalsh.out");
	
	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	printf("Top-k NN Search by QALSH: \n");
	printf("  Top-k\t\tRatio\t\tI/O\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		gettimeofday(&g_start_time, NULL);
		int top_k = TOPK[num];
		MinK_List *list = new MinK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		g_io     = 0;
		for (int i = 0; i < qn; ++i) {
			list->reset();
			g_io += lsh->knn(top_k, query[i], data_folder, list);
			g_recall += calc_recall(top_k, R[i], list);

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
		g_io      = (long long) ceil((float) g_io / (float) qn);

		printf("  %3d\t\t%.4f\t\t%lld\t\t%.2f\t\t%.2f\n", top_k, g_ratio, 
			g_io, g_runtime, g_recall);
		fprintf(fp, "%d\t%f\t%lld\t%f\t%f\n", top_k, g_ratio, g_io, 
			g_runtime, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	assert(g_memory == 0);

	return 0;
}

