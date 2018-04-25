#include "headers.h"

// -----------------------------------------------------------------------------
int ground_truth(					// find ground truth
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float p,							// the p value of Lp norm, p in (0,2]
	const char *data_set,				// address of data  set
	const char *query_set,				// address of query set
	const char *truth_set)				// address of truth set
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read data set and query set
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	float **data = new float*[n];
	for (int i = 0; i < n; ++i) data[i] = new float[d];
	if (read_data(n, d, data_set, data) == 1) {
		printf("Reading Dataset Error!\n");
		return 1;
	}

	float **query = new float*[qn];
	for (int i = 0; i < qn; ++i) query[i] = new float[d];
	if (read_data(qn, d, query_set, query) == 1) {
		printf("Reading Query Set Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Data and Query: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  find ground truth results (using linear scan method)
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	FILE *fp = fopen(truth_set, "w");
	if (!fp) {
		printf("Could not create %s.\n", truth_set);
		return 1;
	}

	MinK_List *list = new MinK_List(MAXK);
	fprintf(fp, "%d %d\n", qn, MAXK);
	for (int i = 0; i < qn; ++i) {
		list->reset();
		for (int j = 0; j < n; ++j) {
			float dist = calc_lp_dist(d, p, data[j], query[i]);
			list->insert(dist, j);
		}

		for (int j = 0; j < MAXK; ++j) {
			fprintf(fp, "%d %f ", list->ith_id(j), list->ith_key(j));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	gettimeofday(&end_time, NULL);
	float truth_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Ground Truth: %f Seconds\n\n", truth_time);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete list; list = NULL;
	for (int i = 0; i < n; ++i) {
		delete[] data[i]; data[i] = NULL;
	}
	delete[] data; data = NULL;

	for (int i = 0; i < qn; ++i) {
		delete[] query[i]; query[i] = NULL;
	}
	delete[] query; query = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int indexing_of_qalsh_plus(			// indexing of qalsh+
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	int   kd_leaf_size,					// leaf size of kd-tree
	int   L,							// number of projection (drusilla)
	int   M,							// number of candidates (drusilla)
	float p,							// the p value of Lp norm, p in (0,2]
	float zeta,							// symmetric factor of p-stable distr.
	float ratio,						// approximation ratio
	const char *data_set,				// address of data set
	const char *data_folder,			// data folder
	const char *output_folder)			// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read dataset
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	float **data = new float*[n];
	for (int i = 0; i < n; ++i) data[i] = new float[d];
	if (read_data(n, d, data_set, data) == 1) {
		printf("Reading Dataset Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Data: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  write the data set in new format to disk
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	write_data_new_form(n, d, B, (const float **) data, data_folder);

	gettimeofday(&end_time, NULL);
	float write_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Write Dataset in New Format: %f Seconds\n\n", write_file_time);

	// -------------------------------------------------------------------------
	//  indexing of QALSH_Plus
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	char index_path[200];
	strcpy(index_path, output_folder);
	strcat(index_path, "qalsh_plus/");

	QALSH_Plus *lsh = new QALSH_Plus();
	lsh->build(n, d, B, kd_leaf_size, L, M, p, zeta, ratio, 
		(const float **) data, index_path);

	gettimeofday(&end_time, NULL);
	float indexing_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  wrtie indexing time to disk
	// -------------------------------------------------------------------------
	char fname[200];
	strcpy(fname, output_folder);
	strcat(fname, "qalsh_plus.index");

	FILE *fp = fopen(fname, "w");
	if (!fp) {
		printf("Could not create %s.\n", fname);
		return 1;
	}
	fprintf(fp, "Indexing Time: %f seconds\n", indexing_time);
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	for (int i = 0; i < n; ++i) {
		delete[] data[i]; data[i] = NULL;
	}
	delete[] data; data = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int knn_of_qalsh_plus(				// k-NN search of qalsh+
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const char *query_set,				// address of query set
	const char *truth_set,				// address of truth set
	const char *data_folder,			// data folder
	const char *output_folder)			// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read query set and truth set
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	float** query = new float*[qn];
	for (int i = 0; i < qn; ++i) query[i] = new float[d];
	if (read_data(qn, d, query_set, query) == 1) {
		printf("Reading Query Set Error!\n");
		return 1;
	}

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Query and Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  load QALSH+ 
	// -------------------------------------------------------------------------
	char index_path[200];
	strcpy(index_path, output_folder);
	strcat(index_path, "qalsh_plus/");

	QALSH_Plus *lsh = new QALSH_Plus();
	if (lsh->load(index_path)) {
		printf("Could not load QALSH_Plus\n");
		return 1;
	}

	// -------------------------------------------------------------------------
	//  c-k-ANN search by QALSH+
	// -------------------------------------------------------------------------
	char output_set[200];
	strcpy(output_set, output_folder);
	strcat(output_set, "qalsh_plus.out");

	FILE *fp = fopen(output_set, "w");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	int kNNs[] = { 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
	int maxRound = 11;
	int top_k = -1;

	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;
	long long io_cost = -1;

	printf("c-k-ANN Search by QALSH+: \n");
	for (int nb = 2; nb <= 10; ++nb) {
		printf("  nb = %d\n", nb);
		fprintf(fp, "nb = %d\n", nb);

		printf("  Top-k\t\tRatio\t\tI/O\t\tTime (ms)\tRecall\n");
		for (int num = 0; num < maxRound; ++num) {
			gettimeofday(&start_time, NULL);
			top_k = kNNs[num];
			MinK_List *list = new MinK_List(top_k);

			overall_ratio = 0.0f;
			recall = 0.0f;
			io_cost = 0;
			for (int i = 0; i < qn; ++i) {
				list->reset();
				io_cost += lsh->knn(top_k, nb, query[i], data_folder, list);
				
				recall += calc_recall(top_k, (const Result *) R[i], list);

				float ratio = 0.0f;
				for (int j = 0; j < top_k; ++j) {
					ratio += list->ith_key(j) / R[i][j].key_;
				}
				overall_ratio += ratio / top_k;
			}
			delete list; list = NULL;
			gettimeofday(&end_time, NULL);
			runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
				start_time.tv_usec) / 1000000.0f;

			overall_ratio = overall_ratio / qn;
			recall = recall / qn;
			runtime = (runtime * 1000.0f) / qn;
			io_cost = (int) ceil((float) io_cost / (float) qn);

			printf("  %3d\t\t%.4f\t\t%lld\t\t%.2f\t\t%.2f\n", 
				top_k, overall_ratio, io_cost, runtime, recall);
			fprintf(fp, "%d\t%f\t%lld\t%f\t%f\n", 
				top_k, overall_ratio, io_cost, runtime, recall);
		}
		printf("\n");
		fprintf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	for (int i = 0; i < qn; ++i) {
		delete[] query[i]; query[i] = NULL;
		delete[] R[i]; R[i] = NULL;
	}
	delete[] query; query = NULL;
	delete[] R; R = NULL;

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
	const char *data_set,				// address of data set
	const char *data_folder,			// data folder
	const char *output_folder)			// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read dataset
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	float **data = new float*[n];
	for (int i = 0; i < n; ++i) data[i] = new float[d];
	if (read_data(n, d, data_set, data) == 1) {
		printf("Reading Dataset Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Data: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  write the data set in new format to disk
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	write_data_new_form(n, d, B, (const float **) data, data_folder);

	gettimeofday(&end_time, NULL);
	float write_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Write Dataset in New Format: %f Seconds\n\n", write_file_time);

	// -------------------------------------------------------------------------
	//  indexing of QALSH
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	char index_path[200];
	strcpy(index_path, output_folder);
	strcat(index_path, "qalsh/");

	QALSH *lsh = new QALSH();
	lsh->build(n, d, B, p, zeta, ratio, (const float **) data, index_path);

	gettimeofday(&end_time, NULL);
	float indexing_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time: %f Seconds\n\n", indexing_time);

	// -------------------------------------------------------------------------
	//  write indexing time to disk
	// -------------------------------------------------------------------------
	char fname[200];
	strcpy(fname, output_folder);
	strcat(fname, "qalsh.index");

	FILE *fp = fopen(fname, "w");
	if (!fp) {
		printf("Could not create %s.\n", fname);
		return 1;
	}
	fprintf(fp, "Indexing Time: %f seconds\n", indexing_time);
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	for (int i = 0; i < n; ++i) {
		delete[] data[i]; data[i] = NULL;
	}
	delete[] data; data = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int knn_of_qalsh(					// k-NN search of qalsh
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const char *query_set,				// address of query set
	const char *truth_set,				// address of truth set
	const char *data_folder,			// data folder
	const char *output_folder)			// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read query set and truth set
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	float** query = new float*[qn];
	for (int i = 0; i < qn; ++i) query[i] = new float[d];
	if (read_data(qn, d, query_set, query) == 1) {
		printf("Reading Query Set Error!\n");
		return 1;
	}

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Query and Ground Truth: %f Seconds\n\n", read_file_time);

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
	
	FILE *fp = fopen(output_set, "w");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	int kNNs[] = { 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
	int maxRound = 11;
	int top_k = -1;

	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;
	long long io_cost = -1;

	printf("c-k-ANN Search by QALSH: \n");
	printf("  Top-k\t\tRatio\t\tI/O\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < maxRound; ++num) {
		gettimeofday(&start_time, NULL);
		top_k = kNNs[num];
		overall_ratio = 0.0f;
		recall = 0.0f;
		io_cost = 0;

		MinK_List *list = new MinK_List(top_k);
		for (int i = 0; i < qn; ++i) {
			list->reset();
			io_cost += lsh->knn(top_k, query[i], data_folder, list);

			recall += calc_recall(top_k, (const Result *) R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			overall_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&end_time, NULL);
		runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
			start_time.tv_usec) / 1000000.0f;

		overall_ratio = overall_ratio / qn;
		recall = recall / qn;
		runtime = (runtime * 1000.0f) / qn;
		io_cost = (int) ceil((float) io_cost / (float) qn);

		printf("  %3d\t\t%.4f\t\t%lld\t\t%.2f\t\t%.2f\n", 
			top_k, overall_ratio, io_cost, runtime, recall);
		fprintf(fp, "%d\t%f\t%lld\t%f\t%f\n", 
			top_k, overall_ratio, io_cost, runtime, recall);
	}
	printf("\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	for (int i = 0; i < qn; ++i) {
		delete[] query[i]; query[i] = NULL;
		delete[] R[i]; R[i] = NULL;
	}
	delete[] query; query = NULL;
	delete[] R; R = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int linear_scan(					// brute-force linear scan (data in disk)
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   B,							// page size
	float p,							// the p value of Lp norm, p in (0,2]
	const char *query_set,				// address of query set
	const char *truth_set,				// address of truth set
	const char *data_folder,			// data folder
	const char *output_folder)			// output folder
{
	timeval start_time, end_time;

	// -------------------------------------------------------------------------
	//  read query set and truth set
	// -------------------------------------------------------------------------
	gettimeofday(&start_time, NULL);
	float** query = new float*[qn];
	for (int i = 0; i < qn; ++i) query[i] = new float[d];
	if (read_data(qn, d, query_set, query) == 1) {
		printf("Reading Query Set Error!\n");
		return 1;
	}

	Result **R = new Result*[qn];
	for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
	if (read_ground_truth(qn, truth_set, R) == 1) {
		printf("Reading Truth Set Error!\n");
		return 1;
	}

	gettimeofday(&end_time, NULL);
	float read_file_time = end_time.tv_sec - start_time.tv_sec + 
		(end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
	printf("Read Query and Ground Truth: %f Seconds\n\n", read_file_time);

	// -------------------------------------------------------------------------
	//  k-NN search by Linear Scan
	// -------------------------------------------------------------------------
	char output_set[200];
	strcpy(output_set, output_folder);
	strcat(output_set, "linear.out");

	FILE *fp = fopen(output_set, "w");
	if (!fp) {
		printf("Could not create %s.\n", output_set);
		return 1;
	}

	int kNNs[] = { 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
	int maxRound = 11;
	int top_k = -1;

	float runtime = -1.0f;
	float overall_ratio = -1.0f;
	float recall = -1.0f;
	long long io_cost = -1;

	printf("k-NN Search by Linear Scan:\n");
	printf("  Top-k\t\tRatio\t\tI/O\t\tTime (ms)\tRecall\n");
	for (int round = 0; round < maxRound; round++) {
		gettimeofday(&start_time, NULL);
		top_k = kNNs[round];
		overall_ratio = 0.0f;
		recall = 0.0f;
		io_cost = 0;

		MinK_List *list = new MinK_List(top_k);
		for (int i = 0; i < qn; ++i) {
			list->reset();
			io_cost += linear(n, d, B, p, top_k, (const float*) query[i], 
				data_folder, list);

			recall += calc_recall(top_k, (const Result *) R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			overall_ratio += ratio / top_k;
		}
		delete list; list = NULL;
		gettimeofday(&end_time, NULL);
		runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - 
			start_time.tv_usec) / 1000000.0f;

		overall_ratio = overall_ratio / qn;
		recall = recall / qn;
		runtime = (runtime * 1000.0f) / qn;
		io_cost = (int) ceil((float) io_cost / (float) qn);

		printf("  %3d\t\t%.4f\t\t%lld\t\t%.2f\t\t%.2f\n", 
			top_k, overall_ratio, io_cost, runtime, recall);
		fprintf(fp, "%d\t%f\t%lld\t%f\t%f\n", 
			top_k, overall_ratio, io_cost, runtime, recall);
	}
	printf("\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	for (int i = 0; i < qn; ++i) {
		delete[] query[i]; query[i] = NULL;
		delete[] R[i]; R[i] = NULL;
	}
	delete[] query; query = NULL;
	delete[] R; R = NULL;
	
	return 0;
}