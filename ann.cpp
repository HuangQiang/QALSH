#include "headers.h"


// -----------------------------------------------------------------------------
int ground_truth(					// output the ground truth results
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	float p,							// the p value of Lp norm, p in (0,2]
	char* data_set,						// address of data set
	char* query_set,					// address of query set
	char* truth_set)					// address of ground truth file
{
	clock_t startTime = (clock_t) -1;
	clock_t endTime   = (clock_t) -1;

	int i, j;
	FILE* fp = NULL;

	// -------------------------------------------------------------------------
	//  Read data set and query set
	// -------------------------------------------------------------------------
	startTime = clock();
	g_memory += SIZEFLOAT * (n + qn) * d;
	if (check_mem()) return 1;

	float** data = new float*[n];
	for (i = 0; i < n; i++) data[i] = new float[d];
	if (read_set(n, d, data_set, data)) {
		error("Reading Dataset Error!\n", true);
	}

	float** query = new float*[qn];
	for (i = 0; i < qn; i++) query[i] = new float[d];
	if (read_set(qn, d, query_set, query) == 1) {
		error("Reading Query Set Error!\n", true);
	}
	endTime = clock();
	printf("Read Dataset and Query Set: %.6f Seconds\n\n", 
		((float) endTime - startTime) / CLOCKS_PER_SEC);

	// -------------------------------------------------------------------------
	//  output ground truth results (using linear scan method)
	// -------------------------------------------------------------------------
	int maxk = MAXK;
	float dist = -1.0F;
	float* knndist = new float[maxk];
	g_memory += SIZEFLOAT * maxk;

	fp = fopen(truth_set, "w");		// open output file
	if (!fp) {
		printf("I could not create %s.\n", truth_set);
		return 1;
	}

	fprintf(fp, "%d %d\n", qn, maxk);
	for (i = 0; i < qn; i++) {
		for (j = 0; j < maxk; j++) {
			knndist[j] = MAXREAL;
		}
									// find k-nn points of query
		for (j = 0; j < n; j++) {
			dist = calc_lp_dist(p, data[j], query[i], d);

			int ii, jj;
			for (jj = 0; jj < maxk; jj++) {
				if (compfloats(dist, knndist[jj]) == -1) {
					break;
				}
			}
			if (jj < maxk) {
				for (ii = maxk - 1; ii >= jj + 1; ii--) {
					knndist[ii] = knndist[ii - 1];
				}
				knndist[jj] = dist;
			}
		}

		fprintf(fp, "%d", i + 1);	// output Lp dist of k-nn points
		for (j = 0; j < maxk; j++) {
			fprintf(fp, " %f", knndist[j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);						// close output file
	endTime = clock();
	printf("Generate Ground Truth: %.6f Seconds\n\n", 
		((float) endTime - startTime) / CLOCKS_PER_SEC);
	
	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	for (i = 0; i < n; i++) {
		delete[] data[i]; data[i] = NULL;
	}
	delete[] data; data = NULL;
	g_memory -= SIZEFLOAT * n * d;
	
	for (i = 0; i < qn; i++) {
		delete[] query[i]; query[i] = NULL;
	}
	delete[] query; query = NULL;
	g_memory -= SIZEFLOAT * qn * d;
	
	delete[] knndist; knndist = NULL;
	g_memory -= SIZEFLOAT * maxk;
	
	//printf("memory = %.2f MB\n", (float) g_memory / (1024.0f * 1024.0f));
	return 0;
}

// -----------------------------------------------------------------------------
int indexing(						// build hash tables for the dataset
	int   n,							// number of data points
	int   d,							// dimension of space
	int   B,							// page size
	float p,							// the p value of Lp norm, p in (0,2]
	float zeta,							// symmetric factor of p-stable distr.
	float ratio,						// approximation ratio
	char* data_set,						// address of data set
	char* data_folder,					// folder to store new format of data
	char* output_folder)				// folder to store info of qalsh
{
	clock_t startTime = (clock_t) -1;
	clock_t endTime   = (clock_t) -1;

	// -------------------------------------------------------------------------
	//  Read data set
	// -------------------------------------------------------------------------
	startTime = clock();
	g_memory += SIZEFLOAT * n * d;	
	if (check_mem()) return 1;

	float** data = new float*[n];
	for (int i = 0; i < n; i++) data[i] = new float[d];
	if (read_set(n, d, data_set, data) == 1) {
		error("Reading Dataset Error!\n", true);
	}
	endTime = clock();
	printf("Read Dataset: %.6f Seconds\n\n", 
		((float) endTime - startTime) / CLOCKS_PER_SEC);

	// -------------------------------------------------------------------------
	//  Write the data set in new format to disk
	// -------------------------------------------------------------------------
	startTime = clock();
	write_data_new_form(n, d, B, data, data_folder);
	endTime = clock();
	printf("Write Dataset in New Format: %.6f Seconds\n\n", 
		((float) endTime - startTime) / CLOCKS_PER_SEC);
	
	// -------------------------------------------------------------------------
	//  Indexing
	// -------------------------------------------------------------------------
	char fname[200];
	strcpy(fname, output_folder);
	get_lp_filename(p, fname);
	strcat(fname, "_index.out");

	FILE* fp = fopen(fname, "w");
	if (!fp) {
		printf("I could not create %s.\n", fname);
		return 1;					// fail to return
	}

	startTime = clock();
	QALSH* lsh = new QALSH();
	lsh->init(n, d, B, p, zeta, ratio, output_folder);
	lsh->bulkload(data);
	endTime = clock();

	float indexing_time = ((float) endTime - startTime) / CLOCKS_PER_SEC;
	printf("\nIndexing Time: %.6f seconds\n\n", indexing_time);
	fprintf(fp, "Indexing Time: %.6f seconds\n", indexing_time);
	fclose(fp);
	
	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	for (int i = 0; i < n; i++) {
		delete[] data[i]; data[i] = NULL;
	}
	delete[] data; data = NULL;
	g_memory -= SIZEFLOAT * n * d;
	
	//printf("memory = %.2f MB\n", (float) g_memory / (1024.0f * 1024.0f));
	return 0;
}

// -----------------------------------------------------------------------------
int lshknn(							// k-nn via qalsh (data in disk)
	int   qn,							// number of query points
	int   d,							// dimensionality
	float p,							// the p value of Lp norm, p in (0,2]
	char* query_set,					// path of query set
	char* truth_set,					// groundtrue file
	char* data_folder,					// folder to store new format of data
	char* output_folder)				// output folder
{
	int ret = 0;
	int maxk = MAXK;
	int i, j;
	FILE* fp = NULL;				// file pointer

	// -------------------------------------------------------------------------
	//  Read query set
	// -------------------------------------------------------------------------
	g_memory += SIZEFLOAT * qn * d;
	float** query = new float*[qn];
	for (i = 0; i < qn; i++) query[i] = new float[d];
	if (read_set(qn, d, query_set, query)) {
		error("Reading Query Set Error!\n", true);
	}

	// -------------------------------------------------------------------------
	//  Read the ground truth file
	// -------------------------------------------------------------------------
	g_memory += SIZEFLOAT * qn * maxk;
	float* R = new float[qn * maxk];

	fp = fopen(truth_set, "r");		// open ground truth file
	if (!fp) {
		printf("Could not open the ground truth file.\n");
		return 1;
	}

	fscanf(fp, "%d %d\n", &qn, &maxk);
	for (int i = 0; i < qn; i++) {
		fscanf(fp, "%d", &j);
		for (j = 0; j < maxk; j ++) {
			fscanf(fp, " %f", &(R[i * maxk + j]));
		}
	}
	fclose(fp);						// close groundtrue file

	// -------------------------------------------------------------------------
	//  K-nearest neighbor (k-nn) search via qalsh
	// -------------------------------------------------------------------------
	int kNNs[] = {1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
	int maxRound = 11;
	int top_k = 0;

	float allTime   = -1.0f;
	float allRatio  = -1.0f;
	float thisRatio = -1.0f;

	int thisIO = 0;
	int allIO  = 0;

	clock_t startTime = (clock_t) -1.0;
	clock_t endTime   = (clock_t) -1.0;
									// init the results
	g_memory += (long) sizeof(ResultItem) * maxk;
	ResultItem* rslt = new ResultItem[maxk];
	for (i = 0; i < maxk; i++) {
		rslt[i].id_ = -1;
		rslt[i].dist_ = MAXREAL;
	}

	QALSH* lsh = new QALSH();		// restore QALSH
	if (lsh->restore(p, output_folder)) {
		error("Could not restore qalsh\n", true);
	}

	char output_set[200];
	strcpy(output_set, output_folder);
	get_lp_filename(p, output_set);
	strcat(output_set, "_qalsh.out");

	fp = fopen(output_set, "w");	// open output file
	if (!fp) {
		printf("Could not create the output file.\n");
		return 1;
	}

	printf("QALSH for c-k-ANN Search: \n");
	printf("  Top-k\tRatio\t\tI/O\t\tTime (ms)\n");
	for (int num = 0; num < maxRound; num++) {
		top_k = kNNs[num];

		allRatio = 0.0f;
		allIO = 0;
		startTime = clock();
		for (i = 0; i < qn; i++) {
			thisIO = lsh->knn(query[i], top_k, rslt, data_folder);

			thisRatio = 0.0f;
			for (j = 0; j < top_k; j++) {
				thisRatio += rslt[j].dist_ / R[i * maxk + j];
			}
			thisRatio /= top_k;
			allRatio += thisRatio;
			allIO += thisIO;
			//printf("        No.%2d: Top-k = %d, IO = %4d, Ratio = %0.6f\n",
			//	i+1, top_k, thisIO, thisRatio);
		}
		endTime = clock();
		allTime = ((float) endTime - startTime) / CLOCKS_PER_SEC;

		allRatio = allRatio / qn;
		allTime  = (allTime * 1000.0f) / qn;
		allIO    = (int) ceil((float) allIO / (float) qn);

		printf("  %3d\t\t%.4f\t\t%d\t%.2f\n", top_k, allRatio, allIO, allTime);
		fprintf(fp, "%d\t%f\t%d\t%f\n", top_k, allRatio, allIO, allTime);
	}
	printf("\n");
	fclose(fp);						// close output file

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;
	for (i = 0; i < qn; i++) {
		delete[] query[i]; query[i] = NULL;
	}
	delete[] query; query = NULL;
	g_memory -= SIZEFLOAT * qn * d;
	
	delete[] R; R = NULL;
	delete[] rslt; rslt = NULL;
	g_memory -= (SIZEFLOAT * qn * maxk + sizeof(ResultItem) * maxk);
	
	//printf("memory = %.2f MB\n", (float) g_memory / (1024.0f * 1024.0f));
	return ret;
}

// -----------------------------------------------------------------------------
int linear_scan(					// brute-force linear scan (data in disk)
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimension of space
	int   B,							// page size
	float p,							// the p value of Lp norm, p in (0,2]
	char* query_set,					// address of query set
	char* truth_set,					// address of ground truth file
	char* data_folder,					// folder to store new format of data
	char* output_folder)				// output folder
{
	clock_t startTime = (clock_t)-1.0f;
	clock_t endTime = (clock_t)-1.0f;

	// -------------------------------------------------------------------------
	//  Allocation and initialzation.
	// -------------------------------------------------------------------------
	int kNNs[] = {1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
	int maxRound = 11;
	
	int i, j, top_k;
	int maxk = MAXK;

	float allTime   = -1.0f;
	float thisRatio = -1.0f;
	float allRatio  = -1.0f;

	g_memory += (SIZEFLOAT * (qn + 1) * maxk);
	float* knn_dist = new float[maxk];
	for (i = 0; i < maxk; i++) knn_dist[i] = MAXREAL;

	float** R = new float*[qn];
	for (i = 0; i < qn; i++) {
		R[i] = new float[maxk];
		for (j = 0; j < maxk; j++) R[i][j] = 0.0f;
	}

	g_memory += (SIZEFLOAT * d * 2 + SIZECHAR * B);
	float* data  = new float[d];	// one data object
	float* query = new float[d];	// one query object
	char* buffer = new char[B];		// every time can read one page

	char fname[200];				// file name for data
	char data_path[200];			// data path
	char out_set[200];				// output file

	// -------------------------------------------------------------------------
	//  Open the output file, and read the ground true results
	// -------------------------------------------------------------------------
	strcpy(out_set, output_folder);	// generate output file
	get_lp_filename(p, out_set);
	strcat(out_set, "_linear.out");

	FILE* ofp = fopen(out_set, "w");
	if (!ofp) {
		printf("I could not create %s.\n", out_set);
		return 1;
	}
									// open ground true file
	FILE* tfp = fopen(truth_set, "r");
	if (!tfp) {
		printf("I could not create %s.\n", truth_set);
		return 1;
	}
									// read top-k nearest distance
	fscanf(tfp, "%d %d\n", &qn, &maxk);
	for (int i = 0; i < qn; i++) {
		fscanf(tfp, "%d", &j);
		for (j = 0; j < maxk; j ++) {
			fscanf(tfp, " %f", &(R[i][j]));
		}
	}
	fclose(tfp);					// close ground true file

	// -------------------------------------------------------------------------
	//  Calc the number of data object in one page and the number of data file.
	//  <num> is the number of data in one data file
	//  <total_file> is the total number of data file
	// -------------------------------------------------------------------------
	int num = (int) floor((float) B / (d * SIZEFLOAT));
	int total_file = (int) ceil((float) n / num);
	if (total_file == 0) return 1;

	// -------------------------------------------------------------------------
	//  Brute-force linear scan method (data in disk)
	//  For each query, we limit that we can ONLY read one page of data.
	// -------------------------------------------------------------------------
	int count = 0;
	float dist = -1.0F;
									// generate the data path
	strcpy(data_path, data_folder);
	strcat(data_path, "data/");

	printf("Linear Scan Search:\n");
	printf("  Top-k\tRatio\t\tI/O\t\tTime (ms)\n");
	for (int round = 0; round < maxRound; round++) {
		top_k = kNNs[round];
		allRatio = 0.0f;

		startTime = clock();
		FILE* qfp = fopen(query_set, "r");
		if (!qfp) error("Could not open the query set.\n", true);

		for (i = 0; i < qn; i++) {
			// -----------------------------------------------------------------
			//  Step 1: read a query from disk and init the k-nn results
			// -----------------------------------------------------------------
			fscanf(qfp, "%d", &j);
			for (j = 0; j < d; j++) {
				fscanf(qfp, " %f", &query[j]);
			}

			for (j = 0; j < maxk; j++) {
				knn_dist[j] = MAXREAL;
			}

			// -----------------------------------------------------------------
			//  Step 2: find k-nn results for the query
			// -----------------------------------------------------------------
			for (j = 0; j < total_file; j++) {
				// -------------------------------------------------------------
				//  Step 2.1: get the file name of current data page
				// -------------------------------------------------------------
				get_data_filename(j, data_path, fname);

				// -------------------------------------------------------------
				//  Step 2.2: read one page of data into buffer
				// -------------------------------------------------------------
				if (read_buffer_from_page(B, fname, buffer) == 1) {
					error("error to read a data page", true);
				}

				// -------------------------------------------------------------
				//  Step 2.3: find the k-nn results in this page. NOTE: the 
				// 	number of data in the last page may be less than <num>
				// -------------------------------------------------------------
				if (j < total_file - 1) count = num;
				else count = n - num * (total_file - 1);

				for (int z = 0; z < count; z++) {
					read_data_from_buffer(z, d, data, buffer);
					dist = calc_lp_dist(p, data, query, d);

					int ii, jj;
					for (jj = 0; jj < top_k; jj++) {
						if (compfloats(dist, knn_dist[jj]) == -1) {
							break;
						}
					}
					if (jj < top_k) {
						for (ii = top_k - 1; ii >= jj + 1; ii--) {
							knn_dist[ii] = knn_dist[ii - 1];
						}
						knn_dist[jj] = dist;
					}
				}
			}

			thisRatio = 0.0f;
			for (j = 0; j < top_k; j++) {
				thisRatio += knn_dist[j] / R[i][j];
			}
			thisRatio /= top_k;
			allRatio += thisRatio;
		}
		// -----------------------------------------------------------------
		//  Step 3: output result of top-k nn points
		// -----------------------------------------------------------------
		fclose(qfp);				// close query file
		endTime  = clock();
		allTime  = ((float) endTime - startTime) / CLOCKS_PER_SEC;
		allTime  = (allTime * 1000.0f) / qn;
		allRatio = allRatio / qn;
									// output results
		printf("  %3d\t\t%.4f\t\t%d\t%.2f\n", top_k, allRatio, 
			total_file, allTime);
		fprintf(ofp, "%d\t%f\t%d\t%f\n", top_k, allRatio, total_file, allTime);
	}
	printf("\n");
	fclose(ofp);					// close output file

	// -------------------------------------------------------------------------
	//  Release space
	// -------------------------------------------------------------------------
	delete[] knn_dist; knn_dist = NULL;
	for (i = 0; i < qn; i++) {
		delete[] R[i]; R[i] = NULL;
	}
	delete[] R; R = NULL;
	g_memory -= (SIZEFLOAT * (qn + 1) * maxk);

	delete[] buffer; buffer = NULL;
	delete[] data; data = NULL;
	delete[] query; query = NULL;
	g_memory -= (SIZEFLOAT * d * 2 + SIZECHAR * B);

	//printf("memory = %.2f MB\n", (float) g_memory / (1024.0f * 1024.0f));
	return 0;
}
