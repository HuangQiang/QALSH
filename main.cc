#include "headers.h"

// -----------------------------------------------------------------------------
void usage() 						// display the usage of QALSH+
{
	printf("\n"
		"--------------------------------------------------------------------\n"
		" Usage of the Package for c-k-Approximate Nearest Neighbor Search   \n"
		"--------------------------------------------------------------------\n"
		"    -alg  (integer)   options of algorithms\n"
		"    -n    (integer)   cardinality of the dataset\n"
		"    -d    (integer)   dimensionality of the dataset\n"
		"    -leaf (integer)   leaf size of kd-tree\n"
		"    -L	   (integer)   drusilla select: l (number of projections)\n"
		"    -M    (integer)   drusilla select: m (number of candidates)\n"
		"    -nb   (integer)   number of blocks to search\n"
		"    -p    (real)      L_{p} Norm, where p in (0, 2]\n"
		"    -z    (real)      symmetric factor of p-stable distr., [-1, 1].\n"
		"    -c    (real)      approximation ratio (c > 1)\n"
		"    -ds   (string)    address of data  set\n"
		"    -qs   (string)    address of query set\n"
		"    -ts   (string)    address of truth set\n"
		"    -op   (string)    output folder\n"
		"\n"
		"--------------------------------------------------------------------\n"
		" The Options of Algorithms (-alg) are:                              \n"
		"--------------------------------------------------------------------\n"
		"    0 - Ground-Truth\n"
		"        Params: -alg 0 -n -qn -d -p -ds -qs -ts\n"
		"\n"
		"    1 - QALSH+\n"
		"        Params: -alg 1 -n -qn -d -leaf -L -M -nb -p -z -c -ds -qs -ts -op\n"
		"\n"
		"    2 - QALSH\n"
		"        Params: -alg 2 -n -qn -d -p -z -c -ds -qs -ts -op\n"
		"\n"
		"    3 - Linear-Scan Method\n"
		"        Params: -alg 3 -n -qn -d -p -ds -qs -ts -op\n"
		"\n"
		"--------------------------------------------------------------------\n"
		" Author: Qiang Huang (huangq2011@gmail.com)                         \n"
		"--------------------------------------------------------------------\n"
		"\n\n\n");
}

// -----------------------------------------------------------------------------
int main(int nargs, char **args)
{
	srand(6);						// srand((unsigned)time(NULL));
	// usage();

	char   data_set[200];			// address of data  set
	char   query_set[200];			// address of query set
	char   truth_set[200];			// address of truth set
	char   out_path[200];			// output path

	int    alg     = -1;			// option of algorithm
	int    n       = -1;			// number of data objects
	int    qn      = -1;			// number of query objects
	int    d       = -1;			// dimensionality
	int    B       = -1;			// page size
	int    leaf    = -1;	 		// leaf size of kd-tree
	int    L       = -1;			// number of projection (drusilla)
	int    M       = -1;			// number of candidates (drusilla)
	int    nb      = -1;			// number of blocks to search
	float  p       = -1.0f;			// L_{p} Norm, where p \in (0,2]
	float  zeta    = -2.0f;			// symmetric factor of p-stable distr. [-1,1]
	float  ratio   = -1.0f;			// approximation ratio

	float  **data  = NULL;			// data set
	float  **query = NULL;			// query set
	Result **R     = NULL;			// k-NN ground truth
	bool   failed  = false;
	int    cnt     = 1;

	while (cnt < nargs && !failed) {
		if (strcmp(args[cnt], "-alg") == 0) {
			alg = atoi(args[++cnt]);
			printf("alg       = %d\n", alg);
			if (alg < 0 || alg > 3) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-n") == 0) {
			n = atoi(args[++cnt]);
			printf("n         = %d\n", n);
			if (n <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-qn") == 0) {
			qn = atoi(args[++cnt]);
			printf("qn        = %d\n", qn);
			if (qn <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-d") == 0) {
			d = atoi(args[++cnt]);
			printf("d         = %d\n", d);
			if (d <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-leaf") == 0) {
			leaf = atoi(args[++cnt]);
			printf("leaf      = %d\n", leaf);
			if (leaf <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-L") == 0) {
			L = atoi(args[++cnt]);
			printf("L         = %d\n", L);
			if (L <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-M") == 0) {
			M = atoi(args[++cnt]);
			printf("M         = %d\n", M);
			if (M <= 0) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-p") == 0) {
			p = (float)atof(args[++cnt]);
			printf("p         = %.1f\n", p);
			if (p <= 0.0f || p > 2.0f) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-z") == 0) {
			zeta = (float)atof(args[++cnt]);
			printf("zeta      = %.1f\n", zeta);
			if (zeta < -1.0f || zeta > 1.0f) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-c") == 0) {
			ratio = (float) atof(args[++cnt]);
			printf("c         = %.1f\n", ratio);
			if (ratio <= 1.0f) {
				failed = true;
				break;
			}
		}
		else if (strcmp(args[cnt], "-ds") == 0) {
			strncpy(data_set, args[++cnt], sizeof(data_set));
			printf("data_set  = %s\n", data_set);
		}
		else if (strcmp(args[cnt], "-qs") == 0) {
			strncpy(query_set, args[++cnt], sizeof(query_set));
			printf("query_set = %s\n", query_set);
		}
		else if (strcmp(args[cnt], "-ts") == 0) {
			strncpy(truth_set, args[++cnt], sizeof(truth_set));
			printf("truth_set = %s\n", truth_set);
		}
		else if (strcmp(args[cnt], "-op") == 0) {
			strncpy(out_path, args[++cnt], sizeof(out_path));
			printf("out_path  = %s\n", out_path);

			int len = (int)strlen(out_path);
			if (out_path[len - 1] != '/') {
				out_path[len] = '/';
				out_path[len + 1] = '\0';
			}
			create_dir(out_path);
		}
		else {
			failed = true;
			usage();
			break;
		}
		cnt++;
	}
	printf("\n");

	// -------------------------------------------------------------------------
	//  read data set, query set, and ground truth file
	// -------------------------------------------------------------------------
	data = new float*[n];
	for (int i = 0; i < n; ++i) data[i] = new float[d];
	if (read_data(n, d, data_set, data) == 1) {
		return 1;
	}

	query = new float*[qn];
	for (int i = 0; i < qn; ++i) query[i] = new float[d];
	if (read_data(qn, d, query_set, query) == 1) {
		return 1;
	}

	if (alg > 0) {
		R = new Result*[qn];
		for (int i = 0; i < qn; ++i) R[i] = new Result[MAXK];
		if (read_ground_truth(qn, truth_set, R) == 1) {
			return 1;
		}
	}
	
	// -------------------------------------------------------------------------
	//  methods
	// -------------------------------------------------------------------------
	switch (alg) {
	case 0:
		ground_truth(n, qn, d, p, (const float **) data, (const float **) query, 
			truth_set);
		break;
	case 1:
		qalsh_plus(n, qn, d, leaf, L, M, p, zeta, ratio, (const float **) data, 
			(const float **) query, (const Result **) R, out_path);
		break;
	case 2:
		qalsh(n, qn, d, p, zeta, ratio, (const float **) data, 
			(const float **) query, (const Result **) R, out_path);
		break;
	case 3:
		linear_scan(n, qn, d, p, (const float **) data, (const float **) query, 
			(const Result **) R, out_path);
		break;
	default:
		printf("Parameters error!\n");
		usage();
		break;
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	for (int i = 0; i < n; ++i) {
		delete[] data[i]; data[i] = NULL;
	}
	delete[] data; data  = NULL;

	for (int i = 0; i < qn; ++i) {
		delete[] query[i]; query[i] = NULL;
	}
	delete[] query; query = NULL;

	if (alg > 0) {
		for (int i = 0; i < qn; ++i) {
			delete[] R[i]; R[i] = NULL;
		}
		delete[] R; R = NULL;
	}

	return 0;
}
