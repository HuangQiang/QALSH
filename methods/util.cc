#include "util.h"

namespace nns {

timeval  g_start_time;				// global parameter: start time
timeval  g_end_time;				// global parameter: end time

float    g_indexing_time = -1.0f;	// global parameter: indexing  time
float    g_estimated_mem = -1.0f;	// global parameter: estimated memory

float    g_runtime       = -1.0f;	// global parameter: running time
float    g_ratio         = -1.0f;	// global parameter: overall ratio
float    g_recall        = -1.0f;	// global parameter: recall
uint64_t g_page_io       = 0;		// global parameter: page i/o

// -----------------------------------------------------------------------------
void create_dir(					// create directory
	char *path)							// input path
{
	int len = (int) strlen(path);
	for (int i = 0; i < len; ++i) {
		if (path[i] != '/') continue;

		char ch = path[i + 1]; path[i+1] = '\0';
		if (access(path, F_OK) != 0) {
			if (mkdir(path, 0755) != 0) printf("Could not create %s\n", path);
		}
		path[i+1] = ch;
	}
}

// -----------------------------------------------------------------------------
int write_buffer_to_page(			// write buffer to one page
	int   B,							// page size
	const char *fname,					// file name of data
	const char *buffer)					// buffer to store data
{
	FILE *fp = fopen(fname, "wb");	// open data file to write
	if (!fp) { printf("Could not open %s\n", fname); return 1; }

	fwrite(buffer, SIZECHAR, B, fp);
	fclose(fp);
	return 0;
}

// -----------------------------------------------------------------------------
int read_buffer_from_page(			// read buffer from page
	int   B,							// page size
	const char *fname,					// file name of data
	char  *buffer)						// buffer to store data (return)
{
	FILE *fp = fopen(fname, "rb");
	if (!fp) { printf("Could not open %s\n", fname); return 1; }

	fread(buffer, SIZECHAR, B, fp);
	fclose(fp);
	return 0;
}

// -----------------------------------------------------------------------------
int write_ground_truth(				// write ground truth to disk
	int   n,							// number of ground truth results
	int   d,			 				// dimension of ground truth results
	float p,							// l_p distance
	const char *prefix,					// prefix of truth set
	const Result *truth)				// ground truth
{
	char fname[200]; sprintf(fname, "%s.gt%3.1f", prefix, p);
	FILE *fp = fopen(fname, "wb");
	if (!fp) { printf("Could not create %s\n", fname); return 1; }
	
	uint64_t size = (uint64_t) n*d;
	uint32_t max_uint = 4294967295U;
	
	if (size < max_uint) {
		// if no longer than size_t, write the whole array to disk directly
		fwrite(truth, sizeof(Result), size, fp);
	}
	else {
		// we store ground truth in linear order, n*d == d*n (save multi times)
		for (int i = 0; i < d; ++i) {
			fwrite(&truth[(uint64_t)i*n], sizeof(Result), n, fp);
		}
	}
	fclose(fp);
	return 0;
}

// -----------------------------------------------------------------------------
float calc_ratio(					// calc overall ratio [1,\infinity)
	int   k,							// top-k value
	const Result *truth,				// ground truth results 
	MinK_List *list)					// top-k approximate results
{
	float ratio = 0.0f;
	for (int i = 0; i < k; ++i) {
		ratio += list->ith_key(i) / truth[i].key_;
	}
	return ratio / k;
}

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int   k,							// top-k value
	const Result *truth,				// ground truth results 
	MinK_List *list)					// top-k approximate results
{
	int i = k - 1;
	int last = k - 1;
	while (i >= 0 && list->ith_key(i) > truth[last].key_) {
		i--;
	}
	return (i + 1) * 100.0f / k;
}

} // end namespace nns
