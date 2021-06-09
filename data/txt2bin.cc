#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstring>
#include <algorithm>

#include <unistd.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>

using namespace std;

// -----------------------------------------------------------------------------
void create_dir(					// create directory
	char *path)							// input path
{
	int len = (int) strlen(path);
	for (int i = 0; i < len; ++i) {
		if (path[i] == '/') {
			char ch = path[i + 1];
			path[i + 1] = '\0';

			int ret = access(path, F_OK);
			if (ret != 0) {
				ret = mkdir(path, 0755);
				if (ret != 0) printf("Could not create %s\n", path);
			}
			path[i + 1] = ch;
		}
	}
}

// -----------------------------------------------------------------------------
template<class DType>
int convert(						// convert text file to binary file
	bool  sign,							// Yes-Data; No-Query
	int   n,							// cardinality
	int   d,							// dimensionality
	const char *ipath,					// input  path
	const char *opath,					// output path
	int   min,							// min coordinate of data
	int   max)							// max coordinate of data
{
	FILE *fp = NULL;
	uint64_t size = (uint64_t) n*d;
	uint32_t max_uint = 4294967295U;

	DType *data = new DType[size];
	char ifname[200], ofname[200];
	if (sign) {
		strcpy(ifname, ipath); strcat(ifname, ".ds");
		strcpy(ofname, opath); strcat(ofname, ".ds");
	} else {
		strcpy(ifname, ipath); strcat(ifname, ".q");
		strcpy(ofname, opath); strcat(ofname, ".q");
	}

	// -------------------------------------------------------------------------
	//  read text file
	// -------------------------------------------------------------------------
	fp = fopen(ifname, "r");
	if (!fp) { printf("Could not open %s\n", ifname); return 1; }

	int i = 0, j = 0; double tmp = -1.0; DType tmp_d = 0;

	while (!feof(fp) && i < n) {
		fscanf(fp, "%d", &j);
		for (j = 0; j < d; ++j) {
			fscanf(fp, " %lf", &tmp); tmp_d = (DType) tmp;
			data[(uint64_t) i*d + j] = tmp_d;
			assert(tmp_d >= min && tmp_d <= max);
		}
		fscanf(fp, "\n");
		++i;
		if (i % 10000 == 0) printf("i = %d\n", i);
	}
	assert(feof(fp) && i == n);
	fclose(fp);

	// -------------------------------------------------------------------------
	//  verification the first data object
	// -------------------------------------------------------------------------
	for (i = 0; i < 30; ++i) cout << (int) data[i] << " "; cout << endl;

	// -------------------------------------------------------------------------
	//  write binary file
	// -------------------------------------------------------------------------
	fp = fopen(ofname, "wb");
	if (!fp) { printf("Could not open %s\n", ofname); return 1; }

	if (size < max_uint) {
		printf("n=%d, d=%d, size=%lu\n", n, d, size);
		fwrite(data, sizeof(DType), (size_t) size, fp);
	} else {
		for (i=0; i<d; ++i) fwrite(&data[(uint64_t)i*n],sizeof(DType),n,fp);
	}
	fclose(fp);
	delete[] data;

	return 0;
}

// -----------------------------------------------------------------------------
int main(int nargs, char** args)
{
	int  n = -1, d = -1, qn = -1, min = -1, max = -1;
	char dtype[20], ipath[200], opath[200];

	// -------------------------------------------------------------------------
	//  read parameters
	// -------------------------------------------------------------------------
	n  = atoi(args[1]);
	d  = atoi(args[2]);
	qn = atoi(args[3]);
	strncpy(dtype, args[4], sizeof(dtype));
	strncpy(ipath, args[5], sizeof(ipath));
	strncpy(opath, args[6], sizeof(opath));
	min = atoi(args[7]);
	max = atoi(args[8]);
	create_dir(opath);

	printf("n         = %d\n", n);
	printf("d         = %d\n", d);
	printf("qn        = %d\n", qn);
	printf("data type = %s\n", dtype);
	printf("in path   = %s\n", ipath);
	printf("out path  = %s\n", opath);

	// -------------------------------------------------------------------------
	//  convert text file to binary
	// -------------------------------------------------------------------------
	if (strcmp(dtype, "int8") == 0) {
		convert<int8_t>(true,  n,  d, ipath, opath, min, max);
		convert<int8_t>(false, qn, d, ipath, opath, min, max);
	}
	else if (strcmp(dtype, "int16") == 0) {
		convert<int16_t>(true,  n,  d, ipath, opath, min, max);
		convert<int16_t>(false, qn, d, ipath, opath, min, max);
	}
	else if (strcmp(dtype, "int32") == 0) {
		convert<int32_t>(true,  n,  d, ipath, opath, min, max);
		convert<int32_t>(false, qn, d, ipath, opath, min, max);
	}
	else if (strcmp(dtype, "int64") == 0) {
		convert<int64_t>(true,  n,  d, ipath, opath, min, max);
		convert<int64_t>(false, qn, d, ipath, opath, min, max);
	}
	else if (strcmp(dtype, "uint8") == 0) {
		convert<uint8_t>(true,  n,  d, ipath, opath, min, max);
		convert<uint8_t>(false, qn, d, ipath, opath, min, max);
	}
	else if (strcmp(dtype, "uint16") == 0) {
		convert<uint16_t>(true,  n,  d, ipath, opath, min, max);
		convert<uint16_t>(false, qn, d, ipath, opath, min, max);
	}
	else if (strcmp(dtype, "uint32") == 0) {
		convert<uint32_t>(true,  n,  d, ipath, opath, min, max);
		convert<uint32_t>(false, qn, d, ipath, opath, min, max);
	}
	else if (strcmp(dtype, "uint64") == 0) {
		convert<uint64_t>(true,  n,  d, ipath, opath, min, max);
		convert<uint64_t>(false, qn, d, ipath, opath, min, max);
	}
	else if (strcmp(dtype, "float32") == 0) {
		convert<float>(true,  n,  d, ipath, opath, min, max);
		convert<float>(false, qn, d, ipath, opath, min, max);
	}
	else if (strcmp(dtype, "double") == 0) {
		convert<double>(true,  n,  d, ipath, opath, min, max);
		convert<double>(false, qn, d, ipath, opath, min, max);
	}

	return 0;
}
