#include "headers.h"

// -----------------------------------------------------------------------------
int ResultComp(						// compare function for qsort (ascending)
	const void *e1,						// 1st element
	const void *e2)						// 2nd element
{
	int ret = 0;
	Result *item1 = (Result*) e1;
	Result *item2 = (Result*) e2;

	if (item1->key_ < item2->key_) {
		ret = -1;
	} 
	else if (item1->key_ > item2->key_) {
		ret = 1;
	} 
	else {
		if (item1->id_ < item2->id_) ret = -1;
		else if (item1->id_ > item2->id_) ret = 1;
	}
	return ret;
}

// -----------------------------------------------------------------------------
int ResultCompDesc(					// compare function for qsort (descending)
	const void *e1,						// 1st element
	const void *e2)						// 2nd element
{
	int ret = 0;
	Result *item1 = (Result*) e1;
	Result *item2 = (Result*) e2;

	if (item1->key_ < item2->key_) {
		ret = 1;
	} 
	else if (item1->key_ > item2->key_) {
		ret = -1;
	} 
	else {
		if (item1->id_ < item2->id_) ret = -1;
		else if (item1->id_ > item2->id_) ret = 1;
	}
	return ret;
}

// -----------------------------------------------------------------------------
void get_lp_filename(				// get file name of L_p norm
	float p,							// the p value of L_p norm
	char  *fname)						// file name (return)
{
	sprintf(fname, "%sL%0.1f", fname, p);
}

// -------------------------------------------------------------------------
int create_dir(						// create directory
	char *path)							// input path
{
	// ---------------------------------------------------------------------
	//  check whether the directory exists. if the directory does not 
	//  exist, we create the directory for each folder)
	// ---------------------------------------------------------------------
	int len = (int)strlen(path);
	for (int i = 0; i < len; ++i) {
		if (path[i] == '/') {
			char ch = path[i + 1];
			path[i + 1] = '\0';

			int ret = access(path, F_OK);
			if (ret != 0) {
				ret = mkdir(path, 0755);
				if (ret != 0) {
					printf("Could not create directory %s\n", path);
					return 1;
				}
			}
			path[i + 1] = ch;
		}
	}
	return 0;
}

// -----------------------------------------------------------------------------
float calc_lp_dist(					// calc L_{p} norm
	int   dim,							// dimension
	float p,							// the p value of Lp norm, p in (0, 2]
	const float *vec1,					// 1st point
	const float *vec2)					// 2nd point
{
	float diff = 0.0f;
	float ret = 0.0f;

	// ---------------------------------------------------------------------
	//  calc L_{0.5} norm
	// ---------------------------------------------------------------------
	if (fabs(p - 0.5f) < FLOATZERO) {
		for (int i = 0; i < dim; ++i) {
			diff = fabs(vec1[i] - vec2[i]);
			ret += sqrt(diff);
		}
		return ret * ret;
	}
	// ---------------------------------------------------------------------
	//  calc L_{1.0} norm
	// ---------------------------------------------------------------------
	else if (fabs(p - 1.0f) < FLOATZERO) {
		for (int i = 0; i < dim; ++i) {
			ret += fabs(vec1[i] - vec2[i]);
		}
		return ret;
	}
	// ---------------------------------------------------------------------
	//  calc L_{2.0} norm
	// ---------------------------------------------------------------------
	else if (fabs(p - 2.0f) < FLOATZERO) {
		for (int i = 0; i < dim; ++i) {
			diff = vec1[i] - vec2[i];
			ret += diff * diff;
		}
		return sqrt(ret);
	}
	// ---------------------------------------------------------------------
	//  calc other L_{p} norm (general way)
	// ---------------------------------------------------------------------
	else {
		for (int i = 0; i < dim; ++i) {
			diff = fabs(vec1[i] - vec2[i]);
			ret += pow(diff, p);
		}
		return pow(ret, 1.0f / p);
	}
}

// -----------------------------------------------------------------------------
int read_set(						// read (data or query) set from disk
	int   n,							// number of data objects
	int   d,							// dimensionality
	const char *set,					// address of dataset
	float **points)						// data or queries (return)
{
	int i = 0;
	int j = 0;
	
	FILE *fp = fopen(set, "r");
	if (!fp) {
		printf("Could not open %s.\n", set);
		return 1;
	}

	i = 0;
	while (!feof(fp) && i < n) {
		fscanf(fp, "%d", &j);
		for (j = 0; j < d; ++j) {
			fscanf(fp, " %f", &points[i][j]);
		}
		fscanf(fp, "\n");

		++i;
	}
	assert(feof(fp) && i == n);
	fclose(fp);

	return 0;
}

// -----------------------------------------------------------------------------
int write_data_new_form(			// write dataset with new format
	int   n,							// cardinality
	int   d,							// dimensionality
	int   B,							// page size
	const float **data,					// data set
	const char *output_path)			// output path
{
	int num = (int) floor((float) B/(d*SIZEFLOAT)); // num of data in one page
	int total_file = (int) ceil((float) n / num); // total number of data file
	if (total_file == 0) return 1;

	// -------------------------------------------------------------------------
	//  check whether the directory exists
	// -------------------------------------------------------------------------
	char data_path[200];
	strcpy(data_path, output_path);
	strcat(data_path, "data/");
	
	create_dir(data_path);
	
	// -------------------------------------------------------------------------
	//  write new format data for qalsh
	// -------------------------------------------------------------------------
	char fname[200];
	char *buffer = new char[B];		// one buffer page size
	for (int i = 0; i < B; ++i) buffer[i] = 0;

	int left  = 0;
	int right = 0;
	for (int i = 0; i < total_file; ++i) {
		get_data_filename(i, data_path, fname);

		// ---------------------------------------------------------------------
		//  write data to buffer
		// ---------------------------------------------------------------------
		left = i * num;
		right = left + num;
		if (right > n) right = n;	
		write_data_to_buffer(d, left, right, data, buffer);

		// ---------------------------------------------------------------------
		//  write one page of data to disk
		// ---------------------------------------------------------------------
		write_buffer_to_page(B, fname, (const char *) buffer);
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] buffer; buffer = NULL;
	
	return 0;
}

// -----------------------------------------------------------------------------
void get_data_filename(				// get file name of data
	int   file_id,						// data file id
	const char *data_path,				// path to store data in new format
	char  *fname)						// file name of data (return)
{
	sprintf(fname, "%s%d.data", data_path, file_id);
}

// -----------------------------------------------------------------------------
void write_data_to_buffer(			// write data to buffer
	int   d,							// dimensionality
	int   left,							// left data id
	int   right,						// right data id
	const float **data,					// data set
	char  *buffer)						// buffer to store data
{
	int c = 0;
	for (int i = left; i < right; ++i) {
		for (int j = 0; j < d; ++j) {
			memcpy(&buffer[c], &data[i][j], SIZEFLOAT);
			c += SIZEFLOAT;
		}
	}
}

// -----------------------------------------------------------------------------
int write_buffer_to_page(			// write buffer to one page
	int   B,							// page size
	const char *fname,					// file name of data
	const char *buffer)					// buffer to store data
{
	assert(fname != NULL && buffer != NULL);

	FILE *fp = fopen(fname, "wb");	// open data file to write
	fwrite(buffer, B, 1, fp);
	fclose(fp);

	return 0;
}

// -----------------------------------------------------------------------------
int read_data_new_format(			// read data with new format from disk
	int   id,							// index of data
	int   d,							// dimensionality
	int   B,							// page size
	const char *output_path, 			// output path
	float *data)						// real data (return)
{
	// -------------------------------------------------------------------------
	//  get file name of data
	// -------------------------------------------------------------------------
	char fname[200];
	char data_path[200];

	strcpy(data_path, output_path);
	strcat(data_path, "data/");
									
	int num = (int) floor((float) B/(d*SIZEFLOAT)); // num of data in one page
	int file_id = (int) floor((float) id / num); // data file id
	get_data_filename(file_id, data_path, fname);

	// -------------------------------------------------------------------------
	//  read buffer (one page of data) in new format from disk
	// -------------------------------------------------------------------------
	char *buffer = new char[B];		// allocate one page size
	for (int i = 0; i < B; ++i) {
		buffer[i] = 0;
	}
	read_buffer_from_page(B, fname, buffer);

	// -------------------------------------------------------------------------
	//  read data from buffer
	// -------------------------------------------------------------------------
	int index = id % num;
	read_data_from_buffer(index, d, (const char *) buffer, data);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] buffer; buffer = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
int read_buffer_from_page(			// read buffer from page
	int   B,								// page size
	const char *fname,					// file name of data
	char  *buffer)						// buffer to store data
{
	assert(fname != NULL && buffer != NULL);

	FILE *fp = fopen(fname, "rb");
	fread(buffer, B, 1, fp);
	fclose(fp);

	return 0;
}

// -----------------------------------------------------------------------------
void read_data_from_buffer(			// read data from buffer
	int   index,						// index of data in buffer
	int   d,							// dimensionality
	const char *buffer,					// buffer to store data
	float *data)						// data set
{
	int c = index * d * SIZEFLOAT;
	for (int i = 0; i < d; ++i) {
		memcpy(&data[i], &buffer[c], SIZEFLOAT);
		c += SIZEFLOAT;
	}
}

// -----------------------------------------------------------------------------
int linear(							// linear scan search
	int   n,							// number of data objects
	int   d,							// dimensionality
	int   B,							// page size
	int   p,							// the p value of L_{p} distance
	int   top_k,						// top-k value
	const float *query,					// query object
	const char *data_folder,			// data folder
	MinK_List *list)					// k-NN results (return)
{
	// -------------------------------------------------------------------------
	//  calc <num> and <total_file>, where <num> is the number of data in one 
	//  data file and <total_file> is the total number of data file
	// -------------------------------------------------------------------------
	int num = (int) floor((float) B / (d * SIZEFLOAT));
	int total_file = (int) ceil((float) n / num);
	if (total_file == 0) return 0;

	// -------------------------------------------------------------------------
	//  linear scan method (data in disk)
	//  For each query, we limit that we can ONLY read one page of data
	// -------------------------------------------------------------------------
	char  data_path[200];
	strcpy(data_path, data_folder);
	strcat(data_path, "data/");

	int   id      = 0;
	int   size    = 0;
	char  *buffer = new char[B];	// one page buffer
	float *data   = new float[d];	// one data object

	for (int i = 0; i < total_file; ++i) {
		// ---------------------------------------------------------------------
		//  read one page of data into buffer
		// ---------------------------------------------------------------------
		char fname[200];
		get_data_filename(i, data_path, fname);
		
		read_buffer_from_page(B, fname, buffer);

		// ---------------------------------------------------------------------
		//  linear scan data objects in this page buffer
		// ---------------------------------------------------------------------
		if (i < total_file - 1) size = num;
		else size = n - num * (total_file - 1);

		for (int j = 0; j < size; ++j) {
			read_data_from_buffer(j, d, (const char *)buffer, data);
			float dist = calc_lp_dist(d, p, (const float *) data, query);

			list->insert(dist, id++);
		}
	}
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] buffer; buffer = NULL;
	delete[] data; data = NULL;
	
	return total_file;
}