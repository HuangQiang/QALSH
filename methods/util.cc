#include "util.h"

namespace nns {

timeval  g_start_time;              // global param: start time
timeval  g_end_time;                // global param: end   time

float    g_indexing_time = -1.0f;   // global param: indexing time
float    g_estimated_mem = -1.0f;   // global param: estimated memory

float    g_runtime       = -1.0f;   // global param: running time
float    g_ratio         = -1.0f;   // global param: overall ratio
float    g_recall        = -1.0f;   // global param: recall
uint64_t g_page_io       = 0;       // global param: page i/o

// -----------------------------------------------------------------------------
void create_dir(                    // create directory
    char *path)                         // input path
{
    int len = (int) strlen(path);
    for (int i = 0; i < len; ++i) {
        if (path[i] != '/') continue;

        char ch = path[i+1]; path[i+1] = '\0';
        if (access(path, F_OK) != 0) {
            if (mkdir(path, 0755) != 0) {
                printf("Could not create %s\n", path); exit(1);
            }
        }
        path[i+1] = ch;
    }
}

// -----------------------------------------------------------------------------
int write_buffer_to_page(           // write buffer to one page
    int   B,                            // page size
    const char *fname,                  // file name of data
    const char *buffer)                 // buffer to store data
{
    FILE *fp = fopen(fname, "wb");
    if (!fp) { printf("Could not open %s\n", fname); return 1; }

    fwrite(buffer, sizeof(char), B, fp);
    fclose(fp);
    return 0;
}

// -----------------------------------------------------------------------------
int read_buffer_from_page(          // read buffer from page
    int   B,                            // page size
    const char *fname,                  // file name of data
    char  *buffer)                      // buffer to store data (return)
{
    FILE *fp = fopen(fname, "rb");
    if (!fp) { printf("Could not open %s\n", fname); return 1; }

    fread(buffer, sizeof(char), B, fp);
    fclose(fp);
    return 0;
}

// -----------------------------------------------------------------------------
int write_ground_truth(             // write ground truth to disk
    int   n,                            // number of ground truth results
    int   d,                            // dimension of ground truth results
    float p,                            // l_p distance
    const char *prefix,                 // prefix of truth set
    const Result *truth)                // ground truth
{
    char fname[200]; sprintf(fname, "%s.gt%3.1f", prefix, p);
    FILE *fp = fopen(fname, "wb");
    if (!fp) { printf("Could not create %s\n", fname); return 1; }
    
    uint64_t size = (uint64_t) n*d;
    fwrite(truth, sizeof(Result), size, fp);
    fclose(fp);
    return 0;
}

// -----------------------------------------------------------------------------
float calc_ratio(                   // calc overall ratio [1,\infinity)
    int   k,                            // top-k value
    const Result *truth,                // ground truth results 
    MinK_List *list)                    // top-k approximate results
{
    float ratio = 0.0f;
    for (int i = 0; i < k; ++i) {
        ratio += list->ith_key(i) / truth[i].key_;
    }
    return ratio / k;
}

// -----------------------------------------------------------------------------
float calc_recall(                  // calc recall (percentage)
    int   k,                            // top-k value
    const Result *truth,                // ground truth results 
    MinK_List *list)                    // top-k approximate results
{
    int i = k - 1;
    int last = k - 1;
    while (i >= 0 && list->ith_key(i) > truth[last].key_) {
        i--;
    }
    return (i+1)*100.0f / k;
}

} // end namespace nns
