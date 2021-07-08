#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "def.h"
#include "pri_queue.h"
#include "util.h"
#include "ann.h"

using namespace nns;

// -----------------------------------------------------------------------------
void usage()                        // usage of the package
{
    printf("\n"
        "--------------------------------------------------------------------\n"
        " Usage of External c-k-Approximate Nearest Neighbor Search (c-ANNS):\n"
        "--------------------------------------------------------------------\n"
        "    -alg  (integer)   options of algorithms\n"
        "    -n    (integer)   number of data points\n"
        "    -qn   (integer)   number of queries\n"
        "    -d    (integer)   dimensionality\n"
        "    -B    (integer)   page size\n"
        "    -p    (real)      l_p distance <==> p-stable distr. (0, 2]\n"
        "    -z    (real)      symmetric factor of p-stable distr. [-1, 1]\n"
        "    -c    (real)      approximation ratio (c > 1)\n"
        "    -lf   (integer)   leaf size of kd-tree\n"
        "    -L       (integer)   number of projections (drusilla)\n"
        "    -M    (integer)   number of candidates  (drusilla)\n"
        "    -dt   (string)    data type\n"
        "    -pf   (string)    prefix folder\n"
        "    -df   (string)    data folder to store new format of data\n"
        "    -of   (string)    output folder\n"
        "\n"
        "--------------------------------------------------------------------\n"
        " The options of algorithms (-alg) are:                              \n"
        "--------------------------------------------------------------------\n"
        "    0 - Ground-Truth\n"
        "        Params: -alg 0 -n -qn -d -p -dt -pf\n"
        "\n"
        "    1 - Two Level Indexing of QALSH+\n"
        "        Params: -alg 1 -n -d -B -lf -L -M -p -z -c -dt -pf -df -of\n"
        "\n"
        "    2 - Two Level c-k-ANNS of QALSH+\n"
        "        Params: -alg 2 -qn -d -p -dt -pf -df -of\n"
        "\n"
        "    3 - Indexing of QALSH\n"
        "        Params: -alg 3 -n -d -B -p -z -c -dt -pf -df -of\n"
        "\n"
        "    4 - c-k-ANN Search of QALSH\n"
        "        Params: -alg 4 -qn -d -p -dt -pf -df -of\n"
        "\n"
        "    5 - Linear Scan Search\n"
        "        Params: -alg 5 -n -qn -d -B -p -dt -pf -df -of\n"
        "\n"
        "--------------------------------------------------------------------\n"
        " Author: HUANG Qiang (huangq@comp.nus.edu.sg)                       \n"
        "--------------------------------------------------------------------\n"
        "\n\n\n");
}

// -----------------------------------------------------------------------------
template<class DType>
void interface(                     // interface for calling function
    int   alg,                          // which algorithm
    int   n,                            // number of data points
    int   qn,                           // number of query points
    int   d,                            // dimensionality
    int   B,                            // page size
    int   leaf,                         // leaf size of kd-tree
    int   L,                            // number of projection (drusilla)
    int   M,                            // number of candidates (drusilla)
    float p,                            // p-stable distr. (0,2]
    float zeta,                         // symmetric factor of p-distr. [-1,1]
    float c,                            // approximation ratio
    const char *prefix,                 // prefix of data, query, and truth
    const char *dfolder,                // data folder
    const char *ofolder)                // output folder
{
    assert(alg >= 0 && alg <= 5);

    // read data set, query set, and ground truth file
    gettimeofday(&g_start_time, NULL);
    DType  *data  = NULL;
    DType  *query = NULL;
    Result *truth = NULL;

    if (alg == 0 || alg == 1 || alg == 3) {
        data = new DType[(uint64_t) n*d];
        if (read_data<DType>(n, d, 0, p, prefix, data)) exit(1);
        if (alg == 1 || alg == 3) {
            assert(B > 0);
            write_data_new_form<DType>(n, d, B, (const DType*) data, dfolder);
        }
    }
    if (alg == 0 || alg == 2 || alg == 4 || alg == 5) {
        query = new DType[(uint64_t) qn*d];
        if (read_data<DType>(qn, d, 1, p, prefix, query)) exit(1);
    }
    if (alg == 2 || alg == 4 || alg == 5) {
        truth = new Result[(uint64_t) qn*MAXK];
        if (read_data<Result>(qn, MAXK, 2, p, prefix, truth)) exit(1);
    }
    gettimeofday(&g_end_time, NULL);

    float running_time = g_end_time.tv_sec - g_start_time.tv_sec + 
        (g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
    printf("Load data and query: %f Seconds\n\n", running_time);

    // methods
    switch (alg) {
    case 0:
        ground_truth<DType>(n, qn, d, p, prefix, (const DType*) data, 
            (const DType*) query);
        break;
    case 1:
        indexing_of_qalsh_plus<DType>(n, d, B, leaf, L, M, p, zeta, c,
            (const DType*) data, ofolder);
        break;
    case 2:
        knn_of_qalsh_plus<DType>(qn, d, (const DType*) query, 
            (const Result*) truth, dfolder, ofolder);
        break;
    case 3:
        indexing_of_qalsh<DType>(n, d, B, p, zeta, c, (const DType*) data, 
            ofolder);
        break;
    case 4:
        knn_of_qalsh<DType>(qn, d, (const DType*) query, (const Result*) truth,
            dfolder, ofolder);
        break;
    case 5:
        linear_scan<DType>(n, qn, d, B, p, (const DType*) query, 
            (const Result*) truth, dfolder, ofolder);
        break;
    default:
        printf("Parameters error!\n");
        usage();
    }
    //  release space
    if (alg == 0 || alg == 1 || alg == 3) delete[] data;
    if (alg == 0 || alg == 2 || alg == 4 || alg == 5)delete[] query;
    if (alg == 2 || alg == 4 || alg == 5) delete[] truth; 
}

// -----------------------------------------------------------------------------
int main(int nargs, char **args)
{
    srand(6);                       // use a fixed seed instead of time(NULL)

    int   cnt  = 1;                 // parameter counter
    int   alg  = -1;                // which algorithm
    int   n    = -1;                // number of data points
    int   qn   = -1;                // number of query points
    int   d    = -1;                // dimensionality
    int   B    = -1;                // page size
    float p    = -1.0f;             // p-stable distr. (0,2]
    float zeta = -2.0f;             // symmetric factor of p-distr. [-1,1]
    float c    = -1.0f;             // approximation ratio
    int   leaf = -1;                // leaf size of kd-tree (QALSH+)
    int   L    = -1;                // #projections for drusilla-select (QALSH+)
    int   M    = -1;                // #candidates  for drusilla-select (QALSH+)
    char  dtype[20];                // data type
    char  prefix[200];              // prefix of data, query, and truth set
    char  dfolder[200];             // data folder
    char  ofolder[200];             // output folder

    while (cnt < nargs) {
        if (strcmp(args[cnt], "-alg") == 0) {
            alg = atoi(args[++cnt]); assert(alg >= 0);
            printf("alg     = %d\n", alg);
        }
        else if (strcmp(args[cnt], "-n") == 0) {
            n = atoi(args[++cnt]); assert(n > 0);
            printf("n       = %d\n", n);
        }
        else if (strcmp(args[cnt], "-qn") == 0) {
            qn = atoi(args[++cnt]); assert(qn > 0);
            printf("qn      = %d\n", qn); 
        }
        else if (strcmp(args[cnt], "-d") == 0) {
            d = atoi(args[++cnt]); assert(d > 0);
            printf("d       = %d\n", d); 
        }
        else if (strcmp(args[cnt], "-B") == 0) {
            B = atoi(args[++cnt]); assert(B > 0);
            printf("B       = %d\n", B); 
        }
        else if (strcmp(args[cnt], "-lf") == 0) {
            leaf = atoi(args[++cnt]); assert(leaf > 0);
            printf("leaf    = %d\n", leaf);
        }
        else if (strcmp(args[cnt], "-L") == 0) {
            L = atoi(args[++cnt]); assert(L > 0);
            printf("L       = %d\n", L);
        }
        else if (strcmp(args[cnt], "-M") == 0) {
            M = atoi(args[++cnt]); assert(M > 0);
            printf("M       = %d\n", M);
        }
        else if (strcmp(args[cnt], "-p") == 0) {
            p = (float) atof(args[++cnt]); assert(p > 0 && p <= 2);
            printf("p       = %.1f\n", p);
        }
        else if (strcmp(args[cnt], "-z") == 0) {
            zeta = (float) atof(args[++cnt]); assert(zeta >= -1 && zeta <= 1);
            printf("zeta    = %.1f\n", zeta);
        }
        else if (strcmp(args[cnt], "-c") == 0) {
            c = (float) atof(args[++cnt]); assert(c > 1);
            printf("c       = %.1f\n", c);
        }
        else if (strcmp(args[cnt], "-dt") == 0) {
            strncpy(dtype, args[++cnt], sizeof(dtype));
            printf("dtype   = %s\n", dtype);
        }
        else if (strcmp(args[cnt], "-pf") == 0) {
            strncpy(prefix, args[++cnt], sizeof(prefix));
            printf("prefix  = %s\n", prefix);
        }
        else if (strcmp(args[cnt], "-df") == 0) {
            strncpy(dfolder, args[++cnt], sizeof(dfolder));
            int len = (int) strlen(dfolder);
            if (dfolder[len-1]!='/') { dfolder[len]='/'; dfolder[len+1]='\0'; }
            printf("dfolder = %s\n", dfolder);
            create_dir(dfolder);
        }
        else if (strcmp(args[cnt], "-of") == 0) {
            strncpy(ofolder, args[++cnt], sizeof(ofolder));
            int len = (int) strlen(ofolder);
            if (ofolder[len-1]!='/') { ofolder[len]='/'; ofolder[len+1]='\0'; }
            printf("ofolder = %s\n", ofolder);
            create_dir(ofolder);
        }
        else {
            printf("Parameters error!\n"); usage(); exit(1);
        }
        ++cnt;
    }
    printf("\n");

    if (strcmp(dtype, "uint8") == 0) {
        interface<uint8_t>(alg, n, qn, d, B, leaf, L, M, p, zeta, c, prefix, 
            dfolder, ofolder);
    }
    else if (strcmp(dtype, "uint16") == 0) {
        interface<uint16_t>(alg, n, qn, d, B, leaf, L, M, p, zeta, c, prefix, 
            dfolder, ofolder);
    }
    else if (strcmp(dtype, "int32") == 0) {
        interface<int>(alg, n, qn, d, B, leaf, L, M, p, zeta, c, prefix, 
            dfolder, ofolder);
    }
    else if (strcmp(dtype, "float32") == 0) {
        interface<float>(alg, n, qn, d, B, leaf, L, M, p, zeta, c, prefix, 
            dfolder, ofolder);
    }
    else {
        printf("Parameters error!\n"); usage();
    }
    return 0;
}
