#ifndef __DEF_H
#define __DEF_H

// -----------------------------------------------------------------------------
//  Typedefs
// -----------------------------------------------------------------------------
typedef char Block[];

// -----------------------------------------------------------------------------
//  Macros
// -----------------------------------------------------------------------------
#define MIN(a, b)	(((a) < (b)) ? (a) : (b))
#define MAX(a, b)	(((a) > (b)) ? (a) : (b))

#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2

// -----------------------------------------------------------------------------
//  Constants
// -----------------------------------------------------------------------------
const float E  = 2.7182818F;		// math constants
const float PI = 3.141592654F;

const float FLOATZERO = 1e-6F;		// accuracy

const long long MAXMEMORY = 1073741823;// max memory, 1 GB

									// max real value
const float MAXREAL = 3.402823466e+38F;
const float MINREAL = -MAXREAL;		// min real value

const int MAXINT = 2147483647;		// max integer value
const int MININT = -MAXINT;			// min integer value

const int SIZEINT   = (int) sizeof(int);
const int SIZECHAR  = (int) sizeof(char);
const int SIZEFLOAT = (int) sizeof(float);
const int SIZEBOOL  = (int) sizeof(bool);

									// file header size
const int BFHEAD_LENGTH = (int) (sizeof(int) * 2);
									// index size of leaf node
const int INDEX_SIZE_LEAF_NODE = 4096;
const int MAXK = 100;				// max top-k value

#endif
