#include "headers.h"

// -----------------------------------------------------------------------------
//	KD_Rect: orthogonal rectangle for bounding rectangle of kd-tree
// -----------------------------------------------------------------------------
KD_Rect::KD_Rect(					// constructor
	int dim,							// dimension
	float l,							// lower bound
	float h)							// higher bound
{
	low_ = new float[dim];
	high_ = new float[dim];

	memset(low_,  l, dim * SIZEFLOAT);
	memset(high_, h, dim * SIZEFLOAT);
}

// -----------------------------------------------------------------------------
KD_Rect::KD_Rect(					// copy constructor
	int dim,							// dimension
	const KD_Rect &rect)				// copy item
{
	low_  = new float[dim];
	high_ = new float[dim];
	for (int i = 0; i < dim; ++i) {
		low_[i]  = rect.low_[i];
		high_[i] = rect.high_[i];
	}
}

// -----------------------------------------------------------------------------
KD_Rect::KD_Rect(					// constrcutor
	int dim,							// dimension
	const float *low,					// lower corner point
	const float *high)					// higher corner point
{
	low_  = new float[dim];
	high_ = new float[dim];
	for (int i = 0; i < dim; ++i) {
		low_[i]  = low[i];
		high_[i] = high[i];
	}
}

// -----------------------------------------------------------------------------
KD_Rect::~KD_Rect()					// destructor
{
	if (low_ != NULL) {
		delete[] low_; low_ = NULL;
	}
	if (high_ != NULL) {
		delete[] high_; high_ = NULL;
	}
}

// -----------------------------------------------------------------------------
bool KD_Rect::inside(				// whether a point inside the rectangle
	int   dim,							// dimension
	const float *point)					// input point
{
	for (int i = 0; i < dim; ++i) {
		if (point[i] < low_[i] || point[i] > high_[i]) return false;
	}
	return true;
}
