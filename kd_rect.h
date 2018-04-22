#ifndef __KD_RECT_H
#define __KD_RECT_H

// -----------------------------------------------------------------------------
//	KD_Rect: orthogonal rectangle for bounding rectangle of kd-tree
// -----------------------------------------------------------------------------
class KD_Rect {
public:
	float *low_;					// rectangle lower bound
	float *high_;					// rectangle upper bound

	// -------------------------------------------------------------------------
	KD_Rect(						// constructor
		int   dim,						// dimension
		float l = 0,					// low  boundary, default is zero
		float h = 0);					// high boundary, default is zero

	// -------------------------------------------------------------------------
	KD_Rect(						// copy constructor
		int   dim,						// dimension
		const KD_Rect &rect);			// rectangle to copy

	// -------------------------------------------------------------------------
	KD_Rect(						// construct from points
		int   dim,						// dimension
		const float *low,				// low point
		const float *high);				// high point

	// -------------------------------------------------------------------------
	~KD_Rect();						// destructor

	// -------------------------------------------------------------------------
	bool inside(					// whether a point inside rectangle
		int   dim,						// dimension
		const float *point);			// one point
};

#endif // __KD_RECT_H