#include "headers.h"

// -----------------------------------------------------------------------------
MinK_List::MinK_List(				// constructor (given max size)
	int max)							// max size
{
	num_ = 0;
	k_ = max;
	list_ = new Result[max + 1];
}

// -----------------------------------------------------------------------------
MinK_List::~MinK_List() 			// destructor
{
	if (list_ != NULL) {
		delete[] list_; list_ = NULL;
	}
}

// -----------------------------------------------------------------------------
bool MinK_List::isFull()			// is full?
{
	if (num_ >= k_) return true;
	else return false;
}

// -----------------------------------------------------------------------------
float MinK_List::insert(			// insert item (inline for speed)
	float key,							// key of item
	int id)								// id of item
{
	int i = 0;
	for (i = num_; i > 0; --i) {
		if (key < list_[i-1].key_) list_[i] = list_[i - 1];
		else break;
	}
	list_[i].key_ = key;			// store new item here
	list_[i].id_ = id;
	if (num_ < k_) num_++;			// increase the number of items

	return max_key();
}
