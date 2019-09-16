#ifndef __QALSH_H
#define __QALSH_H

#include <vector>
using namespace std;

// -----------------------------------------------------------------------------
//  Query-Aware Locality-Sensitive Hashing (QALSH) is used to solve the problem 
//  of c-Approximate Nearest Neighbor (c-ANN) search.
//
//  the idea was introduced by Qiang Huang, Jianlin Feng, Yikai Zhang, Qiong 
//  Fang, and Wilfred Ng in their paper "Query-aware locality-sensitive hashing 
//  for approximate nearest neighbor search", in Proceedings of the VLDB 
//  Endowment (PVLDB), 9(1), pages 1–12, 2015.
// -----------------------------------------------------------------------------
class QALSH {
public:
	QALSH(							// constructor
		int   n,						// cardinality
		int   d,						// dimensionality
		float p,						// l_p distance
		float zeta,						// a parameter of p-stable distr.
		float ratio,					// approximation ratio
		const float **data);			// data objects

	// -------------------------------------------------------------------------
	~QALSH();						// destructor

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	int knn(						// k-NN search
		int top_k,						// top-k value
		const float *query,				// input query object
		MinK_List *list);				// k-NN results (return)

	// -------------------------------------------------------------------------
	int knn(						// k-NN search
		int top_k,						// top-k value
		float R,						// limited search range
		const float *query,				// input query object
		const vector<int> &object_id,	// object id mapping
		MinK_List *list);				// k-NN results (return)

protected:
	// -------------------------------------------------------------------------
	int    n_pts_;					// cardinality
	int    dim_;					// dimensionality
	float  p_;						// l_p distance
	float  zeta_;					// a parameter of p-stable distr.
	float  appr_ratio_;				// approximation ratio
	const  float **data_;			// data objects

	float  w_;						// bucket width
	float  p1_;						// positive probability
	float  p2_;						// negative probability
	float  alpha_;					// collision threshold percentage
	float  beta_;					// false positive percentage
	float  delta_;					// error probability
	int    m_;						// number of hashtables
	int    l_;						// collision threshold
	float  *a_array_;				// hash functions
	Result **tables_;				// hash tables

	int    *freq_;					// frequency		
	int    *lpos_;					// left  position of hash table
	int    *rpos_;					// right position of hash table
	bool   *checked_;				// whether checked
	bool   *bucket_flag_;			// bucket flag
	bool   *range_flag_;			// range flag
	float  *q_val_;					// hash value of query
	
	// -------------------------------------------------------------------------
	float calc_l0_prob(				// calc <p1> and <p2> for L_{0.5} distance
		float x);						// x = w / (2.0 * r)

	float calc_l1_prob(				// calc <p1> and <p2> for L_{1.0} distance
		float x);						// x = w / (2.0 * r)

	float calc_l2_prob(				// calc <p1> and <p2> for L_{2.0} distance
		float x);						// x = w / (2.0 * r)	

	// -------------------------------------------------------------------------
	float calc_hash_value(			// calc hash value
		int   table_id,					// hash table id
		const float *data);				// one data/query object

	// -------------------------------------------------------------------------
	int binary_search_pos(			// find position by binary search
		int   table_id,					// hash table id
		float value);					// hash value
};

#endif // __QALSH_H
