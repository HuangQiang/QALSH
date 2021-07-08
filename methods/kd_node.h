#pragma once

#include <iostream>
#include <algorithm>
#include <vector>

#include "def.h"
#include "util.h"
#include "pri_queue.h"

namespace nns {

// -----------------------------------------------------------------------------
template<class DType>
class KD_Tree;

// -----------------------------------------------------------------------------
//  KD_Node: node of kd-tree (abstract class)
//  There are two types of KD_Node: KD_Leaf and KD_Split
// -----------------------------------------------------------------------------
template<class DType>
class KD_Node {
public:
    virtual ~KD_Node() {}           // virtual destructor

    virtual void search(            // tree search
        float box_dist,                 // box distance to query
        float c,                        // approximation ratio
        const DType *query,             // query point
        MinK_List *list) = 0;           // k-NN results (return)

    virtual void traversal(         // traversal kd-tree
        std::vector<int> &leaf_size) = 0;// leaf size (return)

    friend class KD_Tree<DType>;    // allow kd-tree to access
};


// -----------------------------------------------------------------------------
//    KD_Leaf: leaf node of kd-tree
// -----------------------------------------------------------------------------
template<class DType>
class KD_Leaf: public KD_Node<DType> {
public:
    KD_Leaf(                        // constructor
        int   n,                        // number of data points
        int   d,                        // dimension of data points
        int   *index,                   // data index
        const DType *data);             // data points

    virtual ~KD_Leaf();             // destructor

    virtual void search(            // tree search
        float box_dist,                 // box distance to query
        float c,                        // approximation ratio
        const DType *query,             // query point
        MinK_List *list);               // k-NN results (return)

    virtual void traversal(         // traversal kd-tree
        std::vector<int> &leaf_size);   // leaf size (return)

protected:
    int   n_pts_;                   // number of data points
    int   dim_;                     // dimension of data points
    int   *index_;                  // data index
    const DType *data_;             // data points
};

// -----------------------------------------------------------------------------
template<class DType>
KD_Leaf<DType>::KD_Leaf(            // constructor
    int   n,                            // number of data points
    int   d,                            // dimension of data points
    int   *index,                       // data index
    const DType *data)                  // data points
    : n_pts_(n), dim_(d), index_(index), data_(data)
{
}

// -----------------------------------------------------------------------------
template<class DType>
KD_Leaf<DType>::~KD_Leaf()          // destructor
{
}

// -----------------------------------------------------------------------------
template<class DType>
void KD_Leaf<DType>::search(        // tree search
    float box_dist,                     // box distance to query
    float c,                            // approximation ratio
    const DType *query,                 // query point
    MinK_List *list)                    // k-NN results (return)
{
    for (int i = 0; i < n_pts_; ++i) {
        int id = index_[i];
        const DType *point = &data_[(uint64_t)id*dim_];
        
        // NOTE: calc the l2-dist sqr
        float dist = 0.0F;
        for (int j = 0; j < dim_; ++j) {
            dist += (float) SQR(point[j] - query[j]);
        }
        list->insert(dist, id);
    }
}

// -----------------------------------------------------------------------------
template<class DType>
void KD_Leaf<DType>::traversal(     // traversal kd-tree
    std::vector<int> &leaf_size)        // leaf size (return)
{
    leaf_size.push_back(n_pts_);
}


// -----------------------------------------------------------------------------
//  KD_Split: split node of kd-tree
// -----------------------------------------------------------------------------
template<class DType>
class KD_Split : public KD_Node<DType> {
public:
    KD_Split(                       // constructor
        int   dim,                      // dimension of data points
        int   cut_dim,                  // cutting dimension
        DType cut_val,                  // cutting value
        DType low_val,                  // low  bound in <dim>
        DType high_val,                 // high bound in <dim>
        KD_Node<DType> *left,           // left child
        KD_Node<DType> *right,          // right child
        const DType *data);             // data points

    virtual ~KD_Split();            // destructor

    virtual void search(            // tree search
        float box_dist,                 // box distance to query
        float c,                        // approximation ratio
        const DType *query,             // query point
        MinK_List *list);               // k-NN results (return)

    virtual void traversal(         // traversal kd-tree
        std::vector<int> &leaf_size);   // leaf size (return)

protected:
    int   dim_;                     // dimension of data points
    int   cut_dim_;                 // cutting dimension
    DType cut_val_;                 // cutting value
    DType cd_bnds_[2];              // cutting bounds
    KD_Node<DType> *child_[2];      // children of node
    const DType *data_;             // data points
};

// -----------------------------------------------------------------------------
template<class DType>
KD_Split<DType>::KD_Split(          // constructor
    int   dim,                          // dimension of data points
    int   cut_dim,                      // cutting dimension
    DType cut_val,                      // cutting value
    DType low_val,                      // low  bound in <dim>
    DType high_val,                     // high bound in <dim>
    KD_Node<DType> *left,               // left child
    KD_Node<DType> *right,              // right child
    const DType *data)                  // data points
    : dim_(dim), cut_dim_(cut_dim), cut_val_(cut_val), data_(data)
{
    cd_bnds_[0] = low_val;
    cd_bnds_[1] = high_val;
    child_[0]   = left;
    child_[1]   = right;
}

// -----------------------------------------------------------------------------
template<class DType>
KD_Split<DType>::~KD_Split()        // destructor
{
    if (child_[0] != NULL) { delete child_[0]; child_[0] = NULL; }
    if (child_[1] != NULL) { delete child_[1]; child_[1] = NULL; }
}

// -----------------------------------------------------------------------------
template<class DType>
void KD_Split<DType>::search(       // tree search
    float box_dist,                     // box distance to query
    float c,                            // approximation ratio
    const DType *query,                 // query point
    MinK_List *list)                    // k-NN results (return)
{
    float cut_diff = (float) (query[cut_dim_] - cut_val_);
    if (cut_diff < 0) {
        // query on left of cutting plane
        child_[0]->search(box_dist, c, query, list);

        // visit right child
        float box_diff = (float) (cd_bnds_[0] - query[cut_dim_]);
        if (box_diff < 0) box_diff = 0;

        box_dist += (cut_diff*cut_diff - box_diff*box_diff);
        if (box_dist * c < list->max_key()) {
            child_[1]->search(box_dist, c, query, list);
        }
    }
    else {
        // query on right of cutting plane
        child_[1]->search(box_dist, c, query, list);

        // visit left child
        float box_diff = (float) (query[cut_dim_] - cd_bnds_[1]);
        if (box_diff < 0) box_diff = 0;
        
        box_dist += (cut_diff*cut_diff - box_diff*box_diff);
        if (box_dist * c < list->max_key()) {
            child_[0]->search(box_dist, c, query, list);
        }
    }
}

// -----------------------------------------------------------------------------
template<class DType>
void KD_Split<DType>::traversal(    // traversal kd-tree
    std::vector<int> &leaf_size)        // leaf size (return)
{
    child_[0]->traversal(leaf_size);
    child_[1]->traversal(leaf_size);
}

} // end namespace nns
