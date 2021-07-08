#pragma once

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>

#include "def.h"
#include "util.h"
#include "block_file.h"
#include "b_node.h"

namespace nns {

class BlockFile;
class BNode;

// -----------------------------------------------------------------------------
//  BTree: b-tree to index hash tables produced by qalsh
// -----------------------------------------------------------------------------
class BTree {
public:
    int   root_;                    // disk address of root
    BNode *root_ptr_;               // pointer of root
    BlockFile *file_;               // file in disk to store
    
    // -------------------------------------------------------------------------
    BTree();                        // default constructor
    ~BTree();                       // destructor

    // -------------------------------------------------------------------------
    void init(                      // init a new b-tree
        int   b_length,                 // block length
        const char *fname);             // file name    

    // -------------------------------------------------------------------------
    void init_restore(              // load an exist b-tree
        const char *fname);             // file name

    // -------------------------------------------------------------------------
    int bulkload(                   // bulkload b-tree from hash table in mem
        int   n,                        // number of entries
        const Result *table);           // hash table

protected:
    // -------------------------------------------------------------------------
    inline int read_header(const char *buf) {// read root_ from buffer
        memcpy(&root_, buf, sizeof(int));
        return sizeof(int);
    }

    // -------------------------------------------------------------------------
    inline int write_header(char *buf) {// write root_ into buffer
        memcpy(buf, &root_, sizeof(int));
        return sizeof(int);
    }

    // -------------------------------------------------------------------------
    void load_root();               // load root of b-tree

    // -------------------------------------------------------------------------
    void delete_root();             // delete root of b-tree
};

} // end namespace nns
