#pragma once

#include <iostream>
#include <cassert>
#include <cmath>
#include <cstring>

#include "def.h"

namespace nns {

// -----------------------------------------------------------------------------
//  NOTE: The author of the implementation of class BlockFile is Yufei Tao.
//  Modified by Qiang HUANG
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//  BlockFile: structure of reading and writing file for b-tree
// -----------------------------------------------------------------------------
class BlockFile {
public:
    FILE *fp_;                      // file pointer
    char fname_[200];               // file name
    bool new_flag_;                 // specifies if this is a new file
    
    int  block_length_;             // length of a block
    int  act_block_;                // block num of fp position
    int  num_blocks_;               // total num of blocks

    // -------------------------------------------------------------------------
    BlockFile(                      // constructor
        int  b_length,                  // length of a block
        const char *name);              // file name

    // -------------------------------------------------------------------------
    ~BlockFile();                   // destructor

    // -------------------------------------------------------------------------
    inline void put_bytes(const char *bytes, int num) { // write num bytes
        fwrite(bytes, sizeof(char), num, fp_); 
    }

    // -------------------------------------------------------------------------
    inline void get_bytes(char *bytes, int num) { // read num bytes
        fread(bytes, sizeof(char), num, fp_); 
    }

    // -------------------------------------------------------------------------
    inline void seek_block(int bnum) { // move fp_ to the right with bnum
        fseek(fp_, (bnum-act_block_)*block_length_, SEEK_CUR);
    }

    // -------------------------------------------------------------------------
    inline bool file_new() { return new_flag_; } // is this block modified?

    // -------------------------------------------------------------------------
    inline int get_blocklength() { return block_length_; }

    // -------------------------------------------------------------------------
    inline int get_num_of_blocks() { return num_blocks_; }

    // -------------------------------------------------------------------------
    inline void fwrite_number(int num) { // write a value (type int)
        put_bytes((char*) &num, sizeof(int));
    }

    // -------------------------------------------------------------------------
    inline int fread_number() {     // read a value (type int)
        char ca[sizeof(int)]; get_bytes(ca, sizeof(int)); return *((int *)ca);
    }

    // -------------------------------------------------------------------------
    void read_header(               // read remain bytes excluding header
        char *buffer);                  // contain remain bytes (return)

    // -------------------------------------------------------------------------
    void set_header(                // set remain bytes excluding header
        const char *buffer);            // contain remain bytes

    // -------------------------------------------------------------------------
    bool read_block(                // read a block in the `index` position
        Block block,                    // a block
        int   index);                   // position of the block

    // -------------------------------------------------------------------------
    bool write_block(               // write a block in the `index` position
        Block block,                    // a block
        int   index);                   // pos of the block

    // -------------------------------------------------------------------------
    int append_block(               // append a block at the end of file
        Block block);                   // a block

    // -------------------------------------------------------------------------
    bool delete_last_blocks(        // delete the last `num` blocks
        int num);                       // number of blocks to be deleted
};

} // end namespace nns
