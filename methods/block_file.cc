#include "block_file.h"
#include "def.h"

namespace nns {

// -----------------------------------------------------------------------------
BlockFile::BlockFile(               // constructor
    int   b_length,                     // block length
    const char *name)                   // file name
{
    strcpy(fname_, name);
    block_length_ = b_length;
    num_blocks_   = 0;

    // -------------------------------------------------------------------------
    //  init fp_ and open file_name_. if file_name_ exists, then fp_ != 0
    //  rb+: read binary data from disk
    // -------------------------------------------------------------------------
    if ((fp_ = fopen(fname_, "rb+")) != 0) {
        // since the file exists, new_flag_ is false
        new_flag_ = false;
        block_length_ = fread_number(); // get block_length_ from header
        num_blocks_   = fread_number(); // get num_blocks_   from header
    }
    else {
        // ---------------------------------------------------------------------
        //  wb+: write binary data to disk
        // ---------------------------------------------------------------------
        assert(block_length_ >= BFHEAD_LENGTH);
        fp_ = fopen(fname_, "wb+");     // construct a new file

        // as file is just constructed (new), new_flag_ is true.
        new_flag_ = true;
        fwrite_number(block_length_);   // write block_length_ to header
        fwrite_number(0);               // write num_blocks_ (0) to header

        // ---------------------------------------------------------------------
        //  ftell() return the number of bytes from current position to the 
        //  start position of the file
        //  since block_length_ >= 8 bytes, init 0 for the remain bytes
        // ---------------------------------------------------------------------
        int length = block_length_ - (int) ftell(fp_);
        if (length > 0) {
            char *buffer = new char[length];
            memset(buffer, 0, length*sizeof(char));
            put_bytes(buffer, length);
            delete[] buffer;
        }
    }
    // -------------------------------------------------------------------------
    //  redirect fp_ to the start position & init act_block_ = 0 (no block)
    // -------------------------------------------------------------------------
    fseek(fp_, 0, SEEK_SET);
    act_block_ = 0;
}

// -----------------------------------------------------------------------------
BlockFile::~BlockFile()             // destructor
{
    if (fp_) fclose(fp_);
}

// -----------------------------------------------------------------------------
//  Note: this func does not read the header of block file. it fetches the info 
//  (the root of b+ tree) in the 1st block excluding the header.
// -----------------------------------------------------------------------------
void BlockFile::read_header(        // read remain bytes excluding header
    char *buffer)                       // buffer with remain bytes (return)
{
    fseek(fp_, BFHEAD_LENGTH, SEEK_SET); // jump out of first 8 bytes
    get_bytes(buffer, block_length_ - BFHEAD_LENGTH); // read remaining bytes

    if (num_blocks_ < 1) {          // no remain bytes
        fseek(fp_, 0, SEEK_SET);    // fp_ return to start pos of this file
        act_block_ = 0;             // no act block
    } 
    else {
        // ---------------------------------------------------------------------
        //  since we read the 1st block (header block), act_block_ = 1. 
        //  after get_bytes(), fp_ is at the 2nd block (1st data block)
        // ---------------------------------------------------------------------
        act_block_ = 1;
    }
}

// -----------------------------------------------------------------------------
//  Note: this func does not write the header of block file. it writes the info 
//  (the root of b+ tree) in the 1st block excluding the header.
// -----------------------------------------------------------------------------
void BlockFile::set_header(         // set remain bytes excluding header
    const char *buffer)                 // buffer with remain bytes
{
    fseek(fp_, BFHEAD_LENGTH, SEEK_SET); // jump out of first 8 bytes
    put_bytes(buffer, block_length_ - BFHEAD_LENGTH); // write remain bytes
    
    if (num_blocks_ < 1) {          // no remain bytes
        fseek(fp_, 0, SEEK_SET);    // fp_ return to start pos of this file
        act_block_ = 0;             // no act block
    }
    else {
        // ---------------------------------------------------------------------
        //  since we write the first block (header block), act_block_ = 1. 
        //  after put_bytes(), fp_ is at the 2nd block (1st data block)
        // ---------------------------------------------------------------------
        act_block_ = 1;
    }
}

// -----------------------------------------------------------------------------
//  we point out the difference of counting among the num_blocks_, act_block_, 
//  and index as follows.
//  (1) num_blocks_: record the number of data blocks, exclude the header block
//      start from 1.
//  (2) act_block_: record the number of total blocks currently read or written,
//      include the header block. (when we read or write file, act_block_ = 1.)
//      act_block_ is corresponding to fp_.
//  (3) index: input parameter, record the position of data block we want to 
//      read or write, exclude the header block. start from 0.
//      act_block_ is larger than index by 1. so we should ++index first.
//
//  For example, if num_blocks_ = 3, there are 4 blocks in this block file: 
//  1 header block + 3 data block. 
// 
//  When the block file is opened, act_block_ = 1 and fp_ is at the start 
//  position of 1st data block. 
//
//  If index == 1, it means that we want to read or write the 2nd data block, 
//  thus firstly index++ (index = 2), then call seek_block() to move fp_ to the
//  start position of the 2nd data block (with index-act_block_ = 2-1 = 1).
//
//  After reading or writing the 2nd data block, fp_ is pointed to the start 
//  position of the 3rd data block. As we know it has read or written 3 blocks, 
//  thus currently act_block_ = index + 1 = 2 + 1 = 3.
// -----------------------------------------------------------------------------
bool BlockFile::read_block(         // read a block from index
    Block block,                        // a block (return)
    int   index)                        // position of this block (start from 0)
{
    ++index; assert(index > 0 && index <= num_blocks_);
    seek_block(index);
    get_bytes(block, block_length_);// read this block

    // update act_block_
    if (index + 1 > num_blocks_) {
        // fp_ reaches the end of this block file, so rewinds to start position
        fseek(fp_, 0, SEEK_SET); act_block_ = 0;
    } else {
        act_block_ = index + 1;     // act_block_ to the next position
    }
    return true;
}

// -----------------------------------------------------------------------------
//  Note: this function can ONLY write to an already "allocated" block (in the 
//  range of num_blocks_).
//  If you allocate a new block, please call append_block() instead.
// -----------------------------------------------------------------------------
bool BlockFile::write_block(        // write a block to index
    Block block,                        // a block
    int index)                          // position of this block (start from 0)
{
    ++index; assert(index > 0 && index <= num_blocks_);
    seek_block(index);
    put_bytes(block, block_length_);// write this block

    // update act_block_
    if (index + 1 > num_blocks_) {
        fseek(fp_, 0, SEEK_SET); act_block_ = 0;
    } else {
        act_block_ = index + 1;
    }
    return true;
}

// -----------------------------------------------------------------------------
int BlockFile::append_block(        // append new block at the end of file
    Block block)                        // the new block
{
    fseek(fp_, 0, SEEK_END);        // fp_ points to the end of file
    put_bytes(block, block_length_);// write a block
    ++num_blocks_;                  // add 1 to num_blocks_
    
    // fp_ points to the position to store num_blocks_ & update header
    fseek(fp_, sizeof(int), SEEK_SET); 
    fwrite_number(num_blocks_);

    // -------------------------------------------------------------------------
    //  fp_ points to the start position of new added block
    //  act_block_ = num_blocks_ indicates that fp_ points to new added block.
    //  return the index of new added block
    // -------------------------------------------------------------------------
    fseek(fp_, -block_length_, SEEK_END);
    return (act_block_ = num_blocks_) - 1;
}

// -----------------------------------------------------------------------------
//  NOTE: we just logically (NOT physically) delete the data. The real data is 
//  still stored in file and the size of file is not changed.
// -----------------------------------------------------------------------------
bool BlockFile::delete_last_blocks( // delete the last `num` blocks
    int num)                            // number of blocks to be deleted
{
    if (num > num_blocks_) return false;

    // only update num_blocks_ & re-write it to disk
    num_blocks_ -= num;
    fseek(fp_, sizeof(int), SEEK_SET);
    fwrite_number(num_blocks_);

    // fp_ rewinds to the start position of this block file
    fseek(fp_, 0, SEEK_SET); act_block_ = 0;
    return true;
}

} // end namespace nns
