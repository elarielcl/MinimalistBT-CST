#ifndef BLOCKTREE_PLEAVEBLOCK_H
#define BLOCKTREE_PLEAVEBLOCK_H

#include "Block.h"

class LeafBlock : public Block {
public:

    LeafBlock(Block*, int64_t, int64_t, std::string&);
    ~LeafBlock();

    int64_t size();

    int add_rank_select_support(int);

    int access(int);
    int rank(int, int);
    int select(int, int);

    int add_leaf_rank_select_support();
    int leaf_rank(int);
    int leaf_select(int);

    int add_search_support();
};

#endif //BLOCKTREE_PLEAVEBLOCK_H
