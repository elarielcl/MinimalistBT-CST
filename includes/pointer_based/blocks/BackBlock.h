#ifndef BLOCKTREE_PBACKBLOCK_H
#define BLOCKTREE_PBACKBLOCK_H

#include "Block.h"

class BackBlock : public Block {
public:



    BackBlock(Block*, int64_t, int64_t, std::string&, Block*, Block*, int);
    ~BackBlock();

    int access(int);
    int add_rank_select_support(int);

    int rank(int, int);
    int select(int, int);


    int add_leaf_rank_select_support();

    int leaf_rank(int);
    int leaf_select(int);

    int add_search_support();
};

#endif //BLOCKTREE_PBACKBLOCK_H
