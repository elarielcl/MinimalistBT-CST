#ifndef BLOCKTREE_PLAZYINTERNALBLOCK_H
#define BLOCKTREE_PLAZYINTERNALBLOCK_H

#include "Block.h"

class InternalBlock : public Block {
public:

    InternalBlock(Block*, int64_t, int64_t, std::string&);
    ~InternalBlock();

    std::vector<Block*>& children(int, int);
    void clean_unnecessary_expansions();

    bool is_leaf();
    int access(int);
    int add_rank_select_support(int);

    int rank(int, int);
    int select(int, int);

    int add_leaf_rank_select_support();

    int leaf_rank(int);
    int leaf_select(int);

    int add_search_support();
};

#endif //BLOCKTREE_PLAZYINTERNALBLOCK_H
