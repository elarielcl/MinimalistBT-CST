
#ifndef MINIMALISTBLOCKTREES_BTCT_H
#define MINIMALISTBLOCKTREES_BTCT_H

#include <compressed/CBitBlockTree.h>

class BTCT: public CBitBlockTree {

public:

    sdsl::int_vector<>* bt_first_level_prefix_leaf_ranks_;

    std::vector<sdsl::int_vector<>*> bt_leaf_ranks_;
    std::vector<sdsl::int_vector<>*> bt_second_leaf_ranks_;
    std::vector<sdsl::bit_vector*> bt_starts_with_end_leaf_;
    std::vector<sdsl::bit_vector*> bt_suffix_starts_with_end_leaf_;

    BTCT(BlockTree*, int);
    ~BTCT();

    int leaf_select(int);
    int leaf_rank(int);

};

#endif //MINIMALISTBLOCKTREES_BTCT_H
