
#ifndef MINIMALISTBLOCKTREES_BTCT_H
#define MINIMALISTBLOCKTREES_BTCT_H

#include <compressed/CBitBlockTree.h>

class BTCT: public CBitBlockTree {

public:
    int n_;

    sdsl::int_vector<>* bt_first_level_prefix_leaf_ranks_;

    std::vector<sdsl::int_vector<>*> bt_leaf_ranks_;
    std::vector<sdsl::int_vector<>*> bt_second_leaf_ranks_;
    std::vector<sdsl::bit_vector*> bt_starts_with_end_leaf_;
    std::vector<sdsl::bit_vector*> bt_suffix_starts_with_end_leaf_;


    std::vector<sdsl::int_vector<>*> bt_min_excess_; //Starts with the local min of first level
    std::vector<sdsl::int_vector<>*> bt_min_excess_back_block_;
    std::vector<sdsl::bit_vector*> bt_min_in_first_block_;


    sdsl::int_vector<>* top_excess_;
    sdsl::int_vector<>* top_min_excess_;
    int n_last_level_;
    int n_pre_last_level_;
    int n_internal_nodes_;


    BTCT(BlockTree*, int);
    ~BTCT();

    int first_child(int);
    int tree_depth(int);
    int next_sibling(int);
    int parent(int);
    int level_ancestor(int, int);
    int lca(int, int);
    bool is_leaf(int);
    bool is_leaf_rank(int, int&);
    int lb(int);
    int rb(int);

    int leaf_select(int);
    int leaf_rank(int);

    int fwdsearch(int,int);
    int fwdsearch(int,int,int,int,int, int, int&);
    int next_block_fwdsearch(int,int,int&);

    int bwdsearch(int,int);
    int bwdsearch(int,int,int,int,int, int, int&);
    int next_block_bwdsearch(int,int,int&);

    int min_excess(int,int);
    int min_excess(int,int,int,int, int, int&);
    int next_block_min_excess(int,int,int&);
};

int64_t encoded_excess(int64_t);
int64_t decoded_excess(int64_t);

#endif //MINIMALISTBLOCKTREES_BTCT_H
