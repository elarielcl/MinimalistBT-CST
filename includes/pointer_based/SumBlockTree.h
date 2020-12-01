#ifndef MINIMALISTBT_CST_SUMBLOCKTREE_H
#define MINIMALISTBT_CST_SUMBLOCKTREE_H

#include "blocktree.utils/RabinKarp.h"
#include "blocktree.utils/HashString.h"
#include "pointer_based/blocks/SumBackBlock.h"
#include "pointer_based/blocks/SumBlock.h"

#include <string>
#include <unordered_map>

class SumBlockTree {
    void block_scan(std::vector<SumBlock*>&, int, std::unordered_map<HashString, std::vector<SumBlock*>>&);
public:
    int r_; // Arity
    int leaf_length_;
    std::basic_string<int64_t> input_; // Input sequence of the Tree
    SumBlock* root_block_;

    SumBlockTree(std::basic_string<int64_t >&, int, int);
    ~SumBlockTree();

    void add_access_support();
    int64_t access(int);

    void process_back_pointers_heuristic();
    void process_back_pointers();
    void clean_unnecessary_expansions();

    void process_level_heuristic(std::vector<SumBlock*>&);
    void process_level(std::vector<SumBlock*>&);

    void forward_window_block_scan(std::vector<SumBlock*>& level, int window_size, int N, std::unordered_map<HashString, std::vector<SumBlock*>>& hashtable);
    void forward_pair_window_block_scan(std::vector<SumBlock*>& level, int pair_window_size, int N, std::unordered_map<HashString, std::vector<std::pair<SumBlock*, SumBlock*>>>& pair_hashtable);

    std::vector<SumBlock*> next_level(std::vector<SumBlock*>&);
    // Returns a vector of levels of nodes of the tree where
    // each level is represented by a vector of its nodes (left-to-right).
    //
    // A simple levelwise (left-to-right) traversal of the tree would be:
    //     for (std::vector<SumBlock*> level : bt->levelwise_iterator()) {
    //         for (SumBlock* b : level) {
    //             ...
    std::vector<std::vector<SumBlock*>> levelwise_iterator();
};

#endif //MINIMALISTBT_CST_SUMBLOCKTREE_H
