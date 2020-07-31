#ifndef BLOCKTREE_PBLOCK_H
#define BLOCKTREE_PBLOCK_H

#include <string>
#include <vector>
#include <unordered_map>

class Block {
public:

    Block* parent_;
    int64_t start_index_;
    int64_t end_index_; // In input string represented by the whole BlockTree

    std::string& source_;

    std::unordered_map<int,int> ranks_;
    std::unordered_map<int,int> second_ranks_;

    int leaf_rank_;
    int first_leaf_rank_;
    int second_leaf_rank_;
    bool starts_with_end_leaf_;
    bool suffix_starts_with_end_leaf_;

    Block* first_block_;
    Block* second_block_;
    int offset_;
    bool left_;
    bool right_;
    int pointing_to_me_;
    int level_index_;
    int first_occurrence_level_index_;

    std::vector<Block*> children_;

    Block(Block*, int64_t, int64_t, std::string&);
    virtual ~Block();

    int64_t length();
    std::string represented_string();

    virtual int add_rank_select_support(int);

    virtual int access(int);
    virtual int rank(int, int);
    virtual int select(int, int);

    virtual int add_leaf_rank_select_support();
    virtual int leaf_rank(int);
    virtual int leaf_select(int);

    virtual std::vector<Block*>& children(int, int);
    virtual void clean_unnecessary_expansions();
    void replace_child(Block*, Block*);

    virtual bool is_leaf();
};

#endif //BLOCKTREE_PBLOCK_H
