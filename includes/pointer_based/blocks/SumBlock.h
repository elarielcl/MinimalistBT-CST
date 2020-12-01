#ifndef MINIMALISTBT_CST_SUMBLOCK_H
#define MINIMALISTBT_CST_SUMBLOCK_H

#include <string>
#include <vector>
#include <unordered_map>

class SumBlock {
public:

    SumBlock* parent_;
    int64_t start_index_;
    int64_t end_index_; // In input string represented by the whole SumBlockTree

    std::basic_string<int64_t>& source_;

    int64_t second_sum_;
    int64_t sum_;

    SumBlock* first_block_;
    SumBlock* second_block_;
    int offset_;
    bool left_;
    bool right_;
    int pointing_to_me_;
    int level_index_;
    int first_occurrence_level_index_;

    std::vector<SumBlock*> children_;

    SumBlock(SumBlock*, int64_t, int64_t, std::basic_string<int64_t>&);
    virtual ~SumBlock();

    int64_t length();
    std::basic_string<int64_t> represented_string();

    virtual int64_t add_access_support();
    virtual int64_t access(int);

    virtual std::vector<SumBlock*>& children(int, int);
    virtual void clean_unnecessary_expansions();
    void replace_child(SumBlock*, SumBlock*);

    virtual bool is_leaf();
};

#endif //MINIMALISTBT_CST_SUMBLOCK_H
