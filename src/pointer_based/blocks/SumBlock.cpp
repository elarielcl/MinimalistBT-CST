#include "pointer_based/blocks/SumBlock.h"


SumBlock::SumBlock(SumBlock* parent, int64_t start_index, int64_t end_index, std::basic_string<int64_t >& source):
        parent_(parent), start_index_(start_index), end_index_(end_index), source_(source), left_(false), right_(false), first_block_(
        this), second_block_(nullptr), pointing_to_me_(0), level_index_(0), first_occurrence_level_index_(0), sum_(0), second_sum_(0) {

}


SumBlock::~SumBlock() {

}


int64_t SumBlock::length() {
    return end_index_-start_index_+1;
}


std::basic_string<int64_t> SumBlock::represented_string() {
    return source_.substr(start_index_, length());
}


std::vector<SumBlock*>& SumBlock::children(int leaf_length, int r) {
    return children_;
}


void SumBlock::clean_unnecessary_expansions() {

}


bool SumBlock::is_leaf() {
    return true;
}


int64_t SumBlock::add_access_support() {
    return 0;
}


int64_t SumBlock::access(int i) {
    return -1;
}


void SumBlock::replace_child(SumBlock* old_child, SumBlock* new_child) {
    for (int i = 0; i < children_.size(); ++i) {
        if (children_[i] == old_child) {
            children_[i] = new_child;
            return;
        }
    }
}