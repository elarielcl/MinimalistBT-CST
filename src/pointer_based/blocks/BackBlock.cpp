#include "pointer_based/blocks/BackBlock.h"


BackBlock::BackBlock(Block* parent, int64_t start_index, int64_t end_index, std::string& source, Block* first_block,
                     Block* second_block, int offset) :
        Block(parent, start_index, end_index, source) {
    first_block_ = first_block;
    if (second_block != nullptr)
        if (second_block->start_index_ == start_index && second_block->end_index_ == end_index) second_block_ = this;
        else second_block_ = second_block;
    offset_ = offset;
    if (first_block_ != nullptr) first_block_->pointing_to_me_ = first_block_->pointing_to_me_ + 1;
    if (second_block_ != nullptr) second_block_->pointing_to_me_ = second_block_->pointing_to_me_ + 1;
}

BackBlock::~BackBlock() {
    if (first_block_ != nullptr) {
        first_block_->pointing_to_me_ = first_block_->pointing_to_me_ - 1;
    }
    if (second_block_ != nullptr) {
        second_block_->pointing_to_me_ = second_block_->pointing_to_me_ - 1;
    }
}

int BackBlock::add_rank_select_support(int c) {
    int first_rank = first_block_->rank(c, offset_-1);
    int second_rank = (second_block_ == nullptr) ? first_block_->rank(c, offset_ + length() - 1) - first_rank : first_block_->rank(c, first_block_->length() - 1) - first_rank;
    second_ranks_[c] = second_rank;
    ranks_[c] = (second_block_ == nullptr) ? second_rank : second_rank + second_block_->rank(c, offset_+length()-1-first_block_->length());
    return ranks_[c];
}


int BackBlock::add_leaf_rank_select_support() {
    starts_with_end_leaf_ = source_[start_index_] != source_[0] && source_[start_index_-1] == source_[0];
    suffix_starts_with_end_leaf_ = source_[first_block_->start_index_+offset_] != source_[0] && source_[first_block_->start_index_+offset_-1] == source_[0];
    int first_rank = (offset_ == 0) ? 0 : first_block_->leaf_rank(offset_ - 1);
    int second_rank = (second_block_ == nullptr) ?
                      first_block_->leaf_rank(offset_ + length() - 1) - first_rank :
                      first_block_->leaf_rank(first_block_->length() - 1) - first_rank;
    first_leaf_rank_ = first_rank;
    second_leaf_rank_ = second_rank;

    leaf_rank_ = (second_block_ == nullptr) ? second_leaf_rank_ : second_leaf_rank_ + second_block_->leaf_rank(
            offset_ + length() - 1 - first_block_->length());

    if (starts_with_end_leaf_ && !suffix_starts_with_end_leaf_) ++leaf_rank_;
    if (!starts_with_end_leaf_ && suffix_starts_with_end_leaf_) --leaf_rank_;

    return leaf_rank_;
}


int BackBlock::rank(int c, int i) {
    if (i + offset_ >= first_block_->length()) return second_ranks_[c] + second_block_->rank(c, offset_+i-first_block_->length()); //Loop if it's itself
    return first_block_->rank(c, i+offset_) - (first_block_->ranks_[c] - second_ranks_[c]);
}


int BackBlock::leaf_rank(int i) {
    int different_start_value = 0;
    if (starts_with_end_leaf_ && !suffix_starts_with_end_leaf_) ++different_start_value;
    if (!starts_with_end_leaf_ && suffix_starts_with_end_leaf_) --different_start_value;
    int separated_blocks_value = 0;

    if (!offset_) {
        if (i >= first_block_->length()) {
            return different_start_value + separated_blocks_value + second_leaf_rank_ + second_block_->leaf_rank(i-first_block_->length()); //Loop if it's itself
        }

        return different_start_value + first_block_->leaf_rank(i);
    }
    if (i + offset_ >= first_block_->length()) return separated_blocks_value+different_start_value+ second_leaf_rank_ + second_block_->leaf_rank(offset_+i-first_block_->length()); //Loop if it's itself
    return different_start_value + first_block_->leaf_rank(i+offset_) - (first_block_->leaf_rank_ - second_leaf_rank_);
}


int BackBlock::leaf_select(int j) {
    int different_start_value = 0;
    if (starts_with_end_leaf_ && !suffix_starts_with_end_leaf_) ++different_start_value;
    if (!starts_with_end_leaf_ && suffix_starts_with_end_leaf_) --different_start_value;
    int separated_blocks_value = 0;

    if (offset_ == 0) {
        if (j > different_start_value + second_leaf_rank_){
            if (j - different_start_value - second_leaf_rank_ == 1 && separated_blocks_value == 1) return first_block_->length()-offset_-1;
            return second_block_->leaf_select(j-second_leaf_rank_-different_start_value-separated_blocks_value) + first_block_->length();
        }
        if (j == different_start_value) return -1;
        return first_block_->leaf_select(j+first_leaf_rank_-different_start_value);
    }


    if (j > different_start_value + second_leaf_rank_) {
        if (j - different_start_value - second_leaf_rank_ == 1 && separated_blocks_value == 1) return first_block_->length()-offset_-1;
        return second_block_->leaf_select(j-second_leaf_rank_-different_start_value-separated_blocks_value) + first_block_->length() - offset_;
    }

    if (j == different_start_value) return -1;
    return first_block_->leaf_select(j+(first_block_->leaf_rank_ - second_leaf_rank_)-different_start_value) - offset_;
}


int BackBlock::select(int c, int j) {
    if (j > second_ranks_[c]) return second_block_->select(c, j-second_ranks_[c]) + first_block_->length() - offset_;
    return first_block_->select(c, j+first_block_->ranks_[c] - second_ranks_[c]) - offset_;
}


int BackBlock::access(int i) {
    if (i + offset_ >= first_block_->length()) return second_block_->access(offset_+i-first_block_->length());
    return first_block_->access(i+offset_);
}