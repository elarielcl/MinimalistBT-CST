#include "pointer_based/blocks/SumBackBlock.h"


SumBackBlock::SumBackBlock(SumBlock* parent, int64_t start_index, int64_t end_index, std::basic_string<int64_t >& source, SumBlock* first_block,
                     SumBlock* second_block, int offset) :
        SumBlock(parent, start_index, end_index, source) {
    first_block_ = first_block;
    if (second_block != nullptr)
        if (second_block->start_index_ == start_index && second_block->end_index_ == end_index) second_block_ = this;
        else second_block_ = second_block;
    offset_ = offset;
    if (first_block_ != nullptr) first_block_->pointing_to_me_ = first_block_->pointing_to_me_ + 1;
    if (second_block_ != nullptr) second_block_->pointing_to_me_ = second_block_->pointing_to_me_ + 1;
}

SumBackBlock::~SumBackBlock() {
    if (first_block_ != nullptr) {
        first_block_->pointing_to_me_ = first_block_->pointing_to_me_ - 1;
    }
    if (second_block_ != nullptr) {
        second_block_->pointing_to_me_ = second_block_->pointing_to_me_ - 1;
    }
}



int64_t SumBackBlock::add_access_support() {
    int first_sum = first_block_->access(offset_-1);
    second_sum_ = (second_block_ == nullptr) ? first_block_->access(offset_+length()-1) - first_sum : first_block_->access(first_block_->length()-1) - first_sum;
    sum_ = (second_block_ == nullptr) ? second_sum_ : second_sum_ + second_block_->access(offset_+length()-1-first_block_->length());
    return sum_;
}


int64_t SumBackBlock::access(int i) {
    if (i + offset_ >= first_block_->length()) return second_sum_ + second_block_->access(offset_+i-first_block_->length());
    return first_block_->access(i+offset_) - (first_block_->sum_ - second_sum_);
}