#include "pointer_based/blocks/SumLeafBlock.h"


SumLeafBlock::SumLeafBlock(SumBlock* parent, int64_t start_index, int64_t end_index, std::basic_string<int64_t >& source):
        SumBlock(parent, start_index, end_index, source) {
}


SumLeafBlock::~SumLeafBlock() {

}


int64_t SumLeafBlock::size() {
    int64_t source_end_index = source_.size() - 1;
    return (end_index_ <= source_end_index ? end_index_ : source_end_index)-start_index_+1;
}


int64_t SumLeafBlock::add_access_support() {
    sum_ = access(size()-1);
    return sum_;
}


int64_t SumLeafBlock::access(int i) {
    int s = 0;
    for (int j = 0; j<=i; ++j) {
        s += source_[start_index_+j];
    }
    return s;
}