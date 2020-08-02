#include "pointer_based/blocks/LeafBlock.h"

LeafBlock::LeafBlock(Block* parent, int64_t start_index, int64_t end_index, std::string& source):
        Block(parent, start_index, end_index, source) {
}

LeafBlock::~LeafBlock() {

}

int64_t LeafBlock::size() {
    int64_t source_end_index = source_.size() - 1;
    return (end_index_ <= source_end_index ? end_index_ : source_end_index)-start_index_+1;
}

int LeafBlock::add_rank_select_support(int c) {
    ranks_[c] = rank(c, size()-1);
    return ranks_[c];
}


int LeafBlock::add_leaf_rank_select_support() {
    starts_with_end_leaf_ = source_[start_index_] != source_[0] && source_[start_index_-1] == source_[0];
    leaf_rank_ = (starts_with_end_leaf_) ? 1 : 0;
    std::string rep = represented_string();
    bool one_seen = false;

    for (char c: rep) {
        if (c == source_[0]) {
            one_seen = true;
        }
        else  {
            if (one_seen) {
                ++leaf_rank_;
            }
            one_seen = false;
        }
    }

    return leaf_rank_;
}


int LeafBlock::add_search_support() {
    int excess = 0;
    int min_excess = 1;
    std::string rep = represented_string();
    for (char c : rep) {
        excess += (c == source_[0]) ? 1 : -1;
        if (excess < min_excess) min_excess = excess;
    }
    min_excess_ = min_excess;
    return min_excess_;
}


int LeafBlock::rank(int c, int i) {
    int r = 0;
    for (int j = 0; j<=i; ++j) {
        if (source_[start_index_+j] == c) ++r;
    }
    return r;
}


int LeafBlock::leaf_rank(int i) {
    if (i == -1) return 0;
    int lr = 0;
    bool one_seen = starts_with_end_leaf_;
    for (int j = 0; j<=i; ++j) {
        if (source_[start_index_+j] == source_[0]) {
            one_seen = true;
        } else {
            if (one_seen) {
                ++lr;
            }
            one_seen = false;
        }
    }
    return lr;
}


int LeafBlock::leaf_select(int j) {
    if (starts_with_end_leaf_ && j == 1) return -1;
    bool one_seen = starts_with_end_leaf_;
    for (int i = 0; i < size(); ++i) {
        if (source_[start_index_+i] == source_[0]) {
            one_seen = true;
        } else {
            if (one_seen) {
                --j;
            }
            one_seen = false;
        }
        if (!j) return i-1;
    }
    return -1;
}


int LeafBlock::select(int c, int j) {
    for (int i = 0; i < size(); ++i) {
        if (((int)(source_[start_index_+i])) == c) --j;
        if (!j) return i;
    }
    return -1;
}


int LeafBlock::access(int i) {
    return source_[start_index_+i];
}