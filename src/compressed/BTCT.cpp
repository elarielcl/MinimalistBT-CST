#include <compressed/BTCT.h>

BTCT::BTCT(BlockTree * bt, int one_symbol) : CBitBlockTree(bt, one_symbol) {
    std::vector<Block*> first_level = {bt->root_block_};
    bool is_first_level = false;
    while (!is_first_level) {
        for (Block* b: first_level) {
            is_first_level = is_first_level || b->is_leaf();
        }
        if (is_first_level) break;
        first_level = bt->next_level(first_level);
    }


    int first_level_prefix_ranks_ = 0;
    sdsl::int_vector<>* first_level_prefix_leaf_ranks = new sdsl::int_vector<>(first_level.size());
    sdsl::int_vector<>* first_level_leaf_ranks = new sdsl::int_vector<>(first_level.size());


    int level_length = first_level_length_;

    for (int i = 0; i < first_level.size(); ++i) {
        int excess = 0;
        for (auto pair : first_level[i]->ranks_) {
            if (pair.first == one_symbol) excess += pair.second;
            else excess -= pair.second;
        }
        (*first_level_prefix_leaf_ranks)[i] = first_level_prefix_ranks_;
        (*first_level_leaf_ranks)[i] = first_level[i]->leaf_rank_;
        first_level_prefix_ranks_ += first_level[i]->leaf_rank_;
    }

    sdsl::util::bit_compress(*(first_level_prefix_leaf_ranks));
    bt_first_level_prefix_leaf_ranks_ = first_level_prefix_leaf_ranks;

    sdsl::util::bit_compress(*(first_level_leaf_ranks));
    bt_leaf_ranks_.push_back(first_level_leaf_ranks);

    std::vector<Block*> current_level = first_level;
    std::vector<Block*> next_level = bt->next_level(first_level);
    level_length = first_level_length_;
    while (next_level.size() != 0) {
        sdsl::bit_vector* current_level_starts_with_end_leaf = new sdsl::bit_vector(current_level.size() ,0);

        sdsl::int_vector<>* next_level_leaf_ranks = new sdsl::int_vector<>(next_level.size());

        int number_of_leaves = 0;
        for (int i = 0; i < current_level.size(); ++i) {
            if (current_level[i]->is_leaf()) {
                ++number_of_leaves;
            }
            if (current_level[i]->starts_with_end_leaf_) {
                (* current_level_starts_with_end_leaf)[i] = 1; //HERE! at the end tooo
            } else {
                (* current_level_starts_with_end_leaf)[i] = 0;
            }
        }


        for (int i = 0; i < next_level.size(); ++i) {

            (*next_level_leaf_ranks)[i] = next_level[i]->leaf_rank_;

            int excess = 0;
            for (auto pair : next_level[i]->ranks_) {
                if (pair.first == one_symbol) excess += pair.second;
                else excess -= pair.second;
            }
        }

        sdsl::bit_vector* current_level_suffix_start_with_end_leaf = new sdsl::bit_vector(number_of_leaves ,0);
        sdsl::int_vector<>* current_level_second_leaf_ranks = new sdsl::int_vector<>(number_of_leaves);

        int j = 0;
        for (int i = 0; i < current_level.size(); ++i) {
            if (dynamic_cast<BackBlock*>(current_level[i])) {
                (*current_level_second_leaf_ranks)[j] = current_level[i]->second_leaf_rank_;
                BackBlock* bb  = dynamic_cast<BackBlock*>(current_level[i]);
                if (bb->suffix_starts_with_end_leaf_) {
                    (*current_level_suffix_start_with_end_leaf)[j] = 1;
                } else {
                    (*current_level_suffix_start_with_end_leaf)[j] = 0;
                }
                ++j;
            }
        }

        sdsl::util::bit_compress(*(next_level_leaf_ranks));
        (bt_leaf_ranks_).push_back(next_level_leaf_ranks);

        sdsl::util::bit_compress(*(current_level_second_leaf_ranks));
        (bt_second_leaf_ranks_).push_back(current_level_second_leaf_ranks);

        bt_starts_with_end_leaf_.push_back(current_level_starts_with_end_leaf);

        bt_suffix_starts_with_end_leaf_.push_back(current_level_suffix_start_with_end_leaf);

        current_level = next_level;
        next_level = bt->next_level(current_level);
        level_length /= r_;
    }


    std::vector<Block*> last_level = current_level;
    sdsl::bit_vector* last_level_starts_with_end_leaf = new sdsl::bit_vector(last_level.size() ,0);

    for (int i = 0; i < last_level.size(); ++i) {
        if (last_level[i]->starts_with_end_leaf_) {
            (*last_level_starts_with_end_leaf)[i] = 1;
        } else {
            (*last_level_starts_with_end_leaf)[i] = 0;
        }
    }

    bt_starts_with_end_leaf_.push_back(last_level_starts_with_end_leaf);
}


BTCT::~BTCT() {
    delete bt_first_level_prefix_leaf_ranks_;

    for (auto vector : bt_leaf_ranks_)
        delete vector;
    for (auto vector : bt_second_leaf_ranks_)
        delete vector;
    for (auto vector : bt_starts_with_end_leaf_)
        delete vector;
    for (auto vector : bt_suffix_starts_with_end_leaf_)
        delete vector;
}


int BTCT::leaf_rank(int i) {
    auto& first_level_prefix_leaf_ranks = bt_first_level_prefix_leaf_ranks_;

    auto& leaf_ranks = bt_leaf_ranks_;
    auto& second_leaf_ranks = bt_second_leaf_ranks_;
    auto& starts_with_end_leaf = bt_starts_with_end_leaf_;
    auto& suffix_starts_with_end_leaf = bt_suffix_starts_with_end_leaf_;

    int current_block = i/first_level_length_;
    int current_length = first_level_length_;
    i = i-current_block*current_length;
    int level = 0;

    int r = (*first_level_prefix_leaf_ranks)[current_block];
    while (level < number_of_levels_-1) {
        if ((*bt_bv_[level])[current_block]) { // Case InternalBlock
            current_length /= r_;
            int child_number = i/current_length;
            i -= child_number*current_length;

            int firstChild = (*bt_bv_rank_[level])(current_block)*r_;
            for (int child = firstChild; child < firstChild + child_number; ++child)
                r += (*leaf_ranks[level+1])[child];
            current_block = firstChild + child_number;
            ++level;
        } else { // Case BackBlock
            int index = current_block - (*bt_bv_rank_[level])(current_block+1);
            int encoded_offset = (*bt_offsets_[level])[index];
            int f_condition = (*starts_with_end_leaf[level])[current_block];
            int s_condition = (*suffix_starts_with_end_leaf[level])[index];
            if (f_condition && !s_condition) ++r;
            if (!f_condition && s_condition) --r;

            current_block = encoded_offset/current_length;
            i += encoded_offset%current_length;
            r += (*second_leaf_ranks[level])[index];
            if (i >= current_length) {
                ++current_block;
                i -= current_length;
            } else {
                r -= (*leaf_ranks[level])[current_block];
            }

        }
    }

    auto& leaf_bv = *leaf_bv_;
    int chunk = (current_block*current_length)/64;
    uint64_t chunk_info = *(leaf_bv.m_data + chunk);

    i  += current_block*current_length;

    bool one_seen = (*starts_with_end_leaf[level])[current_block];

    for (int j = current_block*current_length; j <= i; ++j) {
        int value = (chunk_info >> (j%64))%2;
        if (value == 1) {
            one_seen = true;
        } else {
            if (one_seen) ++r;
            one_seen = false;
        }
        if ((j + 1)%64 == 0) {
            ++chunk;
            chunk_info = *(leaf_bv.m_data + chunk);
        }
    }

    return r;
}


int BTCT::leaf_select(int k) {
    auto& first_level_prefix_ranks = bt_first_level_prefix_ranks_;
    auto& first_level_prefix_leaf_ranks = bt_first_level_prefix_leaf_ranks_;

    auto& leaf_ranks = bt_leaf_ranks_;
    auto& second_leaf_ranks = bt_second_leaf_ranks_;
    auto& starts_with_end_leaf = bt_starts_with_end_leaf_;
    auto& suffix_starts_with_end_leaf = bt_suffix_starts_with_end_leaf_;

    int current_block = (k-1)/first_level_length_;

    int end_block = (*first_level_prefix_ranks).size()-1;
    while (current_block != end_block) {
        int m = current_block + (end_block-current_block)/2;
        int f = (*first_level_prefix_leaf_ranks)[m];
        if (f < k) {
            if (end_block - current_block == 1) {
                if ((*first_level_prefix_leaf_ranks)[m+1] < k) {
                    current_block = m+1;
                }
                break;
            }
            current_block = m;
        } else {
            end_block = m-1;
        }
    }

    int current_length = first_level_length_;
    int s = current_block*current_length;
    k -= (*first_level_prefix_leaf_ranks)[current_block];
    int level = 0;
    while (level < number_of_levels_-1) {
        if ((*bt_bv_[level])[current_block]) { // Case InternalBlock
            int firstChild = (*bt_bv_rank_[level])(current_block)*r_;
            int child = firstChild;
            int r = (*leaf_ranks[level+1])[child];
            int last_possible_child = (firstChild + r_-1 > (*leaf_ranks[level+1]).size()-1) ?  (*leaf_ranks[level+1]).size()-1 : firstChild + r_-1;
            while ( child < last_possible_child && k > r) { //Border conditions?
                ++child;
                r+= (*leaf_ranks[level+1])[child];
            }
            k -= r - (*leaf_ranks[level+1])[child];
            current_length /= r_;
            s += (child-firstChild)*current_length;
            current_block = child;
            ++level;
        } else { // Case BackBlock
            int index = current_block -  (*bt_bv_rank_[level])(current_block+1);
            int encoded_offset = (*bt_offsets_[level])[index];
            int f_condition = (*starts_with_end_leaf[level])[current_block];
            int s_condition = (*suffix_starts_with_end_leaf[level])[index];
            if (f_condition && !s_condition) --k;
            if (!f_condition && s_condition) ++k;
            current_block = encoded_offset/current_length;

            int second_rank = (*second_leaf_ranks[level])[index];

            if (k > second_rank) {
                k -= second_rank;
                s += current_length - encoded_offset%current_length;
                ++current_block;
                if (k == 0) return s-1;
            } else {
                if (k == 0) return s-1;
                s -= encoded_offset%current_length;
                k += (*leaf_ranks[level])[current_block] - second_rank;;
            }

        }
    }

    bool one_seen = (*starts_with_end_leaf[level])[current_block];
    if (one_seen && k == 1) return s-1;

    auto& leaf_bv = *leaf_bv_;
    int chunk = (current_block*current_length)/64;
    uint64_t chunk_info = *(leaf_bv.m_data + chunk);

    for (int j = current_block*current_length; ; ++j) {
        int value = (chunk_info >> (j%64))%2;
        if (value == 1) {
            one_seen = true;
        } else {
            if (one_seen) {
                --k;
            }
            one_seen = false;
        }
        if (!k) return s + j - current_block*current_length - 1;
        if ((j + 1)%64 == 0) {
            ++chunk;
            chunk_info = *(leaf_bv.m_data + chunk);
        }
    }

    return -1;
}
