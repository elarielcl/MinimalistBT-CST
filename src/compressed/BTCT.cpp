#include <compressed/BTCT.h>

BTCT::BTCT(BlockTree * bt, int one_symbol) : CBitBlockTree(bt, one_symbol) {
    n_ = bt->input_.length();
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

    sdsl::int_vector<>* local_first_level_min_excess = new sdsl::int_vector<>(first_level.size());

    int level_length = first_level_length_;

    for (int i = 0; i < first_level.size(); ++i) {
        (*local_first_level_min_excess)[i] = 1- first_level[i]->min_excess_;

        int excess = 0;
        for (auto pair : first_level[i]->ranks_) {
            if (pair.first == one_symbol) excess += pair.second;
            else excess -= pair.second;
        }
        (*first_level_prefix_leaf_ranks)[i] = first_level_prefix_ranks_;
        (*first_level_leaf_ranks)[i] = first_level[i]->leaf_rank_;
        first_level_prefix_ranks_ += first_level[i]->leaf_rank_;
    }

    sdsl::util::bit_compress(*(local_first_level_min_excess));
    (bt_min_excess_).push_back(local_first_level_min_excess);

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

        sdsl::int_vector<>* next_level_min_excess = new sdsl::int_vector<>(next_level.size());

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

            (*next_level_min_excess)[i] = 1-next_level[i]->min_excess_;

            int excess = 0;
            for (auto pair : next_level[i]->ranks_) {
                if (pair.first == one_symbol) excess += pair.second;
                else excess -= pair.second;
            }
        }

        sdsl::bit_vector* current_level_suffix_start_with_end_leaf = new sdsl::bit_vector(number_of_leaves ,0);
        sdsl::int_vector<>* current_level_second_leaf_ranks = new sdsl::int_vector<>(number_of_leaves);

        sdsl::bit_vector* current_level_min_in_first_block = new sdsl::bit_vector(number_of_leaves ,0);

        int j = 0;
        sdsl::int_vector<>* current_level_min_excess_back_block = new sdsl::int_vector<>(number_of_leaves);
        for (int i = 0; i < current_level.size(); ++i) {
            if (dynamic_cast<BackBlock*>(current_level[i])) {
                (*current_level_second_leaf_ranks)[j] = current_level[i]->second_leaf_rank_;
                BackBlock* bb  = dynamic_cast<BackBlock*>(current_level[i]);
                if (bb->suffix_starts_with_end_leaf_) {
                    (*current_level_suffix_start_with_end_leaf)[j] = 1;
                } else {
                    (*current_level_suffix_start_with_end_leaf)[j] = 0;
                }

                if (bb->min_in_first_block_) {
                    (*current_level_min_in_first_block)[j] = 1;
                    (*current_level_min_excess_back_block)[j] = 1-(bb->second_block_min_excess_);
                } else {
                    (*current_level_min_in_first_block)[j] = 0;
                    (*current_level_min_excess_back_block)[j] = 1-(bb->first_block_min_excess_);
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

        sdsl::util::bit_compress(*(next_level_min_excess));
        (bt_min_excess_).push_back(next_level_min_excess);

        sdsl::util::bit_compress(*current_level_min_excess_back_block);
        bt_min_excess_back_block_.push_back(current_level_min_excess_back_block);

        bt_min_in_first_block_.push_back(current_level_min_in_first_block);

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

    int number_of_blocks_first_level = (*bt_first_level_prefix_ranks_).size();
    int n = number_of_blocks_first_level;
    int i = 0;
    int r_power_i = 1;
    int last_block_index = 1;
    while (n/r_ > 0) {
        r_power_i *= r_;
        last_block_index += r_power_i;
        n /= r_;
        ++i;
    }
    --last_block_index;

    int p = (number_of_blocks_first_level-r_power_i)/(r_-1);
    int m = (number_of_blocks_first_level-r_power_i)% (r_-1);
    if (m>0) ++p;

    int n_pre_last_level = r_power_i - p;
    int n_last_level = number_of_blocks_first_level+p-r_power_i;
    n_last_level_ = n_last_level;
    n_pre_last_level_ = n_pre_last_level;
    int nodes = n_last_level + (r_power_i*r_-1)/(r_-1);
    n_internal_nodes_ = nodes-n_last_level-n_pre_last_level;

    int* top_excess = new int[nodes];
    int* top_min_excess = new int[nodes];

    for (int i = 0; i < n_last_level; ++i) {
        top_excess[nodes-n_last_level+i] = 2*(*bt_ranks_[0])[i] - first_level_length_;
        top_min_excess[nodes-n_last_level+i] = 1-(*bt_min_excess_[0])[i];
    }

    for (int i = 0; i < n_pre_last_level; ++i) {
        top_excess[nodes-n_last_level-n_pre_last_level+i] = 2*(*bt_ranks_[0])[n_last_level+i] - first_level_length_;
        top_min_excess[nodes-n_last_level-n_pre_last_level+i] = 1-(*bt_min_excess_[0])[n_last_level+i];
    }

    if (n_ % first_level_length_) { // If the last first level block is not complete, fix its excess
        top_excess[last_block_index] += first_level_length_ - (n_ % first_level_length_);
    }

    for (int i = nodes-n_last_level-n_pre_last_level-1; i >=0; --i) {
        top_excess[i] = 0;
        int min_excess = 1;
        for (int j = 0; j < r_; ++j) {
            if (i*r_+j+1 < nodes) {
                if (top_excess[i] + top_min_excess[i*r_+j+1] < min_excess) {
                    min_excess = top_excess[i] + top_min_excess[i*r_+j+1];
                }

                top_excess[i] +=  top_excess[i*r_+j+1];
            }
        }
        top_min_excess[i] = min_excess;
    }

    top_excess_ = new sdsl::int_vector<>(nodes-n_last_level-n_pre_last_level);
    top_min_excess_ = new sdsl::int_vector<>(nodes-n_last_level-n_pre_last_level);
    for (int i = 0; i < nodes-n_last_level-n_pre_last_level; ++i) {
        (*top_excess_)[i] = encoded_excess(top_excess[i]);
        (*top_min_excess_)[i] = 1-top_min_excess[i];
    }

    delete[] top_excess;
    delete[] top_min_excess;

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

    for (auto vector : bt_min_excess_)
        delete vector;
    for (auto vector : bt_min_excess_back_block_)
        delete vector;
    for (auto vector : bt_min_in_first_block_)
        delete vector;

    delete top_excess_;
    delete top_min_excess_;
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
                k += (*leaf_ranks[level])[current_block] - second_rank;
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


int BTCT::fwdsearch(int i, int d) {
    int e = 0;
    int a = fwdsearch(((i+1)%first_level_length_) - 1,first_level_length_-1,d,0, (i+1)/first_level_length_, first_level_length_, e);
    if (a != -1) return a + ((i+1)/first_level_length_)*first_level_length_;
    int index = next_block_fwdsearch((i+1)/first_level_length_, d , e);
    if (index == -1) return n_;
    a = fwdsearch(-1,first_level_length_-1,d,0, index, first_level_length_, e);
    if (a != -1) return a + index*first_level_length_;
    return n_;
}


int BTCT::next_block_fwdsearch(int initial_block, int d, int& e) {
    int number_of_blocks_first_level = (*bt_first_level_prefix_ranks_).size();
    int r_initial_block = 0;
    int index_last = (top_excess_->size()+n_pre_last_level_)-1;
    if (initial_block < n_last_level_) {
        r_initial_block = initial_block + n_pre_last_level_ + n_internal_nodes_;
        index_last = r_*(index_last+1);
    }
    else
        r_initial_block = initial_block - n_last_level_ + n_internal_nodes_;
    bool going_up = true;
    while (true) {
        if (going_up) {
            if (r_initial_block == 0) break;
            while (r_initial_block%r_ != 0) {
                ++r_initial_block;
                if (r_initial_block >= number_of_blocks_first_level+n_internal_nodes_) break;
                if (r_initial_block >= n_internal_nodes_) {
                    if (r_initial_block-n_internal_nodes_ < n_pre_last_level_) initial_block = n_last_level_ + r_initial_block-n_internal_nodes_;
                    else initial_block = r_initial_block - n_pre_last_level_-n_internal_nodes_;
                    if (initial_block >= number_of_blocks_first_level) break;
                    if (e + 1-(int)((*(bt_min_excess_[0]))[initial_block]) <= d) {
                        int a = 1-(int)((*(bt_min_excess_[0]))[initial_block]);
                        return initial_block;
                    } else {
                        int rank_ones = ((int) ((*bt_ranks_[0])[initial_block]));
                        e += 2*rank_ones - first_level_length_;
                    }
                }
                else {
                    if (e + 1-(int)((*(top_min_excess_))[r_initial_block]) <= d) {
                        going_up = false;
                        break;
                    } else {
                        e += decoded_excess((*top_excess_)[r_initial_block]);
                    }
                }
            }

            // LOOKUP
            if (going_up && r_initial_block != index_last && r_initial_block+1 < number_of_blocks_first_level+n_internal_nodes_) {
                if (r_initial_block+1 >= n_internal_nodes_) {
                    if (r_initial_block+1-n_internal_nodes_ < n_pre_last_level_) initial_block = n_last_level_ + r_initial_block+1-n_internal_nodes_;
                    else initial_block = r_initial_block+1 - n_pre_last_level_-n_internal_nodes_;
                    if (initial_block >= number_of_blocks_first_level) break;
                    if (e + 1-(int)((*(bt_min_excess_[0]))[initial_block]) <= d) {
                        return initial_block;
                    }
                }
                else {
                    if (e + 1-(int)((*(top_min_excess_))[r_initial_block+1]) <= d) {
                        going_up = false;
                        ++r_initial_block;
                    }
                }
            }

            if (going_up) {
                r_initial_block = (r_initial_block-1)/r_;
                index_last = index_last/r_ - 1;
            }
        } else {
            if (r_initial_block >= n_internal_nodes_) {
                r_initial_block -= n_internal_nodes_;
                if (r_initial_block < n_pre_last_level_) return  n_last_level_ + r_initial_block;
                return  r_initial_block - n_pre_last_level_;
            } else {
                r_initial_block = r_initial_block*r_ +1;
                while (true) {
                    if (r_initial_block >= n_internal_nodes_) {
                        if (r_initial_block-n_internal_nodes_ < n_pre_last_level_) initial_block = n_last_level_ + r_initial_block-n_internal_nodes_;
                        else initial_block = r_initial_block - n_pre_last_level_-n_internal_nodes_;
                        if (initial_block >= number_of_blocks_first_level) break;
                        if (e + 1-(int) ((*(bt_min_excess_[0]))[initial_block]) <= d) {
                            return initial_block;
                        } else {
                            int rank_ones = ((int) ((*bt_ranks_[0])[initial_block]));
                            e += 2*rank_ones - first_level_length_;
                        }
                        ++r_initial_block;
                    } else {
                        if (e + 1-(int) ((*(top_min_excess_))[r_initial_block]) <= d) {
                            break;
                        } else {
                            e += decoded_excess((*top_excess_)[r_initial_block]);
                        }
                        ++r_initial_block;
                    }
                }
            }
        }
    }
    return -1;

}


int BTCT::fwdsearch(int i, int j, int d, int level, int level_index, int level_length, int& e) {
    if (i == -1 && j == level_length-1 && e + 1-(int)((*(bt_min_excess_[level]))[level_index]) > d) {
        int rank_ones = (*(bt_ranks_[level]))[level_index];
        e += 2*rank_ones - level_length;
        return -1;
    }

    int a = -1;
    if (level == number_of_levels_-1) { // Case LeafBlock
        auto& leaf_string = *leaf_bv_;
        int leaf_string_size = leaf_string.size();
        int k = i+1;
        int chunk = (level_length*level_index + k)/64;
        uint64_t chunk_info = *(leaf_string.m_data+chunk);

        for (; k <= j && level_length*level_index + k < leaf_string_size; ++k) {
            e += ((chunk_info >> ((level_length*level_index + k)%64))%2)?1:-1;

            if (e == d) return k;

            if ((level_length*level_index + k + 1)%64 == 0) {
                chunk++;
                chunk_info = *(leaf_string.m_data+chunk);
            }
        }
        return -1;
    }
    if ((*bt_bv_[level])[level_index]) { // Case InternalBlock
        int child_length = level_length/r_;
        int initial_block = (i+1)/child_length;
        int end_block = j/child_length;
        int rank = (*bt_bv_rank_[level])(level_index);
        int next_level_size = (*bt_min_excess_[level+1]).size();
        if (rank*r_ + end_block >= next_level_size) {
            end_block = (next_level_size-1)%r_;
            j = (end_block+1)*child_length-1;
        }
        if (rank*r_ + initial_block >= next_level_size) return -1;
        if (initial_block == end_block) {
            a = fwdsearch(((i+1)%child_length) - 1, j%child_length, d, level+1, rank*r_ + initial_block,child_length, e);
            if (a != -1) return a + initial_block*child_length;
            return -1;
        } else {
            a = fwdsearch(((i+1)%child_length) - 1, child_length-1, d, level+1, rank*r_ + initial_block,child_length, e);
            if (a != -1) return a + initial_block*child_length;
            for (int block = initial_block+1; block<end_block; ++block) {
                a = fwdsearch(-1, child_length - 1, d, level+1, rank*r_ + block, child_length, e);
                if (a != -1) return a + block*child_length;
            }
            a = fwdsearch(-1, j%child_length, d, level+1, rank*r_ + end_block, child_length, e);
            if (a != -1) return a + end_block*child_length;
            return -1;
        }
    } else { // Case BackBlock
        int back_index = level_index - (*bt_bv_rank_[level])(level_index + 1);
        int encoded_offset = (*bt_offsets_[level])[back_index];

        int first_block = encoded_offset / level_length;
        int offset = encoded_offset % level_length;

        if (i + 1 + offset < level_length && j + offset < level_length) {
            int min_first_block = 1-((*(bt_min_excess_back_block_[level]))[back_index]);
            if ((*(bt_min_in_first_block_[level]))[back_index]) {
                min_first_block = (1-((*(bt_min_excess_[level]))[level_index]));
            }
            if (i + 1 == 0 && j + offset == level_length - 1 &&
                e + min_first_block > d) {
                int rank_ones = (*bt_second_ranks_[level])[back_index];
                e += 2*rank_ones - (level_length-offset);
                return -1;
            }
            a = fwdsearch(i + offset, j + offset, d, level, first_block, level_length, e);
            if (a != -1) return a - offset;
            return -1;
        } else if (i + 1 + offset >= level_length && j + offset >= level_length) {
            int min_second_block = 1-((*(bt_min_excess_back_block_[level]))[back_index]);
            if (!((*(bt_min_in_first_block_[level]))[back_index])) {
                int first_block_excess = 0;
                int rank_ones = (*bt_second_ranks_[level])[back_index];
                first_block_excess += 2*rank_ones - (level_length-offset);
                min_second_block = (1-((*(bt_min_excess_[level]))[level_index])) - first_block_excess;
            }
            if (i + 1 + offset == level_length && j == level_length - 1 &&
                e + min_second_block > d) {
                int rank_ones = (*(bt_ranks_[level]))[level_index];
                e += 2*rank_ones - level_length;

                rank_ones = (*bt_second_ranks_[level])[back_index];
                e -= 2*rank_ones - (level_length-offset);

                return -1;
            }
            a = fwdsearch(offset + i - level_length, offset + j - level_length, d, level, first_block + 1, level_length,
                                   e);
            if (a != -1) return a + level_length - offset;
            return -1;
        } else {
            int min_first_block = 1-((*(bt_min_excess_back_block_[level]))[back_index]);
            if ((*(bt_min_in_first_block_[level]))[back_index]) {
                min_first_block = (1-((*(bt_min_excess_[level]))[level_index]));
            }
            if (i + 1 == 0 && e + min_first_block > d) {
                int rank_ones = (*bt_second_ranks_[level])[back_index];
                e += 2*rank_ones - (level_length-offset);
            } else {
                a = fwdsearch(i + offset, level_length - 1, d, level, first_block, level_length, e);
                if (a != -1) return a - offset;
            }

            int min_second_block = 1-((*(bt_min_excess_back_block_[level]))[back_index]);
            if (!((*(bt_min_in_first_block_[level]))[back_index])) {
                int first_block_excess = 0;
                int rank_ones = (*bt_second_ranks_[level])[back_index];
                first_block_excess += 2*rank_ones - (level_length-offset);
                min_second_block = (1-((*(bt_min_excess_[level]))[level_index])) - first_block_excess;
            }
            if (j == level_length - 1 && e + min_second_block > d) {
                int rank_ones = (*(bt_ranks_[level]))[level_index];
                e += 2*rank_ones - level_length;
                rank_ones = (*bt_second_ranks_[level])[back_index];
                e -= 2*rank_ones - (level_length-offset);
            } else {
                a = fwdsearch(-1, offset + j - level_length, d, level, first_block + 1, level_length, e);
                if (a != -1) return a + level_length - offset;
            }
            return a;
        }
    }
}




int64_t encoded_excess(int64_t e) {
    bool sgn = e < 0;
    if (sgn) e = -e;
    return 2*e+((sgn)?1:0);
}

int64_t decoded_excess(int64_t e) {
    return (e % 2) ? -(e / 2) : (e / 2);
}
