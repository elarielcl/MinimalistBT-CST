#include <compressed/DABT.h>
#include <unordered_set>
#include <algorithm>



DABT::DABT(SumBlockTree * bt) : r_(bt->r_) {
    std::vector<SumBlock*> first_level = {bt->root_block_};
    bool is_first_level = false;
    while (!is_first_level) {
        for (SumBlock* b: first_level) {
            is_first_level = is_first_level || b->is_leaf();
        }
        if (is_first_level) break;
        first_level = bt->next_level(first_level);
    }

    int64_t first_level_prefix_sum_ = 0;
    sdsl::int_vector<>* first_level_prefix_sum = new sdsl::int_vector<>(first_level.size());
    sdsl::int_vector<>* first_level_sum = new sdsl::int_vector<>(first_level.size());

    for (int i = 0; i < first_level.size(); ++i) {
        (*first_level_prefix_sum)[i] = encoded_sum(first_level_prefix_sum_);

        (*first_level_sum)[i] = encoded_sum(first_level[i]->sum_);
        first_level_prefix_sum_ += first_level[i]->sum_;

    }

    sdsl::util::bit_compress(*first_level_prefix_sum);
    bt_first_level_prefix_sum_ = first_level_prefix_sum;

    sdsl::util::bit_compress(*first_level_sum);
    bt_sum_.push_back(first_level_sum);


    first_level_length_ = first_level[0]->length();
    number_of_levels_ = 0;

    std::vector<SumBlock*> current_level = first_level;
    std::vector<SumBlock*> next_level = bt->next_level(first_level);

    while (next_level.size() != 0) {
        sdsl::bit_vector* current_level_bv = new sdsl::bit_vector(current_level.size() ,0);

        sdsl::int_vector<>* next_level_sum = new sdsl::int_vector<>(next_level.size());

        int number_of_leaves = 0;
        int current_length = current_level.front()->length();
        for (int i = 0; i < current_level.size(); ++i) {
            current_level[i]->level_index_ = i;

            if (current_level[i]->is_leaf()) {
                (*current_level_bv)[i] = 0;
                ++number_of_leaves;
            }
            else {
                (*current_level_bv)[i] = 1;
            }
        }

        for (int i = 0; i < next_level.size(); ++i) {
            (*next_level_sum)[i] = encoded_sum(next_level[i]->sum_);
        }



        sdsl::int_vector<>* current_level_offsets = new sdsl::int_vector<>(number_of_leaves);

        sdsl::int_vector<>* current_level_second_sum;
        current_level_second_sum = new sdsl::int_vector<>(number_of_leaves);


        int j = 0;
        for (int i = 0; i < current_level.size(); ++i) {
            if (!(*current_level_bv)[i]) {
                (*current_level_second_sum)[j] = encoded_sum(current_level[i]->second_sum_);
                (*current_level_offsets)[j++] = current_level[i]->first_block_->level_index_ * current_length + current_level[i]->offset_;
            }
        }

        sdsl::util::bit_compress(*current_level_offsets);
        bt_offsets_.push_back(current_level_offsets);

        sdsl::util::bit_compress(*next_level_sum);
        bt_sum_.push_back(next_level_sum);

        sdsl::util::bit_compress(*current_level_second_sum);
        bt_second_sum_.push_back(current_level_second_sum);

        bt_bv_.push_back(current_level_bv);
        sdsl::rank_support_v<1>* current_level_bv_rank = new sdsl::rank_support_v<1>(current_level_bv);
        bt_bv_rank_.push_back(current_level_bv_rank);

        current_level = next_level;
        next_level = bt->next_level(current_level);
        ++number_of_levels_;
    }

    ++number_of_levels_;

    std::vector<SumBlock*> last_level = current_level;

    std::basic_string<int64_t> leaf_string;
    std::basic_string<int64_t> sum_leaf_string;
    for (SumBlock* b: last_level) {
        std::basic_string<int64_t> wstring = b->represented_string();
        leaf_string += wstring;
        int64_t sum = 0;
        for (int64_t s : wstring) {
            sum += s;
            sum_leaf_string += sum;
        }
    }
    

    std::set<uint64_t> sum_leaf_string_set;
    for (int i = 0; i<sum_leaf_string.size(); ++i) {
        sum_leaf_string_set.insert(encoded_sum(sum_leaf_string[i]));
    }
    std::vector<uint64_t> sum_leaf_string_array;
    for (int64_t e : sum_leaf_string_set) {
        sum_leaf_string_array.push_back(e);
    }

    std::sort(sum_leaf_string_array.begin(), sum_leaf_string_array.end());

    sum_leaf_string_map_ = new sdsl::sd_vector<>(sum_leaf_string_array.begin(), sum_leaf_string_array.end());
    sum_leaf_string_map_m_rank_support_ = new sdsl::rank_support_sd<1>();
    sum_leaf_string_map_m_select_support_ = new sdsl::select_support_sd<1>();
    sdsl::util::init_support(*sum_leaf_string_map_m_rank_support_, sum_leaf_string_map_);
    sdsl::util::init_support(*sum_leaf_string_map_m_select_support_, sum_leaf_string_map_);

    sum_mapped_leaf_string_ = new sdsl::int_vector<>(sum_leaf_string.size());
    for (int i = 0; i<sum_leaf_string.size(); ++i) {
        (*sum_mapped_leaf_string_)[i] = (*sum_leaf_string_map_m_rank_support_)(encoded_sum(sum_leaf_string[i])+1);
    }
    sdsl::util::bit_compress(*sum_mapped_leaf_string_);

}



DABT::DABT(std::istream& in) {
    in.read((char *) &r_, sizeof(int));
    in.read((char *) &first_level_length_, sizeof(int));
    in.read((char *) &number_of_levels_, sizeof(int));

    for (int i = 0; i < number_of_levels_-1; ++i) {
        bt_bv_.push_back(new sdsl::bit_vector());
        (*bt_bv_[i]).load(in);
    }

    for (sdsl::bit_vector* bv : bt_bv_) {
        bt_bv_rank_.push_back(new sdsl::rank_support_v<1>(bv));
    }

    for (int i = 0; i < number_of_levels_-1; ++i) {
        bt_offsets_.push_back(new sdsl::int_vector<>());
        (*bt_offsets_[i]).load(in);
    }

    sum_leaf_string_map_ = new sdsl::sd_vector<>();
    (*sum_leaf_string_map_).load(in);

    sum_leaf_string_map_m_rank_support_ = new sdsl::rank_support_sd<1>();
    sum_leaf_string_map_m_select_support_ = new sdsl::select_support_sd<1>();
    sdsl::util::init_support(*sum_leaf_string_map_m_rank_support_, sum_leaf_string_map_);
    sdsl::util::init_support(*sum_leaf_string_map_m_select_support_, sum_leaf_string_map_);


    sum_mapped_leaf_string_ = new sdsl::int_vector<>();
    (*sum_mapped_leaf_string_).load(in);

    bt_first_level_prefix_sum_ = new sdsl::int_vector<>();
    (*bt_first_level_prefix_sum_).load(in);

    for (int i = 0; i < number_of_levels_; ++i) {
            bt_sum_.push_back(new sdsl::int_vector<>());
            (*bt_sum_[i]).load(in);
    }

    for (int i = 0; i < number_of_levels_-1; ++i) {
        bt_second_sum_.push_back(new sdsl::int_vector<>());
        (*bt_second_sum_[i]).load(in);
    }

}



DABT::~DABT() {

    for (sdsl::bit_vector* bv : bt_bv_)
        delete bv;
    for (sdsl::rank_support_v<1>* rank : bt_bv_rank_)
        delete rank;


    for (sdsl::int_vector<>* offsets : bt_offsets_)
        delete offsets;

    delete bt_first_level_prefix_sum_;

    for (sdsl::int_vector<>* sum : bt_sum_) {
        delete sum;
    }
    for (sdsl::int_vector<>* sum : bt_second_sum_) {
        delete sum;
    }

    delete sum_leaf_string_map_;
    delete sum_leaf_string_map_m_rank_support_;
    delete sum_leaf_string_map_m_select_support_;
    delete sum_mapped_leaf_string_;

}



int64_t DABT::access(int i) {

    int current_block = i/first_level_length_;
    int current_length = first_level_length_;
    i = i-current_block*current_length;
    int level = 0;

    int64_t r = decoded_sum((*bt_first_level_prefix_sum_)[current_block]);
    while (level < number_of_levels_-1) {
        if ((*bt_bv_[level])[current_block]) { // Case SumInternalBlock
            current_length /= r_;
            int child_number = i/current_length;
            i -= child_number*current_length;

            int firstChild = (*bt_bv_rank_[level])(current_block)*r_;
            for (int child = firstChild; child < firstChild + child_number; ++child)
                r += decoded_sum((*bt_sum_[level+1])[child]);
            current_block = firstChild + child_number;
            ++level;
        } else { // Case SumBackBlock
            int index = current_block - (*bt_bv_rank_[level])(current_block+1);
            int encoded_offset = (*bt_offsets_[level])[index];
            current_block = encoded_offset/current_length;
            i += encoded_offset%current_length;
            r += decoded_sum((*bt_second_sum_[level])[index]);
            if (i >= current_length) {
                ++current_block;
                i -= current_length;
            } else {
                r -= decoded_sum((*bt_sum_[level])[current_block]);
            }

        }
    }

    i  += current_block*current_length;
    r += decoded_sum((*sum_leaf_string_map_m_select_support_)((*sum_mapped_leaf_string_)[i]));

    return r;
}



int DABT::size() {

    int bt_bv_size = sizeof(void*);
    for (sdsl::bit_vector* bv : bt_bv_) {
        bt_bv_size += sdsl::size_in_bytes(*bv);
    }

    int bt_bv_rank_size = sizeof(void*);
    for (sdsl::rank_support_v<1>* bvr : bt_bv_rank_) {
        bt_bv_rank_size += sdsl::size_in_bytes(*bvr);
    }

    int bt_offsets_size = sizeof(void*);
    for (sdsl::int_vector<>* offsets: bt_offsets_) {
        bt_offsets_size += sdsl::size_in_bytes(*offsets);
    }

    int sum_mapped_leaf_string_size = sdsl::size_in_bytes(*sum_mapped_leaf_string_);
    int map_sum_size = sdsl::size_in_bytes(*sum_leaf_string_map_);

    int bt_prefix_sum_first_level_size = 0;
    bt_prefix_sum_first_level_size += sdsl::size_in_bytes(*bt_first_level_prefix_sum_);

    int bt_sum_total_size = (bt_sum_.size()+1) * sizeof(void*);
    int size = 0;
    for (sdsl::int_vector<>* sum: bt_sum_) {
        size += sdsl::size_in_bytes(*sum);
    }
    bt_sum_total_size += size;



    int bt_second_sum_total_size = (bt_second_sum_.size()+1) * sizeof(void*);
    size = 0;
    for (sdsl::int_vector<>* sum: bt_second_sum_) {
        size += sdsl::size_in_bytes(*sum);
    }
    bt_second_sum_total_size += size;

    int partial_total_size = bt_bv_size+ bt_bv_rank_size+ bt_offsets_size + sum_mapped_leaf_string_size + map_sum_size;
    int sum_alternative_version_size = bt_second_sum_total_size + bt_sum_total_size + bt_prefix_sum_first_level_size;

    return partial_total_size + sum_alternative_version_size;
}


void DABT::serialize(std::ostream& out) {

    out.write((char *) &r_, sizeof(int));
    out.write((char *) &first_level_length_, sizeof(int));
    out.write((char *) &number_of_levels_, sizeof(int));

    for (sdsl::bit_vector* bv : bt_bv_) {
        (*bv).serialize(out);
    }

    for (sdsl::int_vector<>* offsets : bt_offsets_) {
        (*offsets).serialize(out);
    }

    (*sum_leaf_string_map_).serialize(out);

    (*sum_mapped_leaf_string_).serialize(out);

    (*bt_first_level_prefix_sum_).serialize(out);


    for (sdsl::int_vector<>* ranks: bt_sum_) {
            (*ranks).serialize(out);
    }

    for (sdsl::int_vector<>* second_ranks: bt_second_sum_) {
            (*second_ranks).serialize(out);
    }
}


int64_t encoded_sum(int64_t sum) {
    bool sgn = sum < 0;
    if (sgn) sum = -sum;
    return 2*sum+((sgn)?1:0);
}

int64_t decoded_sum(int64_t sum) {
    return (sum%2) ? -(sum/2) : (sum/2);
}