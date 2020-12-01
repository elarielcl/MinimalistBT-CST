#ifndef MINIMALISTBT_CST_DABT_H
#define MINIMALISTBT_CST_DABT_H

#include <string>
#include <unordered_map>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>
#include <pointer_based/SumBlockTree.h>

class DABT {
public:
    int r_; // Arity
    int first_level_length_;
    int number_of_levels_;

    std::vector<sdsl::bit_vector*> bt_bv_; // 1 when is Internal SumBlock
    std::vector<sdsl::rank_support_v<1>*> bt_bv_rank_;
    std::vector<sdsl::int_vector<>*> bt_offsets_;

    sdsl::sd_vector<>* sum_leaf_string_map_;
    sdsl::rank_support_sd<1>* sum_leaf_string_map_m_rank_support_;
    sdsl::select_support_sd<1>* sum_leaf_string_map_m_select_support_;
    sdsl::int_vector<>* sum_mapped_leaf_string_;

    sdsl::int_vector<>* bt_first_level_prefix_sum_;
    std::vector<sdsl::int_vector<>*> bt_sum_;
    std::vector<sdsl::int_vector<>*> bt_second_sum_;

    DABT(SumBlockTree*);
    DABT(std::istream&);
    virtual ~DABT();

    int64_t access(int);

    int size();
    void serialize(std::ostream&);
};

int64_t encoded_sum(int64_t);
int64_t decoded_sum(int64_t);


#endif //MINIMALISTBT_CST_DABT_H
