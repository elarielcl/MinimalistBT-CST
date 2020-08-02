#include <pointer_based/blocks/LeafBlock.h>
#include <unordered_set>
#include <compressed/BTCT.h>
#include "gtest/gtest.h"

#include "pointer_based/BlockTree.h"

using ::testing::Combine;
using ::testing::Values;

typedef BlockTree* CreateBlockTreeFunc(int, int, std::string);

BlockTree* block_tree(int r, int max_leaf_length, std::string input) {

    BlockTree* block_tree_ = new BlockTree(input, r, max_leaf_length);
    block_tree_->process_back_pointers();
    block_tree_->clean_unnecessary_expansions();
    return block_tree_;
}

BlockTree* block_tree_without_cleanning(int r, int max_leaf_length, std::string input) {
    BlockTree* block_tree_ = new BlockTree(input, r, max_leaf_length);
    block_tree_->process_back_pointers();
    return block_tree_;
}


BlockTree* heuristic_bit_block_tree(int r, int max_leaf_length, std::string input) {
    BlockTree* block_tree_ = new BlockTree(input, r, max_leaf_length);
    block_tree_->process_back_pointers_heuristic();
    return block_tree_;
}


class BTCTFixture : public ::testing::TestWithParam<::testing::tuple<int, int, std::string, CreateBlockTreeFunc*>> {
protected:
    virtual void TearDown() {
        delete block_tree_;
        delete btct_;
    }

    virtual void SetUp() {
        CreateBlockTreeFunc* create_blocktree = ::testing::get<3>(GetParam());
        r_ = ::testing::get<0>(GetParam());
        max_leaf_length_ = ::testing::get<1>(GetParam());

        std::ifstream t(::testing::get<2>(GetParam()));
        std::stringstream buffer;
        buffer << t.rdbuf();
        input_= buffer.str();
        one_symbol = input_[0];
        block_tree_ = (*create_blocktree)(r_ , max_leaf_length_, input_);

        std::unordered_set<int> characters;
        for (char c: input_)
            characters.insert(c);
        for (int c: characters) {
            block_tree_->add_rank_select_support(c);
        }
        block_tree_->add_leaf_rank_select_support();
        block_tree_->add_search_support();

        btct_ = new BTCT(block_tree_, one_symbol);
    }

public:
    BlockTree* block_tree_;

    BTCT* btct_;

    std::string input_;
    int one_symbol;
    int r_;
    int max_leaf_length_;


    BTCTFixture() : ::testing::TestWithParam<::testing::tuple<int, int, std::string, CreateBlockTreeFunc*>>() {
    }

    virtual ~BTCTFixture() {
    }
};

INSTANTIATE_TEST_CASE_P(PBTCTTest,
                        BTCTFixture,
                        Combine(Values(2),
                                Values(4),
                                Values("../../../tests/data/dna.par"),
                                Values(&block_tree, &block_tree_without_cleanning, &heuristic_bit_block_tree)));


// This test checks if the fields and
// number_of_levels_ are correct
TEST_P(BTCTFixture, general_fields_check) {
    EXPECT_EQ(btct_->r_, r_);
    auto iterator = block_tree_->levelwise_iterator();
    std::vector<Block*> level;
    bool contains_back_block = false;
    int i;
    for (i = 0; i < iterator.size(); ++i) {
        level = iterator[i];
        for (Block *b : level) {
            if (dynamic_cast<BackBlock*>(b)) contains_back_block = true;
        }
        if (contains_back_block) break;
    }
    EXPECT_EQ(iterator.size()-i, btct_->number_of_levels_);

}

// This test checks if the BTCT has the same
// structure that its corresponding BlockTree
// in particular the bt_bv field is checked
TEST_P(BTCTFixture, bt_bv_field_check) {
    auto iterator = block_tree_->levelwise_iterator();
    std::vector<Block*> level;
    bool contains_back_block = false;
    int i;
    for (i = 0; i < iterator.size(); ++i) {
        level = iterator[i];
        for (Block *b : level) {
            if (dynamic_cast<BackBlock*>(b)) contains_back_block = true;
        }
        if (contains_back_block) break;
    }

    for (int j = 0; j < btct_->number_of_levels_-1; ++j) {
        level = iterator[i + j];
        auto level_bt_bv = *(btct_->bt_bv_[j]);
        EXPECT_EQ(level.size(), level_bt_bv.size());
        for (int k = 0; k < level.size(); ++k) {
            Block *b = level[k];
            if (dynamic_cast<BackBlock *>(b)) EXPECT_FALSE(level_bt_bv[k]);
            else
                EXPECT_TRUE(level_bt_bv[k]);
        }
    }
}


// This test checks if the BTCT has the same
// structure that its corresponding BlockTree
// in particular the bt_offsets field is checked
TEST_P(BTCTFixture, bt_offsets_field_check) {
    auto iterator = block_tree_->levelwise_iterator();
    std::vector<Block*> level;
    bool contains_back_block = false;
    int i = 0;
    for (; i < iterator.size(); ++i) {
        level = iterator[i];
        for (Block *b : level) {
            if (dynamic_cast<BackBlock*>(b)) contains_back_block = true;
        }
        if (contains_back_block) break;
    }

    for (int j = 0; j < btct_->number_of_levels_-1; ++j) {
        level = iterator[i+j];
        auto level_bt_offsets = *(btct_->bt_offsets_[j]);

        int max_size_level = level.front()->length();
        int l = 0;
        for (Block* b: level) {
            if (dynamic_cast<BackBlock*>(b)) {
                EXPECT_EQ(level_bt_offsets[l], max_size_level*b->first_block_->level_index_ + b->offset_);
                ++l;
            }
        }
        EXPECT_EQ(l, level_bt_offsets.size());
    }
}



// This test checks if the BTCT has the same
// structure that its corresponding BlockTree
// in particular the bt_leaf_string field is checked
TEST_P(BTCTFixture, bt_leaf_string_field_check) {
    auto iterator = block_tree_->levelwise_iterator();
    std::string leaf_bv = "";
    for (Block* b : iterator.back()) {
        for (char c : b->represented_string()) {
            if (c == one_symbol) {
                leaf_bv += "1";
            } else {
                leaf_bv += "0";
            }
        }
    }

    std::string leaf_c_bv = "";
    for (int i : (*btct_->leaf_bv_)) {
        leaf_c_bv += i ? "1" : "0";
    }

    EXPECT_EQ(leaf_c_bv, leaf_bv);
}




// This test checks if the BTCT has the same
// structure that its corresponding BlockTree
// in particular the bt_second_ranks field are checked
TEST_P(BTCTFixture, bt_second_ranks_field_check) {
    auto iterator = block_tree_->levelwise_iterator();
    std::vector<Block*> level;
    bool contains_back_block = false;
    int i = 0;
    for (; i < iterator.size(); ++i) {
        level = iterator[i];
        for (Block *b : level) {
            if (dynamic_cast<BackBlock*>(b)) contains_back_block = true;
        }
        if (contains_back_block) break;
    }

    for (int j = 0; j < btct_->number_of_levels_-1; ++j) {
        level = iterator[i+j];

        auto level_bt_second_ranks = *(btct_->bt_second_ranks_[j]);

        int l = 0;
        for (Block *b: level) {
            if (dynamic_cast<BackBlock *>(b)) {
                EXPECT_EQ(level_bt_second_ranks[l], b->second_ranks_[one_symbol]) ;
                ++l;
            }
        }
        EXPECT_EQ(l, level_bt_second_ranks.size());
    }
}


// This test checks if the BTCT has the same
// structure that its corresponding BlockTree
// in particular the bt_ranks_ is checked
TEST_P(BTCTFixture, bt_bv_ranks_prefix_check) {
    auto iterator = block_tree_->levelwise_iterator();
    std::vector<Block*> level;
    bool contains_back_block = false;
    int i = 0;
    for (; i < iterator.size(); ++i) {
        level = iterator[i];
        for (Block *b : level) {
            if (dynamic_cast<BackBlock*>(b)) contains_back_block = true;
        }
        if (contains_back_block) break;
    }



    level = iterator[i];
    auto level_bt_ranks = *(btct_->bt_ranks_[0]);
    EXPECT_EQ(level.size(), level_bt_ranks.size());

    for (int k = 0; k < level.size(); ++k) {
        Block* b = level[k];
        EXPECT_EQ(b->ranks_[one_symbol], level_bt_ranks[k]);
    }


    for (int j = 1; j < btct_->number_of_levels_; ++j) {
        level = iterator[i + j];
        auto level_bt_ranks = *(btct_->bt_ranks_[j]);
        EXPECT_EQ(level.size(), level_bt_ranks.size());

        for (int k = 0; k < level.size(); ++k) {
            Block* b = level[k];
            EXPECT_EQ(b->ranks_[one_symbol], level_bt_ranks[k]);
        }
    }
}


// This test checks if the BTCT has the same
// structure that its corresponding BlockTree
// in particular the first level for bt_prefix_ranks_,
// is checked
TEST_P(BTCTFixture, bt_bv_first_level_prefix_ranks_check) {
    auto iterator = block_tree_->levelwise_iterator();
    std::vector<Block*> level;
    bool contains_back_block = false;
    int i = 0;
    for (; i < iterator.size(); ++i) {
        level = iterator[i];
        for (Block *b : level) {
            if (dynamic_cast<BackBlock*>(b)) contains_back_block = true;
        }
        if (contains_back_block) break;
    }

    level = iterator[i];
    auto first_level_bt_prefix_ranks = *(btct_->bt_first_level_prefix_ranks_);
    int r  = 0;

    EXPECT_EQ(first_level_bt_prefix_ranks.size(), level.size());
    for (int k = 0; k < level.size(); ++k) {
        EXPECT_EQ(r, first_level_bt_prefix_ranks[k]);
        r += level[k]->ranks_[one_symbol];
    }

}


// This test checks the leaf_rank method for every
// position in the input
TEST_P(BTCTFixture, leaf_rank_check) {
    int r = 0;
    bool one_seen = false;
    for (int i = 0; i < input_.length(); ++i) {
        if (input_[i] == input_[0]) {
            one_seen = true;
        } else {
            if (one_seen) {
                ++r;
            }
            one_seen = false;
        }
        EXPECT_EQ(btct_->leaf_rank(i), r);
    }
}


// This test checks the leaf_select for every
// position in the input
TEST_P(BTCTFixture, leaf_select_check) {
    int r = 0;
    bool one_seen = false;
    for (int i = 0; i < input_.length(); ++i) {
        if (input_[i] == input_[0]) {
            one_seen = true;
        } else {
            if (one_seen) {
                ++r;
                EXPECT_EQ(btct_->leaf_select(r), i-1);
            }
            one_seen = false;
        }

    }
}


// This test checks the fwdsearch method for every character
// in the input and d = {-1, -2} works correctly
TEST_P(BTCTFixture, fwdsearch_check) {
    std::unordered_set<int> ds = {-1, -2};
    for (int d : ds) {
        for (int i = 0; i < input_.length(); ++i) {
            int search = btct_->fwdsearch(i, d);
            int excess = 0;
            int j = i;
            while (true) {
                ++j;
                if (j == input_.length()) break;
                excess += (input_[0] == input_[j]) ? 1 : -1;
                if (excess == d) break;
            }

            EXPECT_EQ(search, j);
        }
    }
}


// This test checks if the BTCT has the same
// structure that its corresponding BlockTree
// in particular the bt_second_leaf_rank field is checked
TEST_P(BTCTFixture, bt_second_leaf_ranks_field_check) {
    auto iterator = block_tree_->levelwise_iterator();
    std::vector<Block*> level;
    bool contains_back_block = false;
    int i = 0;
    for (; i < iterator.size(); ++i) {
        level = iterator[i];
        for (Block *b : level) {
            if (dynamic_cast<BackBlock*>(b)) contains_back_block = true;
        }
        if (contains_back_block) break;
    }

    for (int j = 0; j < btct_->number_of_levels_-1; ++j) {
        level = iterator[i+j];
        auto level_bt_second_ranks = *(btct_->bt_second_leaf_ranks_[j]);

        int l = 0;
        for (Block *b: level) {
            if (dynamic_cast<BackBlock *>(b)) {
                EXPECT_EQ(level_bt_second_ranks[l], b->second_leaf_rank_) ;
                ++l;
            }
        }
        EXPECT_EQ(l, level_bt_second_ranks.size());
    }
}


// This test checks if the BTCT has the same
// structure that its corresponding BlockTree
// in particular the bt_leaf_ranks_field is checked
TEST_P(BTCTFixture, bt_bv_leaf_ranks_check) {
    auto iterator = block_tree_->levelwise_iterator();
    std::vector<Block*> level;
    bool contains_back_block = false;
    int i = 0;
    for (; i < iterator.size(); ++i) {
        level = iterator[i];
        for (Block *b : level) {
            if (dynamic_cast<BackBlock*>(b)) contains_back_block = true;
        }
        if (contains_back_block) break;
    }

    level = iterator[i];
    auto level_bt_ranks = *(btct_->bt_leaf_ranks_[0]);
    EXPECT_EQ(level.size(), level_bt_ranks.size());

    for (int k = 0; k < level.size(); ++k) {
        Block* b = level[k];
        EXPECT_EQ(b->leaf_rank_, level_bt_ranks[k]);
    }

    for (int j = 1; j < btct_->number_of_levels_; ++j) {
        level = iterator[i + j];
        auto level_bt_ranks = *(btct_->bt_leaf_ranks_[j]);
        EXPECT_EQ(level.size(), level_bt_ranks.size());

        for (int k = 0; k < level.size(); ++k) {
            Block* b = level[k];
            EXPECT_EQ(b->leaf_rank_, level_bt_ranks[k]);
        }
    }
}


// This test checks if the BTCT has the same
// structure that its corresponding BlockTree
// in particular the first level for bt_prefix_leaf_ranks_
// is checked
TEST_P(BTCTFixture, bt_bv_first_level_prefix_leaf_ranks_check) {
    auto iterator = block_tree_->levelwise_iterator();
    std::vector<Block*> level;
    bool contains_back_block = false;
    int i = 0;
    for (; i < iterator.size(); ++i) {
        level = iterator[i];
        for (Block *b : level) {
            if (dynamic_cast<BackBlock*>(b)) contains_back_block = true;
        }
        if (contains_back_block) break;
    }

    level = iterator[i];
    auto first_level_bt_prefix_ranks = *(btct_->bt_first_level_prefix_leaf_ranks_);
    int r  = 0;

    EXPECT_EQ(first_level_bt_prefix_ranks.size(), level.size());
    for (int k = 0; k < level.size(); ++k) {
        EXPECT_EQ(r, first_level_bt_prefix_ranks[k]);
        r += level[k]->leaf_rank_;
    }
}


// This test checks if the BTCT has the same
// structure that its corresponding BlockTree
// in particular the bt_starts_with_end_leaf_ and bt_suffix_starts_with_end_leaf_
// fields are checked
TEST_P(BTCTFixture, bt_bv_border_fields_check) {
    auto iterator = block_tree_->levelwise_iterator();
    std::vector<Block*> level;
    bool contains_back_block = false;
    int i = 0;
    for (; i < iterator.size(); ++i) {
        level = iterator[i];
        for (Block *b : level) {
            if (dynamic_cast<BackBlock*>(b)) contains_back_block = true;
        }
        if (contains_back_block) break;
    }

    for (int j = 0; j < btct_->number_of_levels_-1; ++j) {
        level = iterator[i+j];
        auto level_bt_starts_with_end_leaf = *(btct_->bt_starts_with_end_leaf_[j]);
        auto level_bt_suffix_starts_with_end_leaf = *(btct_->bt_suffix_starts_with_end_leaf_[j]);

        EXPECT_EQ(level_bt_starts_with_end_leaf.size(), level.size());
        int l = 0;
        int k = 0;
        for (Block *b: level) {
            EXPECT_EQ((b->starts_with_end_leaf_)?1:0, level_bt_starts_with_end_leaf[k]);
            if (dynamic_cast<BackBlock *>(b)) {
                BackBlock* bb = dynamic_cast<BackBlock*>(b);
                EXPECT_EQ(level_bt_suffix_starts_with_end_leaf[l], bb->suffix_starts_with_end_leaf_) ;
                ++l;
            }
            ++k;
        }

        EXPECT_EQ(l, level_bt_suffix_starts_with_end_leaf.size());
    }
}
