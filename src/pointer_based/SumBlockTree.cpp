#include "pointer_based/SumBlockTree.h"

#include <pointer_based/blocks/SumInternalBlock.h>
#include <pointer_based/blocks/SumLeafBlock.h>

#include <stack>


SumBlockTree::SumBlockTree(std::basic_string<int64_t >& input, int r, int leaf_length): r_(r), input_(input), leaf_length_(leaf_length) {

    if (input_.size() <= leaf_length_ || input_.size()<r)
        root_block_ = new SumLeafBlock(nullptr, 0, input_.size() - 1, input_);
    else {
        int number_of_leaves = (input_.size()%leaf_length_ == 0) ? input_.size()/leaf_length_ : input_.size()/leaf_length_+1;
        int height =  0;

        int nl = number_of_leaves-1;
        int64_t block_length = leaf_length_;
        while (nl){
            height++;
            block_length*=r_;
            nl/=r_;
        }

        root_block_ = new SumInternalBlock(nullptr, 0, block_length - 1, input_);
    }

}


SumBlockTree::~SumBlockTree() {
    delete root_block_;
}


std::vector<std::vector<SumBlock*>> SumBlockTree::levelwise_iterator() {
    std::vector<std::vector<SumBlock*>> result = {{root_block_}};
    while (!dynamic_cast<SumLeafBlock*>(result.back()[0])) {
        std::vector<SumBlock*> next_level = {};
        for (SumBlock *b : result.back())
            for (SumBlock *child : b->children(leaf_length_, r_))
                next_level.push_back(child);
        result.push_back(next_level);
    }

    return result;
}


void SumBlockTree::clean_unnecessary_expansions() {
    root_block_->clean_unnecessary_expansions();
    for (std::vector<SumBlock*> level : levelwise_iterator()) {
        for (int i = 0; i<level.size(); ++i) {
            level[i]->level_index_ = i;
            level[i]->first_occurrence_level_index_ = level[i]->first_block_->level_index_;
        }
    }
}


void SumBlockTree::add_access_support() {
    root_block_->add_access_support();
}


int64_t SumBlockTree::access(int i) {
    return root_block_->access(i);
}


std::vector<SumBlock*> SumBlockTree::next_level(std::vector<SumBlock*>& level) {
    std::vector<SumBlock*> next_level;
    for (int i = 0; i < level.size(); ++i) {
        SumBlock* b = level[i];
        for (SumBlock *child : b->children(leaf_length_, r_)) { // Do it in order
            child->level_index_ = next_level.size();
            child->first_occurrence_level_index_ = next_level.size();
            next_level.push_back(child);
        }
    }
    return next_level;
}


void SumBlockTree::forward_pair_window_block_scan(std::vector<SumBlock*>& level, int pair_window_size, int N, std::unordered_map<HashString, std::vector<std::pair<SumBlock*, SumBlock*>>>& pair_hashtable) {
    for (std::vector<SumBlock *>::iterator it = level.begin(); it != level.end();) {
        SumBlock *b = (*it);
        b->right_ = true;
        int offset = 0;
        RabinKarp rk(input_, (*it)->start_index_ + offset, pair_window_size, N); // offset is always 0 here
        for (; it != level.end() && ((*it) == b || (*(it-1))->end_index_ == (*it)->start_index_ - 1); it++) {
            SumBlock* current = *(it);
            bool last_block = ((it+1) == level.end() ||  current->end_index_ != (*(it+1))->start_index_ - 1);
            for (offset = 0; offset < current->length(); ++offset) {
                if (last_block && current->length() - offset < pair_window_size)  break;
                HashString hS(rk.hash(), input_, current->start_index_ + offset, current->start_index_ + offset + pair_window_size - 1);
                std::unordered_map<HashString, std::vector<std::pair<SumBlock*, SumBlock*>>>::const_iterator result = pair_hashtable.find(hS);
                if (result != pair_hashtable.end()) { // Here, It could be that the scanning should have finished with the penultimate, but it never should enter this ''if''
                    // when We're on the penultimate block and the window exceeds the last block because if that is a first occurrence should have been occured before in a pair of blocks
                    // maybe use a condition more like rk's condition below could work fine too
                    // Same logic: for when passing a window of size 2l + 2 over 2 block of length l
                    for (std::pair<SumBlock*,SumBlock*> p: result->second) {
                        if (current->start_index_ + offset < p.first->start_index_) {
                            p.first->left_ = true;
                            p.second->right_ = true;
                        }
                    }
                    pair_hashtable.erase(hS);
                }
                if (current->start_index_+offset+pair_window_size < input_.size()) rk.wnext();
            }
        }
        (*(it-1))->left_ = true;
    }
}


void SumBlockTree::forward_window_block_scan(std::vector<SumBlock*>& level, int window_size, int N, std::unordered_map<HashString, std::vector<SumBlock*>>& hashtable) {
    int i = 0;
    for (std::vector<SumBlock *>::iterator it = level.begin(); it != level.end();) {
        SumBlock *b = (*it);
        int offset = 0;
        RabinKarp rk(input_, (*it)->start_index_ + offset, window_size, N);
        for (; it != level.end() && ((*it) == b || (*(it-1))->end_index_ == (*it)->start_index_ - 1); it++, i++) {
            SumBlock* current = *(it);
            bool last_block = ((it+1) == level.end() ||  current->end_index_ != (*(it+1))->start_index_ - 1);
            for (offset = 0; offset < current->length(); ++offset) {
                if (last_block && current->length() - offset < window_size)  break;
                HashString hS(rk.hash(), input_, current->start_index_ + offset, current->start_index_ + offset + window_size - 1);
                std::unordered_map<HashString, std::vector<SumBlock *>>::const_iterator result = hashtable.find(hS);
                if (result != hashtable.end()) {
                    std::vector<SumBlock*> blocks = result->second;
                    for (SumBlock* b : blocks) {
                        b->first_occurrence_level_index_ = i;
                        b->first_block_ = current;
                        b->offset_ = offset;
                        if (offset + window_size > b->first_block_->length()) b->second_block_ = (*(it+1));
                        else b->second_block_ = nullptr;
                    }
                    hashtable.erase(hS);
                }
                if (current->start_index_+offset+window_size < input_.size()) rk.wnext();
            }
        }
    }
}


void SumBlockTree::block_scan(std::vector<SumBlock *>& level, int N , std::unordered_map<HashString, std::vector<SumBlock*>>& hashtable) {
    for (SumBlock* b : level) {
        RabinKarp rk(input_, b->start_index_, b->length(), N);
        HashString hS(rk.hash(),  input_, b->start_index_, b->end_index_);

        std::unordered_map<HashString, std::vector<SumBlock*>>::const_iterator result = hashtable.find(hS);

        if (result == hashtable.end())
            hashtable[hS] = {b};
        else
            hashtable[hS].push_back(b);
    }
}


void SumBlockTree::process_level(std::vector<SumBlock*>& level) {

    int N = 6700417; //Large prime
    int level_length = level.front()->length();

    // SumBlock scan
    std::unordered_map<HashString, std::vector<SumBlock*>> hashtable;
    block_scan(level, N, hashtable);

    // Pairs of blocks scan
    std::unordered_map<HashString, std::vector<std::pair<SumBlock *,SumBlock*>>> pair_hashtable;
    for (std::vector<SumBlock *>::iterator it = level.begin(); it != level.end();) {
        for (++it; (it != level.end() && (*(it-1))->end_index_ == (*it)->start_index_ - 1); ++it) {
            SumBlock* current = (*(it - 1));
            SumBlock* next = (*it);
            RabinKarp rk(input_, current->start_index_, current->length() + next->length(), N);
            HashString hS(rk.hash(), input_, current->start_index_, current->start_index_ + current->length() + next->length()-1); // Second parameter is next->end_index

            std::unordered_map<HashString, std::vector<std::pair<SumBlock *,SumBlock*>>>::const_iterator result = pair_hashtable.find(hS);

            if (result == pair_hashtable.end())
                pair_hashtable[hS] = {{current, next}};
            else
                pair_hashtable[hS].push_back({current, next});
        }
    }


    // Window block scan
    //Establishes first occurrences of blocks
    forward_window_block_scan(level, level_length, N, hashtable);



    // Window Pair of blocks scans
    if (level.size() > 1)
        forward_pair_window_block_scan(level, level_length*2, N, pair_hashtable);




    // SumBackBlock creation
    for (int i = 0; i < level.size(); ++i) {
        SumBlock* b = level[i];
        if (b->left_ && b->right_ && b->first_occurrence_level_index_ < b->level_index_) {
            // This doesn't have the bug of the dangling reference fixed with first_occurrence_level_index, because it shouldn't happen that
            // A block points back to a SumBackBlock
            SumBackBlock* bb = new SumBackBlock(b->parent_, b->start_index_, b->end_index_, input_,
                                          level[b->first_occurrence_level_index_], (b->second_block_ ==
                                                                                    nullptr) ? nullptr : level[b->first_occurrence_level_index_ +1], b->offset_);
            bb->level_index_ = b->level_index_;
            bb->first_occurrence_level_index_ = b->first_occurrence_level_index_;
            bb->left_ = true;
            bb->right_ = true;
            b->parent_->replace_child(b, bb);
            delete b;
            level[i] = bb;
        }
    }

}


void SumBlockTree::process_back_pointers() {
    std::vector<SumBlock*> current_level = {root_block_};
    std::stack<SumBlock*> none_blocks;
    while ((current_level = next_level(current_level)).size() != 0) {
        if (current_level[0]->length() < r_ ||  current_level[0]->length() <= leaf_length_) break;
        while (current_level.size() != 0 && current_level.back()->end_index_ >= input_.size()) {
            none_blocks.push(current_level.back());
            current_level.pop_back();
        }
        process_level(current_level);
        while (!none_blocks.empty()) {
            current_level.push_back(none_blocks.top());
            none_blocks.pop();
        }
    }
}


void SumBlockTree::process_level_heuristic(std::vector<SumBlock*>& level) {

    int N = 6700417; // Large prime
    int level_length = level.front()->length();

    // SumBlock scan
    std::unordered_map<HashString, std::vector<SumBlock*>> hashtable;
    block_scan(level, N, hashtable);

    // Window block scan
    // This is almost the same as forward_window_block_scan, as well as the SumBackBlock creation
    forward_window_block_scan(level, level_length, N, hashtable);



    // SumBackBlock creation
    for (int i = 0; i < level.size(); ++i) {
        SumBlock* b = level[i];
        if (b->first_occurrence_level_index_ < b->level_index_) {

            SumBackBlock* bb = new SumBackBlock(b->parent_, b->start_index_, b->end_index_, input_,
                                          level[b->first_occurrence_level_index_], (b->second_block_ ==
                                                                                    nullptr) ? nullptr : level[b->first_occurrence_level_index_ +1], b->offset_);
            bb->level_index_ = b->level_index_;
            bb->first_occurrence_level_index_ = b->first_occurrence_level_index_;
            b->parent_->replace_child(b, bb);
            delete b;
            level[i] = bb;
        }
    }

}


void SumBlockTree::process_back_pointers_heuristic() {
    std::vector<SumBlock *> current_level = {root_block_};
    std::stack<SumBlock*> none_blocks;
    while ((current_level = next_level(current_level)).size() != 0) {
        if (current_level[0]->length() < r_ ||
            current_level[0]->length() <= leaf_length_)
            break;
        while (current_level.size() != 0 && current_level.back()->end_index_ >= input_.size()) {
            none_blocks.push(current_level.back());
            current_level.pop_back();
        }
        process_level_heuristic(current_level);
        while (!none_blocks.empty()) {
            current_level.push_back(none_blocks.top());
            none_blocks.pop();
        }
    }
}