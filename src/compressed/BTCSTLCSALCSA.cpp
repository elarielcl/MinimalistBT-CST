#include <compressed/BTCSTLCSALCSA.h>
#include <sdsl/suffix_trees.hpp>

BTCSTLCSALCSA::BTCSTLCSALCSA(std::string& input_string, int block_tree_version, int r, int leaf_length,
                     int block_size) {

    sdsl::cst_sada<> cst;
    construct_im(cst, input_string, 1);
    std::string topology;
    for (int i = 0; i < cst.bp.size(); ++i) {
        if (cst.bp[i]) topology += "(";
        else topology += ")";
    }

    sdsl::csa_bitcompressed<> csa;
    sdsl::construct_im(csa, input_string, 1);
    dsa_ = new LCSALENGTHS::LCSALENGTHS(csa.sa_sample, 128); // Fix Sweet point
    sdsl::int_vector<> isa(csa.size(),0);
    for (int i = 0; i < isa.size(); ++i) {
        isa[csa[i]] = i;
    }
    disa_ = new LCSALENGTHS::LCSALENGTHS(isa, 128);


    std::set<int> alphabet;
    for (char c : input_string) {
        alphabet.insert(c);
    }
    sigma_ = alphabet.size();
    children = new int64_t[sigma_+1];

    BlockTree* bt_top = new BlockTree(topology, r, leaf_length);
    switch (block_tree_version) {
        case PAPER:
            bt_top->process_back_pointers();
            bt_top->clean_unnecessary_expansions();
            break;
        case PAPER_NO_CLEAN:
            bt_top->process_back_pointers();
            break;
        case HEURISTIC:
            bt_top->process_back_pointers_heuristic();
            break;
    }
    bt_top->add_rank_select_support('(');
    bt_top->add_rank_select_support(')');
    bt_top->add_leaf_rank_select_support();
    bt_top->add_search_support();


    btct_ = new BTCT(bt_top, '(');
    delete bt_top;

    unsigned char* data = new unsigned char[input_string.size()+1];
    for (int i = 0; i < input_string.size(); ++i) {
        data[i] = input_string[i];
    }
    data[input_string.size()] = 0;

    TextIndexes::RLCSA* temp_rlcsa = new TextIndexes::RLCSA(data, input_string.size()+1, block_size, 2, false, false);;
    lcp_rlcsa_ = new LCP_FMN_RLCSA(temp_rlcsa, (char*)data, input_string.size()+1);

    temp_rlcsa->support_display = false;
    temp_rlcsa->support_locate = false;
    temp_rlcsa->sample_rate = 1;

    std::ofstream temp_file_out("save.tf");
    temp_rlcsa->save(temp_file_out);
    temp_file_out.close();

    delete temp_rlcsa;
    delete[] data;

    std::ifstream temp_file_in("save.tf");
    rlcsa_ = new TextIndexes::RLCSA(temp_file_in);
    temp_file_in.close();
    remove("save.tf");

    index_ = new TextIndexRLCSA();

    index_->rlcsa = rlcsa_;
    n_ = input_string.size()+1;
}

BTCSTLCSALCSA::BTCSTLCSALCSA(std::ifstream& in) {
    in.read((char *) &n_, sizeof(int));
    in.read((char *) &sigma_, sizeof(int));
    children = new int64_t[sigma_+1];
    btct_ = new BTCT(in);
    rlcsa_ = RLCSA::load(in);
    index_ = new TextIndexRLCSA();
    index_->rlcsa = rlcsa_;
    lcp_rlcsa_ = LCP_FMN_RLCSA::load(in);
    dsa_ = new LCSALENGTHS::LCSALENGTHS(in);
    disa_ = new LCSALENGTHS::LCSALENGTHS(in);
}


BTCSTLCSALCSA::~BTCSTLCSALCSA() {
    delete lcp_rlcsa_;
    delete index_;
    delete btct_;
    delete dsa_;
    delete disa_;
    delete[] children;
}


int BTCSTLCSALCSA::first_child(int node) {
    return btct_->first_child(node);
}


int BTCSTLCSALCSA::tree_depth(int node) {
    return btct_->tree_depth(node);
}


int BTCSTLCSALCSA::next_sibling(int node) {
    return btct_->next_sibling(node);
}


int BTCSTLCSALCSA::parent(int node) {
    return btct_->parent(node);
}


int BTCSTLCSALCSA::level_ancestor(int node, int d) {
    return btct_->level_ancestor(node, d);
}


int BTCSTLCSALCSA::lca(int node1, int node2) {
    return btct_->lca(node1, node2);
}


int BTCSTLCSALCSA::suffix_link(int node) {
    if (node == 0 || node == 1) return 0;

    if (btct_->is_leaf(node)) {
        int lml = btct_->lb(node);
        int psi_lml = rlcsa_->getPsi(lml,1);
        return btct_->leaf_select(psi_lml+1);
    }

    int lml = btct_->lb(node);
    int rml = btct_->rb(node);
    if (lml == 0 && rml == 1) return 0;
    int psi_lml = rlcsa_->getPsi(lml,1);
    int psi_rml = rlcsa_->getPsi(rml-1,1);


    return lca(btct_->leaf_select(psi_lml+1), btct_->leaf_select(psi_rml+1));
}


int BTCSTLCSALCSA::string_depth(int node) {
    if (node == 0) return 0;

    int lr = 0;
    if (btct_->is_leaf_rank(node, lr)) {
        return n_ - (*dsa_)[lr-1];
    }

    int sc = btct_->fwdsearch(node+1, -1) + 1;
    int i = btct_->leaf_rank(sc);
    return lcp_rlcsa_->get_LCP(i, (*dsa_)[i]);
}


int BTCSTLCSALCSA::string_ancestor(int node, int d) {

    int initial_bottom_depth = btct_->tree_depth(node);
    int bottom_depth = initial_bottom_depth;
    int top_depth = 0;
    int top = 0;
    int bottom = node;

    while (true) {
        if (bottom == top) return bottom;
        int guess_depth = top_depth + (bottom_depth-top_depth)/2;
        int middle_node = btct_->level_ancestor(node, initial_bottom_depth  - guess_depth);
        int s_depth = string_depth(middle_node);
        if (s_depth == d) return middle_node;
        if (s_depth < d) {
            if (top_depth == bottom_depth-1) return bottom;
            top = middle_node;
            top_depth = guess_depth;
        } else {
            if (top_depth == bottom_depth-1) return top;
            bottom = middle_node;
            bottom_depth = guess_depth;
        }
    }

}


int BTCSTLCSALCSA::child(int node, int c) {
    if (c == 0) c = '$';
    int sdepth = string_depth(node);
    int child = node+1;
    while (true) {
        int ch = string(child, sdepth);
        if (ch > c) return -1;
        if (ch == c) return child;
        child = btct_->next_sibling(child);
        if (child == -1) return -1;
    }
}


int BTCSTLCSALCSA::bin_search_child(int node, int c) {
    if (c == 0) c = '$';
    int sdepth = string_depth(node);
    int child = node+1;
    int j = -1;
    while (true) {
        children[++j] = child;
        child = btct_->next_sibling(child);
        if (child == -1) break;
    }

    int min = 0;
    int max = j;
    while (min <= max) {
        // Key is in a[lo..hi] or not present.
        int mid = min + (max - min) / 2;
        int ch = string(children[mid], sdepth);
        if      (c < ch) max = mid - 1;
        else if (c > ch) min = mid + 1;
        else return children[mid];
    }
    return -1;
}


int BTCSTLCSALCSA::string(int node, int i) {
    if (node == node_ && i == i_+1) {
        ++i_;
        psi_ = rlcsa_->getPsi(psi_,1);
        int r = rlcsa_->getT(psi_);
        if (r == 0) return '$';
        return r;
    } else {
        node_ = node;
        i_ = i;
        int lr = btct_->leaf_rank(node);
        if (i > disa_constant) psi_ = (*disa_)[(*dsa_)[lr]+i];
        else psi_ = rlcsa_->getPsi(lr,i);

        int r = rlcsa_->getT(psi_);
        if (r == 0) return '$';
        return r;
    }
}


int BTCSTLCSALCSA::size() {
    return btct_->size() + rlcsa_->getSize() + lcp_rlcsa_->getSize() + dsa_->size_in_bytes() + disa_->size_in_bytes();
}


void BTCSTLCSALCSA::serialize(std::ofstream& out) {
    out.write((char *) &n_, sizeof(int));
    out.write((char *) &sigma_, sizeof(int));
    btct_->serialize(out);
    rlcsa_->save(out);
    lcp_rlcsa_->save(out);
    dsa_->serialize(out);
    disa_->serialize(out);
}