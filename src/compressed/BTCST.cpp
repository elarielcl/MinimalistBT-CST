#include <compressed/BTCST.h>
#include <sdsl/suffix_trees.hpp>

BTCST::BTCST(std::string& input_string, int block_tree_version, int r, int leaf_length,
                         int block_size, int sa_sampling_rate) {

    sdsl::cst_sada<> cst;
    construct_im(cst, input_string, 1);
    std::string topology;
    for (int i = 0; i < cst.bp.size(); ++i) {
        if (cst.bp[i]) topology += "(";
        else topology += ")";
    }

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

    rlcsa_ = new TextIndexes::RLCSA(data, input_string.size()+1, block_size, sa_sampling_rate, false, false);
    index_ = new TextIndexRLCSA();

    index_->rlcsa = rlcsa_;

    lcp_rlcsa_ = new LCP_FMN_RLCSA(rlcsa_, (char*)data, input_string.size()+1);
    delete[] data;

    n_ = input_string.size()+1;
}


BTCST::BTCST(std::ifstream& in) {
    in.read((char *) &n_, sizeof(int));
    in.read((char *) &sigma_, sizeof(int));
    children = new int64_t[sigma_+1];
    btct_ = new BTCT(in);
    rlcsa_ = RLCSA::load(in);
    index_ = new TextIndexRLCSA();
    index_->rlcsa = rlcsa_;
    lcp_rlcsa_ = LCP_FMN_RLCSA::load(in);
}


BTCST::~BTCST() {
    delete lcp_rlcsa_;
    delete index_;
    delete btct_;
    delete[] children;
}


int BTCST::first_child(int node) {
    return btct_->first_child(node);
}


int BTCST::tree_depth(int node) {
    return btct_->tree_depth(node);
}


int BTCST::next_sibling(int node) {
    return btct_->next_sibling(node);
}


int BTCST::parent(int node) {
    return btct_->parent(node);
}


int BTCST::level_ancestor(int node, int d) {
    return btct_->level_ancestor(node, d);
}


int BTCST::lca(int node1, int node2) {
    return btct_->lca(node1, node2);
}


int BTCST::suffix_link(int node) {
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


int BTCST::string_depth(int node) {
    if (node == 0) return 0;

    int lr = 0;
    if (btct_->is_leaf_rank(node, lr)) {
        return n_ - rlcsa_->locate(lr-1);
    }

    int sc = btct_->fwdsearch(node+1, -1) + 1;
    return lcp_rlcsa_->get_LCP(btct_->leaf_rank(sc), index_);
}


int BTCST::string_ancestor(int node, int d) {

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


int BTCST::child(int node, int c) {
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


int BTCST::bin_search_child(int node, int c) {
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


int BTCST::string(int node, int i) {
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
        psi_ = rlcsa_->getPsi(lr,i);

        int r = rlcsa_->getT(psi_);
        if (r == 0) return '$';
        return r;
    }
}


int BTCST::size() {
    return btct_->size() + rlcsa_->getSize() + lcp_rlcsa_->getSize();
}


void BTCST::serialize(std::ofstream& out) {
    out.write((char *) &n_, sizeof(int));
    out.write((char *) &sigma_, sizeof(int));
    btct_->serialize(out);
    rlcsa_->save(out);
    lcp_rlcsa_->save(out);
}