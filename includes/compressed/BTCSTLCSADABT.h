#ifndef MINIMALISTBT_CST_BTCSTLCSADABT_H
#define MINIMALISTBT_CST_BTCSTLCSADABT_H

#include <compressed/BTCT.h>
#include <compressed/DABT.h>
#include <LCP_RLCSA.h>
#include <rlcsa/rlcsa.h>
#include <lcsa/LCSALENGTHS.hpp>


class BTCSTLCSADABT {
public:

    int64_t node_ = -2;
    int64_t i_ = -2;
    int64_t psi_ = -2;
    int64_t* children;
    static const int PAPER = 0;
    static const int PAPER_NO_CLEAN = 1;
    static const int HEURISTIC = 2;
    static const int disa_constant = 32;

    TextIndexes::RLCSA* rlcsa_;
    TextIndexRLCSA* index_;
    LCSALENGTHS::LCSALENGTHS* dsa_;
    DABT* disa_;
    cds_static::LCP_RLCSA* lcp_rlcsa_;
    BTCT* btct_;

    int n_;
    int sigma_;

    BTCSTLCSADABT(std::string&, int = PAPER, int = 2, int = 16, int = 32);
    BTCSTLCSADABT(std::ifstream&);
    ~BTCSTLCSADABT();


    int first_child(int);
    int tree_depth(int);
    int next_sibling(int);
    int parent(int);
    int level_ancestor(int, int);
    int lca(int, int);
    int suffix_link(int);
    int string_depth(int);
    int string_ancestor(int, int);
    int child(int, int);
    int bin_search_child(int, int);
    int string(int, int);

    int size();
    void serialize(std::ofstream&);
};


#endif //MINIMALISTBT_CST_BTCSTLCSADABT_H
