#ifndef MINIMALISTBT_CST_BTCST_H
#define MINIMALISTBT_CST_BTCST_H

#include <compressed/BTCT.h>
#include <LCP_RLCSA.h>
#include <rlcsa/rlcsa.h>


class BTCST {
public:

    int64_t node_ = -2;
    int64_t i_ = -2;
    int64_t psi_ = -2;
    int64_t* children;
    static const int PAPER = 0;
    static const int PAPER_NO_CLEAN = 1;
    static const int HEURISTIC = 2;

    TextIndexes::RLCSA* rlcsa_;
    TextIndexRLCSA* index_;
    cds_static::LCP_RLCSA* lcp_rlcsa_;
    BTCT* btct_;

    int n_;

    BTCST(std::string&, int = PAPER, int = 2, int = 16, int = 32, int = 128);
    //BTCST(std::string&, std::string&, int = PAPER, int = 2, int = 16, int = 32, int = 128);
    ~BTCST();


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
    int string(int, int);
};


#endif //MINIMALISTBT_CST_BTCST_H
