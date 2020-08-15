#include <pointer_based/BlockTree.h>
#include <compressed/CBitBlockTree.h>

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <compressed/BTCT.h>
#include <compressed/BTCST.h>
#include <compressed/BTCSTLCSA.h>


int main() {

    std::string input;
    std::ifstream t("../../tests/data/dna.par");
    std::stringstream buffer;
    buffer << t.rdbuf();
    input = buffer.str();
    std::cout << input.length() << std::endl;

    std::unordered_set<int> characters;
    for (char c: input) {
        characters.insert(c);
    }


    BlockTree* bt = new BlockTree(input, 2, 32);
    bt->process_back_pointers();
    bt->clean_unnecessary_expansions();
    for (char c: characters)
        bt->add_rank_select_support(c);
    bt->add_leaf_rank_select_support();
    bt->add_search_support();


    CBitBlockTree* cbt = new CBitBlockTree(bt, input[0]);
    cbt->access(0);
    std::cout << cbt->rank_1(10) << std::endl;

    std::ofstream ot("dna.par.bt");
    cbt->serialize(ot);
    ot.close();
    std::ifstream it("dna.par.bt");
    CBitBlockTree* lcbt = new CBitBlockTree(it);
    lcbt->access(0);
    std::cout << lcbt->rank_1(10) << std::endl;

    BTCT* btct = new BTCT(bt, input[0]);
    btct->access(0);
    std::cout << btct->next_sibling(1) << std::endl;

    std::ofstream ot2("dna.par.btct");
    btct->serialize(ot2);
    ot2.close();
    std::ifstream it2("dna.par.btct");
    BTCT* lbtct = new BTCT(it2);
    lbtct->access(0);
    std::cout << btct->next_sibling(1) << std::endl;

    std::string banana = "banananananananana";
    BTCST* btcst = new BTCST(banana);
    std::cout << btcst->string_depth(6) << std::endl;
    std::ofstream ott("banana.btcst");
    btcst->serialize(ott);
    ott.close();
    std::ifstream itt("banana.btcst");
    BTCST* lbtcst = new BTCST(itt);
    std::cout << lbtcst->string_depth(6) << std::endl;


    BTCSTLCSA* btcstlcsa = new BTCSTLCSA(banana);
    std::cout << btcstlcsa->suffix_link(4) << std::endl;
    std::ofstream ott2("banana.btcstlcsa");
    btcstlcsa->serialize(ott2);
    ott2.close();
    std::ifstream itt2("banana.btcstlcsa");
    BTCSTLCSA* lbtcstlcsa = new BTCSTLCSA(itt2);
    std::cout << lbtcstlcsa->suffix_link(4) << std::endl;
    
    delete bt;
    delete cbt;
    delete lcbt;
    delete btct;
    delete lbtct;
    delete btcst;
    delete lbtcst;
    delete btcstlcsa;
    delete lbtcstlcsa;
    return 0;
}

