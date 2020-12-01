#ifndef MINIMALISTBT_CST_SUMINTERNALBLOCK_H
#define MINIMALISTBT_CST_SUMINTERNALBLOCK_H

#include "SumBlock.h"

class SumInternalBlock : public SumBlock {
public:

    SumInternalBlock(SumBlock*, int64_t, int64_t, std::basic_string<int64_t >&);
    ~SumInternalBlock();

    std::vector<SumBlock*>& children(int, int);
    void clean_unnecessary_expansions();

    bool is_leaf();

    int64_t add_access_support();
    int64_t access(int);
};

#endif //MINIMALISTBT_CST_SUMINTERNALBLOCK_H
