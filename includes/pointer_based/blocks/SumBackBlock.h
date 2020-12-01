#ifndef MINIMALISTBT_CST_SUMBACKBLOCK_H
#define MINIMALISTBT_CST_SUMBACKBLOCK_H

#include "SumBlock.h"

class SumBackBlock : public SumBlock {
public:

    SumBackBlock(SumBlock*, int64_t, int64_t, std::basic_string<int64_t >&, SumBlock*, SumBlock*, int);
    ~SumBackBlock();

    int64_t add_access_support();
    int64_t access(int);
};

#endif //MINIMALISTBT_CST_SUMBACKBLOCK_H
