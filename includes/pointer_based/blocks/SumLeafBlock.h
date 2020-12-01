#ifndef MINIMALISTBT_CST_SUMLEAFBLOCK_H
#define MINIMALISTBT_CST_SUMLEAFBLOCK_H

#include "SumBlock.h"

class SumLeafBlock : public SumBlock {
public:

    SumLeafBlock(SumBlock*, int64_t, int64_t, std::basic_string<int64_t >&);
    ~SumLeafBlock();

    int64_t size();

    int64_t add_access_support();
    int64_t access(int);
};

#endif //MINIMALISTBT_CST_SUMLEAFBLOCK_H
