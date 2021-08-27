#ifndef JULY_TZ_SORTER_H
#define JULY_TZ_SORTER_H

#include "common/Space.h"
#include <numeric>
#include <utility>

class Tz_Sorter {
public:
    explicit Tz_Sorter(const Spaces::Indexes& indexes);

    Space& operator()(Space& space) const;

private:
    const Spaces::Indexes& indexes_;
    int max_ntz_proj;
};

#endif //JULY_TZ_SORTER_H
