#ifndef JULY_TZ_SORTER_H
#define JULY_TZ_SORTER_H

#include "entities/Space.h"
#include <numeric>
#include <utility>

class Tz_Sorter {
public:
    explicit Tz_Sorter(const spaces::LexicographicIndexWorker& indexes);

    Space& operator()(Space& space) const;

private:
    const spaces::LexicographicIndexWorker& indexes_;
    uint32_t max_ntz_proj;
};

#endif //JULY_TZ_SORTER_H
