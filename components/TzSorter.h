#ifndef JULY_TZSORTER_H
#define JULY_TZSORTER_H

#include "entities/space/Space.h"
#include <numeric>
#include <utility>

class TzSorter {
public:
    explicit TzSorter(spaces::LexicographicIndexConverter converter);

    Space apply(Space&& space) const;

private:
    const spaces::LexicographicIndexConverter converter_;
    uint32_t max_ntz_proj;
};

#endif //JULY_TZSORTER_H
