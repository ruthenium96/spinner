#ifndef JULY_TZSORTER_H
#define JULY_TZSORTER_H

#include <numeric>
#include <utility>

#include "src/entities/space/Space.h"

class TzSorter {
  public:
    explicit TzSorter(lexicographic::IndexConverter converter);

    Space apply(Space&& space) const;

  private:
    const lexicographic::IndexConverter converter_;
    uint32_t max_ntz_proj;
};

#endif  //JULY_TZSORTER_H
