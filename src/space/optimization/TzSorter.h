#ifndef SPINNER_TZSORTER_H
#define SPINNER_TZSORTER_H

#include <numeric>
#include <utility>

#include "src/entities/data_structures/FactoriesList.h"
#include "src/space/Space.h"

namespace space::optimization {
class TzSorter {
  public:
    explicit TzSorter(
        lexicographic::IndexConverter converter,
        quantum::linear_algebra::FactoriesList factories);

    Space apply(Space&& space) const;

  private:
    const lexicographic::IndexConverter converter_;
    const quantum::linear_algebra::FactoriesList factories_;
    uint32_t max_ntz_proj;
};
}  // namespace space::optimization

#endif  //SPINNER_TZSORTER_H
