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
        std::function<uint8_t(uint32_t)> index_to_tz_projection_functor,
        uint32_t max_ntz_proj,
        quantum::linear_algebra::FactoriesList factories);

    Space apply(Space&& space) const;

  private:
    std::function<uint8_t(uint32_t)> index_to_tz_projection_functor_;
    std::shared_ptr<const lexicographic::IndexConverter> converter_;
    const quantum::linear_algebra::FactoriesList factories_;
    uint32_t max_ntz_proj_;
};
}  // namespace space::optimization

#endif  //SPINNER_TZSORTER_H
