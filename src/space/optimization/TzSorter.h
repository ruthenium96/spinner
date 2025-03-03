#ifndef SPINNER_TZSORTER_H
#define SPINNER_TZSORTER_H

#include <numeric>
#include <utility>

#include "src/entities/data_structures/FactoriesList.h"
#include "src/space/Space.h"
#include "src/common/index_converter/AbstractIndexConverter.h"

namespace space::optimization {
class TzSorter {
  public:
    explicit TzSorter(
        std::shared_ptr<const index_converter::AbstractIndexConverter> indexConverter,
        quantum::linear_algebra::FactoriesList factories);

    Space apply(Space&& space) const;

  private:
    std::shared_ptr<const index_converter::AbstractIndexConverter> converter_;
    const quantum::linear_algebra::FactoriesList factories_;
    uint32_t max_ntz_proj_;
};
}  // namespace space::optimization

#endif  //SPINNER_TZSORTER_H
