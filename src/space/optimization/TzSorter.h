#ifndef SPINNER_TZSORTER_H
#define SPINNER_TZSORTER_H

#include <numeric>
#include <utility>

#include "src/space/Space.h"

namespace space::optimization {
class TzSorter {
  public:
    explicit TzSorter(lexicographic::IndexConverter converter);

    Space apply(Space&& space) const;

  private:
    const lexicographic::IndexConverter converter_;
    uint32_t max_ntz_proj;
};
}  // namespace space::optimization

#endif  //SPINNER_TZSORTER_H
