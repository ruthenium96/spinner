#ifndef SPINNER_NONMINIMALPROJECTIONSELIMINATOR_H
#define SPINNER_NONMINIMALPROJECTIONSELIMINATOR_H

#include "src/space/Space.h"

namespace space::optimization {
class NonMinimalProjectionsEliminator {
  public:
    explicit NonMinimalProjectionsEliminator(uint32_t max_total_mult);
    Space apply(Space&& space) const;

  private:
    uint32_t max_total_mult_;
};
} // namespace space::optimization

#endif // SPINNER_NONMINIMALPROJECTIONSELIMINATOR_H