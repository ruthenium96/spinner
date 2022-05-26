#ifndef SPINNER_POSITIVEPROJECTIONSELIMINATOR_H
#define SPINNER_POSITIVEPROJECTIONSELIMINATOR_H

#include "src/space/Space.h"

namespace space::optimization {

class PositiveProjectionsEliminator {
  public:
    explicit PositiveProjectionsEliminator(uint32_t max_ntz_proj);
    Space apply(Space&& space) const;

  private:
    uint32_t max_ntz_proj_;
};
}  // namespace space::optimization

#endif  //SPINNER_POSITIVEPROJECTIONSELIMINATOR_H
