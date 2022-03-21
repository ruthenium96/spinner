#ifndef JULY_POSITIVEPROJECTIONSELIMINATOR_H
#define JULY_POSITIVEPROJECTIONSELIMINATOR_H

#include "src/entities/space/Space.h"

class PositiveProjectionsEliminator {
  public:
    explicit PositiveProjectionsEliminator(uint32_t max_ntz_proj);
    Space apply(Space&& space) const;

  private:
    uint32_t max_ntz_proj_;
};

#endif  //JULY_POSITIVEPROJECTIONSELIMINATOR_H
