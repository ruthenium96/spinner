#ifndef JULY_NON_ABELIAN_SIMPLIFIER_H
#define JULY_NON_ABELIAN_SIMPLIFIER_H

#include "src/entities/space/Space.h"

class NonAbelianSimplifier {
  public:
    Space apply(Space&& space) const;
};

#endif  //JULY_NON_ABELIAN_SIMPLIFIER_H
