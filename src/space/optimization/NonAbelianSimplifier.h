#ifndef JULY_NON_ABELIAN_SIMPLIFIER_H
#define JULY_NON_ABELIAN_SIMPLIFIER_H

#include "src/space/Space.h"

namespace space::optimization {
class NonAbelianSimplifier {
  public:
    Space apply(Space&& space) const;
};
}  // namespace space::optimization

#endif  //JULY_NON_ABELIAN_SIMPLIFIER_H
