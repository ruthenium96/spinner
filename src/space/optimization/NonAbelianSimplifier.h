#ifndef SPINNER_NON_ABELIAN_SIMPLIFIER_H
#define SPINNER_NON_ABELIAN_SIMPLIFIER_H

#include "src/linalg_structures/FactoriesList.h"
#include "src/space/Space.h"

namespace spinner::space::optimization {
class NonAbelianSimplifier {
  public:
    explicit NonAbelianSimplifier(linalg_structures::FactoriesList factories);
    Space apply(Space&& space) const;

  private:
    const linalg_structures::FactoriesList factories_;
};
} // namespace spinner::space::optimization

#endif  //SPINNER_NON_ABELIAN_SIMPLIFIER_H
