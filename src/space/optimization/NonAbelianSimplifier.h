#ifndef SPINNER_NON_ABELIAN_SIMPLIFIER_H
#define SPINNER_NON_ABELIAN_SIMPLIFIER_H

#include "src/entities/data_structures/FactoriesList.h"
#include "src/space/Space.h"

namespace spinner::space::optimization {
class NonAbelianSimplifier {
  public:
    explicit NonAbelianSimplifier(linear_algebra::FactoriesList factories);
    Space apply(Space&& space) const;

  private:
    const linear_algebra::FactoriesList factories_;
};
} // namespace spinner::space::optimization

#endif  //SPINNER_NON_ABELIAN_SIMPLIFIER_H
