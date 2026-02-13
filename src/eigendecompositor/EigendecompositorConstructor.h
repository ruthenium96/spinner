#ifndef SPINNER_EIGENDECOMPOSITORCONSTRUCTOR_H
#define SPINNER_EIGENDECOMPOSITORCONSTRUCTOR_H

#include "AbstractEigendecompositor.h"
#include "src/common/runner/ConsistentModelOptimizationList.h"

namespace spinner::eigendecompositor {

class EigendecompositorConstructor {
  public:
    static std::unique_ptr<AbstractEigendecompositor> construct(
        const runner::ConsistentModelOptimizationList& consistentModelOptimizationList,
        const linalg_structures::FactoriesList& factories);

};

} // namespace spinner::eigendecompositor

#endif  //SPINNER_EIGENDECOMPOSITORCONSTRUCTOR_H
