#ifndef SPINNER_WORKERCONSTRUCTOR_H
#define SPINNER_WORKERCONSTRUCTOR_H

#include <memory>

#include "AbstractWorker.h"
#include "src/common/runner/ConsistentModelOptimizationList.h"
#include "src/eigendecompositor/FlattenedSpectra.h"

namespace magnetic_susceptibility::worker {

class WorkerConstructor {
  public:
    static std::unique_ptr<AbstractWorker> construct(
        const runner::ConsistentModelOptimizationList& consistentModelOptimizationList,
        std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra,
        const quantum::linear_algebra::FactoriesList& factories);
};

}  // namespace magnetic_susceptibility::worker

#endif  //SPINNER_WORKERCONSTRUCTOR_H
