#ifndef SPINNER_OPTIMIZEDSPACECONSTRUCTOR_H
#define SPINNER_OPTIMIZEDSPACECONSTRUCTOR_H

#include "src/common/runner/ConsistentModelOptimizationList.h"
#include "src/space/Space.h"
namespace space::optimization {
// This class constructs Space using ConsistentModelOptimizationList.
// TODO: Can it be rewritten as function?
class OptimizedSpaceConstructor {
  public:
    static Space construct(const runner::ConsistentModelOptimizationList&);
};
}  // namespace space::optimization

#endif  //SPINNER_OPTIMIZEDSPACECONSTRUCTOR_H
