#include "ConsistentModelOptimizationList.h"

#include <utility>

#include "ModelOptimizationListConsistence.h"

namespace runner {

const model::Model& ConsistentModelOptimizationList::getModel() const {
    return model_;
}

model::Model& ConsistentModelOptimizationList::getModel() {
    return model_;
}

const common::physical_optimization::OptimizationList&
ConsistentModelOptimizationList::getOptimizationList() const {
    return optimizationList_;
}

ConsistentModelOptimizationList::ConsistentModelOptimizationList(
    model::Model model,
    common::physical_optimization::OptimizationList optimizationList) :
    model_(std::move(model)),
    optimizationList_(std::move(optimizationList)) {
    // this function throw an exception if Model and OptimizationList are inconsistent
    ModelOptimizationListConsistence::check(model_, optimizationList_);
}
}  // namespace runner