#include "AllDenseFactories.h"

std::vector<std::shared_ptr<quantum::linear_algebra::AbstractDenseFactory>>
constructAllDenseFactories() {
    return {
        std::make_shared<quantum::linear_algebra::ArmaDenseFactory>(),
        std::make_shared<quantum::linear_algebra::EigenDenseFactory>()};
}
