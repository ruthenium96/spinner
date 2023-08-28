#include "AllSymmetricMatrixFactories.h"

std::vector<std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory>>
constructAllSymmetricMatrixFactories() {
    return {
        std::make_shared<quantum::linear_algebra::ArmaDenseTransformAndDiagonalizeFactory>(),
        std::make_shared<quantum::linear_algebra::EigenSymmetricMatrixFactory>()};
}
