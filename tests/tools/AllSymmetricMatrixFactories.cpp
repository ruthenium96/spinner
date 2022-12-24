#include "AllSymmetricMatrixFactories.h"

std::vector<std::shared_ptr<quantum::linear_algebra::AbstractSymmetricMatrixFactory>>
constructAllSymmetricMatrixFactories() {
    return {
        std::make_shared<quantum::linear_algebra::ArmaSymmetricMatrixFactory>(),
        std::make_shared<quantum::linear_algebra::EigenSymmetricMatrixFactory>()};
}
