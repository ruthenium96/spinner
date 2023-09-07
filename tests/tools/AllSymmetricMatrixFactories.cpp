#include "AllSymmetricMatrixFactories.h"

#include "src/entities/data_structures/arma/ArmaFactories.h"
#include "src/entities/data_structures/eigen/EigenFactories.h"

std::vector<std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory>>
constructAllDenseTransformAndDiagonalizeFactories() {
    return {
        std::make_shared<quantum::linear_algebra::ArmaDenseTransformAndDiagonalizeFactory>(),
        std::make_shared<quantum::linear_algebra::EigenDenseTransformAndDiagonalizeFactory>()};
}
