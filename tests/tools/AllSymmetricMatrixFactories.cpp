#include "AllSymmetricMatrixFactories.h"

#include "src/entities/data_structures/arma/ArmaFactories.h"
#include "src/entities/data_structures/eigen/EigenFactories.h"

std::vector<std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory>>
constructAllDenseTransformAndDiagonalizeFactories() {
    auto armaDouble = std::make_shared<quantum::linear_algebra::ArmaDenseTransformAndDiagonalizeFactory>();
    armaDouble->setPrecision(quantum::linear_algebra::DOUBLE);
    auto armaFloat = std::make_shared<quantum::linear_algebra::ArmaDenseTransformAndDiagonalizeFactory>();
    armaFloat->setPrecision(quantum::linear_algebra::SINGLE);
    auto eigen = std::make_shared<quantum::linear_algebra::EigenDenseTransformAndDiagonalizeFactory>();

    return {armaDouble, armaFloat, eigen};
}
