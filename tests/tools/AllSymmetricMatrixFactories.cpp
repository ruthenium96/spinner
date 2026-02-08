#include "AllSymmetricMatrixFactories.h"

#include "src/entities/data_structures/arma/ArmaFactories.h"
#include "src/entities/data_structures/eigen/EigenFactories.h"

using namespace spinner::linear_algebra;

std::vector<std::shared_ptr<AbstractDenseTransformAndDiagonalizeFactory>>
constructAllDenseTransformAndDiagonalizeFactories() {
    auto armaDouble = std::make_shared<ArmaDenseTransformAndDiagonalizeFactory>();
    armaDouble->setPrecision(DOUBLE);
    auto armaFloat = std::make_shared<ArmaDenseTransformAndDiagonalizeFactory>();
    armaFloat->setPrecision(SINGLE);
    auto eigenDouble = std::make_shared<EigenDenseTransformAndDiagonalizeFactory>();
    eigenDouble->setPrecision(DOUBLE);
    auto eigenSingle = std::make_shared<EigenDenseTransformAndDiagonalizeFactory>();
    eigenSingle->setPrecision(SINGLE);
    return {armaDouble, armaFloat, eigenDouble, eigenSingle};
}
