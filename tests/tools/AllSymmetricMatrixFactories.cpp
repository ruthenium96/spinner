#include "AllSymmetricMatrixFactories.h"

#include "src/linalg_structures/arma/ArmaFactories.h"
#include "src/linalg_structures/eigen/EigenFactories.h"

using namespace spinner::linalg_structures;

std::vector<std::shared_ptr<AbstractDenseTransformAndDiagonalizeFactory>>
constructAllDenseTransformAndDiagonalizeFactories() {
    auto armaDouble = std::make_shared<armadillo::ArmaDenseTransformAndDiagonalizeFactory>();
    armaDouble->setPrecision(DOUBLE);
    auto armaFloat = std::make_shared<armadillo::ArmaDenseTransformAndDiagonalizeFactory>();
    armaFloat->setPrecision(SINGLE);
    auto eigenDouble = std::make_shared<eigen::EigenDenseTransformAndDiagonalizeFactory>();
    eigenDouble->setPrecision(DOUBLE);
    auto eigenSingle = std::make_shared<eigen::EigenDenseTransformAndDiagonalizeFactory>();
    eigenSingle->setPrecision(SINGLE);
    return {armaDouble, armaFloat, eigenDouble, eigenSingle};
}
