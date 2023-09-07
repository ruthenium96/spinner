#include "EigenFactories.h"

#include "EigenDenseDiagonalizableMatrix.h"
#include "EigenDenseSemiunitaryMatrix.h"
#include "EigenDenseVector.h"
#include "EigenSparseDiagonalizableMatrix.h"

namespace quantum::linear_algebra {

std::unique_ptr<AbstractDiagonalizableMatrix>
EigenDenseTransformAndDiagonalizeFactory::createDenseDiagonalizableMatrix(uint32_t size) {
    auto matrix = std::make_unique<EigenDenseDiagonalizableMatrix>();
    matrix->resize(size);
    return matrix;
}

std::unique_ptr<AbstractDiagonalizableMatrix>
EigenDenseTransformAndDiagonalizeFactory::createSparseDiagonalizableMatrix(uint32_t size) {
    auto matrix = std::make_unique<EigenSparseDiagonalizableMatrix>();
    matrix->resize(size);
    return matrix;
}

std::unique_ptr<AbstractDenseSemiunitaryMatrix>
EigenDenseTransformAndDiagonalizeFactory::createDenseSemiunitaryMatrix(
    uint32_t cols,
    uint32_t rows) {
    auto matrix = std::make_unique<EigenDenseSemiunitaryMatrix>();
    matrix->resize(rows, cols);
    return matrix;
}

std::unique_ptr<AbstractDenseVector> EigenDenseTransformAndDiagonalizeFactory::createVector() {
    auto vector = std::make_unique<EigenDenseVector>();
    return vector;
}

}  // namespace quantum::linear_algebra