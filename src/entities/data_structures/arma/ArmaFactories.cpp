#include "ArmaFactories.h"

#include "ArmaDenseDiagonalizableMatrix.h"
#include "ArmaDenseSemiunitaryMatrix.h"
#include "ArmaDenseVector.h"
#include "ArmaSparseDiagonalizableMatrix.h"
#include "ArmaSparseSemiunitaryMatrix.h"

namespace quantum::linear_algebra {

std::unique_ptr<AbstractDiagonalizableMatrix>
ArmaDenseTransformAndDiagonalizeFactory::createDenseDiagonalizableMatrix(uint32_t size) {
    auto matrix = std::make_unique<ArmaDenseDiagonalizableMatrix>();
    matrix->resize(size);
    return matrix;
}

std::unique_ptr<AbstractDiagonalizableMatrix>
ArmaDenseTransformAndDiagonalizeFactory::createSparseDiagonalizableMatrix(uint32_t size) {
    auto matrix = std::make_unique<ArmaSparseDiagonalizableMatrix>();
    matrix->resize(size);
    return matrix;
}

std::unique_ptr<AbstractDenseSemiunitaryMatrix>
ArmaDenseTransformAndDiagonalizeFactory::createDenseSemiunitaryMatrix(
    uint32_t cols,
    uint32_t rows) {
    auto matrix = std::make_unique<ArmaDenseSemiunitaryMatrix>();
    matrix->resize(rows, cols);
    return matrix;
}

std::unique_ptr<AbstractDenseVector> ArmaDenseTransformAndDiagonalizeFactory::createVector() {
    auto vector = std::make_unique<ArmaDenseVector>();
    return vector;
}

std::unique_ptr<AbstractSparseSemiunitaryMatrix>
ArmaSparseTransformFactory::createSparseSemiunitaryMatrix(uint32_t cols, uint32_t rows) {
    auto matrix = std::make_unique<ArmaSparseSemiunitaryMatrix>();
    matrix->resize(cols, rows);
    return matrix;
}

std::unique_ptr<AbstractSymmetricMatrix>
ArmaSparseTransformFactory::createSparseSymmetricMatrix(uint32_t size) {
    auto matrix = std::make_unique<ArmaSparseSymmetricMatrix>();
    matrix->resize(size);
    return matrix;
}
}  // namespace quantum::linear_algebra