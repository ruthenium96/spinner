#include "ArmaFactories.h"

#include "ArmaDenseSemiunitaryMatrix.h"
#include "ArmaDenseSymmetricMatrix.h"
#include "ArmaDenseVector.h"
#include "ArmaSparseSemiunitaryMatrix.h"
#include "ArmaSparseSymmetricMatrix.h"

namespace quantum::linear_algebra {

std::unique_ptr<AbstractSymmetricMatrix>
ArmaSymmetricMatrixFactory::createSymmetricMatrix(uint32_t size) {
    auto matrix = std::make_unique<ArmaDenseSymmetricMatrix>();
    matrix->resize(size);
    return matrix;
}

std::unique_ptr<AbstractSymmetricMatrix>
ArmaSymmetricMatrixFactory::createSparseSymmetricMatrix(uint32_t size) {
    auto matrix = std::make_unique<ArmaSparseSymmetricMatrix>();
    matrix->resize(size);
    return matrix;
}

std::unique_ptr<AbstractDenseVector> ArmaDenseVectorFactory::createVector() {
    auto vector = std::make_unique<ArmaDenseVector>();
    return vector;
}
std::unique_ptr<AbstractSparseSemiunitaryMatrix>
ArmaSparseSemiunitaryFactory::createSparseSemiunitaryMatrix(uint32_t cols, uint32_t rows) {
    auto matrix = std::make_unique<ArmaSparseSemiunitaryMatrix>();
    matrix->resize(cols, rows);
    return matrix;
}
}  // namespace quantum::linear_algebra