#include "EigenFactories.h"

#include "EigenDenseSemiunitaryMatrix.h"
#include "EigenDenseSymmetricMatrix.h"
#include "EigenDenseVector.h"
#include "EigenSparseSymmetricMatrix.h"

namespace quantum::linear_algebra {

std::unique_ptr<AbstractSymmetricMatrix>
EigenSymmetricMatrixFactory::createDenseSymmetricMatrix(uint32_t size) {
    auto matrix = std::make_unique<EigenDenseSymmetricMatrix>();
    matrix->resize(size);
    return matrix;
}

std::unique_ptr<AbstractSymmetricMatrix>
EigenSymmetricMatrixFactory::createSparseSymmetricMatrix(uint32_t size) {
    auto matrix = std::make_unique<EigenSparseSymmetricMatrix>();
    matrix->resize(size);
    return matrix;
}

std::unique_ptr<AbstractDenseSemiunitaryMatrix>
EigenSymmetricMatrixFactory::createDenseSemiunitaryMatrix(uint32_t cols, uint32_t rows) {
    auto matrix = std::make_unique<EigenDenseSemiunitaryMatrix>();
    matrix->resize(rows, cols);
    return matrix;
}

std::unique_ptr<AbstractDenseVector> EigenDenseVectorFactory::createVector() {
    auto vector = std::make_unique<EigenDenseVector>();
    return vector;
}

}  // namespace quantum::linear_algebra