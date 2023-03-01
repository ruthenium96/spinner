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

std::unique_ptr<AbstractDenseVector> EigenDenseVectorFactory::createVector() {
    auto vector = std::make_unique<EigenDenseVector>();
    return vector;
}

}  // namespace quantum::linear_algebra