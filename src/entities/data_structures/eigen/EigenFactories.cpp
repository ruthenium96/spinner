#include "EigenFactories.h"

#include "EigenDenseSemiunitaryMatrix.h"
#include "EigenDenseSymmetricMatrix.h"
#include "EigenDenseVector.h"

namespace quantum::linear_algebra {

std::unique_ptr<AbstractSymmetricMatrix>
EigenSymmetricMatrixFactory::createDenseSymmetricMatrix(uint32_t size) {
    auto matrix = std::make_unique<EigenDenseSymmetricMatrix>();
    matrix->resize(size);
    return matrix;
}

std::unique_ptr<AbstractSymmetricMatrix>
EigenSymmetricMatrixFactory::createSparseSymmetricMatrix(uint32_t size) {
    // TODO: implement it ASAP
    throw std::logic_error("Not implemented");
}

std::unique_ptr<AbstractDenseVector> EigenDenseVectorFactory::createVector() {
    auto vector = std::make_unique<EigenDenseVector>();
    return vector;
}

}  // namespace quantum::linear_algebra