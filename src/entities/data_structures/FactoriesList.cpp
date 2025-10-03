#include "FactoriesList.h"
#include <vector>

namespace quantum::linear_algebra {

std::unique_ptr<AbstractDiagonalizableMatrix>
FactoriesList::createDenseDiagonalizableMatrix(uint32_t size) const {
    return denseFactory_->createDenseDiagonalizableMatrix(size);
}

std::unique_ptr<AbstractDiagonalizableMatrix>
FactoriesList::createSparseDiagonalizableMatrix(uint32_t size) const {
    return denseFactory_->createSparseDiagonalizableMatrix(size);
}

std::unique_ptr<AbstractDenseSemiunitaryMatrix>
FactoriesList::createDenseSemiunitaryMatrix(uint32_t cols, uint32_t rows) const {
    return denseFactory_->createDenseSemiunitaryMatrix(cols, rows);
}

std::vector<std::unique_ptr<AbstractDenseVector>> 
FactoriesList::createRandomUnitVectors(uint32_t size_of_vector, uint32_t number_of_vectors) const {
    return denseFactory_->createRandomUnitVectors(size_of_vector, number_of_vectors);
}

std::unique_ptr<AbstractDenseVector> FactoriesList::createVector() const {
    return denseFactory_->createVector();
}

std::unique_ptr<AbstractSparseSemiunitaryMatrix>
FactoriesList::createSparseSemiunitaryMatrix(uint32_t cols, uint32_t rows) const {
    return sparseFactory_->createSparseSemiunitaryMatrix(cols, rows);
}

std::unique_ptr<AbstractSymmetricMatrix>
FactoriesList::createSparseSymmetricMatrix(uint32_t size) const {
    return sparseFactory_->createSparseSymmetricMatrix(size);
}

FactoriesList::FactoriesList(
    std::shared_ptr<AbstractDenseTransformAndDiagonalizeFactory> symmetricMatrixFactory,
    std::shared_ptr<AbstractSparseTransformFactory> sparseMatrix) {
    denseFactory_ = symmetricMatrixFactory;
    sparseFactory_ = sparseMatrix;
}
}  // namespace quantum::linear_algebra