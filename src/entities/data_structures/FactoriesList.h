#ifndef SPINNER_FACTORIESLIST_H
#define SPINNER_FACTORIESLIST_H

#include "AbstractFactories.h"

// TODO: refactor the whole project and pass FactoriesList only to some controlling class.
//  in case of space/optimization it can be done by means of replacing move_vector_from to
//  move_vector_to function. In other cases it can be done through refactoring of constructors.

namespace quantum::linear_algebra {
class FactoriesList {
  public:
    explicit FactoriesList(
        std::shared_ptr<AbstractDenseTransformAndDiagonalizeFactory> symmetricMatrixFactory =
            AbstractDenseTransformAndDiagonalizeFactory::defaultFactory(),
        std::shared_ptr<AbstractSparseTransformFactory> sparseMatrix =
            AbstractSparseTransformFactory::defaultSparseFactory());

    std::unique_ptr<AbstractDiagonalizableMatrix>
    createDenseDiagonalizableMatrix(uint32_t size) const;
    std::unique_ptr<AbstractDiagonalizableMatrix>
    createSparseDiagonalizableMatrix(uint32_t size) const;
    std::unique_ptr<AbstractDenseSemiunitaryMatrix>
    createDenseSemiunitaryMatrix(uint32_t cols, uint32_t rows) const;
    std::unique_ptr<AbstractDenseVector> createRandomUnitVector(uint32_t size) const;
    std::unique_ptr<AbstractDenseVector> createVector() const;
    std::unique_ptr<AbstractSparseSemiunitaryMatrix>
    createSparseSemiunitaryMatrix(uint32_t cols, uint32_t rows) const;
    std::unique_ptr<AbstractSymmetricMatrix> createSparseSymmetricMatrix(uint32_t size) const;

  private:
    std::shared_ptr<AbstractDenseTransformAndDiagonalizeFactory> denseFactory_;
    std::shared_ptr<AbstractSparseTransformFactory> sparseFactory_;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_FACTORIESLIST_H
