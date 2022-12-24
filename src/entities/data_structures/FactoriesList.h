#ifndef SPINNER_ABSTRACTFACTORY__H
#define SPINNER_ABSTRACTFACTORY__H

#include "AbstractFactories.h"

// TODO: refactor the whole project and pass FactoriesList only to some controlling class.
//  in case of space/optimization it can be done by means of replacing move_vector_from to
//  move_vector_to function. In other cases it can be done through refactoring of constructors.

namespace quantum::linear_algebra {
class FactoriesList {
  public:
    explicit FactoriesList(
        std::shared_ptr<AbstractSymmetricMatrixFactory> symmetricMatrixFactory =
            AbstractSymmetricMatrixFactory::defaultFactory(),
        std::shared_ptr<AbstractSparseSemiunitaryFactory> sparseMatrix =
            AbstractSparseSemiunitaryFactory::defaultSparseFactory(),
        std::shared_ptr<AbstractDenseVectorFactory> denseVectorMatrix =
            AbstractDenseVectorFactory::defaultFactory());

    std::unique_ptr<AbstractSymmetricMatrix> createSymmetricMatrix(uint32_t size) const;
    std::unique_ptr<AbstractDenseVector> createVector() const;
    std::unique_ptr<AbstractSparseSemiunitaryMatrix>
    createSparseSemiunitaryMatrix(uint32_t cols, uint32_t rows) const;
    std::unique_ptr<AbstractSymmetricMatrix> createSparseSymmetricMatrix(uint32_t size) const;

  private:
    std::shared_ptr<AbstractSymmetricMatrixFactory> symmetricMatrixFactory_;
    std::shared_ptr<AbstractSparseSemiunitaryFactory> sparseSemiunitaryMatrixFactory_;
    std::shared_ptr<AbstractDenseVectorFactory> denseVectorFactory_;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_ABSTRACTFACTORY__H
