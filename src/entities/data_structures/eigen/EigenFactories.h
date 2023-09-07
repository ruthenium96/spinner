#ifndef SPINNER_EIGENFACTORIES_H
#define SPINNER_EIGENFACTORIES_H

#include "src/entities/data_structures/AbstractFactories.h"

namespace quantum::linear_algebra {

class EigenDenseTransformAndDiagonalizeFactory: public AbstractDenseTransformAndDiagonalizeFactory {
  public:
    std::unique_ptr<AbstractDiagonalizableMatrix>
    createDenseDiagonalizableMatrix(uint32_t matrix_in_space_basis_size_i) override;
    std::unique_ptr<AbstractDiagonalizableMatrix>
    createSparseDiagonalizableMatrix(uint32_t size) override;
    std::unique_ptr<AbstractDenseSemiunitaryMatrix>
    createDenseSemiunitaryMatrix(uint32_t cols, uint32_t rows) override;
    std::unique_ptr<AbstractDenseVector> createVector() override;
};

}  // namespace quantum::linear_algebra
#endif  //SPINNER_EIGENFACTORIES_H
