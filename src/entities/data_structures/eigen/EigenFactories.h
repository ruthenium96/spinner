#ifndef SPINNER_EIGENFACTORIES_H
#define SPINNER_EIGENFACTORIES_H

#include "src/entities/data_structures/AbstractFactories.h"

namespace quantum::linear_algebra {

class EigenSymmetricMatrixFactory: public AbstractSymmetricMatrixFactory {
  public:
    std::unique_ptr<AbstractSymmetricMatrix>
    createSymmetricMatrix(uint32_t matrix_in_space_basis_size_i) override;
    std::unique_ptr<AbstractSymmetricMatrix> createSparseSymmetricMatrix(uint32_t size) override;
};

class EigenDenseVectorFactory: public AbstractDenseVectorFactory {
  public:
    std::unique_ptr<AbstractDenseVector> createVector() override;
};

}  // namespace quantum::linear_algebra
#endif  //SPINNER_EIGENFACTORIES_H
