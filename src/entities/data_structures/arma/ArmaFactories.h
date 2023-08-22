#ifndef SPINNER_ARMAFACTORIES_H
#define SPINNER_ARMAFACTORIES_H

#include "src/entities/data_structures/AbstractFactories.h"

namespace quantum::linear_algebra {
class ArmaSymmetricMatrixFactory: public AbstractSymmetricMatrixFactory {
  public:
    std::unique_ptr<AbstractSymmetricMatrix> createDenseSymmetricMatrix(uint32_t size) override;
    std::unique_ptr<AbstractSymmetricMatrix> createSparseSymmetricMatrix(uint32_t size) override;
    std::unique_ptr<AbstractDenseSemiunitaryMatrix>
    createDenseSemiunitaryMatrix(uint32_t cols, uint32_t rows) override;
};

class ArmaDenseVectorFactory: public AbstractDenseVectorFactory {
  public:
    std::unique_ptr<AbstractDenseVector> createVector() override;
};

class ArmaSparseSemiunitaryFactory: public AbstractSparseSemiunitaryFactory {
    std::unique_ptr<AbstractSparseSemiunitaryMatrix>
    createSparseSemiunitaryMatrix(uint32_t rows, uint32_t cols) override;
};

}  // namespace quantum::linear_algebra
#endif  //SPINNER_ARMAFACTORIES_H
