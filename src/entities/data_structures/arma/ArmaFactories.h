#ifndef SPINNER_ARMAFACTORIES_H
#define SPINNER_ARMAFACTORIES_H

#include "src/entities/data_structures/AbstractFactories.h"

namespace quantum::linear_algebra {
class ArmaDenseTransformAndDiagonalizeFactory: public AbstractDenseTransformAndDiagonalizeFactory {
  public:
    std::unique_ptr<AbstractDiagonalizableMatrix>
    createDenseDiagonalizableMatrix(uint32_t size) override;
    std::unique_ptr<AbstractDiagonalizableMatrix>
    createSparseDiagonalizableMatrix(uint32_t size) override;
    std::unique_ptr<AbstractDenseSemiunitaryMatrix>
    createDenseSemiunitaryMatrix(uint32_t cols, uint32_t rows) override;
    std::unique_ptr<AbstractDenseVector> createRandomUnitVector(uint32_t size) override;
    std::unique_ptr<AbstractDenseVector> createVector() override;
};

class ArmaSparseTransformFactory: public AbstractSparseTransformFactory {
  public:
    std::unique_ptr<AbstractSparseSemiunitaryMatrix>
    createSparseSemiunitaryMatrix(uint32_t rows, uint32_t cols) override;
    std::unique_ptr<AbstractSymmetricMatrix> createSparseSymmetricMatrix(uint32_t size) override;
};

}  // namespace quantum::linear_algebra
#endif  //SPINNER_ARMAFACTORIES_H
