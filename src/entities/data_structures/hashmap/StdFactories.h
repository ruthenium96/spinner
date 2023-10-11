#ifndef SPINNER_STDFACTORIES_H
#define SPINNER_STDFACTORIES_H

#include "src/entities/data_structures/AbstractFactories.h"

namespace quantum::linear_algebra {
class StdSparseTransformFactory: public AbstractSparseTransformFactory {
  public:
    std::unique_ptr<AbstractSparseSemiunitaryMatrix>
    createSparseSemiunitaryMatrix(uint32_t cols, uint32_t rows) override;
    std::unique_ptr<AbstractSymmetricMatrix> createSparseSymmetricMatrix(uint32_t size) override;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_STDFACTORIES_H
