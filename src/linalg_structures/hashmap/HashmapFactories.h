#ifndef SPINNER_HASHMAPFACTORIES_H
#define SPINNER_HASHMAPFACTORIES_H

#include "src/linalg_structures/AbstractFactories.h"

namespace spinner::linalg_structures {
class EmhashSparseTransformFactory: public AbstractSparseTransformFactory {
  public:
    std::unique_ptr<AbstractSparseSemiunitaryMatrix>
    createSparseSemiunitaryMatrix(uint32_t cols, uint32_t rows) override;
    std::unique_ptr<AbstractSymmetricMatrix> createSparseSymmetricMatrix(uint32_t size) override;
};
} // namespace spinner::linalg_structures
#endif  //SPINNER_HASHMAPFACTORIES_H
