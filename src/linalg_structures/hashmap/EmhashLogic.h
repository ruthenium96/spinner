#ifndef SPINNER_EMHASHLOGIC_H
#define SPINNER_EMHASHLOGIC_H

#include "EmhashSparseSemiunitaryMatrix.h"
#include "EmhashSparseSymmetricMatrix.h"

namespace spinner::linalg_structures::hashmap {

class EmhashLogic {
  public:
    void unitaryTransform(
        const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrixToTransform,
        std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrixToAdd,
        const EmhashSparseSemiunitaryMatrix& unitaryMatrix) const;
};

} // namespace spinner::linalg_structures::hashmap

#endif  //SPINNER_EMHASHLOGIC_H
