#ifndef SPINNER_EMHASHLOGIC_H
#define SPINNER_EMHASHLOGIC_H

#include "EmhashSparseSymmetricMatrix.h"
#include "StdSparseSemiunitaryMatrix.h"

namespace quantum::linear_algebra {

class EmhashLogic {
  public:
    void unitaryTransform(
        const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrixToTransform,
        std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrixToAdd,
        const AbstractSparseSemiunitaryMatrix& unitaryMatrix) const;
};

}  // namespace quantum::linear_algebra

#endif  //SPINNER_EMHASHLOGIC_H
