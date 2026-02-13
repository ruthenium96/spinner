#ifndef SPINNER_EIGENLOGIC_H
#define SPINNER_EIGENLOGIC_H

#include <memory>

#include "src/linalg_structures/AbstractDenseSemiunitaryMatrix.h"
#include "src/linalg_structures/AbstractDenseVector.h"
#include "src/linalg_structures/AbstractDiagonalizableMatrix.h"

namespace spinner::linalg_structures {
template <typename T>
class EigenLogic {
  public:
    std::unique_ptr<AbstractDenseVector>
    diagonalizeValues(const AbstractDiagonalizableMatrix& symmetricMatrix) const;

    EigenCouple diagonalizeValuesVectors(const AbstractDiagonalizableMatrix& symmetricMatrix) const;

    KrylovCouple krylovDiagonalizeValues(
      const AbstractDiagonalizableMatrix& diagonalizableMatrix,
      const AbstractDenseVector& seed_vector,
      size_t krylov_subspace_size) const;
    KrylovTriple krylovDiagonalizeValuesVectors(
      const AbstractDiagonalizableMatrix& diagonalizableMatrix,
      const AbstractDenseVector& seed_vector,
      size_t krylov_subspace_size) const;
};

} // namespace spinner::linalg_structures

#endif  //SPINNER_EIGENLOGIC_H
