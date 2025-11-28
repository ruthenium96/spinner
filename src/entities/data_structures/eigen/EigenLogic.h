#ifndef SPINNER_EIGENLOGIC_H
#define SPINNER_EIGENLOGIC_H

#include <memory>

#include "src/entities/data_structures/AbstractDenseSemiunitaryMatrix.h"
#include "src/entities/data_structures/AbstractDenseVector.h"
#include "src/entities/data_structures/AbstractDiagonalizableMatrix.h"

namespace quantum::linear_algebra {
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

}  // namespace quantum::linear_algebra

#endif  //SPINNER_EIGENLOGIC_H
