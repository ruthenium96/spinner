#ifndef SPINNER_ARMALOGIC_H
#define SPINNER_ARMALOGIC_H

#include <memory>

#include "src/linalg_structures/AbstractDenseSemiunitaryMatrix.h"
#include "src/linalg_structures/AbstractDenseVector.h"
#include "src/linalg_structures/AbstractDiagonalizableMatrix.h"

namespace spinner::linalg_structures::armadillo {
template <typename T>
class ArmaLogic {
  public:
    std::unique_ptr<AbstractDenseVector>
    diagonalizeValues(const AbstractDiagonalizableMatrix& diagonalizableMatrix) const;

    EigenCouple
    diagonalizeValuesVectors(const AbstractDiagonalizableMatrix& diagonalizableMatrix) const;

    KrylovCouple krylovDiagonalizeValues(
        const AbstractDiagonalizableMatrix& diagonalizableMatrix,
        const AbstractDenseVector& seed_vector,
        size_t krylov_subspace_size) const;

    KrylovTriple krylovDiagonalizeValuesVectors(
        const AbstractDiagonalizableMatrix& diagonalizableMatrix,
        const AbstractDenseVector& seed_vector,
        size_t krylov_subspace_size) const;
};
} // namespace spinner::linalg_structures::armadillo

#endif  //SPINNER_ARMALOGIC_H
