#ifndef SPINNER_ARMALOGIC_H
#define SPINNER_ARMALOGIC_H

#include <memory>

#include "src/entities/data_structures/AbstractDenseSemiunitaryMatrix.h"
#include "src/entities/data_structures/AbstractDenseVector.h"
#include "src/entities/data_structures/AbstractDiagonalizableMatrix.h"

namespace quantum::linear_algebra {
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
}  // namespace quantum::linear_algebra

#endif  //SPINNER_ARMALOGIC_H
