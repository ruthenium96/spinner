#ifndef SPINNER_ARMALOGIC_H
#define SPINNER_ARMALOGIC_H

#include <memory>

#include "src/entities/data_structures/AbstractDenseSemiunitaryMatrix.h"
#include "src/entities/data_structures/AbstractDenseVector.h"
#include "src/entities/data_structures/AbstractSparseSemiunitaryMatrix.h"
#include "src/entities/data_structures/AbstractSymmetricMatrix.h"
#include "src/entities/data_structures/arma/ArmaDenseSemiunitaryMatrix.h"
#include "src/entities/data_structures/arma/ArmaKrylovDenseSemiunitaryMatrix.h"

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
    
    std::unique_ptr<AbstractDiagonalizableMatrix> unitaryTransform(
        const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
        const ArmaDenseSemiunitaryMatrix<T>& denseSemiunitaryMatrix) const;
};
}  // namespace quantum::linear_algebra

#endif  //SPINNER_ARMALOGIC_H
