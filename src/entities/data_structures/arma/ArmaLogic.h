#ifndef SPINNER_ARMALOGIC_H
#define SPINNER_ARMALOGIC_H

#include <memory>

#include "src/entities/data_structures/AbstractDenseSemiunitaryMatrix.h"
#include "src/entities/data_structures/AbstractDenseVector.h"
#include "src/entities/data_structures/AbstractSparseSemiunitaryMatrix.h"
#include "src/entities/data_structures/AbstractSymmetricMatrix.h"

namespace quantum::linear_algebra {
class ArmaLogic {
  public:
    std::unique_ptr<AbstractDenseVector>
    diagonalizeValues(const AbstractSymmetricMatrix& symmetricMatrix) const;
    std::unique_ptr<AbstractDenseVector>
    diagonalizeValues(const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrix) const;

    EigenCouple diagonalizeValuesVectors(const AbstractSymmetricMatrix& symmetricMatrix) const;
    EigenCouple
    diagonalizeValuesVectors(const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrix) const;

    std::unique_ptr<AbstractDenseVector> unitaryTransformAndReturnMainDiagonal(
        const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrix,
        const AbstractDenseSemiunitaryMatrix& denseSemiunitaryMatrix) const;
    std::unique_ptr<AbstractDenseVector> unitaryTransformAndReturnMainDiagonal(
        const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrix,
        const std::unique_ptr<AbstractDenseSemiunitaryMatrix>& denseSemiunitaryMatrix) const;
};
}  // namespace quantum::linear_algebra

#endif  //SPINNER_ARMALOGIC_H
