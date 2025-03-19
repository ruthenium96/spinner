#ifndef SPINNER_ABSTRACTDIAGONALIZABLESYMMETRICMATRIX_H
#define SPINNER_ABSTRACTDIAGONALIZABLESYMMETRICMATRIX_H

#include "AbstractDenseVector.h"
#include "AbstractSymmetricMatrix.h"

namespace quantum::linear_algebra {
class AbstractDenseSemiunitaryMatrix;
struct EigenCouple {
    std::unique_ptr<AbstractDenseVector> eigenvalues;
    std::unique_ptr<AbstractDenseSemiunitaryMatrix> eigenvectors;
};

class AbstractDiagonalizableMatrix: public AbstractSymmetricMatrix {
  public:
    virtual EigenCouple diagonalizeValuesVectors() const = 0;
    virtual std::unique_ptr<AbstractDenseVector> diagonalizeValues() const = 0;
    virtual std::unique_ptr<AbstractDiagonalizableMatrix> multiply_by(double multiplier) const = 0;
    ~AbstractDiagonalizableMatrix() override = default;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_ABSTRACTDIAGONALIZABLESYMMETRICMATRIX_H
