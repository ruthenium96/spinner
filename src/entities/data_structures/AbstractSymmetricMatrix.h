#ifndef SPINNER_ABSTRACTSYMMETRICMATRIX_H
#define SPINNER_ABSTRACTSYMMETRICMATRIX_H

#include <cstdint>
#include <memory>

#include "AbstractDenseVector.h"

namespace quantum::linear_algebra {
class AbstractDenseSemiunitaryMatrix;

struct EigenCouple {
    std::unique_ptr<AbstractDenseVector> eigenvalues;
    std::unique_ptr<AbstractDenseSemiunitaryMatrix> eigenvectors;
};

class AbstractSymmetricMatrix {
  public:
    virtual void add_to_position(double value, uint32_t i, uint32_t j) = 0;
    virtual EigenCouple diagonalizeValuesVectors() const = 0;
    virtual std::unique_ptr<AbstractDenseVector> diagonalizeValues() const = 0;

    virtual uint32_t size() const = 0;
    virtual double at(uint32_t i, uint32_t j) const = 0;

    virtual void print(std::ostream& os) const = 0;

    virtual ~AbstractSymmetricMatrix() = default;
};
}  // namespace quantum::linear_algebra

#endif  //SPINNER_ABSTRACTSYMMETRICMATRIX_H
