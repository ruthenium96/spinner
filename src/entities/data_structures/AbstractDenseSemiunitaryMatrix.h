#ifndef SPINNER_ABSTRACTDENSESEMIUNITARYMATRIX_H
#define SPINNER_ABSTRACTDENSESEMIUNITARYMATRIX_H

#include <memory>

#include "AbstractDenseVector.h"
#include "AbstractDenseSemiunitaryTransformer.h"

namespace quantum::linear_algebra {
class AbstractDiagonalizableMatrix;
class AbstractDenseSemiunitaryMatrix {
  public:
    virtual uint32_t size_rows() const = 0;
    virtual uint32_t size_cols() const = 0;

    virtual double at(uint32_t i, uint32_t j) const = 0;
    virtual void add_to_position(double value, uint32_t i, uint32_t j) = 0;

    virtual void print(std::ostream& os) const = 0;

    virtual std::unique_ptr<AbstractDiagonalizableMatrix> unitaryTransform(
        const std::unique_ptr<AbstractDiagonalizableMatrix>& matrix_to_transform) const = 0;

    virtual const std::unique_ptr<AbstractDenseSemiunitaryTransformer>& getUnitaryTransformer() const = 0;

    virtual void normalize() = 0;

    virtual ~AbstractDenseSemiunitaryMatrix() = default;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_ABSTRACTDENSESEMIUNITARYMATRIX_H
