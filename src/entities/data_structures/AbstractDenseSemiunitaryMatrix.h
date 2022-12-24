#ifndef SPINNER_ABSTRACTDENSESEMIUNITARYMATRIX_H
#define SPINNER_ABSTRACTDENSESEMIUNITARYMATRIX_H

#include <memory>

#include "AbstractDenseVector.h"

namespace quantum::linear_algebra {
class AbstractSymmetricMatrix;
class AbstractDenseSemiunitaryMatrix {
  public:
    virtual uint32_t size_rows() const = 0;
    virtual uint32_t size_cols() const = 0;

    virtual double at(uint32_t i, uint32_t j) const = 0;

    virtual void print(std::ostream& os) const = 0;

    // TODO: I guess, we should do it faster. If B = U * A * U^T,
    //  we need only B_{ii} = \sum_{k} \sum_{l} U_{ik}*U_{il}*A_{kl}
    //  We can multiply A * U^T directly for no worse than O(N^3),
    //  then calculate all B_{ii} for O(N^2).
    virtual std::unique_ptr<AbstractDenseVector> unitaryTransformAndReturnMainDiagonal(
        const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrix) const = 0;

    virtual ~AbstractDenseSemiunitaryMatrix() = default;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_ABSTRACTDENSESEMIUNITARYMATRIX_H