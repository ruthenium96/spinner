#ifndef SPINNER_ABSTRACTDENSESEMIUNITARYTRANSFORMER_H
#define SPINNER_ABSTRACTDENSESEMIUNITARYTRANSFORMER_H

#include <functional>
#include <memory>

namespace quantum::linear_algebra {
class AbstractDenseSemiunitaryMatrix;
class AbstractDenseVector;
class AbstractDiagonalizableMatrix;

// This class calculates UAU* and UABU* main diagonals:
class AbstractDenseSemiunitaryTransformer {
public:
    virtual std::unique_ptr<AbstractDenseVector> calculateUnitaryTransformationOfMatrix(
        std::reference_wrapper<const std::unique_ptr<AbstractDiagonalizableMatrix>> matrix) const = 0;
};
} // namespace quantum::linear_algebra

#endif  //SPINNER_ABSTRACTDENSESEMIUNITARYTRANSFORMER_H