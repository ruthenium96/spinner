#ifndef SPINNER_ABSTRACTDENSESEMIUNITARYTRANSFORMER_H
#define SPINNER_ABSTRACTDENSESEMIUNITARYTRANSFORMER_H

#include <functional>
#include <memory>

namespace spinner::linalg_structures {
class AbstractDenseVector;
class AbstractDiagonalizableMatrix;

/**
 * @class AbstractDenseSemiunitaryTransformer
 * @brief Performs unitary transformations of the form @f$U^* * A * U@f$.
 *
 * This class is constructed in AbstractDenseSemiunitaryMatrix and
 * is responsible for applying the transformation to a diagonalizable matrix,
 * returning the main diagonal of the result @f$U^* * A * U@f$ 
 * or @f$<n(r)|A|r><r|n(r)>@f$ in the case of Krylov's approximation.
 */
class AbstractDenseSemiunitaryTransformer {
public:
    /**
     * @brief Computes the main diagonal of @f$U^* * A * U@f$.
     *
     * @param matrix A reference wrapper around a unique_ptr to the digonalizable matrix @f$A@f$ to be transformed.
     * @return Dense vector containing the diagonal entries of the transformed matrix.
     */
    virtual std::unique_ptr<AbstractDenseVector> calculateUnitaryTransformationOfMatrix(
        std::reference_wrapper<const std::unique_ptr<AbstractDiagonalizableMatrix>> matrix) const = 0;
    /// Virtual destructor.
    virtual ~AbstractDenseSemiunitaryTransformer() = default;
};
} // namespace spinner::linalg_structures

#endif  //SPINNER_ABSTRACTDENSESEMIUNITARYTRANSFORMER_H