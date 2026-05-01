#ifndef SPINNER_ABSTRACTDENSESEMIUNITARYMATRIX_H
#define SPINNER_ABSTRACTDENSESEMIUNITARYMATRIX_H

#include <cstdint>
#include <memory>

#include "AbstractDenseSemiunitaryTransformer.h"

namespace spinner::linalg_structures {
/**
 * @class AbstractDenseSemiunitaryMatrix
 * @brief Abstract interface for a dense rectangle matrix that is semi‑unitary, @f$U^* U = I@f$.
 *
 * Such matrices typically represent collections of orthonormal vectors (columns).
 * They can be used to store eigenvectors or their approximation from Krylov's procedure.
 *
 */
class AbstractDenseSemiunitaryMatrix {
  public:
    /**
     * @brief Returns the number of rows of the matrix.
     * @return Row count.
     */
    virtual uint32_t size_rows() const = 0;
    /**
     * @brief Returns the number of columns of the matrix.
     * @return Column count.
     */
    virtual uint32_t size_cols() const = 0;

    /**
     * @brief Returns the element at position (i, j).
     * @param i Row index.
     * @param j Column index.
     * @return Matrix entry.
     */
    virtual double at(uint32_t i, uint32_t j) const = 0;

    /**
     * @brief Prints a human‑readable representation of the matrix to an output stream.
     * @param os Output stream.
     */
    virtual void print(std::ostream& os) const = 0;

    /**
     * @brief Returns a pointer to the transformer object associated with this matrix.
     *
     * The transformer is responsible for applying the semi‑unitary transformation
     * to other matrices (e.g., computing the main diagonal of  @f$U^* * A * U@f$).
     *
     * @return Const reference to the unique_ptr holding the transformer.
     */
    virtual const std::unique_ptr<AbstractDenseSemiunitaryTransformer>& getUnitaryTransformer() const = 0;

    /// Virtual destructor.
    virtual ~AbstractDenseSemiunitaryMatrix() = default;
};
} // namespace spinner::linalg_structures
#endif  //SPINNER_ABSTRACTDENSESEMIUNITARYMATRIX_H
