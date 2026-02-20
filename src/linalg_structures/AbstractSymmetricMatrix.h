#ifndef SPINNER_ABSTRACTSYMMETRICMATRIX_H
#define SPINNER_ABSTRACTSYMMETRICMATRIX_H

#include <cstdint>
#include <memory>

namespace spinner::linalg_structures {
/**
 * @class AbstractSymmetricMatrix
 * @brief Abstract base class for a symmetric matrix (dense or sparse).
 *
 * Provides a common interface for symmetric matrices, supporting element
 * access and modification. The matrix is always square; its size is
 * obtained via `size()`.
 */
class AbstractSymmetricMatrix {
  public:
    /**
     * @brief Adds a value to the element at position (i, j).
     *
     * Since the matrix is symmetric, the element at (j, i) is also affected
     * to maintain symmetry. Implementations must ensure symmetry.
     *
     * @param value Value to add.
     * @param i     Row/column index (0-based).
     * @param j     Column/row index (0-based).
     */
    virtual void add_to_position(double value, uint32_t i, uint32_t j) = 0;

    /**
     * @brief Returns the number of rows/columns of the square matrix.
     * @return Matrix dimension.
     */
    virtual uint32_t size() const = 0;
    /**
     * @brief Returns the element at position (i, j).
     *
     * Because of symmetry, `at(i,j)` should equal `at(j,i)`.
     *
     * @param i Row/column index.
     * @param j Column/row index.
     * @return Value at (i,j).
     */
    virtual double at(uint32_t i, uint32_t j) const = 0;

    /**
     * @brief Prints a human‑readable representation of the matrix to an output stream.
     * @param os Output stream.
     */
    virtual void print(std::ostream& os) const = 0;

    /// Virtual destructor for proper cleanup of derived objects.
    virtual ~AbstractSymmetricMatrix() = default;
};
} // namespace spinner::linalg_structures

#endif  //SPINNER_ABSTRACTSYMMETRICMATRIX_H
