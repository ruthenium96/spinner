#ifndef SPINNER_ABSTRACTSPARSESEMIUNITARYMATRIX_H
#define SPINNER_ABSTRACTSPARSESEMIUNITARYMATRIX_H

#include <memory>

#include "AbstractDiagonalizableMatrix.h"
#include "AbstractSymmetricMatrix.h"

namespace spinner::linalg_structures {

/**
 * @class AbstractSparseSemiunitaryMatrix
 * @brief Abstract interface for a sparse matrix whose columns are orthonormal.
 *
 * Provides an iterator over the non‑zero entries of a given column,
 * allowing efficient sparse operations. The matrix is semi‑unitary:
 * its columns are orthonormal, @f$U^* U = I@f$.
 * There are no internal checks for matrix orthonormality, 
 * since this is computationally expensive.
 */
class AbstractSparseSemiunitaryMatrix {
  public:
    /**
     * @struct Iterator
     * @brief Iterator over the non‑zero elements of a specific column.
     *
     * Provides sequential access to index‑value pairs. The iterator is
     * designed for manual iteration (hasNext()/getNext()) because
     * range‑based loops would require a different interface.
     */
    struct Iterator {
        /**
         * @struct IndexValueItem
         * @brief A single non‑zero element: (row index, value).
         */
        struct IndexValueItem {
            uint32_t index; ///< Row index (0‑based) of the element.
            double value; ///< Numerical value.
        };
        /**
         * @brief Checks if there are more elements to iterate.
         * @return true if another element exists, false otherwise.
         */
        virtual bool hasNext() const = 0;
        /**
         * @brief Returns the next non‑zero element and advances the iterator.
         * @return IndexValueItem containing the row index and value.
         * @pre hasNext() == true.
         */
        virtual IndexValueItem getNext() = 0;
        /**
         * @brief Returns the total number of non‑zero elements in the column.
         * @return Column non‑zero count.
         */
        virtual size_t size() const = 0;
        /// Virtual destructor.
        virtual ~Iterator() = default;
    };

    /**
     * @brief Creates a new iterator over the non‑zero elements of a given column.
     * @param col Column index (0‑based).
     * @return std::unique_ptr<Iterator> to iterate over that column.
     */
    virtual std::unique_ptr<Iterator> GetNewIterator(size_t col) const = 0;

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
     * @brief Checks if the entire matrix has no stored entries.
     * @return true if the matrix is completely empty, false otherwise.
     */
    virtual bool empty() const = 0;
    /**
     * @brief Checks if a specific column has no stored entries.
     * @param col Column index.
     * @return true if the column is empty, false otherwise.
     */
    virtual bool vempty(uint32_t col) const = 0;
    /**
     * @brief Removes all entries from the matrix, making it empty.
     */
    virtual void clear() = 0;

    /**
     * @brief Removes entries that are explicitly stored but numerically zero.
     *
     * Useful for cleaning up after arithmetic operations that may introduce zeros.
     */
    virtual void eraseExplicitZeros() = 0;

    /**
     * @brief Checks whether the element at (col, row) is implicit zero (i.e., not stored).
     * @param col Column index.
     * @param row Row index.
     * @return true if element is implicit zero, false if stored value.
     */
    virtual bool is_zero(uint32_t col, uint32_t row) const = 0;
    /**
     * @brief Moves the entire column `col` from another sparse semi‑unitary matrix.
     *
     * After the operation, the source column becomes empty (ownership transferred).
     *
     * @param col            Target column index in this matrix.
     * @param subspace_from  Source matrix; its column col is moved into this matrix.
     */
    virtual void move_vector_from(
        uint32_t col,
        std::unique_ptr<AbstractSparseSemiunitaryMatrix>& subspace_from) = 0;

    /**
     * @brief Adds a value to the element at (col, row).
     *
     * Could break orthonormality of the matrix.
     *
     * @param value Value to add.
     * @param col   Column index.
     * @param row   Row index.
     */
    virtual void add_to_position(double value, uint32_t col, uint32_t row) = 0;
    /**
     * @brief Returns the value at position (col, row) (zero if not stored).
     * @param col Column index.
     * @param row Row index.
     * @return Stored value or 0.0.
     */
    virtual double at(uint32_t col, uint32_t row) const = 0;

    /**
     * @brief Normalizes each column to unit length (Euclidean norm).
     */
    virtual void normalize() = 0;
    /**
     * @brief Applies the semi‑unitary transformation to a symmetric matrix,
     *        adding the result to another diagonalizable matrix.
     *
     * Computes `U^* * symmetricMatrixToTransform * U` and adds it to
     * `symmetricMatrixToAdd`.
     *
     * @param symmetricMatrixToTransform The symmetric matrix to be transformed.
     * @param symmetricMatrixToAdd       Destination matrix (updated in‑place).
     */
    virtual void unitaryTransform(
        const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrixToTransform,
        std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrixToAdd) const = 0;

    /**
     * @brief Prints a human‑readable representation of the matrix.
     * @param os Output stream.
     */
    virtual void print(std::ostream& os) const = 0;

    /// Virtual destructor.
    virtual ~AbstractSparseSemiunitaryMatrix() = default;
};
} // namespace spinner::linalg_structures
#endif  //SPINNER_ABSTRACTSPARSESEMIUNITARYMATRIX_H
