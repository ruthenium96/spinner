#ifndef SPINNER_ABSTRACTDIAGONALIZABLESYMMETRICMATRIX_H
#define SPINNER_ABSTRACTDIAGONALIZABLESYMMETRICMATRIX_H

#include "AbstractDenseVector.h"
#include "AbstractSymmetricMatrix.h"

namespace spinner::linalg_structures {
// Forward declaration (defined in AbstractDenseSemiunitaryMatrix.h)
class AbstractDenseSemiunitaryMatrix;
/**
 * @struct EigenCouple
 * @brief Holds the result of a full diagonalization (eigenvalues + eigenvectors).
 * Eigenvectors @f$U@f$, such that @f$AU=EU@f$.
 */
struct EigenCouple {
    std::unique_ptr<AbstractDenseVector> eigenvalues; ///< Eigenvalues as a dense vector.
    std::unique_ptr<AbstractDenseSemiunitaryMatrix> eigenvectors; ///< Eigenvectors as an unitary matrix.
};

/**
 * @struct KrylovCouple
 * @brief Result of a Krylov‑based diagonalization returning only eigenvalues and weights.
 */
struct KrylovCouple {
    std::unique_ptr<AbstractDenseVector> eigenvalues; ///< Approximate eigenvalues.
    std::unique_ptr<AbstractDenseVector> ftlm_weights_of_states; ///< Weights for the states.
};

/**
 * @struct KrylovTriple
 * @brief Result of a Krylov‑based diagonalization returning eigenvalues, eigenvectors, and weights.
 */
struct KrylovTriple {
  std::unique_ptr<AbstractDenseVector> eigenvalues; ///< Approximate eigenvalues.
  std::unique_ptr<AbstractDenseSemiunitaryMatrix> eigenvectors; ///< Approximate eigenvectors as an unitary matrix.
  std::unique_ptr<AbstractDenseVector> ftlm_weights_of_states; ///< Weights for the states.
};

/**
 * @class AbstractDiagonalizableMatrix
 * @brief Abstract interface for a symmetric matrix that can be diagonalized.
 *
 * Inherits from AbstractSymmetricMatrix and adds diagonalization capabilities:
 * full diagonalization (all eigenvalues/vectors) and Krylov subspace methods.
 */
class AbstractDiagonalizableMatrix: public AbstractSymmetricMatrix {
  public:
    /**
     * @brief Computes all eigenvalues and eigenvectors (full diagonalization).
     * @return EigenCouple containing eigenvalues and eigenvectors.
     */
    virtual EigenCouple diagonalizeValuesVectors() const = 0;
    /**
     * @brief Computes all eigenvalues only (no eigenvectors).
     * @return Dense vector of eigenvalues.
     */
    virtual std::unique_ptr<AbstractDenseVector> diagonalizeValues() const = 0;

    /**
     * @brief Performs Krylov subspace diagonalization to obtain eigenvalues and weights.
     *
     * @param seed_vector        Starting vector for the Krylov subspace.
     * @param krylov_subspace_size Dimension of the Krylov subspace.
     * @return KrylovCouple containing approximate eigenvalues and FTLM weights.
     */
    virtual KrylovCouple krylovDiagonalizeValues(
      const std::unique_ptr<AbstractDenseVector>& seed_vector,
      size_t krylov_subspace_size) const = 0;
    /**
     * @brief Performs Krylov subspace diagonalization to obtain eigenvalues, eigenvectors, and weights.
     *
     * @param seed_vector        Starting vector for the Krylov subspace.
     * @param krylov_subspace_size Dimension of the Krylov subspace.
     * @return KrylovTriple containing approximate eigenvalues, eigenvectors, and FTLM weights.
     */
    virtual KrylovTriple krylovDiagonalizeValuesVectors(
      const std::unique_ptr<AbstractDenseVector>& seed_vector,
      size_t krylov_subspace_size) const = 0;

    /**
     * @brief Returns a new matrix that is a scalar multiple of the current one.
     * @param multiplier Scalar factor.
     * @return std::unique_ptr<AbstractDiagonalizableMatrix> New matrix = multiplier * this.
     */
    virtual std::unique_ptr<AbstractDiagonalizableMatrix> multiply_by(double multiplier) const = 0;

    /// Virtual destructor.
    ~AbstractDiagonalizableMatrix() override = default;
};
} // namespace spinner::linalg_structures
#endif  //SPINNER_ABSTRACTDIAGONALIZABLESYMMETRICMATRIX_H
