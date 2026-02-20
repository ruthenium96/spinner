#ifndef SPINNER_ABSTRACTDENSEVECTOR_H
#define SPINNER_ABSTRACTDENSEVECTOR_H

#include <cstdint>
#include <memory>

namespace spinner::linalg_structures {
/**
 * @class AbstractDenseVector
 * @brief Abstract interface for a dense vector of real numbers.
 *
 * Provides arithmetic operations, element access, and utilities like
 * concatenation, exponentiation, and different products.
 */
class AbstractDenseVector {
  public:
    /**
     * @brief Appends the contents of another vector to the end of this vector.
     * @param rhs Const vector to concatenate with.
     */
    virtual void concatenate_with(const std::unique_ptr<AbstractDenseVector>& rhs) = 0;
    /**
     * @brief Appends a given number of identical values to the end of the vector.
     * @param number How many entries to add.
     * @param value  The value to insert.
     */
    virtual void add_identical_values(size_t number, double value) = 0;
    /**
     * @brief Shifts the vector so that its minimum element becomes zero.
     *
     * After the call, `new_value = old_value - min(old_value)`.
     */
    virtual void subtract_minimum() = 0;
    /**
     * @brief Applies the exponential function element‑wise to the vector.
     *
     * Each element `x` is replaced by `exp(x)`.
     */
    virtual void wise_exp() = 0;

    /**
     * @brief Returns a new vector that is the current vector multiplied by a scalar.
     * @param multiplier Scalar factor.
     * @return New vector = multiplier * this.
     */
    virtual std::unique_ptr<AbstractDenseVector> multiply_by(double multiplier) const = 0;
    /**
     * @brief Computes the dot product with another vector (einsum "ii->").
     * Both vectors must have the same size.
     * @param rhs The other vector.
     * @return `this · rhs`.
     */
    virtual double dot(const std::unique_ptr<AbstractDenseVector>& rhs) const = 0;
    /**
     * @brief Computes a triple dot product (einsum "iii->").
     *
     * Equivalent to sum over i of this[i] * second[i] * third[i].
     * All vectors must have the same size.
     *
     * @param second Second vector.
     * @param third  Third vector.
     * @return Scalar result.
     */    
    virtual double triple_dot(const std::unique_ptr<AbstractDenseVector>& second, 
      const std::unique_ptr<AbstractDenseVector>& third) const = 0;
    /**
     * @brief Returns a new vector whose elements are the element‑wise product of this and rhs (einsum "ii->i").
     * Both vectors must have the same size.
     * @param rhs Other vector.
     * @return New vector with elements `this[i] * rhs[i]`.
     */
    virtual std::unique_ptr<AbstractDenseVector>
    element_wise_multiplication(const std::unique_ptr<AbstractDenseVector>& rhs) const = 0;

    /**
     * @brief Returns the number of elements in the vector.
     * @return Vector size.
     */
    virtual uint32_t size() const = 0;
    /**
     * @brief Returns the element at index i.
     * @param i Index (0‑based).
     * @return Value.
     */
    virtual double at(uint32_t i) const = 0;

    /**
     * @brief Prints a human‑readable representation of the vector to an output stream.
     * @param os Output stream.
     */
    virtual void print(std::ostream& os) const = 0;

    /// Virtual destructor.
    virtual ~AbstractDenseVector() = default;
};
} // namespace spinner::linalg_structures

#endif  //SPINNER_ABSTRACTDENSEVECTOR_H