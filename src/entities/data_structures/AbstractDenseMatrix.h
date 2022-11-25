#ifndef SPINNER_ABSTRACTDENSEMATRIX_H
#define SPINNER_ABSTRACTDENSEMATRIX_H

#include <cstdint>
#include <memory>

#include "AbstractDenseVector.h"

// TODO: rename namespace
namespace quantum::linear_algebra {

struct EigenCouple {
    std::unique_ptr<AbstractDenseVector> eigenvalues;
    std::unique_ptr<AbstractDenseMatrix> eigenvectors;
};

class AbstractDenseMatrix {
  public:
    virtual void add_to_position(double value, uint32_t i, uint32_t j) = 0;
    virtual void assign_to_position(double value, uint32_t i, uint32_t j) = 0;
    virtual EigenCouple diagonalizeValuesVectors() const = 0;
    virtual std::unique_ptr<AbstractDenseVector> diagonalizeValues() const = 0;
    virtual std::unique_ptr<AbstractDenseVector> return_main_diagonal() const = 0;

    // TODO: should we change input matrix instead?
    virtual std::unique_ptr<AbstractDenseMatrix>
    unitary_transform(const std::unique_ptr<AbstractDenseMatrix>& matrix_to_transform) const = 0;

    virtual uint32_t size() const = 0;
    virtual uint32_t size_rows() const = 0;
    virtual uint32_t size_cols() const = 0;
    virtual double at(uint32_t i, uint32_t j) const = 0;

    virtual void print(std::ostream& os) const = 0;

    virtual ~AbstractDenseMatrix() = default;

  protected:
    // TODO: can we remove it from this abstract class?
    virtual void
    resize(uint32_t matrix_in_space_basis_size_i, uint32_t matrix_in_space_basis_size_j) = 0;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_ABSTRACTDENSEMATRIX_H
