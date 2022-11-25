#ifndef SPINNER_ARMADENSEMATRIX_H
#define SPINNER_ARMADENSEMATRIX_H

#include <armadillo>
#include <cstdint>

#include "ArmaDenseVector.h"
#include "src/entities/data_structures/AbstractDenseMatrix.h"

namespace quantum::linear_algebra {
class ArmaDenseFactory;
class ArmaDenseMatrix: public AbstractDenseMatrix {
    friend ArmaDenseFactory;

  public:
    void add_to_position(double value, uint32_t i, uint32_t j) override;
    void assign_to_position(double value, uint32_t i, uint32_t j) override;
    void
    resize(uint32_t matrix_in_space_basis_size_i, uint32_t matrix_in_space_basis_size_j) override;
    std::unique_ptr<AbstractDenseMatrix> unitary_transform(
        const std::unique_ptr<AbstractDenseMatrix>& matrix_to_transform) const override;

    EigenCouple diagonalizeValuesVectors() const override;
    std::unique_ptr<AbstractDenseVector> diagonalizeValues() const override;
    std::unique_ptr<AbstractDenseVector> return_main_diagonal() const override;

    uint32_t size() const override;
    uint32_t size_rows() const override;
    uint32_t size_cols() const override;
    double at(uint32_t i, uint32_t j) const override;
    void print(std::ostream& os) const override;

  protected:
    void
    resize(uint32_t matrix_in_space_basis_size_i, uint32_t matrix_in_space_basis_size_j) override;

  private:
    arma::dmat matrix_;
    // c-like pointers are necessary to avoid double-free error
    static const ArmaDenseMatrix* downcast_ptr(const std::unique_ptr<AbstractDenseMatrix>& ptr);
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_ARMADENSEMATRIX_H
