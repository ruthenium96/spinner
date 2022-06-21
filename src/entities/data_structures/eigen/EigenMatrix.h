#ifndef SPINNER_EIGENMATRIX_H
#define SPINNER_EIGENMATRIX_H

#include <Eigen/Dense>

#include "EigenVector.h"
#include "src/entities/data_structures/AbstractMatrix.h"

namespace quantum::linear_algebra {
class EigenFactory;
class EigenMatrix: public AbstractMatrix {
  public:
    void add_to_position(double value, uint32_t i, uint32_t j) override;
    void assign_to_position(double value, uint32_t i, uint32_t j) override;
    void
    resize(uint32_t matrix_in_space_basis_size_i, uint32_t matrix_in_space_basis_size_j) override;
    void resize_with_nans(
        uint32_t matrix_in_space_basis_size_i,
        uint32_t matrix_in_space_basis_size_j) override;
    EigenCouple diagonalizeValuesVectors() const override;
    std::unique_ptr<AbstractVector> diagonalizeValues() const override;
    std::unique_ptr<AbstractVector> return_main_diagonal() const override;
    std::unique_ptr<AbstractMatrix>
    unitary_transform(const std::unique_ptr<AbstractMatrix>& matrix_to_transform) const override;
    uint32_t size() const override;
    uint32_t size_rows() const override;
    uint32_t size_cols() const override;
    double at(uint32_t i, uint32_t j) const override;
    void print(std::ostream& os) const override;

  private:
    friend EigenFactory;
    Eigen::MatrixXd matrix_;

    static const EigenMatrix* downcast_ptr(const std::unique_ptr<AbstractMatrix>& ptr);
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_EIGENMATRIX_H
