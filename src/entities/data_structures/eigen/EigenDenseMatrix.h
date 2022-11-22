#ifndef SPINNER_EIGENDENSEMATRIX_H
#define SPINNER_EIGENDENSEMATRIX_H

#include <Eigen/Dense>

#include "EigenDenseVector.h"
#include "src/entities/data_structures/AbstractDenseMatrix.h"

namespace quantum::linear_algebra {
class EigenDenseFactory;
class EigenDenseMatrix: public AbstractDenseMatrix {
  public:
    void add_to_position(double value, uint32_t i, uint32_t j) override;
    void assign_to_position(double value, uint32_t i, uint32_t j) override;
    void
    resize(uint32_t matrix_in_space_basis_size_i, uint32_t matrix_in_space_basis_size_j) override;
    EigenCouple diagonalizeValuesVectors() const override;
    std::unique_ptr<AbstractDenseVector> diagonalizeValues() const override;
    std::unique_ptr<AbstractDenseVector> return_main_diagonal() const override;
    std::unique_ptr<AbstractDenseMatrix> unitary_transform(
        const std::unique_ptr<AbstractDenseMatrix>& matrix_to_transform) const override;
    uint32_t size() const override;
    uint32_t size_rows() const override;
    uint32_t size_cols() const override;
    double at(uint32_t i, uint32_t j) const override;
    void print(std::ostream& os) const override;

  private:
    friend EigenDenseFactory;
    Eigen::MatrixXd matrix_;

    static const EigenDenseMatrix* downcast_ptr(const std::unique_ptr<AbstractDenseMatrix>& ptr);
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_EIGENDENSEMATRIX_H
