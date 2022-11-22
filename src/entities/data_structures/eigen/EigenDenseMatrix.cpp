#include "EigenDenseMatrix.h"

namespace quantum::linear_algebra {
void EigenDenseMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    matrix_(i, j) += value;
}

void EigenDenseMatrix::assign_to_position(double value, uint32_t i, uint32_t j) {
    matrix_(i, j) = value;
}

void EigenDenseMatrix::resize(
    uint32_t matrix_in_space_basis_size_i,
    uint32_t matrix_in_space_basis_size_j) {
    matrix_.resize(matrix_in_space_basis_size_i, matrix_in_space_basis_size_j);
    matrix_.fill(0);
}

EigenCouple EigenDenseMatrix::diagonalizeValuesVectors() const {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(matrix_, Eigen::ComputeEigenvectors);

    auto eigenvalues_ = std::make_unique<EigenDenseVector>();
    eigenvalues_->vector_ = es.eigenvalues();
    auto eigenvectors_ = std::make_unique<EigenDenseMatrix>();
    eigenvectors_->matrix_ = es.eigenvectors();

    EigenCouple answer = {std::move(eigenvalues_), std::move(eigenvectors_)};

    return answer;
}

std::unique_ptr<AbstractDenseVector> EigenDenseMatrix::diagonalizeValues() const {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(matrix_, Eigen::EigenvaluesOnly);

    auto eigenvalues_ = std::make_unique<EigenDenseVector>();
    eigenvalues_->vector_ = es.eigenvalues();

    return eigenvalues_;
}

std::unique_ptr<AbstractDenseVector> EigenDenseMatrix::return_main_diagonal() const {
    auto answer = std::make_unique<EigenDenseVector>();

    answer->vector_ = matrix_.diagonal();
    return answer;
}

std::unique_ptr<AbstractDenseMatrix> EigenDenseMatrix::unitary_transform(
    const std::unique_ptr<AbstractDenseMatrix>& matrix_to_transform) const {
    auto matrix_to_transform_ = downcast_ptr(matrix_to_transform);

    auto transformed_matrix = std::make_unique<EigenDenseMatrix>();
    transformed_matrix->matrix_ =
        matrix_.transpose() * matrix_to_transform_->matrix_.transpose() * matrix_;

    return std::move(transformed_matrix);
}

uint32_t EigenDenseMatrix::size() const {
    // TODO: is it good idea?
    return matrix_.rows();
}

uint32_t EigenDenseMatrix::size_rows() const {
    return matrix_.rows();
}

uint32_t EigenDenseMatrix::size_cols() const {
    return matrix_.cols();
}

double EigenDenseMatrix::at(uint32_t i, uint32_t j) const {
    return matrix_(i, j);
}

void EigenDenseMatrix::print(std::ostream& os) const {
    os << matrix_ << std::endl;
}

const EigenDenseMatrix*
EigenDenseMatrix::downcast_ptr(const std::unique_ptr<AbstractDenseMatrix>& ptr) {
    auto answer = dynamic_cast<const EigenDenseMatrix*>(ptr.get());
    if (answer == nullptr) {
        throw std::bad_cast();
    }
    return answer;
}

}  // namespace quantum::linear_algebra