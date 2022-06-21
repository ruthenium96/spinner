#include "EigenMatrix.h"

namespace quantum::linear_algebra {
void EigenMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    matrix_(i, j) += value;
}

void EigenMatrix::assign_to_position(double value, uint32_t i, uint32_t j) {
    matrix_(i, j) = value;
}

void EigenMatrix::resize(
    uint32_t matrix_in_space_basis_size_i,
    uint32_t matrix_in_space_basis_size_j) {
    matrix_.resize(matrix_in_space_basis_size_i, matrix_in_space_basis_size_j);
    matrix_.fill(0);
}

void EigenMatrix::resize_with_nans(
    uint32_t matrix_in_space_basis_size_i,
    uint32_t matrix_in_space_basis_size_j) {
    matrix_.resize(matrix_in_space_basis_size_i, matrix_in_space_basis_size_j);
    matrix_.fill(NAN);
}

EigenCouple EigenMatrix::diagonalizeValuesVectors() const {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(matrix_, Eigen::ComputeEigenvectors);

    auto eigenvalues_ = std::make_unique<EigenVector>();
    eigenvalues_->vector_ = es.eigenvalues();
    auto eigenvectors_ = std::make_unique<EigenMatrix>();
    eigenvectors_->matrix_ = es.eigenvectors();

    EigenCouple answer = {std::move(eigenvalues_), std::move(eigenvectors_)};

    return answer;
}

std::unique_ptr<AbstractVector> EigenMatrix::diagonalizeValues() const {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(matrix_, Eigen::EigenvaluesOnly);

    auto eigenvalues_ = std::make_unique<EigenVector>();
    eigenvalues_->vector_ = es.eigenvalues();

    return eigenvalues_;
}

std::unique_ptr<AbstractVector> EigenMatrix::return_main_diagonal() const {
    auto answer = std::make_unique<EigenVector>();

    answer->vector_ = matrix_.diagonal();
    return answer;
}

std::unique_ptr<AbstractMatrix>
EigenMatrix::unitary_transform(const std::unique_ptr<AbstractMatrix>& matrix_to_transform) const {
    auto matrix_to_transform_ = downcast_ptr(matrix_to_transform);

    auto transformed_matrix = std::make_unique<EigenMatrix>();
    transformed_matrix->matrix_ =
        matrix_.transpose() * matrix_to_transform_->matrix_.transpose() * matrix_;

    return std::move(transformed_matrix);
}

uint32_t EigenMatrix::size() const {
    // TODO: is it good idea?
    return matrix_.rows();
}

uint32_t EigenMatrix::size_rows() const {
    return matrix_.rows();
}

uint32_t EigenMatrix::size_cols() const {
    return matrix_.cols();
}

double EigenMatrix::at(uint32_t i, uint32_t j) const {
    return matrix_(i, j);
}

void EigenMatrix::print(std::ostream& os) const {
    os << matrix_ << std::endl;
}

const EigenMatrix* EigenMatrix::downcast_ptr(const std::unique_ptr<AbstractMatrix>& ptr) {
    return dynamic_cast<const EigenMatrix*>(ptr.get());
}

}  // namespace quantum::linear_algebra