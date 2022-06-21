#include "ArmaMatrix.h"

#include "src/entities/data_structures/AbstractVector.h"

namespace quantum::linear_algebra {

void ArmaMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    matrix_(i, j) += value;
}

void ArmaMatrix::assign_to_position(double value, uint32_t i, uint32_t j) {
    matrix_(i, j) = value;
}

void ArmaMatrix::resize(
    uint32_t matrix_in_space_basis_size_i,
    uint32_t matrix_in_space_basis_size_j) {
    // TODO: is it the fastest way to initialize zero matrix?
    matrix_.resize(matrix_in_space_basis_size_i, matrix_in_space_basis_size_j);
    matrix_.fill(arma::fill::zeros);
}

void ArmaMatrix::resize_with_nans(
    uint32_t matrix_in_space_basis_size_i,
    uint32_t matrix_in_space_basis_size_j) {
    // TODO: is it the fastest way to initialize NaNs matrix?
    matrix_.resize(matrix_in_space_basis_size_i, matrix_in_space_basis_size_j);
    matrix_.fill(arma::datum::nan);
}

uint32_t ArmaMatrix::size() const {
    return matrix_.n_rows;
}

uint32_t ArmaMatrix::size_rows() const {
    return matrix_.n_rows;
}

uint32_t ArmaMatrix::size_cols() const {
    return matrix_.n_rows;
}

double ArmaMatrix::at(uint32_t i, uint32_t j) const {
    return matrix_(i, j);
}

std::unique_ptr<AbstractMatrix>
ArmaMatrix::unitary_transform(const std::unique_ptr<AbstractMatrix>& matrix_to_transform) const {
    auto matrix_to_transform_ = downcast_ptr(matrix_to_transform);

    auto transformed_matrix = std::make_unique<ArmaMatrix>();
    transformed_matrix->matrix_ = matrix_.t() * matrix_to_transform_->matrix_.t() * matrix_;

    return std::move(transformed_matrix);
}

EigenCouple ArmaMatrix::diagonalizeValuesVectors() const {
    auto eigenvalues_ = std::make_unique<ArmaVector>();
    auto eigenvectors_ = std::make_unique<ArmaMatrix>();

    arma::eig_sym(eigenvalues_->vector_, eigenvectors_->matrix_, matrix_);

    EigenCouple answer;
    answer.eigenvalues = std::move(eigenvalues_);
    answer.eigenvectors = std::move(eigenvectors_);

    return answer;
}

std::unique_ptr<AbstractVector> ArmaMatrix::diagonalizeValues() const {
    auto eigenvalues = std::make_unique<ArmaVector>();

    arma::eig_sym(eigenvalues->vector_, matrix_);

    return eigenvalues;
}

std::unique_ptr<AbstractVector> ArmaMatrix::return_main_diagonal() const {
    auto main_diagonal = std::make_unique<ArmaVector>();
    main_diagonal->vector_ = matrix_.diag();
    return std::move(main_diagonal);
}

void ArmaMatrix::print(std::ostream& os) const {
    os << matrix_ << std::endl;
}

const ArmaMatrix* ArmaMatrix::downcast_ptr(const std::unique_ptr<AbstractMatrix>& ptr) {
    return dynamic_cast<const ArmaMatrix*>(ptr.get());
}
}  // namespace quantum::linear_algebra