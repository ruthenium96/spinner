#include "ArmaDenseMatrix.h"

#include "src/entities/data_structures/AbstractDenseVector.h"

namespace quantum::linear_algebra {

void ArmaDenseMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    matrix_(i, j) += value;
}

void ArmaDenseMatrix::assign_to_position(double value, uint32_t i, uint32_t j) {
    matrix_(i, j) = value;
}

void ArmaDenseMatrix::resize(
    uint32_t matrix_in_space_basis_size_i,
    uint32_t matrix_in_space_basis_size_j) {
    // TODO: is it the fastest way to initialize zero matrix?
    matrix_.resize(matrix_in_space_basis_size_i, matrix_in_space_basis_size_j);
    matrix_.fill(arma::fill::zeros);
}

uint32_t ArmaDenseMatrix::size() const {
    return matrix_.n_rows;
}

uint32_t ArmaDenseMatrix::size_rows() const {
    return matrix_.n_rows;
}

uint32_t ArmaDenseMatrix::size_cols() const {
    return matrix_.n_rows;
}

double ArmaDenseMatrix::at(uint32_t i, uint32_t j) const {
    return matrix_(i, j);
}

std::unique_ptr<AbstractDenseMatrix> ArmaDenseMatrix::unitary_transform(
    const std::unique_ptr<AbstractDenseMatrix>& matrix_to_transform) const {
    auto matrix_to_transform_ = downcast_ptr(matrix_to_transform);

    auto transformed_matrix = std::make_unique<ArmaDenseMatrix>();
    transformed_matrix->matrix_ = matrix_.t() * matrix_to_transform_->matrix_.t() * matrix_;

    return std::move(transformed_matrix);
}

EigenCouple ArmaDenseMatrix::diagonalizeValuesVectors() const {
    auto eigenvalues_ = std::make_unique<ArmaDenseVector>();
    auto eigenvectors_ = std::make_unique<ArmaDenseMatrix>();

    arma::eig_sym(eigenvalues_->vector_, eigenvectors_->matrix_, matrix_);

    EigenCouple answer;
    answer.eigenvalues = std::move(eigenvalues_);
    answer.eigenvectors = std::move(eigenvectors_);

    return answer;
}

std::unique_ptr<AbstractDenseVector> ArmaDenseMatrix::diagonalizeValues() const {
    auto eigenvalues = std::make_unique<ArmaDenseVector>();

    arma::eig_sym(eigenvalues->vector_, matrix_);

    return eigenvalues;
}

std::unique_ptr<AbstractDenseVector> ArmaDenseMatrix::return_main_diagonal() const {
    auto main_diagonal = std::make_unique<ArmaDenseVector>();
    main_diagonal->vector_ = matrix_.diag();
    return std::move(main_diagonal);
}

void ArmaDenseMatrix::print(std::ostream& os) const {
    os << matrix_ << std::endl;
}

const ArmaDenseMatrix*
ArmaDenseMatrix::downcast_ptr(const std::unique_ptr<AbstractDenseMatrix>& ptr) {
    auto answer = dynamic_cast<const ArmaDenseMatrix*>(ptr.get());
    if (answer == nullptr) {
        throw std::bad_cast();
    }
    return answer;
}
}  // namespace quantum::linear_algebra