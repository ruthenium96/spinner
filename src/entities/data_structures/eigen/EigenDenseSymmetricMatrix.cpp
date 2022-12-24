#include "EigenDenseSymmetricMatrix.h"

#include <Eigen/Dense>

#include "EigenDenseSemiunitaryMatrix.h"
#include "EigenDenseVector.h"
namespace quantum::linear_algebra {

void EigenDenseSymmetricMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    denseSymmetricMatrix_(i, j) += value;
}

EigenCouple EigenDenseSymmetricMatrix::diagonalizeValuesVectors() const {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(denseSymmetricMatrix_, Eigen::ComputeEigenvectors);

    auto eigenvalues_ = std::make_unique<EigenDenseVector>();
    eigenvalues_->vector_ = es.eigenvalues();
    auto eigenvectors_ = std::make_unique<EigenDenseSemiunitaryMatrix>();
    eigenvectors_->denseSemiunitaryMatrix_ = es.eigenvectors();

    EigenCouple answer = {std::move(eigenvalues_), std::move(eigenvectors_)};

    return answer;
}

std::unique_ptr<AbstractDenseVector> EigenDenseSymmetricMatrix::diagonalizeValues() const {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(denseSymmetricMatrix_, Eigen::EigenvaluesOnly);

    auto eigenvalues_ = std::make_unique<EigenDenseVector>();
    eigenvalues_->vector_ = es.eigenvalues();

    return eigenvalues_;
}

uint32_t EigenDenseSymmetricMatrix::size() const {
    return denseSymmetricMatrix_.cols();
}

double EigenDenseSymmetricMatrix::at(uint32_t i, uint32_t j) const {
    return denseSymmetricMatrix_(i, j);
}

void EigenDenseSymmetricMatrix::print(std::ostream& os) const {
    os << denseSymmetricMatrix_ << std::endl;
}
void EigenDenseSymmetricMatrix::resize(uint32_t matrix_in_space_basis_size_i) {
    denseSymmetricMatrix_.resize(matrix_in_space_basis_size_i, matrix_in_space_basis_size_i);
    denseSymmetricMatrix_.fill(0);
}
}  // namespace quantum::linear_algebra