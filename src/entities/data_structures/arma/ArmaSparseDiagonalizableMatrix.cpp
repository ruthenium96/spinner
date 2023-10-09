#include "ArmaSparseDiagonalizableMatrix.h"

#include "ArmaLogic.h"

namespace quantum::linear_algebra {

void ArmaSparseDiagonalizableMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    sparseDiagonalizableMatrix_.at(i, j) += value;
    if (i != j) {
        sparseDiagonalizableMatrix_(j, i) += value;
    }
}

void ArmaSparseDiagonalizableMatrix::resize(uint32_t matrix_in_space_basis_size_i) {
    sparseDiagonalizableMatrix_.resize(matrix_in_space_basis_size_i, matrix_in_space_basis_size_i);
}

EigenCouple ArmaSparseDiagonalizableMatrix::diagonalizeValuesVectors() const {
    ArmaLogic logic;

    return logic.diagonalizeValuesVectors(*this);
}

std::unique_ptr<AbstractDenseVector> ArmaSparseDiagonalizableMatrix::diagonalizeValues() const {
    ArmaLogic logic;

    return logic.diagonalizeValues(*this);
}

std::unique_ptr<AbstractDiagonalizableMatrix>
ArmaSparseDiagonalizableMatrix::multiply_by(double multiplier) const {
    auto answer = std::make_unique<ArmaSparseDiagonalizableMatrix>();
    answer->resize(size());
    answer->sparseDiagonalizableMatrix_ = sparseDiagonalizableMatrix_ * multiplier;
    return answer;
}

uint32_t ArmaSparseDiagonalizableMatrix::size() const {
    return sparseDiagonalizableMatrix_.n_rows;
}

double ArmaSparseDiagonalizableMatrix::at(uint32_t i, uint32_t j) const {
    return sparseDiagonalizableMatrix_.at(i, j);
}

void ArmaSparseDiagonalizableMatrix::print(std::ostream& os) const {
    os << sparseDiagonalizableMatrix_ << std::endl;
}

const arma::sp_mat& ArmaSparseDiagonalizableMatrix::getSparseSymmetricMatrix() const {
    return sparseDiagonalizableMatrix_;
}

arma::sp_mat& ArmaSparseDiagonalizableMatrix::modifySparseSymmetricMatrix() {
    return sparseDiagonalizableMatrix_;
}
}  // namespace quantum::linear_algebra