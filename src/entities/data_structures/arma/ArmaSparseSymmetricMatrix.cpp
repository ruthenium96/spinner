#include "ArmaSparseSymmetricMatrix.h"

#include "ArmaLogic.h"

namespace quantum::linear_algebra {

void ArmaSparseSymmetricMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    sparseSymmetricMatrix_.at(i, j) += value;
}

void ArmaSparseSymmetricMatrix::resize(uint32_t matrix_in_space_basis_size_i) {
    sparseSymmetricMatrix_.resize(matrix_in_space_basis_size_i, matrix_in_space_basis_size_i);
}

EigenCouple ArmaSparseSymmetricMatrix::diagonalizeValuesVectors() const {
    ArmaLogic logic;

    return logic.diagonalizeValuesVectors(*this);
}

std::unique_ptr<AbstractDenseVector> ArmaSparseSymmetricMatrix::diagonalizeValues() const {
    ArmaLogic logic;

    return logic.diagonalizeValues(*this);
}

uint32_t ArmaSparseSymmetricMatrix::size() const {
    return sparseSymmetricMatrix_.n_rows;
}

double ArmaSparseSymmetricMatrix::at(uint32_t i, uint32_t j) const {
    return sparseSymmetricMatrix_.at(i, j);
}

void ArmaSparseSymmetricMatrix::print(std::ostream& os) const {
    os << sparseSymmetricMatrix_ << std::endl;
}

const arma::sp_mat& ArmaSparseSymmetricMatrix::getSparseSymmetricMatrix() const {
    return sparseSymmetricMatrix_;
}

arma::sp_mat& ArmaSparseSymmetricMatrix::modifySparseSymmetricMatrix() {
    return sparseSymmetricMatrix_;
}
}  // namespace quantum::linear_algebra