#include "ArmaDenseSymmetricMatrix.h"

#include "ArmaLogic.h"

namespace quantum::linear_algebra {

void ArmaDenseSymmetricMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    denseSymmetricMatrix_(i, j) += value;
}

EigenCouple ArmaDenseSymmetricMatrix::diagonalizeValuesVectors() const {
    ArmaLogic logic;

    return logic.diagonalizeValuesVectors(*this);
}

std::unique_ptr<AbstractDenseVector> ArmaDenseSymmetricMatrix::diagonalizeValues() const {
    ArmaLogic logic;

    return logic.diagonalizeValues(*this);
}

uint32_t ArmaDenseSymmetricMatrix::size() const {
    return denseSymmetricMatrix_.n_rows;
}

double ArmaDenseSymmetricMatrix::at(uint32_t i, uint32_t j) const {
    return denseSymmetricMatrix_.at(i, j);
}

void ArmaDenseSymmetricMatrix::print(std::ostream& os) const {
    os << denseSymmetricMatrix_;
}

void ArmaDenseSymmetricMatrix::resize(uint32_t matrix_in_space_basis_size_i) {
    denseSymmetricMatrix_.resize(matrix_in_space_basis_size_i, matrix_in_space_basis_size_i);
}

const arma::dmat& ArmaDenseSymmetricMatrix::getDenseSymmetricMatrix() const {
    return denseSymmetricMatrix_;
}

arma::dmat& ArmaDenseSymmetricMatrix::modifyDenseSymmetricMatrix() {
    return denseSymmetricMatrix_;
}
}  // namespace quantum::linear_algebra