#include "EigenDenseSymmetricMatrix.h"

#include "EigenLogic.h"

namespace quantum::linear_algebra {

void EigenDenseSymmetricMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    denseSymmetricMatrix_(i, j) += value;
    if (i != j) {
        denseSymmetricMatrix_(j, i) += value;
    }
}

EigenCouple EigenDenseSymmetricMatrix::diagonalizeValuesVectors() const {
    EigenLogic eigenLogic;

    return eigenLogic.diagonalizeValuesVectors(*this);
}

std::unique_ptr<AbstractDenseVector> EigenDenseSymmetricMatrix::diagonalizeValues() const {
    EigenLogic eigenLogic;

    return eigenLogic.diagonalizeValues(*this);
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

const Eigen::MatrixXd& EigenDenseSymmetricMatrix::getDenseSymmetricMatrix() const {
    return denseSymmetricMatrix_;
}

Eigen::MatrixXd& EigenDenseSymmetricMatrix::modifyDenseSymmetricMatrix() {
    return denseSymmetricMatrix_;
}
}  // namespace quantum::linear_algebra