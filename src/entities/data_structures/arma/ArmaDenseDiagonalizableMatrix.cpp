#include "ArmaDenseDiagonalizableMatrix.h"

#include "ArmaLogic.h"

namespace quantum::linear_algebra {

void ArmaDenseDiagonalizableMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    denseDiagonalizableMatrix_(i, j) += value;
    if (i != j) {
        denseDiagonalizableMatrix_(j, i) += value;
    }
}

EigenCouple ArmaDenseDiagonalizableMatrix::diagonalizeValuesVectors() const {
    ArmaLogic logic;

    return logic.diagonalizeValuesVectors(*this);
}

std::unique_ptr<AbstractDenseVector> ArmaDenseDiagonalizableMatrix::diagonalizeValues() const {
    ArmaLogic logic;

    return logic.diagonalizeValues(*this);
}

std::unique_ptr<AbstractDiagonalizableMatrix>
ArmaDenseDiagonalizableMatrix::multiply_by(double multiplier) const {
    auto answer = std::make_unique<ArmaDenseDiagonalizableMatrix>();
    answer->resize(size());
    answer->denseDiagonalizableMatrix_ = denseDiagonalizableMatrix_ * multiplier;
    return answer;
}

uint32_t ArmaDenseDiagonalizableMatrix::size() const {
    return denseDiagonalizableMatrix_.n_rows;
}

double ArmaDenseDiagonalizableMatrix::at(uint32_t i, uint32_t j) const {
    return denseDiagonalizableMatrix_.at(i, j);
}

void ArmaDenseDiagonalizableMatrix::print(std::ostream& os) const {
    os << denseDiagonalizableMatrix_ << std::endl;
}

void ArmaDenseDiagonalizableMatrix::resize(uint32_t size) {
    denseDiagonalizableMatrix_.resize(size, size);
}

const arma::dmat& ArmaDenseDiagonalizableMatrix::getDenseDiagonalizableMatrix() const {
    return denseDiagonalizableMatrix_;
}

arma::dmat& ArmaDenseDiagonalizableMatrix::modifyDenseDiagonalizableMatrix() {
    return denseDiagonalizableMatrix_;
}
}  // namespace quantum::linear_algebra