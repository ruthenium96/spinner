#include "EigenDenseDiagonalizableMatrix.h"

#include "EigenLogic.h"

namespace quantum::linear_algebra {

void EigenDenseDiagonalizableMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    denseDiagonalizableMatrix_(i, j) += value;
    if (i != j) {
        denseDiagonalizableMatrix_(j, i) += value;
    }
}

EigenCouple EigenDenseDiagonalizableMatrix::diagonalizeValuesVectors() const {
    EigenLogic eigenLogic;

    return eigenLogic.diagonalizeValuesVectors(*this);
}

std::unique_ptr<AbstractDenseVector> EigenDenseDiagonalizableMatrix::diagonalizeValues() const {
    EigenLogic eigenLogic;

    return eigenLogic.diagonalizeValues(*this);
}

std::unique_ptr<AbstractDiagonalizableMatrix>
EigenDenseDiagonalizableMatrix::multiply_by(double multiplier) const {
    auto answer = std::make_unique<EigenDenseDiagonalizableMatrix>();
    answer->resize(size());

    answer->denseDiagonalizableMatrix_ = denseDiagonalizableMatrix_ * multiplier;

    return answer;
}

uint32_t EigenDenseDiagonalizableMatrix::size() const {
    return denseDiagonalizableMatrix_.cols();
}

double EigenDenseDiagonalizableMatrix::at(uint32_t i, uint32_t j) const {
    return denseDiagonalizableMatrix_(i, j);
}

void EigenDenseDiagonalizableMatrix::print(std::ostream& os) const {
    os << denseDiagonalizableMatrix_ << std::endl;
}
void EigenDenseDiagonalizableMatrix::resize(uint32_t size) {
    denseDiagonalizableMatrix_.resize(size, size);
    denseDiagonalizableMatrix_.fill(0);
}

const Eigen::Matrix<double, -1, -1>& EigenDenseDiagonalizableMatrix::getDenseDiagonalizableMatrix() const {
    return denseDiagonalizableMatrix_;
}

Eigen::Matrix<double, -1, -1>& EigenDenseDiagonalizableMatrix::modifyDenseDiagonalizableMatrix() {
    return denseDiagonalizableMatrix_;
}
}  // namespace quantum::linear_algebra