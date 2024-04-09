#include "EigenDenseSemiunitaryMatrix.h"

#include "EigenDenseVector.h"
#include "EigenLogic.h"

namespace quantum::linear_algebra {

uint32_t EigenDenseSemiunitaryMatrix::size_rows() const {
    return denseSemiunitaryMatrix_.rows();
}

uint32_t EigenDenseSemiunitaryMatrix::size_cols() const {
    return denseSemiunitaryMatrix_.cols();
}

void EigenDenseSemiunitaryMatrix::resize(size_t size_rows, size_t size_cols) {
    denseSemiunitaryMatrix_.resize(size_rows, size_cols);
    denseSemiunitaryMatrix_.fill(0);
}

double EigenDenseSemiunitaryMatrix::at(uint32_t i, uint32_t j) const {
    return denseSemiunitaryMatrix_(j, i);
}

void EigenDenseSemiunitaryMatrix::print(std::ostream& os) const {
    os << denseSemiunitaryMatrix_ << std::endl;
}

std::unique_ptr<AbstractDenseVector>
EigenDenseSemiunitaryMatrix::unitaryTransformAndReturnMainDiagonal(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& matrix_to_transform) const {
    EigenLogic eigenLogic;

    return eigenLogic.unitaryTransformAndReturnMainDiagonal(matrix_to_transform, *this);
}

std::unique_ptr<AbstractDiagonalizableMatrix> EigenDenseSemiunitaryMatrix::unitaryTransform(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& matrix_to_transform) const {
    EigenLogic eigenLogic;

    return eigenLogic.unitaryTransform(matrix_to_transform, *this);
}

const Eigen::Matrix<double, -1, -1>& EigenDenseSemiunitaryMatrix::getDenseSemiunitaryMatrix() const {
    return denseSemiunitaryMatrix_;
}

Eigen::Matrix<double, -1, -1>& EigenDenseSemiunitaryMatrix::modifyDenseSemiunitaryMatrix() {
    return denseSemiunitaryMatrix_;
}

void EigenDenseSemiunitaryMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    denseSemiunitaryMatrix_(j, i) += value;
}

void EigenDenseSemiunitaryMatrix::normalize() {
    for (int i = 0; i < denseSemiunitaryMatrix_.cols(); i++) {
        denseSemiunitaryMatrix_.col(i).normalize();
    }
}
}  // namespace quantum::linear_algebra