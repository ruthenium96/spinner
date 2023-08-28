#include "EigenSparseDiagonalizableMatrix.h"

#include "EigenLogic.h"

namespace quantum::linear_algebra {

void EigenSparseDiagonalizableMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    sparseDiagonalizableMatrix_.coeffRef(i, j) += value;
    if (i != j) {
        sparseDiagonalizableMatrix_.coeffRef(j, i) += value;
    }
}

EigenCouple EigenSparseDiagonalizableMatrix::diagonalizeValuesVectors() const {
    EigenLogic eigenLogic;

    return eigenLogic.diagonalizeValuesVectors(*this);
}

std::unique_ptr<AbstractDenseVector> EigenSparseDiagonalizableMatrix::diagonalizeValues() const {
    EigenLogic eigenLogic;

    return eigenLogic.diagonalizeValues(*this);
}

uint32_t EigenSparseDiagonalizableMatrix::size() const {
    return sparseDiagonalizableMatrix_.innerSize();
}

double EigenSparseDiagonalizableMatrix::at(uint32_t i, uint32_t j) const {
    return sparseDiagonalizableMatrix_.coeff(i, j);
}

void EigenSparseDiagonalizableMatrix::resize(size_t size) {
    sparseDiagonalizableMatrix_.resize(size, size);
}

const Eigen::SparseMatrix<double>&
EigenSparseDiagonalizableMatrix::getSparseDiagonalizableMatrix() const {
    return sparseDiagonalizableMatrix_;
}

Eigen::SparseMatrix<double>& EigenSparseDiagonalizableMatrix::modifySparseDiagonalizableMatrix() {
    return sparseDiagonalizableMatrix_;
}

void EigenSparseDiagonalizableMatrix::print(std::ostream& os) const {
    os << sparseDiagonalizableMatrix_ << std::endl;
}
}  // namespace quantum::linear_algebra