#include "EigenSparseSymmetricMatrix.h"

#include "EigenLogic.h"

namespace quantum::linear_algebra {

void EigenSparseSymmetricMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    sparseSymmetricMatrix_.coeffRef(i, j) += value;
    //    sparseSymmetricMatrix_.insert(i, j) = value;
}

EigenCouple EigenSparseSymmetricMatrix::diagonalizeValuesVectors() const {
    EigenLogic eigenLogic;

    return eigenLogic.diagonalizeValuesVectors(*this);
}

std::unique_ptr<AbstractDenseVector> EigenSparseSymmetricMatrix::diagonalizeValues() const {
    EigenLogic eigenLogic;

    eigenLogic.diagonalizeValues(*this);
}

uint32_t EigenSparseSymmetricMatrix::size() const {
    return sparseSymmetricMatrix_.innerSize();
}

double EigenSparseSymmetricMatrix::at(uint32_t i, uint32_t j) const {
    return sparseSymmetricMatrix_.coeff(i, j);
}

void EigenSparseSymmetricMatrix::resize(size_t size) {
    sparseSymmetricMatrix_.resize(size, size);
}

const Eigen::SparseMatrix<double>& EigenSparseSymmetricMatrix::getSparseSymmetricMatrix() const {
    return sparseSymmetricMatrix_;
}

Eigen::SparseMatrix<double>& EigenSparseSymmetricMatrix::modifySparseSymmetricMatrix() {
    return sparseSymmetricMatrix_;
}

void EigenSparseSymmetricMatrix::print(std::ostream& os) const {
    os << sparseSymmetricMatrix_ << std::endl;
}
}  // namespace quantum::linear_algebra