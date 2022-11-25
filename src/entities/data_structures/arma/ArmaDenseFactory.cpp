#include "ArmaDenseFactory.h"

#include "ArmaDenseMatrix.h"
#include "ArmaDenseVector.h"

namespace quantum::linear_algebra {
std::unique_ptr<AbstractDenseMatrix> ArmaDenseFactory::createMatrix(
    uint32_t matrix_in_space_basis_size_i,
    uint32_t matrix_in_space_basis_size_j) {
    auto matrix = std::make_unique<ArmaDenseMatrix>();
    matrix->resize(matrix_in_space_basis_size_i, matrix_in_space_basis_size_j);
    return matrix;
}

std::unique_ptr<AbstractDenseVector> ArmaDenseFactory::createVector() {
    auto vector = std::make_unique<ArmaDenseVector>();
    return vector;
}

std::vector<double>
ArmaDenseFactory::concatenate(const std::vector<std::unique_ptr<AbstractDenseVector>>& vectors) {
    size_t size = 0;
    for (const auto& dense_vector : vectors) {
        size += dense_vector->size();
    }
    std::vector<double> result_vector;
    result_vector.reserve(size);
    for (auto& vector : vectors) {
        auto vector_ = ArmaDenseVector::downcast_ptr(vector);
        result_vector.insert(result_vector.end(), vector_->vector_.begin(), vector_->vector_.end());
    }
    return result_vector;
}
}  // namespace quantum::linear_algebra