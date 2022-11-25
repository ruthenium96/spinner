#include "EigenDenseFactory.h"

#include "EigenDenseMatrix.h"
#include "EigenDenseVector.h"

namespace quantum::linear_algebra {
std::unique_ptr<AbstractDenseMatrix> EigenDenseFactory::createMatrix(
    uint32_t matrix_in_space_basis_size_i,
    uint32_t matrix_in_space_basis_size_j) {
    auto matrix = std::make_unique<EigenDenseMatrix>();
    matrix->resize(matrix_in_space_basis_size_i, matrix_in_space_basis_size_j);
    return matrix;
}

std::unique_ptr<AbstractDenseVector> EigenDenseFactory::createVector() {
    auto vector = std::make_unique<EigenDenseVector>();
    return vector;
}

std::vector<double>
EigenDenseFactory::concatenate(const std::vector<std::unique_ptr<AbstractDenseVector>>& vectors) {
    size_t size = 0;
    for (const auto& dense_vector : vectors) {
        size += dense_vector->size();
    }
    std::vector<double> result_vector;
    result_vector.reserve(size);
    for (auto& vector : vectors) {
        auto vector_ = EigenDenseVector::downcast_ptr(vector);
        // Until Eigen 3.4:
        result_vector.insert(
            result_vector.end(),
            vector_->vector_.data(),
            vector_->vector_.data() + vector_->vector_.size());
    }
    return result_vector;
}

}  // namespace quantum::linear_algebra