#include "EigenFactory.h"

#include "EigenMatrix.h"
#include "EigenVector.h"

namespace quantum::linear_algebra {
std::unique_ptr<AbstractMatrix> EigenFactory::createMatrix() {
    auto matrix = std::make_unique<EigenMatrix>();
    return matrix;
}

std::unique_ptr<AbstractVector> EigenFactory::createVector() {
    auto vector = std::make_unique<EigenVector>();
    return vector;
}

std::vector<double>
EigenFactory::concatenate(const std::vector<std::unique_ptr<AbstractVector>>& vectors) {
    size_t size = 0;
    for (const auto& dense_vector : vectors) {
        size += dense_vector->size();
    }
    std::vector<double> result_vector;
    result_vector.reserve(size);
    for (auto& vector : vectors) {
        auto vector_ = EigenVector::downcast_ptr(vector);
        // Until Eigen 3.4:
        result_vector.insert(
            result_vector.end(),
            vector_->vector_.data(),
            vector_->vector_.data() + vector_->vector_.size());
    }
    return result_vector;
}

}  // namespace quantum::linear_algebra