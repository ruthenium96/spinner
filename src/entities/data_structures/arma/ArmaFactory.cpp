#include "ArmaFactory.h"

#include "ArmaMatrix.h"
#include "ArmaVector.h"

namespace quantum::linear_algebra {
std::unique_ptr<AbstractMatrix> ArmaFactory::createMatrix() {
    auto matrix = std::make_unique<ArmaMatrix>();
    return matrix;
}

std::unique_ptr<AbstractVector> ArmaFactory::createVector() {
    auto vector = std::make_unique<ArmaVector>();
    return vector;
}

std::vector<double>
ArmaFactory::concatenate(const std::vector<std::unique_ptr<AbstractVector>>& vectors) {
    size_t size = 0;
    for (const auto& dense_vector : vectors) {
        size += dense_vector->size();
    }
    std::vector<double> result_vector;
    result_vector.reserve(size);
    for (auto& vector : vectors) {
        auto vector_ = ArmaVector::downcast_ptr(vector);
        result_vector.insert(result_vector.end(), vector_->vector_.begin(), vector_->vector_.end());
    }
    return result_vector;
}
}  // namespace quantum::linear_algebra