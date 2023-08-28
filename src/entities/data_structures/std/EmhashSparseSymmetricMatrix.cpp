#include "EmhashSparseSymmetricMatrix.h"

namespace quantum::linear_algebra {

uint64_t SzudzikPair::operator()(const std::pair<uint32_t, uint32_t>& x) const noexcept {
    uint64_t result;
    if (x.first > x.second) {
        result = x.second * x.second + x.first;
    } else {
        result = x.first * x.first + x.first + x.second;
    }
    return result;
}

void EmhashSparseSymmetricMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    hashmap_[{i, j}] += value;
    if (i != j) {
        hashmap_[{j, i}] += value;
    }
}

uint32_t EmhashSparseSymmetricMatrix::size() const {
    return size_;
}

double EmhashSparseSymmetricMatrix::at(uint32_t i, uint32_t j) const noexcept {
    std::pair<uint32_t, uint32_t> pair = {i, j};
    auto value = hashmap_.at(pair);
    return value;
}

void EmhashSparseSymmetricMatrix::resize(size_t size) {
    size_ = size;
}

void EmhashSparseSymmetricMatrix::print(std::ostream& os) const {
    for (const auto& p : hashmap_) {
        auto [i, j] = p.first;
        os << "(" << i << ", " << j << ") :" << p.second << std::endl;
    }
}

}  // namespace quantum::linear_algebra