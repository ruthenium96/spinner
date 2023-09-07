#include "EmhashSparseSymmetricMatrix.h"

namespace quantum::linear_algebra {

void EmhashSparseSymmetricMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    if (!hashmap_.contains(i)) {
        hashmap_[i] = emhash7::HashMap<uint32_t, double>();
    }
    hashmap_[i][j] += value;
    if (i != j) {
        if (!hashmap_.contains(j)) {
            hashmap_[j] = emhash7::HashMap<uint32_t, double>();
        }
        hashmap_[j][i] += value;
    }
}

uint32_t EmhashSparseSymmetricMatrix::size() const {
    return size_;
}

double EmhashSparseSymmetricMatrix::at(uint32_t i, uint32_t j) const noexcept {
    if (!hashmap_.contains(i)) {
        return 0;
    }
    return hashmap_.at(i).at(j);
}

void EmhashSparseSymmetricMatrix::resize(size_t size) {
    size_ = size;
}

void EmhashSparseSymmetricMatrix::print(std::ostream& os) const {
    for (const auto& p : hashmap_) {
        auto i = p.first;
        for (const auto& pp : p.second) {
            auto j = pp.first;
            auto value = pp.second;
            os << "(" << i << ", " << j << ") :" << value << std::endl;
        }
    }
}

}  // namespace quantum::linear_algebra