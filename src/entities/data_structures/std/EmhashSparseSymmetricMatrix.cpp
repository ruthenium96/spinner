#include "EmhashSparseSymmetricMatrix.h"

namespace quantum::linear_algebra {

void EmhashSparseSymmetricMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    if (!hashmap_.contains(i)) {
        hashmap_.emplace(i, Map());
    }
    if (!hashmap_[i].contains(j)) {
        hashmap_[i].emplace(j, value);
    } else {
        hashmap_[i][j] += value;
    }
    if (i != j) {
        if (!hashmap_.contains(j)) {
            hashmap_.emplace(j, Map());
        }
        if (!hashmap_[j].contains(i)) {
            hashmap_[j].emplace(i, value);
        } else {
            hashmap_[j][i] += value;
        }
    }
}

uint32_t EmhashSparseSymmetricMatrix::size() const {
    return size_;
}

double EmhashSparseSymmetricMatrix::at(uint32_t i, uint32_t j) const noexcept {
    if (!hashmap_.contains(i) || !hashmap_.at(i).contains(j)) {
        return 0;
    }
    return hashmap_.at(i).at(j);
}

void EmhashSparseSymmetricMatrix::resize(size_t size) {
    size_ = size;
}

const emhash8::HashMap<uint32_t, EmhashSparseSymmetricMatrix::Map>&
EmhashSparseSymmetricMatrix::getSparseSymmetricMatrix() const {
    return hashmap_;
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