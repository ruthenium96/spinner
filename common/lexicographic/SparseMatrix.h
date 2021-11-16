#ifndef JULY_SPARSEMATRIX_H
#define JULY_SPARSEMATRIX_H

#include "IndexConverter.h"

#include <cstdint>
#include <memory>
#include <ostream>

namespace lexicographic {

class SparseMatrix {
public:
    explicit SparseMatrix(lexicographic::IndexConverter converter);
    SparseMatrix(const SparseMatrix&) = delete;
    SparseMatrix& operator=(const SparseMatrix&) = delete;
    SparseMatrix(SparseMatrix&&) noexcept;
    SparseMatrix& operator=(SparseMatrix&&) noexcept;
    ~SparseMatrix();

    [[nodiscard]] uint32_t size() const;

    void add_to_position(double value, uint32_t i, uint32_t j);
    double operator()(uint32_t i, uint32_t j) const;

    void resize(uint32_t new_size);

    void add_scalar_product(uint32_t index_of_vector, uint32_t center_a, uint32_t center_b, double factor);

    friend std::ostream &operator<<(std::ostream &os, const SparseMatrix &data);

private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
    lexicographic::IndexConverter converter_;

    void add_scalar_product_nondiagonal_part(uint32_t index_of_vector, uint32_t plus_center, uint32_t minus_center,
                                             uint32_t projection_of_plus_center, uint32_t projection_of_minus_center,
                                             double factor);
};
}

#endif //JULY_SPARSEMATRIX_H
