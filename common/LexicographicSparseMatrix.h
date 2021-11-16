#ifndef JULY_LEXICOGRAPHICSPARSEMATRIX_H
#define JULY_LEXICOGRAPHICSPARSEMATRIX_H

#include "LexicographicIndexConverter.h"

#include <cstdint>
#include <memory>
#include <ostream>


class LexicographicSparseMatrix {
public:
    explicit LexicographicSparseMatrix(spaces::LexicographicIndexConverter converter);
    LexicographicSparseMatrix(const LexicographicSparseMatrix&) = delete;
    LexicographicSparseMatrix& operator=(const LexicographicSparseMatrix&) = delete;
    LexicographicSparseMatrix(LexicographicSparseMatrix&&) noexcept;
    LexicographicSparseMatrix& operator=(LexicographicSparseMatrix&&) noexcept;
    ~LexicographicSparseMatrix();

    [[nodiscard]] uint32_t size() const;

    void add_to_position(double value, uint32_t i, uint32_t j);
    double operator()(uint32_t i, uint32_t j) const;

    void resize(uint32_t new_size);

    void add_scalar_product(uint32_t index_of_vector, uint32_t center_a, uint32_t center_b, double factor);

    friend std::ostream &operator<<(std::ostream &os, const LexicographicSparseMatrix &data);

private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
    spaces::LexicographicIndexConverter converter_;

    void add_scalar_product_nondiagonal_part(uint32_t index_of_vector, uint32_t plus_center, uint32_t minus_center,
                                             uint32_t projection_of_plus_center, uint32_t projection_of_minus_center,
                                             double factor);
};


#endif //JULY_LEXICOGRAPHICSPARSEMATRIX_H
