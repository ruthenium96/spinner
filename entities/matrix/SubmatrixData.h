#ifndef JULY_SUBMATRIXDATA_H
#define JULY_SUBMATRIXDATA_H

#include <iostream>
#include <memory>

class SubmatrixData {
public:
    void add_to_position(double value, uint32_t i, uint32_t j);
    void resize(uint32_t matrix_in_space_basis_size_i, uint32_t matrix_in_space_basis_size_j);
    // TODO: is it possible to implement diagonalize() in the other place?
    void diagonalize();

    SubmatrixData();
    SubmatrixData(const SubmatrixData&) = delete;
    SubmatrixData& operator=(const SubmatrixData&) = delete;
    SubmatrixData(SubmatrixData&&) noexcept;
    SubmatrixData& operator=(SubmatrixData&&) noexcept;
    ~SubmatrixData();

    friend std::ostream &operator<<(std::ostream &os, const SubmatrixData &raw_data);

private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
};

#endif //JULY_SUBMATRIXDATA_H
