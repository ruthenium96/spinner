#ifndef JULY_RAWSUBMATRIXDATA_H
#define JULY_RAWSUBMATRIXDATA_H

#include <memory>

class RawSubmatrixData {
public:
    void add_to_position(double value, uint32_t i, uint32_t j);
    void resize(uint32_t matrix_in_space_basis_size_i, uint32_t matrix_in_space_basis_size_j);

    RawSubmatrixData();
    RawSubmatrixData(const RawSubmatrixData&) = delete;
    RawSubmatrixData& operator=(const RawSubmatrixData&) = delete;
    RawSubmatrixData(RawSubmatrixData&&) noexcept;
    RawSubmatrixData& operator=(RawSubmatrixData&&) noexcept;
    ~RawSubmatrixData();

    friend std::ostream &operator<<(std::ostream &os, const RawSubmatrixData &raw_data);

private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
};


#endif //JULY_RAWSUBMATRIXDATA_H
