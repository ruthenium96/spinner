#ifndef SPINNER_EMHASHSPARSESEMIUNITARYMATRIX_H
#define SPINNER_EMHASHSPARSESEMIUNITARYMATRIX_H

#include <cstdint>
#include <hash_table8.hpp>
#include <map>
#include <memory>
#include <ostream>
#include <vector>

#include "src/entities/data_structures/AbstractSparseSemiunitaryMatrix.h"

namespace quantum::linear_algebra {
class EmhashSparseSemiunitaryMatrix: public AbstractSparseSemiunitaryMatrix {
  public:
    using Map = emhash8::HashMap<uint32_t, double>;

    EmhashSparseSemiunitaryMatrix() = default;
    EmhashSparseSemiunitaryMatrix(const EmhashSparseSemiunitaryMatrix&) = delete;
    EmhashSparseSemiunitaryMatrix& operator=(const EmhashSparseSemiunitaryMatrix&) = delete;
    EmhashSparseSemiunitaryMatrix(EmhashSparseSemiunitaryMatrix&&) noexcept = default;
    EmhashSparseSemiunitaryMatrix& operator=(EmhashSparseSemiunitaryMatrix&&) noexcept = default;
    ~EmhashSparseSemiunitaryMatrix() override = default;

    std::unique_ptr<Iterator> GetNewIterator(size_t index_of_vector) const override;

    uint32_t size_rows() const override;
    uint32_t size_cols() const override;
    bool empty() const override;
    bool vempty(uint32_t index_of_vector) const override;
    void clear() override;

    void eraseExplicitZeros() override;

    bool is_zero(uint32_t i, uint32_t j) const override;
    void move_vector_from(
        uint32_t i,
        std::unique_ptr<AbstractSparseSemiunitaryMatrix>& subspace_from) override;
    void resize(uint32_t cols, uint32_t rows);

    void add_to_position(double value, uint32_t i, uint32_t j) override;
    double at(uint32_t i, uint32_t j) const override;

    void normalize() override;

    void print(std::ostream& os) const override;
    void unitaryTransform(
        const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrixToTransform,
        std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrixToAdd) const override;

    const std::vector<Map>& getSparseSemiunitaryMatrix() const;

  private:
    std::vector<Map> basis_;
    uint32_t rows_;
    // c-like pointers are necessary to avoid double-free error
    static const EmhashSparseSemiunitaryMatrix*
    downcast_ptr(const std::unique_ptr<AbstractSparseSemiunitaryMatrix>& ptr);
    static EmhashSparseSemiunitaryMatrix*
    downcast_ptr(std::unique_ptr<AbstractSparseSemiunitaryMatrix>& ptr);
};
}  // namespace quantum::linear_algebra

#endif  //SPINNER_EMHASHSPARSESEMIUNITARYMATRIX_H
