#ifndef SPINNER_STDSPARSEMATRIX_H
#define SPINNER_STDSPARSEMATRIX_H

#include <cstdint>
#include <map>
#include <memory>
#include <ostream>
#include <vector>

#include "src/entities/data_structures/AbstractSparseMatrix.h"

namespace quantum::linear_algebra {
class StdSparseMatrix: public AbstractSparseMatrix {
  public:
    StdSparseMatrix() = default;
    StdSparseMatrix(const StdSparseMatrix&) = delete;
    StdSparseMatrix& operator=(const StdSparseMatrix&) = delete;
    StdSparseMatrix(StdSparseMatrix&&) noexcept = default;
    StdSparseMatrix& operator=(StdSparseMatrix&&) noexcept = default;
    ~StdSparseMatrix() = default;

    std::unique_ptr<Iterator> GetNewIterator(size_t index_of_vector) const override;

    uint32_t size() const override;
    bool empty() const override;
    bool vempty(uint32_t index_of_vector) const override;
    void clear() override;

    void erase_if_zero() override;

    bool is_zero(uint32_t i, uint32_t j) const override;
    void
    move_vector_from(uint32_t i, std::unique_ptr<AbstractSparseMatrix>& subspace_from) override;
    void move_all_from(std::unique_ptr<AbstractSparseMatrix>& subspace_from) override;
    void copy_vector_from(uint32_t i, const std::unique_ptr<AbstractSparseMatrix>& subspace_from)
        override;
    void copy_all_from(const std::unique_ptr<AbstractSparseMatrix>& subspace_from) override;
    void resize(uint32_t new_size) override;

    void add_to_position(double value, uint32_t i, uint32_t j) override;
    double at(uint32_t i, uint32_t j) const override;

    void normalize() override;

    bool
    is_equal_up_to_vector_order(const std::unique_ptr<AbstractSparseMatrix>& rhs) const override;
    void print(std::ostream& os) const override;

  private:
    std::vector<std::map<uint32_t, double>> basis_;
    // c-like pointers are necessary to avoid double-free error
    static const StdSparseMatrix* downcast_ptr(const std::unique_ptr<AbstractSparseMatrix>& ptr);
    static StdSparseMatrix* downcast_ptr(std::unique_ptr<AbstractSparseMatrix>& ptr);
};
}  // namespace quantum::linear_algebra

#endif  //SPINNER_STDSPARSEMATRIX_H
