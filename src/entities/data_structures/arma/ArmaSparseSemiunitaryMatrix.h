#ifndef SPINNER_ARMASPARSESEMIUNITARYMATRIX_H
#define SPINNER_ARMASPARSESEMIUNITARYMATRIX_H

#include <armadillo>

#include "src/entities/data_structures/AbstractSparseSemiunitaryMatrix.h"

namespace quantum::linear_algebra {
class ArmaSparseSemiunitaryMatrix: public AbstractSparseSemiunitaryMatrix {
  public:
    std::unique_ptr<Iterator> GetNewIterator(size_t index_of_vector) const override;
    uint32_t size_cols() const override;
    uint32_t size_rows() const override;
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

  private:
    arma::sp_mat sparseSemiunitaryMatrix_;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_ARMASPARSESEMIUNITARYMATRIX_H
