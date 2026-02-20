#ifndef SPINNER_ARMASPARSESEMIUNITARYMATRIX_H
#define SPINNER_ARMASPARSESEMIUNITARYMATRIX_H

#include <armadillo>

#include "src/linalg_structures/AbstractSparseSemiunitaryMatrix.h"

namespace spinner::linalg_structures::armadillo {
class ArmaSparseSemiunitaryMatrix: public AbstractSparseSemiunitaryMatrix {
  public:
    std::unique_ptr<Iterator> GetNewIterator(size_t col) const override;
    uint32_t size_cols() const override;
    uint32_t size_rows() const override;
    bool empty() const override;
    bool vempty(uint32_t col) const override;
    void clear() override;
    void eraseExplicitZeros() override;
    bool is_zero(uint32_t col, uint32_t row) const override;
    void move_vector_from(
        uint32_t col,
        std::unique_ptr<AbstractSparseSemiunitaryMatrix>& subspace_from) override;
    void resize(uint32_t cols, uint32_t rows);
    void add_to_position(double value, uint32_t col, uint32_t row) override;
    double at(uint32_t col, uint32_t row) const override;
    void normalize() override;
    void print(std::ostream& os) const override;

    const arma::sp_mat& getSparseSemiunitaryMatrix() const;
    void unitaryTransform(
        const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrixToTransform,
        std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrixToAdd) const override;

  private:
    arma::sp_mat sparseSemiunitaryMatrix_;
};
} // namespace spinner::linalg_structures::armadillo
#endif  //SPINNER_ARMASPARSESEMIUNITARYMATRIX_H
