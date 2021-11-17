#ifndef JULY_SCALARPRODUCT_H
#define JULY_SCALARPRODUCT_H

#include "components/operator/Interaction.h"

class ScalarProduct: public TwoCenterTerm {
  public:
    explicit ScalarProduct(size_t number_of_spins);
    explicit ScalarProduct(std::shared_ptr<const DenseMatrix> isotropic_exchange_parameters);

    void construct(
        lexicographic::SparseMatrix& matrix_in_lexicografical_basis,
        uint32_t index_of_vector,
        uint32_t center_a,
        uint32_t center_b) const override;

    [[nodiscard]] std::shared_ptr<const DenseMatrix> get_parameters() const override;

  private:
    std::shared_ptr<const DenseMatrix> coefficients;
    // TODO: abstract from arma::dmat
    //arma::dmat coefficients;
};

#endif  //JULY_SCALARPRODUCT_H
