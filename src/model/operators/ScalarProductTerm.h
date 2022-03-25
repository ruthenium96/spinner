#ifndef JULY_SCALARPRODUCTTERM_H
#define JULY_SCALARPRODUCTTERM_H

#include "Term.h"

namespace model::operators {
class ScalarProductTerm: public TwoCenterTerm {
  public:
    explicit ScalarProductTerm(size_t number_of_spins);
    explicit ScalarProductTerm(std::shared_ptr<const DenseMatrix> isotropic_exchange_parameters);

    std::unique_ptr<TwoCenterTerm> clone() const override;

    void construct(
        lexicographic::SparseMatrix& matrix_in_lexicografical_basis,
        uint32_t index_of_vector,
        uint32_t center_a,
        uint32_t center_b) const override;

    std::shared_ptr<const DenseMatrix> get_parameters() const override;

  private:
    std::shared_ptr<const DenseMatrix> coefficients;
};
}  // namespace model::operators

#endif  //JULY_SCALARPRODUCTTERM_H
