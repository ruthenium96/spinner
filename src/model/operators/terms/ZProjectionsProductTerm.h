#ifndef SPINNER_ZPROJECTIONSPRODUCTTERM_H
#define SPINNER_ZPROJECTIONSPRODUCTTERM_H

#include "Term.h"
namespace model::operators {
class ZProjectionsProductTerm: public TwoCenterTerm {
  public:
    ZProjectionsProductTerm(
        lexicographic::IndexConverter converter,
        std::shared_ptr<const TwoDNumericalParameters<double>> parameters);
    std::unique_ptr<TwoCenterTerm> clone() const override;
    void construct(
        UnitarySparseMatrix& matrix_in_lexicografical_basis,
        uint32_t index_of_vector,
        uint32_t center_a,
        uint32_t center_b) const override;
    std::shared_ptr<const TwoDNumericalParameters<double>> get_parameters() const override;

  private:
    const lexicographic::IndexConverter converter_;
    std::shared_ptr<const TwoDNumericalParameters<double>> coefficients_;
};
}  // namespace model::operators
#endif  //SPINNER_ZPROJECTIONSPRODUCTTERM_H
