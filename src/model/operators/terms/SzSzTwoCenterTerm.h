#ifndef SPINNER_SZSZTWOCENTERTERM_H
#define SPINNER_SZSZTWOCENTERTERM_H

#include "Term.h"

#include "src/common/index_converter/lexicographic/IndexConverter.h"
#include "src/model/NumericalParameters.h"

namespace model::operators {
class SzSzTwoCenterTerm: public TwoCenterTerm {
  public:
    SzSzTwoCenterTerm(
        std::shared_ptr<const index_converter::lexicographic::IndexConverter> converter,
        std::shared_ptr<const TwoDNumericalParameters<double>> parameters,
        double prefactor);
    SzSzTwoCenterTerm(
        std::shared_ptr<const index_converter::lexicographic::IndexConverter> converter,
        std::shared_ptr<const TwoDNumericalParameters<double>> parameters);
    std::unique_ptr<Term> clone() const override;
    void construct(
        quantum::linear_algebra::AbstractSymmetricMatrix&
            matrix_in_lexicografical_basis,
        const std::set<unsigned int>& indexes_of_vectors,
        uint32_t center_a,
        uint32_t center_b) const override;

  private:
    std::shared_ptr<const index_converter::lexicographic::IndexConverter> converter_;
    std::shared_ptr<const TwoDNumericalParameters<double>> coefficients_;
    double prefactor_ = 1;
};
}  // namespace model::operators
#endif  //SPINNER_SZSZTWOCENTERTERM_H
