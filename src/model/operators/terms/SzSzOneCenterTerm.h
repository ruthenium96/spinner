#ifndef SPINNER_SZSZONECENTERTERM_H
#define SPINNER_SZSZONECENTERTERM_H
#include "Term.h"

#include "src/common/index_converter/lexicographic/IndexConverter.h"
#include "src/model/NumericalParameters.h"

namespace model::operators {

class SzSzOneCenterTerm: public OneCenterTerm {
  public:
    SzSzOneCenterTerm(
        std::shared_ptr<const index_converter::lexicographic::IndexConverter> converter,
        std::shared_ptr<const OneDNumericalParameters<double>> coefficients);
    std::unique_ptr<Term> clone() const override;
    void construct(
        quantum::linear_algebra::AbstractSymmetricMatrix&
            matrix_in_lexicografical_basis,
        const std::set<unsigned int>& indexes_of_vectors,
        uint32_t center_a) const override;

  private:
    std::shared_ptr<const index_converter::lexicographic::IndexConverter> converter_;
    std::shared_ptr<const OneDNumericalParameters<double>> coefficients_;
};
}  // namespace model::operators
#endif  //SPINNER_SZSZONECENTERTERM_H
