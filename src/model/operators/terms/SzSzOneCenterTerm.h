#ifndef SPINNER_SZSZONECENTERTERM_H
#define SPINNER_SZSZONECENTERTERM_H
#include "Term.h"

namespace model::operators {

class SzSzOneCenterTerm: public OneCenterTerm {
  public:
    SzSzOneCenterTerm(
        lexicographic::IndexConverter converter,
        std::shared_ptr<const OneDNumericalParameters<double>> coefficients);
    std::unique_ptr<OneCenterTerm> clone() const override;
    void construct(
        quantum::linear_algebra::AbstractSymmetricMatrix&
            matrix_in_lexicografical_basis,
        uint32_t index_of_vector,
        uint32_t center_a) const override;

  private:
    const lexicographic::IndexConverter converter_;
    std::shared_ptr<const OneDNumericalParameters<double>> coefficients_;
};
}  // namespace model::operators
#endif  //SPINNER_SZSZONECENTERTERM_H
