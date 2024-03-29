#ifndef SPINNER_LOCALSSQUAREDONECENTERTERM_H
#define SPINNER_LOCALSSQUAREDONECENTERTERM_H

#include "Term.h"

namespace model::operators {
class LocalSSquaredOneCenterTerm: public OneCenterTerm {
  public:
    LocalSSquaredOneCenterTerm(
        lexicographic::IndexConverter converter,
        std::shared_ptr<const OneDNumericalParameters<double>> coefficients);
    LocalSSquaredOneCenterTerm(
        lexicographic::IndexConverter converter,
        std::shared_ptr<const OneDNumericalParameters<double>> coefficients,
        double prefactor);
    std::unique_ptr<OneCenterTerm> clone() const override;
    void construct(
        std::unique_ptr<quantum::linear_algebra::AbstractSymmetricMatrix>&
            matrix_in_lexicografical_basis,
        uint32_t index_of_vector,
        uint32_t center_a) const override;

  private:
    const lexicographic::IndexConverter converter_;
    std::shared_ptr<const OneDNumericalParameters<double>> coefficients_;
    const double prefactor_ = 1;
};

}  // namespace model::operators

#endif  //SPINNER_LOCALSSQUAREDONECENTERTERM_H
