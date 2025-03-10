#ifndef SPINNER_LOCALSSQUAREDONECENTERTERM_H
#define SPINNER_LOCALSSQUAREDONECENTERTERM_H

#include "Term.h"

#include "src/common/index_converter/AbstractIndexConverter.h"
#include "src/model/NumericalParameters.h"

namespace model::operators {
class LocalSSquaredOneCenterTerm: public OneCenterTerm {
  public:
    LocalSSquaredOneCenterTerm(
        std::shared_ptr<const index_converter::AbstractIndexConverter> converter,
        std::shared_ptr<const OneDNumericalParameters<double>> coefficients);
    LocalSSquaredOneCenterTerm(
        std::shared_ptr<const index_converter::AbstractIndexConverter> converter,
        std::shared_ptr<const OneDNumericalParameters<double>> coefficients,
        double prefactor);
    std::unique_ptr<Term> clone() const override;
    void construct(
        quantum::linear_algebra::AbstractSymmetricMatrix&
            matrix_in_lexicografical_basis,
        const std::set<unsigned int>& indexes_of_vectors,
        uint32_t center_a) const override;

  private:
    std::shared_ptr<const index_converter::AbstractIndexConverter> converter_;
    std::shared_ptr<const OneDNumericalParameters<double>> coefficients_;
    const double prefactor_ = 1;
};

}  // namespace model::operators

#endif  //SPINNER_LOCALSSQUAREDONECENTERTERM_H
