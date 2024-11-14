#ifndef SPINNER_SCALARPRODUCTITOTERM_H
#define SPINNER_SCALARPRODUCTITOTERM_H

#include "Term.h"
#include "src/spin_algebra/SSquaredConverter.h"

namespace model::operators {

class ScalarProductITOTerm : public TwoCenterTerm {
  public:
    ScalarProductITOTerm(
        std::shared_ptr<const spin_algebra::SSquaredConverter> ssquared_converter,
        std::shared_ptr<const TwoDNumericalParameters<double>> isotropic_exchange_parameters);

    std::unique_ptr<Term> clone() const override;

    void construct(
        quantum::linear_algebra::AbstractSymmetricMatrix&
            matrix,
        uint32_t index_of_vector,
        uint32_t center_a,
        uint32_t center_b) const override;

  private:
    std::shared_ptr<const spin_algebra::SSquaredConverter> ssquared_converter_;
    std::shared_ptr<const TwoDNumericalParameters<double>> coefficients_;

    void add_scalar_product(
    quantum::linear_algebra::AbstractSymmetricMatrix& matrix,
    uint32_t index_of_vector,
    uint32_t center_a,
    uint32_t center_b,
    double factor) const;
};

}  // namespace model::operators

#endif  //SPINNER_SCALARPRODUCTITOTERM_H
