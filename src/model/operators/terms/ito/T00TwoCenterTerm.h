#ifndef SPINNER_ITOSCALARPRODUCTTERM_H
#define SPINNER_ITOSCALARPRODUCTTERM_H

#include "src/model/operators/terms/Term.h"
#include "src/common/index_converter/s_squared/IndexConverter.h"
#include "WignerEckartHelper.h"
#include "src/model/NumericalParameters.h"

namespace model::operators::ito {

class T00TwoCenterTerm : public TwoCenterTerm {
  public:
    T00TwoCenterTerm(
        std::shared_ptr<const index_converter::s_squared::IndexConverter> converter,
        std::shared_ptr<const TwoDNumericalParameters<double>> isotropic_exchange_parameters,
        double prefactor = 1);

    std::unique_ptr<Term> clone() const override;

    void construct(
        quantum::linear_algebra::AbstractSymmetricMatrix&
            matrix,
        const std::set<unsigned int>& indexes_of_vectors,
        uint32_t center_a,
        uint32_t center_b) const override;

  private:
    std::shared_ptr<const index_converter::s_squared::IndexConverter> converter_;
    std::shared_ptr<const TwoDNumericalParameters<double>> coefficients_;
    const double prefactor_;
    WignerEckartHelper wigner_eckart_helper_;

    void add_scalar_product(
    quantum::linear_algebra::AbstractSymmetricMatrix& matrix,
    const std::set<unsigned int>& indexes_of_vectors,
    uint32_t center_a,
    uint32_t center_b,
    double factor) const;

    std::vector<uint8_t> constructRanksOfTZero(uint32_t center_a, uint32_t center_b) const;
};

}  // namespace model::operators::ito

#endif  //SPINNER_SCALARPRODUCTTERM_H
