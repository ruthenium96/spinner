#ifndef SPINNER_T20ONECENTERTERM_H
#define SPINNER_T20ONECENTERTERM_H

#include "src/model/operators/terms/Term.h"

#include "src/common/index_converter/s_squared/IndexConverter.h"
#include "src/model/NumericalParameters.h"
#include "WignerEckartHelper.h"

namespace model::operators::ito {
class T20OneCenterTerm : public OneCenterTerm {
  public:
    T20OneCenterTerm(
        std::shared_ptr<const index_converter::s_squared::IndexConverter> converter,
        std::shared_ptr<const OneDNumericalParameters<double>> coefficients,
        double prefactor = 1);
    std::unique_ptr<Term> clone() const override;
    void construct(
        quantum::linear_algebra::AbstractSymmetricMatrix& matrix,
        const std::set<unsigned int>& indexes_of_vectors,
        uint32_t center_a) const override;

  private:
    std::shared_ptr<const index_converter::s_squared::IndexConverter> converter_;
    std::shared_ptr<const OneDNumericalParameters<double>> coefficients_;
    const double prefactor_;
    WignerEckartHelper wigner_eckart_helper_;

    void add_ttwo_term(
        quantum::linear_algebra::AbstractSymmetricMatrix& matrix,
        const std::set<unsigned int>& indexes_of_vectors,
        uint32_t center_a,
        double factor) const;
    std::vector<uint8_t> constructRanksOfTTwo(uint32_t center_a) const;
};

} // namespace model::operators::ito

#endif // SPINNER_T20ONECENTERTERM_H