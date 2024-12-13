#ifndef SPINNER_ITOSCALARPRODUCTTERM_H
#define SPINNER_ITOSCALARPRODUCTTERM_H

#include "src/model/operators/terms/Term.h"
#include "src/common/index_converter/s_squared/IndexConverter.h"
#include "src/spin_algebra/ClebshGordanCalculator.h"
#include "src/model/NumericalParameters.h"

namespace model::operators::ito {

class ScalarProductTerm : public TwoCenterTerm {
  public:
    ScalarProductTerm(
        std::shared_ptr<const index_converter::s_squared::IndexConverter> converter,
        std::shared_ptr<const TwoDNumericalParameters<double>> isotropic_exchange_parameters);

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
    spin_algebra::ClebshGordanCalculator clebshGordanCalculator_;

    void add_scalar_product(
    quantum::linear_algebra::AbstractSymmetricMatrix& matrix,
    const std::set<unsigned int>& indexes_of_vectors,
    uint32_t center_a,
    uint32_t center_b,
    double factor) const;

    std::vector<uint8_t> constructRanksOfTZero(uint32_t center_a, uint32_t center_b) const;

    double total_9j_coefficient(
        const index_converter::s_squared::Level& left,
        const index_converter::s_squared::Level& right,
        const std::vector<uint8_t>& ranks,
        double local_prod) const;

    void construct_overlapping_levels(const index_converter::s_squared::Level& level, 
        const std::vector<uint8_t>& ranks,
        std::vector<index_converter::s_squared::Level>& answer) const;
};

}  // namespace model::operators::ito

#endif  //SPINNER_SCALARPRODUCTTERM_H
