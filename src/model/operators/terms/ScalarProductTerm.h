#ifndef SPINNER_SCALARPRODUCTTERM_H
#define SPINNER_SCALARPRODUCTTERM_H

#include "Term.h"

#include "src/common/index_converter/lexicographic/IndexConverter.h"
#include "src/model/NumericalParameters.h"

namespace model::operators {
class ScalarProductTerm: public TwoCenterTerm {
  public:
    explicit ScalarProductTerm(std::shared_ptr<const index_converter::lexicographic::IndexConverter> converter);
    ScalarProductTerm(
        std::shared_ptr<const index_converter::lexicographic::IndexConverter> converter,
        std::shared_ptr<const TwoDNumericalParameters<double>> isotropic_exchange_parameters);

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
    void add_scalar_product(
        quantum::linear_algebra::AbstractSymmetricMatrix& matrix,
        uint32_t index_of_vector,
        uint32_t center_a,
        uint32_t center_b,
        double factor) const;
    void add_scalar_product_nondiagonal_part(
        quantum::linear_algebra::AbstractSymmetricMatrix& matrix,
        uint32_t index_of_vector,
        uint32_t plus_center,
        uint32_t minus_center,
        uint32_t projection_of_plus_center,
        uint32_t projection_of_minus_center,
        double factor) const;
};
}  // namespace model::operators

#endif  //SPINNER_SCALARPRODUCTTERM_H
