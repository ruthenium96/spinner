#include "ScalarProductTerm.h"

#include <cmath>
#include <utility>

namespace model::operators {
ScalarProductTerm::ScalarProductTerm(lexicographic::IndexConverter converter) :
    converter_(std::move(converter)) {
    size_t number_of_spins = converter_.get_mults().size();
    auto mutable_coefficients =
        std::make_shared<TwoDNumericalParameters<double>>(number_of_spins, NAN);
    for (size_t i = 0; i < number_of_spins; ++i) {
        for (size_t j = i + 1; j < number_of_spins; ++j) {
            // TODO: it does not looks good
            mutable_coefficients->at(i, j) = -1;
            mutable_coefficients->at(j, i) = -1;
        }
    }
    coefficients_ = mutable_coefficients;
}

ScalarProductTerm::ScalarProductTerm(
    lexicographic::IndexConverter converter,
    std::shared_ptr<const TwoDNumericalParameters<double>> parameters) :
    converter_(std::move(converter)),
    coefficients_(std::move(parameters)) {}

void ScalarProductTerm::construct(
    lexicographic::SparseMatrix& matrix_in_lexicografical_basis,
    uint32_t index_of_vector,
    uint32_t center_a,
    uint32_t center_b) const {
    if (!std::isnan(coefficients_->at(center_a, center_b))) {
        double factor = -2 * coefficients_->at(center_a, center_b);
        matrix_in_lexicografical_basis
            .add_scalar_product(index_of_vector, center_a, center_b, factor);
    }
}

std::shared_ptr<const TwoDNumericalParameters<double>> ScalarProductTerm::get_parameters() const {
    return coefficients_;
}

std::unique_ptr<TwoCenterTerm> ScalarProductTerm::clone() const {
    return std::make_unique<ScalarProductTerm>(coefficients_);
}
}  // namespace model::operators