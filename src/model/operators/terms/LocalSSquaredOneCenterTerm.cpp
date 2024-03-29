#include "LocalSSquaredOneCenterTerm.h"

#include <cmath>
#include <utility>

namespace model::operators {

std::unique_ptr<OneCenterTerm> LocalSSquaredOneCenterTerm::clone() const {
    return std::make_unique<LocalSSquaredOneCenterTerm>(converter_, coefficients_, prefactor_);
}
void LocalSSquaredOneCenterTerm::construct(
    std::unique_ptr<quantum::linear_algebra::AbstractSymmetricMatrix>&
        matrix_in_lexicografical_basis,
    uint32_t index_of_vector,
    uint32_t center_a) const {
    double factor = coefficients_->at(center_a) * prefactor_;

    if (!std::isnan(factor)) {
        double spin = converter_.get_spins()[center_a];

        double diagonal_value = factor * spin * (spin + 1);
        matrix_in_lexicografical_basis->add_to_position(
            diagonal_value,
            index_of_vector,
            index_of_vector);
    }
}

LocalSSquaredOneCenterTerm::LocalSSquaredOneCenterTerm(
    lexicographic::IndexConverter converter,
    std::shared_ptr<const OneDNumericalParameters<double>> coefficients,
    double prefactor) :
    converter_(std::move(converter)),
    coefficients_(std::move(coefficients)),
    prefactor_(prefactor) {}
LocalSSquaredOneCenterTerm::LocalSSquaredOneCenterTerm(
    lexicographic::IndexConverter converter,
    std::shared_ptr<const OneDNumericalParameters<double>> coefficients) :
    converter_(std::move(converter)),
    coefficients_(std::move(coefficients)) {}
}  // namespace model::operators