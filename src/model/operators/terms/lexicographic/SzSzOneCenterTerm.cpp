#include "SzSzOneCenterTerm.h"

#include <cmath>
#include <utility>

namespace model::operators::lexicographic {
std::unique_ptr<Term> SzSzOneCenterTerm::clone() const {
    return std::make_unique<SzSzOneCenterTerm>(converter_, coefficients_);
}

void SzSzOneCenterTerm::construct(
    quantum::linear_algebra::AbstractSymmetricMatrix&
        matrix_in_lexicografical_basis,
    const std::set<unsigned int>& indexes_of_vectors,
    uint32_t center_a) const {
    double factor = coefficients_->at(center_a);

    if (!std::isnan(factor)) {
        for (const auto& index_of_vector : indexes_of_vectors) {
            uint32_t projection_of_center_a =
                converter_->convert_lex_index_to_one_sz_projection(index_of_vector, center_a);

            // Saz Saz
            double diagonal_value = (projection_of_center_a - converter_->get_spins()[center_a])
                * (projection_of_center_a - converter_->get_spins()[center_a]) * factor;
            matrix_in_lexicografical_basis.add_to_position(
                diagonal_value,
                index_of_vector,
                index_of_vector);
        }
    }
}

SzSzOneCenterTerm::SzSzOneCenterTerm(
    std::shared_ptr<const index_converter::lexicographic::IndexConverter> converter,
    std::shared_ptr<const OneDNumericalParameters<double>> coefficients) :
    OneCenterTerm(converter->get_mults().size()),
    converter_(std::move(converter)),
    coefficients_(std::move(coefficients)) {}
}  // namespace model::operators::lexicographic