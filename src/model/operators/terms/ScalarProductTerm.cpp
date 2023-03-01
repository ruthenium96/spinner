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
    std::unique_ptr<quantum::linear_algebra::AbstractSymmetricMatrix>&
        matrix_in_lexicografical_basis,
    uint32_t index_of_vector,
    uint32_t center_a,
    uint32_t center_b) const {
    if (!std::isnan(coefficients_->at(center_a, center_b))) {
        double factor = -2 * coefficients_->at(center_a, center_b);
        add_scalar_product(
            matrix_in_lexicografical_basis,
            index_of_vector,
            center_a,
            center_b,
            factor);
    }
}

std::unique_ptr<TwoCenterTerm> ScalarProductTerm::clone() const {
    return std::make_unique<ScalarProductTerm>(converter_, coefficients_);
}

void ScalarProductTerm::add_scalar_product(
    std::unique_ptr<quantum::linear_algebra::AbstractSymmetricMatrix>& matrix,
    uint32_t index_of_vector,
    uint32_t center_a,
    uint32_t center_b,
    double factor) const {
    uint32_t projection_of_center_a =
        converter_.convert_lex_index_to_one_sz_projection(index_of_vector, center_a);
    uint32_t projection_of_center_b =
        converter_.convert_lex_index_to_one_sz_projection(index_of_vector, center_b);

    // (Sa, Sb) = Sax*Sbx + Say*Sby + Saz*Sbz = 0.5 * (Sa+Sb- + Sa-Sb+) + Saz*Sbz

    // Saz Sbz
    double diagonal_value = (projection_of_center_a - converter_.get_spins()[center_a])
        * (projection_of_center_b - converter_.get_spins()[center_b]) * factor;
    matrix->add_to_position(diagonal_value, index_of_vector, index_of_vector);

    // Sa+ Sb-
    add_scalar_product_nondiagonal_part(
        matrix,
        index_of_vector,
        center_a,
        center_b,
        projection_of_center_a,
        projection_of_center_b,
        factor);

    // Sa- Sb+
    add_scalar_product_nondiagonal_part(
        matrix,
        index_of_vector,
        center_b,
        center_a,
        projection_of_center_b,
        projection_of_center_a,
        factor);
}

void ScalarProductTerm::add_scalar_product_nondiagonal_part(
    std::unique_ptr<quantum::linear_algebra::AbstractSymmetricMatrix>& matrix,
    uint32_t index_of_vector,
    uint32_t plus_center,
    uint32_t minus_center,
    uint32_t projection_of_plus_center,
    uint32_t projection_of_minus_center,
    double factor) const {
    if (projection_of_plus_center == converter_.get_mults()[plus_center] - 1
        || projection_of_minus_center == 0) {
        return;
    }
    uint32_t index_of_new_vector = converter_.ladder_projection(
        converter_.ladder_projection(index_of_vector, plus_center, +1),
        minus_center,
        -1);
    if (index_of_new_vector < index_of_vector) {
        // For index_of_vector we add to (index_of_vector, index_of_new_vector).
        // If there will be no breaking return, lately index_of_new_vector
        // will add to (index_of_new_vector, index_of_vector), and it will lead to double adding
        // because of the properties of assign_to_position (it adds to both (i,j) and (j,i)).
        return;
    }

    // projection m = number n - spin S
    // so S(S+1)-m(m+1) = (2S-n)(n+1)
    // so S(S+1)-m(m-1) = n(2S+1-n)
    double factor_a = (2 * converter_.get_spins()[plus_center] - projection_of_plus_center)
        * (projection_of_plus_center + 1);
    double factor_b = projection_of_minus_center
        * (2 * converter_.get_spins()[minus_center] + 1 - projection_of_minus_center);

    // TODO: I do not know, why we need here 0.25 instead of 0.5.
    //  I guess, there is some double-addition.
    matrix->add_to_position(
        0.5 * sqrt(factor_a * factor_b) * factor,
        index_of_vector,
        index_of_new_vector);
}
}  // namespace model::operators