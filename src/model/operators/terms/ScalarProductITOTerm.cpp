#include "ScalarProductITOTerm.h"

#include <cmath>
#include <utility>

#include "src/common/Logger.h"

namespace model::operators {

ScalarProductITOTerm::ScalarProductITOTerm(
    std::shared_ptr<const spin_algebra::SSquaredConverter> ssquared_converter,
    std::shared_ptr<const TwoDNumericalParameters<double>> isotropic_exchange_parameters) :
    TwoCenterTerm(isotropic_exchange_parameters->size()),
    ssquared_converter_(std::move(ssquared_converter)),
    coefficients_(std::move(isotropic_exchange_parameters)) {}

void ScalarProductITOTerm::construct(
    quantum::linear_algebra::AbstractSymmetricMatrix& matrix,
    uint32_t index_of_vector,
    uint32_t center_a,
    uint32_t center_b) const {
    if (!std::isnan(coefficients_->at(center_a, center_b))) {
        double factor = 2 * sqrt(3) * coefficients_->at(center_a, center_b);
        add_scalar_product(
            matrix,
            index_of_vector,
            center_a,
            center_b,
            factor);
    }
//    matrix.print(std::cout);
}

void ScalarProductITOTerm::add_scalar_product(
    quantum::linear_algebra::AbstractSymmetricMatrix& matrix,
    uint32_t index_of_vector,
    uint32_t center_a,
    uint32_t center_b,
    double factor) const {
    const auto& ssquared_converter = *ssquared_converter_;

    auto ranks = ssquared_converter.constructRanksOfTZero(center_a, center_b);

    auto row = ssquared_converter.number_in_block(index_of_vector);
    const auto& ssquared_state_left = ssquared_converter.at(index_of_vector);

    const auto& ssquared_states_right = ssquared_converter.block_with_number(index_of_vector);

    for (size_t col = row; col < ssquared_states_right.size(); ++col) {
        const auto& ssquared_state_right = ssquared_states_right.at(col);
        double total_9j = ssquared_converter.total_9j_coefficient(ssquared_state_left, ssquared_state_right, ranks);
        double value = factor * total_9j;

        matrix.add_to_position(value, row, col);
    }
}

std::unique_ptr<Term> ScalarProductITOTerm::clone() const {
    return std::make_unique<ScalarProductITOTerm>(ssquared_converter_, coefficients_);
}

}  // namespace model::operators