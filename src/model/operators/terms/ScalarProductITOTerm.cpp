#include "ScalarProductITOTerm.h"

#include <cmath>
#include <string>
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
    const std::set<unsigned int>& indexes_of_vectors,
    uint32_t center_a,
    uint32_t center_b) const {
    if (!std::isnan(coefficients_->at(center_a, center_b))) {
        double factor = 2 * sqrt(3) * coefficients_->at(center_a, center_b);
        add_scalar_product(
            matrix,
            indexes_of_vectors,
            center_a,
            center_b,
            factor);
    }
//    matrix.print(std::cout);
}

void ScalarProductITOTerm::add_scalar_product(
    quantum::linear_algebra::AbstractSymmetricMatrix& matrix,
    const std::set<unsigned int>& indexes_of_vectors,
    uint32_t center_a,
    uint32_t center_b,
    double factor) const {
    const auto& ssquared_converter = *ssquared_converter_;

    auto ranks = ssquared_converter.constructRanksOfTZero(center_a, center_b);

    for (const auto& index_of_vector_row : indexes_of_vectors) {
        const auto& ssquared_state_left = ssquared_converter.at(index_of_vector_row);
        for (const auto& index_of_vector_col : indexes_of_vectors) {
            if (index_of_vector_col < index_of_vector_row) {
                continue;
            }
            const auto& ssquared_state_right = ssquared_converter.at(index_of_vector_col);
            double total_9j = ssquared_converter.total_9j_coefficient(ssquared_state_left, ssquared_state_right, ranks);
            double value = factor * total_9j;

            matrix.add_to_position(value, index_of_vector_row, index_of_vector_col);
        }
    }
}

std::unique_ptr<Term> ScalarProductITOTerm::clone() const {
    return std::make_unique<ScalarProductITOTerm>(ssquared_converter_, coefficients_);
}

}  // namespace model::operators