#include "T00TwoCenterTerm.h"

#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#include "src/common/index_converter/s_squared/OrderOfSummation.h"
#include "src/common/index_converter/s_squared/Level.h"
#include "src/model/operators/terms/ito/WignerEckartHelper.h"

namespace model::operators::ito {

T00TwoCenterTerm::T00TwoCenterTerm(
    std::shared_ptr<const index_converter::s_squared::IndexConverter> converter,
    std::shared_ptr<const TwoDNumericalParameters<double>> isotropic_exchange_parameters,
    double prefactor) :
    TwoCenterTerm(isotropic_exchange_parameters->size()),
    converter_(converter),
    coefficients_(std::move(isotropic_exchange_parameters)),
    prefactor_(prefactor),
    wigner_eckart_helper_(converter->getOrderOfSummation()) {}

void T00TwoCenterTerm::construct(
    quantum::linear_algebra::AbstractSymmetricMatrix& matrix,
    const std::set<unsigned int>& indexes_of_vectors,
    uint32_t center_a,
    uint32_t center_b) const {
    if (!std::isnan(coefficients_->at(center_a, center_b))) {
        double factor = prefactor_ * coefficients_->at(center_a, center_b);
        add_scalar_product(
            matrix,
            indexes_of_vectors,
            center_a,
            center_b,
            factor);
    }
}

void T00TwoCenterTerm::add_scalar_product(
    quantum::linear_algebra::AbstractSymmetricMatrix& matrix,
    const std::set<unsigned int>& indexes_of_vectors,
    uint32_t center_a,
    uint32_t center_b,
    double factor) const {

    auto ranks = constructRanksOfTZero(center_a, center_b);

    std::vector<index_converter::s_squared::Level> levels_right;

    double local_prod = wigner_eckart_helper_.local_product(converter_->get_mults(), ranks);

    for (const auto& index_of_vector_row : indexes_of_vectors) {
        const auto& state_left = converter_->convert_index_to_state(index_of_vector_row);
        const auto& level_left = state_left.first;
        auto projection = state_left.second;

        wigner_eckart_helper_.construct_overlapping_levels(level_left, ranks, levels_right);
        for (const auto& level_right : levels_right) {
            auto mb_index_of_vector_col = converter_->convert_state_to_index(level_right, projection);
            if (!mb_index_of_vector_col.has_value()) {
                continue;
            }
            auto index_of_vector_col = mb_index_of_vector_col.value();
            if (index_of_vector_col < index_of_vector_row) {
                continue;
            }
            double total_9j = wigner_eckart_helper_.total_9j_coefficient(level_left, level_right, ranks, local_prod);
            // due to the Wigher-Echart theorem, we also need to multiply by CG-coefficient,
            // but it is (S M 0 0; S M) = 1
            double value = factor * total_9j;
            if (value != 0.0) {
                matrix.add_to_position(value, index_of_vector_row, index_of_vector_col);
            }
        }
    }
}

std::unique_ptr<Term> T00TwoCenterTerm::clone() const {
    return std::make_unique<T00TwoCenterTerm>(converter_, coefficients_, prefactor_);
}

std::vector<uint8_t>
T00TwoCenterTerm::constructRanksOfTZero(uint32_t center_a, uint32_t center_b) const {
    // todo: cache results, there is just O(N^3) memory required
    auto number_of_initial_mults = converter_->get_mults().size();
    auto number_of_all_mults = number_of_initial_mults + 
    converter_->getOrderOfSummation()->getInstructions().size();

    std::vector<uint8_t> answer(number_of_all_mults, 255); // 255 is for "not initialized yet".

    for (size_t i = 0; i < number_of_initial_mults; ++i) {
        if (i == center_a || i == center_b) {
            answer[i] = 1;
        } else {
            answer[i] = 0;
        }
    }

    for (const auto& instruction : *converter_->getOrderOfSummation()) {
        size_t first_spin = instruction.positions_of_summands.at(0);
        size_t second_spin = instruction.positions_of_summands.at(1);
        size_t result_spin = instruction.position_of_sum;

        uint8_t result_rank = (answer[first_spin] + answer[second_spin]) % 2;
        answer[result_spin] = result_rank;
    }

    return answer;
}

}  // namespace model::operators::ito