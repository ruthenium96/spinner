#include "ScalarProductTerm.h"

#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#include "src/common/index_converter/s_squared/OrderOfSummation.h"
#include "src/common/index_converter/s_squared/Level.h"
#include "src/spin_algebra/Multiplicity.h"

namespace {
// todo: this value can be calculated once per ranks
inline double local_product(const std::vector<spin_algebra::Multiplicity>& mults, const std::vector<uint8_t>& ranks) {
    double answer = 1;

    for (size_t i = 0; i < mults.size(); ++i) {
        double spin = ((double)mults[i] - 1) / 2;
        if (ranks[i] == 0) {
            answer *= sqrt(2 * spin + 1);
        } else if (ranks[i] == 1) {
            answer *= sqrt(spin * (spin + 1) * (2 * spin + 1));
        } else {
            throw std::invalid_argument("We can calculate local products only for ranks 0 and 1");
        }
    }

    return answer;
}
} // namespace

namespace model::operators::ito {

ScalarProductTerm::ScalarProductTerm(
    std::shared_ptr<const index_converter::s_squared::IndexConverter> converter,
    std::shared_ptr<const TwoDNumericalParameters<double>> isotropic_exchange_parameters) :
    TwoCenterTerm(isotropic_exchange_parameters->size()),
    converter_(converter),
    coefficients_(std::move(isotropic_exchange_parameters)) {}

void ScalarProductTerm::construct(
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
}

void ScalarProductTerm::add_scalar_product(
    quantum::linear_algebra::AbstractSymmetricMatrix& matrix,
    const std::set<unsigned int>& indexes_of_vectors,
    uint32_t center_a,
    uint32_t center_b,
    double factor) const {

    auto ranks = constructRanksOfTZero(center_a, center_b);

    std::vector<index_converter::s_squared::Level> levels_right;

    for (const auto& index_of_vector_row : indexes_of_vectors) {
        const auto& state_left = converter_->convert_index_to_state(index_of_vector_row);
        const auto& level_left = state_left.first;
        auto projection = state_left.second;

        construct_overlapping_levels(level_left, ranks, levels_right);
        for (const auto& level_right : levels_right) {
            auto mb_index_of_vector_col = converter_->convert_state_to_index(level_right, projection);
            if (!mb_index_of_vector_col.has_value()) {
                continue;
            }
            auto index_of_vector_col = mb_index_of_vector_col.value();
            if (index_of_vector_col < index_of_vector_row) {
                continue;
            }
            double total_9j = total_9j_coefficient(level_left, level_right, ranks);
            double value = factor * total_9j;
            if (value != 0.0) {
                matrix.add_to_position(value, index_of_vector_row, index_of_vector_col);
            }
        }
    }
}

std::unique_ptr<Term> ScalarProductTerm::clone() const {
    return std::make_unique<ScalarProductTerm>(converter_, coefficients_);
}

std::vector<uint8_t>
ScalarProductTerm::constructRanksOfTZero(uint32_t center_a, uint32_t center_b) const {
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

double ScalarProductTerm::total_9j_coefficient(
    const index_converter::s_squared::Level& left,
    const index_converter::s_squared::Level& right,
    const std::vector<uint8_t>& ranks) const {
    // product of ninejs
    double ninejs = 1;
    // product of square roots
    double square_roots_prod = 1;

    for (const auto& instruction : converter_->getOrderOfSummation()->getInstructions()) {
        size_t pos_one = instruction.positions_of_summands[0];
        size_t pos_two = instruction.positions_of_summands[1];
        size_t pos_fin = instruction.position_of_sum;

        double left_one = left.getSpin(pos_one);
        double left_two = left.getSpin(pos_two);
        double left_fin = left.getSpin(pos_fin);

        double right_one = right.getSpin(pos_one);
        double right_two = right.getSpin(pos_two);
        double right_fin = right.getSpin(pos_fin);

        uint8_t rank_one = ranks.at(pos_one);
        uint8_t rank_two = ranks.at(pos_two);
        uint8_t rank_fin = ranks.at(pos_fin);

        ninejs *= clebshGordanCalculator_.ninej_element(left_one, left_two, left_fin,
                                                        right_one, right_two, right_fin,
                                                        rank_one, rank_two, rank_fin);
        if (ninejs == 0) {
            return 0;
        }

        double mult_left_fin = left.getMultiplicity(pos_fin);
        double mult_right_fin = right.getMultiplicity(pos_fin);

        square_roots_prod *= (2 * rank_fin + 1) * (mult_left_fin) * (mult_right_fin);
    }
    square_roots_prod = sqrt(square_roots_prod);

    // product of strange things:
    // TODO: can be calculated only once if we swap left and right
    double local_prod = local_product(*right.getInitialMultiplicities(), ranks);

    double final_mult = left.getMultiplicity(left.getSize() - 1);

    return ninejs * square_roots_prod * local_prod / sqrt(final_mult);
}

void ScalarProductTerm::construct_overlapping_levels(const index_converter::s_squared::Level& level, 
    const std::vector<uint8_t>& ranks,
    std::vector<index_converter::s_squared::Level>& answer) const {
    auto number_of_summations = level.getInitialMultiplicities()->size();

    auto empty_level = index_converter::s_squared::Level(
        level.getInitialMultiplicities(),
        converter_->getOrderOfSummation()->size());

    answer.clear();
    answer.push_back(empty_level);

    for (const auto& instruction : converter_->getOrderOfSummation()->getInstructions()) {
        auto pos_one = instruction.positions_of_summands[0];
        auto pos_two = instruction.positions_of_summands[1];
        auto pos_fin = instruction.position_of_sum;
        std::vector<index_converter::s_squared::Level> temp_result;
        temp_result.reserve(answer.capacity());

        for (const auto& incomplete_state : answer) {
            if (ranks[pos_one] == 0) {
                if (incomplete_state.getMultiplicity(pos_one) != level.getMultiplicity(pos_one)) {
                    continue;
                }
            }
            if (ranks[pos_two] == 0) {
                if (incomplete_state.getMultiplicity(pos_two) != level.getMultiplicity(pos_two)) {
                    continue;
                }
            }

            std::vector<spin_algebra::Multiplicity> to_add_multiplicities = {level.getMultiplicity(pos_fin)};
            if (ranks[pos_fin] == 1) {
                auto mult_one = incomplete_state.getMultiplicity(pos_one);
                auto mult_two = incomplete_state.getMultiplicity(pos_two);
                spin_algebra::Multiplicity min_multiplicity = std::min(mult_one, mult_two);
                spin_algebra::Multiplicity max_multiplicity = std::max(mult_one, mult_two);

                auto min_multiplicity_sum = max_multiplicity - min_multiplicity + 1;
                auto max_multiplicity_sum = max_multiplicity + min_multiplicity - 1;

                if (min_multiplicity_sum + 2 <= level.getMultiplicity(pos_fin)) {
                    to_add_multiplicities.push_back(level.getMultiplicity(pos_fin) - 2);
                }
                if (max_multiplicity_sum >= level.getMultiplicity(pos_fin) + 2) {
                    to_add_multiplicities.push_back(level.getMultiplicity(pos_fin) + 2);
                }
            }
            for (const auto& to_add_multiplicity : to_add_multiplicities) {
                temp_result.push_back(incomplete_state);
                temp_result.back().setMultiplicity(pos_fin, to_add_multiplicity);
            }
        }

        std::swap(answer, temp_result);
        temp_result.clear();
    }
}

}  // namespace model::operators::ito