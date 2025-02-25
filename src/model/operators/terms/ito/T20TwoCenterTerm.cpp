#include "T20TwoCenterTerm.h"

#include <cmath>

namespace model::operators::ito {

T20TwoCenterTerm::T20TwoCenterTerm(
	std::shared_ptr<const index_converter::s_squared::IndexConverter> converter,
	std::shared_ptr<const TwoDNumericalParameters<double>> coefficients,
	double prefactor) :
	TwoCenterTerm(coefficients->size()),
	converter_(converter), 
	coefficients_(coefficients),
	prefector_(prefactor),
    wigner_eckart_helper_(converter->getOrderOfSummation()) {}

std::unique_ptr<Term> T20TwoCenterTerm::clone() const {
	return std::make_unique<T20TwoCenterTerm>(converter_, coefficients_, prefector_);
};

void T20TwoCenterTerm::construct(
	quantum::linear_algebra::AbstractSymmetricMatrix& matrix,
	const std::set<unsigned int>& indexes_of_vectors,
	uint32_t center_a, uint32_t center_b) const {
	
    if (!std::isnan(coefficients_->at(center_a, center_b))) {
        double factor_ttwo = coefficients_->at(center_a, center_b) * prefector_; 
        add_ttwo_term(
            matrix,
            indexes_of_vectors,
            center_a, center_b,
            factor_ttwo);
    }
};

void T20TwoCenterTerm::add_ttwo_term(
	quantum::linear_algebra::AbstractSymmetricMatrix& matrix,
	const std::set<unsigned int>& indexes_of_vectors,
	uint32_t center_a, uint32_t center_b,
	double factor) const {

	auto ranks = constructRanksOfTTwo(center_a, center_b);

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
            double value = factor * total_9j;

            double real_projection = - ((double)level_left.total() - 1.0) / 2.0 + projection;
            double left_spin = ((double)level_left.total() - 1.0) / 2.0;
            double right_spin = ((double)level_right.total() - 1.0) / 2.0;

            value *= wigner_eckart_helper_.clebsh_gordan_coefficient(left_spin, 2, right_spin, real_projection, 0);
			if (value != 0.0) {
				matrix.add_to_position(value, index_of_vector_row, index_of_vector_col);
			}
		}
	}	
}

std::vector<uint8_t> T20TwoCenterTerm::constructRanksOfTTwo(uint32_t center_a, uint32_t center_b) const {
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

        uint8_t result_rank = answer[first_spin] + answer[second_spin];
        answer[result_spin] = result_rank;
    }

    return answer;
}

} // namespace model::operators::ito