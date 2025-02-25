#include "WignerEckartHelper.h"
#include <cmath>

namespace model::operators::ito {

WignerEckartHelper::WignerEckartHelper(
	std::shared_ptr<const index_converter::s_squared::OrderOfSummation> order_of_summation
) : order_of_summation_(order_of_summation) {}

double WignerEckartHelper::clebsh_gordan_coefficient(
	double l1,
	double l2,
	double l3,
	double m1,
	double m2) const {
	return clebshGordanCalculator_.clebsh_gordan_coefficient(l1, l2, l3, m1, m2);
}

double WignerEckartHelper::local_product(
	const std::vector<spin_algebra::Multiplicity>& mults, 
	const std::vector<uint8_t>& ranks) const {

	double answer = 1;

	for (size_t i = 0; i < mults.size(); ++i) {
		double spin = ((double)mults[i] - 1) / 2;
		if (ranks[i] == 0) {
			answer *= 2 * spin + 1;
		} else if (ranks[i] == 1) {
			answer *= spin * (spin + 1) * (2 * spin + 1);
		} else if (ranks[i] == 2) {
			answer *= (2 * spin + 3) * (2 * spin + 1) * (spin + 1) * spin * (2 * spin - 1) / 6;
		} else {
			throw std::invalid_argument("We can calculate local products only for ranks 0, 1 and 2");
		}
	}

	answer = sqrt(answer);

	return answer;
}


double WignerEckartHelper::total_9j_coefficient(
	const index_converter::s_squared::Level& left,
	const index_converter::s_squared::Level& right,
	const std::vector<uint8_t>& ranks,
	double local_prod) const {
	// product of ninejs
	double ninejs = 1;
	// product of square roots
	double square_roots_prod = 1;

	for (const auto& instruction : order_of_summation_->getInstructions()) {
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

		double mult_left_fin = 2.0 * left_fin + 1.0;
		double mult_right_fin = 2.0 * right_fin + 1.0;

		square_roots_prod *= (2 * rank_fin + 1) * (mult_left_fin) * (mult_right_fin);
	}
	square_roots_prod = sqrt(square_roots_prod);

	double final_mult = left.getMultiplicity(left.getSize() - 1);

	return ninejs * square_roots_prod * local_prod / sqrt(final_mult);
}


void WignerEckartHelper::construct_overlapping_levels(
	const index_converter::s_squared::Level& level, 
	const std::vector<uint8_t>& ranks,
	std::vector<index_converter::s_squared::Level>& answer) const {
	    auto number_of_summations = level.getInitialMultiplicities()->size();

    auto empty_level = index_converter::s_squared::Level(
        level.getInitialMultiplicities(),
        order_of_summation_->size());

    answer.clear();
    answer.push_back(empty_level);
    std::vector<index_converter::s_squared::Level> temp_result;

    for (const auto& instruction : order_of_summation_->getInstructions()) {
        auto pos_one = instruction.positions_of_summands[0];
        auto pos_two = instruction.positions_of_summands[1];
        auto pos_fin = instruction.position_of_sum;
        temp_result.reserve(answer.capacity());

        for (int i = 0; i < answer.size(); ++i) {
            auto incomplete_state = std::move(answer[i]);
            auto mult_fin_level = level.getMultiplicity(pos_fin);

            std::vector<spin_algebra::Multiplicity> to_add_multiplicities = {mult_fin_level};
            if (ranks[pos_fin] > 0) {
                auto mult_one = incomplete_state.getMultiplicity(pos_one);
                auto mult_two = incomplete_state.getMultiplicity(pos_two);
                spin_algebra::Multiplicity min_multiplicity = std::min(mult_one, mult_two);
                spin_algebra::Multiplicity max_multiplicity = std::max(mult_one, mult_two);

                auto min_multiplicity_sum = max_multiplicity - min_multiplicity + 1;
                auto max_multiplicity_sum = max_multiplicity + min_multiplicity - 1;

                if (min_multiplicity_sum + 2 <= mult_fin_level) {
                    to_add_multiplicities.push_back(mult_fin_level - 2);
                }
                if (max_multiplicity_sum >= mult_fin_level + 2) {
                    to_add_multiplicities.push_back(mult_fin_level + 2);
                }
				if (ranks[pos_fin] == 2) {
					if (min_multiplicity_sum + 4 <= mult_fin_level) {
						to_add_multiplicities.push_back(mult_fin_level - 4);
					}	
					if (max_multiplicity_sum >= mult_fin_level + 4) {
						to_add_multiplicities.push_back(mult_fin_level + 4);
					}	
				}
            }
            // Reuse incomplete_state at least once.
            // Surprisingly, it significantly speeds up the execution of this function.
            for (int i = 0; i + 1 < to_add_multiplicities.size(); ++i) {
                auto to_add_multiplicity = to_add_multiplicities[i];
                temp_result.push_back(incomplete_state);
                temp_result.back().setMultiplicity(pos_fin, to_add_multiplicity);
            }
            auto to_add_multiplicity = to_add_multiplicities.back();
            temp_result.push_back(std::move(incomplete_state));
            temp_result.back().setMultiplicity(pos_fin, to_add_multiplicity);
        }

        std::swap(answer, temp_result);
        temp_result.clear();
    }
}

} // namespace model::operators::ito