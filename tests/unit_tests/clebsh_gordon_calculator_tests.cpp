#include <gtest/gtest.h>
#include <cmath>

#include "src/spin_algebra/ClebshGordanCalculator.h"

const double EPSILON = 1e-9;

TEST(ClebshGordanCalculatorTest, Wigner9j) {
	auto cg_calculator = spin_algebra::ClebshGordanCalculator();

	std::vector<std::pair<std::array<double, 9>, double>> test_cases = {
		{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0}, 1.0}, 
		{{0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0, 0, 0}, sqrt(3.0)/6.0},
		{{1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0, 0, 0}, 1.0/3.0},
		{{2.0, 2.0, 0.0, 2.0, 2.0, 0.0, 0, 0, 0}, 1.0/5.0},
		{{2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0, 0, 0}, sqrt(5.0)/25.0},
		{{3.0, 0.0, 3.0, 3.0, 0.0, 3.0, 0, 0, 0}, 1.0/7.0},

		{{1.0, 4.0, 4.0, 1.0, 3.0, 3.0, 0, 1, 1}, sqrt(105.0)/252.0},
		{{2.0, 0.5, 1.5, 2.0, 1.5, 1.5, 0, 1, 1}, -sqrt(3.0)/30.0},
		{{2.0, 2.0, 2.0, 2.0, 1.0, 1.0, 0, 1, 1}, sqrt(35.0)/150.0},
		{{3.0, 3.0, 0.0, 3.0, 2.0, 1.0, 0, 1, 1}, 1.0/21.0},
		{{4.0, 3.0, 2.0, 4.0, 4.0, 2.0, 0, 1, 1}, -sqrt(770.0)/1260.0},

		{{1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1, 1, 0}, 1.0/9.0},
		{{1.5, 1.5, 1.0, 1.5, 1.5, 1.0, 1, 1, 0}, 11.0/180.0},
		{{2.0, 1.5, 0.5, 1.0, 0.5, 0.5, 1, 1, 0}, sqrt(2.0)/12.0},

		{{2.0, 1.5, 1.5, 2.0, 1.5, 0.5, 1, 0, 1}, sqrt(30.0)/120.0},
		{{2.0, 1.5, 1.5, 2.0, 1.5, 1.5, 1, 0, 1}, sqrt(6.0)/60.0 },
		{{2.0, 2.0, 0.0, 1.0, 2.0, 1.0, 1, 0, 1}, 1.0/15.0},

		{{0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 1, 1, 2}, 1.0/9.0},
		{{1.0, 1.0, 0.0, 1.0, 1.0, 2.0, 1, 1, 2}, -0.04303314829119352},
		{{3.0, 0.0, 3.0, 2.0, 1.0, 1.0, 1, 1, 2}, sqrt(21.0)/105.0},
		{{3.0, 0.0, 3.0, 2.0, 1.0, 2.0, 1, 1, 2}, -sqrt(14.0)/105.0},
		{{3.0, 0.0, 3.0, 2.0, 1.0, 3.0, 1, 1, 2}, sqrt(6.0)/105.0},
		{{3.5, 3.0, 0.5, 2.5, 2.0, 1.5, 1, 1, 2}, -sqrt(2.0)/210.0},
		{{4.0, 3.0, 1.0, 3.0, 4.0, 1.0, 1, 1, 2}, 19.0/756.0},
		{{4.0, 3.0, 1.0, 3.0, 4.0, 2.0, 1, 1, 2}, -2.0*sqrt(5.0)/315.0},
		{{4.0, 3.0, 1.0, 3.0, 4.0, 3.0, 1, 1, 2}, sqrt(77.0)/1764.0},
		{{4.0, 4.0, 4.0, 3.0, 4.0, 2.0, 1, 1, 2}, -sqrt(858.0)/3780.0},
		{{4.0, 4.0, 4.0, 3.0, 4.0, 3.0, 1, 1, 2}, sqrt(390.0)/1890.0},
		{{4.0, 4.0, 4.0, 3.0, 4.0, 4.0, 1, 1, 2}, -4*sqrt(858.0)/31185.0},

		{{1.0, 3.0, 3.0, 1.0, 2.0, 3.0, 0, 2, 2}, sqrt(2.0)/70.0},
		{{1.0, 3.0, 3.0, 1.0, 3.0, 2.0, 0, 2, 2}, -sqrt(2.0)/70.0},
		{{1.0, 3.0, 3.0, 1.0, 3.0, 3.0, 0, 2, 2}, sqrt(15.0)/140.0},
		{{2.5, 1.0, 1.5, 2.5, 3.0, 0.5, 0, 2, 2}, -1.0/30.0},
		{{2.5, 1.0, 1.5, 2.5, 3.0, 1.5, 0, 2, 2}, sqrt(5.0)/75.0},

		{{3.0, 3.5, 3.5, 2.0, 3.5, 3.5, 2, 0, 2}, -sqrt(66.0)/840.0},
		{{3.0, 3.5, 3.5, 3.0, 3.5, 1.5, 2, 0, 2}, sqrt(165.0)/840.0},
		{{3.5, 3.5, 1.0, 2.5, 3.5, 1.0, 2, 0, 2}, sqrt(70.0)/560.0},
		{{3.5, 3.5, 1.0, 2.5, 3.5, 2.0, 2, 0, 2}, 13*sqrt(42.0)/5040.0},
		{{3.5, 3.5, 1.0, 2.5, 3.5, 3.0, 2, 0, 2}, -sqrt(15.0)/315.0},
	};

	for (const auto& pair : test_cases) {
		const auto& v = pair.first;
		double value = cg_calculator.ninej_element(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
		EXPECT_NEAR(value, pair.second, EPSILON);
	}
}