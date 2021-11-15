#include "gtest/gtest.h"
#include "common/Runner.h"

TEST(add_isotropic_exchange, actually_added_new_term_22_333_4444_23456) {
    std::vector<std::vector<int>> vector_of_mults = {{2, 2},
                                                     {3, 3, 3},
                                                     {4, 4, 4, 4},
                                                     {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);
        spaces::LexicographicIndexConverter converter(mults);

        arma::dmat js(mults.size(), mults.size());
        for (int i = 0; i < mults.size(); ++i) {
            for (int j = 0; j < mults.size(); ++j) {
                js(i, j) = NAN;
            }
        }

        runner.AddIsotropicExchange(js);

        EXPECT_EQ(runner.getOperator(common::QuantityEnum::Energy).two_center_terms.size(), 1);
    }
}
