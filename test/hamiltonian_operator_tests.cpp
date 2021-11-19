#include "common/runner/Runner.h"
#include "gtest/gtest.h"

TEST(finalize_isotropic_exchange, forget_to_finalize) {
    std::vector<std::vector<int>> vector_of_mults =
        {{2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        runner.AddSymbol("J", 10);
        runner.AddIsotropicExchange("J", 0, 1);
        //runner.FinalizeIsotropicInteraction();

        EXPECT_EQ(runner.getOperator(common::QuantityEnum::Energy).two_center_terms.size(), 0);
    }
}

TEST(finalize_isotropic_exchange, actually_added_new_term_22_333_4444_23456) {
    std::vector<std::vector<int>> vector_of_mults =
        {{2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        runner.AddSymbol("J", 10);
        runner.AddIsotropicExchange("J", 0, 1);
        runner.FinalizeIsotropicInteraction();

        EXPECT_EQ(runner.getOperator(common::QuantityEnum::Energy).two_center_terms.size(), 1);
    }
}
