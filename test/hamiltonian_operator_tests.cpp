#include "common/runner/Runner.h"
#include "gtest/gtest.h"

TEST(hamiltonian_operator, throw_isotropic_exchange_same_center_22_333_4444_23456) {
    std::vector<std::vector<int>> vector_of_mults =
        {{2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        runner.AddSymbol("J", 10);
        EXPECT_THROW(runner.AddIsotropicExchange("J", 0, 0), std::invalid_argument);
    }
}

TEST(hamiltonian_operator, exchange_interaction_with_and_without_finalization_22_333_4444_23456) {
    std::vector<std::vector<int>> vector_of_mults =
        {{2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        runner.AddSymbol("J", 10);
        runner.AddIsotropicExchange("J", 0, 1);
        EXPECT_EQ(runner.getOperator(common::QuantityEnum::Energy).two_center_terms.size(), 0);
        runner.FinalizeIsotropicInteraction();
        EXPECT_EQ(runner.getOperator(common::QuantityEnum::Energy).two_center_terms.size(), 1);
    }
}

TEST(hamiltonian_operator, throw_adding_isotropic_exchange_after_finalization_333_4444_23456) {
    std::vector<std::vector<int>> vector_of_mults = {{3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        runner.AddSymbol("J", 10);
        runner.AddIsotropicExchange("J", 0, 1);
        runner.FinalizeIsotropicInteraction();
        EXPECT_THROW(runner.AddIsotropicExchange("J", 1, 2), std::invalid_argument);
    }
}

TEST(hamiltonian_operator, throw_finalization_without_inizialization_333_4444_23456) {
    std::vector<std::vector<int>> vector_of_mults = {{3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        runner.AddSymbol("J", 10);
        EXPECT_THROW(runner.FinalizeIsotropicInteraction(), std::length_error);
    }
}