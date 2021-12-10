#include "common/runner/Runner.h"
#include "gtest/gtest.h"

TEST(hamiltonian_operator, throw_isotropic_exchange_same_center_22_333_4444_23456) {
    std::vector<std::vector<int>> vector_of_mults =
        {{2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        auto J = runner.AddSymbol("J", 10);
        EXPECT_THROW(runner.AssignSymbolToIsotropicExchange(J, 0, 0), std::invalid_argument);
    }
}

TEST(hamiltonian_operator, exchange_interaction_22_333_4444_23456) {
    std::vector<std::vector<int>> vector_of_mults =
        {{2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);

        auto J = runner.AddSymbol("J", 10);
        runner.AssignSymbolToIsotropicExchange(J, 0, 1);
        EXPECT_EQ(runner.getOperator(common::QuantityEnum::Energy).two_center_terms.size(), 1);
    }
}

TEST(
    hamiltonian_operator,
    throw_adding_isotropic_exchange_after_build_matrices_build_spectra_333_4444_23456) {
    std::vector<std::vector<int>> vector_of_mults = {{3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        {
            runner::Runner runner(mults);

            auto J = runner.AddSymbol("J", 10);
            runner.AssignSymbolToIsotropicExchange(J, 0, 1);
            runner.BuildMatrices();
            EXPECT_THROW(runner.AssignSymbolToIsotropicExchange(J, 1, 2), std::invalid_argument);
        }
        {
            runner::Runner runner(mults);

            auto J = runner.AddSymbol("J", 10);
            runner.AssignSymbolToIsotropicExchange(J, 0, 1);
            runner.BuildSpectra();
            EXPECT_THROW(runner.AssignSymbolToIsotropicExchange(J, 1, 2), std::invalid_argument);
        }
    }
}
