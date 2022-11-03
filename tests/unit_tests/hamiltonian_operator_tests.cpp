#include "gtest/gtest.h"
#include "src/common/runner/Runner.h"

TEST(hamiltonian_operator, throw_isotropic_exchange_same_center_22_333_4444_23456) {
    std::vector<std::vector<int>> vector_of_mults =
        {{2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        model::ModelInput model(mults);

        auto J = model.getSymbols().addSymbol("J", 10);
        EXPECT_THROW(
            model.getSymbols().assignSymbolToIsotropicExchange(J, 0, 0),
            std::invalid_argument);
    }
}

TEST(hamiltonian_operator, exchange_interaction_22_333_4444_23456) {
    std::vector<std::vector<int>> vector_of_mults =
        {{2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        model::ModelInput modelInput(mults);

        auto J = modelInput.getSymbols().addSymbol("J", 10);
        modelInput.getSymbols().assignSymbolToIsotropicExchange(J, 0, 1);
        // explicitly initialize isotropic exchange:

        model::Model model(modelInput);

        EXPECT_EQ(model.getOperator(common::QuantityEnum::Energy).getTwoCenterTerms().size(), 1);
    }
}
