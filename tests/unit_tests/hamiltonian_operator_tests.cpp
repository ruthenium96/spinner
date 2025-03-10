#include "gtest/gtest.h"
#include "src/common/runner/ConsistentModelOptimizationList.h"

TEST(hamiltonian_operator, throw_isotropic_exchange_same_center_22_333_4444_23456) {
    std::vector<std::vector<spin_algebra::Multiplicity>> vector_of_mults =
        {{2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        model::ModelInput model(mults);

        auto J = model.addSymbol("J", 10);
        EXPECT_THROW(
            model.assignSymbolToIsotropicExchange(J, 0, 0),
            std::invalid_argument);
    }
}

TEST(hamiltonian_operator, exchange_interaction_22_333_4444_23456) {
    std::vector<std::vector<spin_algebra::Multiplicity>> vector_of_mults =
        {{2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        model::ModelInput modelInput(mults);

        // explicitly initialize isotropic exchange:
        auto J = modelInput.addSymbol("J", 10);
        modelInput.assignSymbolToIsotropicExchange(J, 0, 1);

        runner::ConsistentModelOptimizationList model_and_optimization(modelInput, common::physical_optimization::OptimizationList());

        EXPECT_EQ(model_and_optimization.getModel().getOperator(common::QuantityEnum::Energy).value()->getTerms().size(), 1);
    }
}
