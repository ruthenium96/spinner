#include "gtest/gtest.h"
#include "src/common/runner/ConsistentModelOptimizationList.h"

TEST(consistentModelOptimizationList_tests, throw_wrong_size_of_pemutation) {
    std::vector<int> mults = {4, 4, 4};
    model::ModelInput model(mults);
    common::physical_optimization::OptimizationList optimizationList;
    optimizationList.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});

    EXPECT_THROW(
        runner::ConsistentModelOptimizationList(model, optimizationList),
        std::length_error);
}

TEST(consistentModelOptimizationList_tests, throw_permutes_different_multiplicities) {
    std::vector<int> mults = {4, 4, 4, 3};
    model::ModelInput model(mults);
    common::physical_optimization::OptimizationList optimizationList;
    optimizationList.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});

    EXPECT_THROW(
        runner::ConsistentModelOptimizationList(model, optimizationList),
        std::invalid_argument);
}

TEST(consistentModelOptimizationList_tests, throw_2222_isotropic_inconsistent_symmetry) {
    std::vector<int> mults = {2, 2, 2, 2};
    model::ModelInput model(mults);
    double J_value = 10;
    auto tripledJ = model.getSymbols().addSymbol("3J", 3 * J_value);
    auto J = model.getSymbols().addSymbol("J", J_value);
    auto doubledJ = model.getSymbols().addSymbol("2J", 2 * J_value);
    model.getSymbols()
        .assignSymbolToIsotropicExchange(tripledJ, 0, 1)
        .assignSymbolToIsotropicExchange(J, 1, 2)
        .assignSymbolToIsotropicExchange(doubledJ, 2, 3)
        .assignSymbolToIsotropicExchange(J, 3, 0);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList.TzSort()
        .Symmetrize(group::Group::S2, {{1, 0, 3, 2}})
        .Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

    EXPECT_THROW(
        runner::ConsistentModelOptimizationList(model, optimizationList),
        std::invalid_argument);
}

TEST(consistentModelOptimizationList_tests, throw_2222_isotropic_accidental_symmetry) {
    std::vector<int> mults = {2, 2, 2, 2};
    model::ModelInput model(mults);
    double J = 10;
    auto firstJ = model.getSymbols().addSymbol("J1", J);
    auto secondJ = model.getSymbols().addSymbol("J2", J);
    model.getSymbols()
        .assignSymbolToIsotropicExchange(firstJ, 0, 1)
        .assignSymbolToIsotropicExchange(firstJ, 1, 2)
        .assignSymbolToIsotropicExchange(firstJ, 2, 3)
        .assignSymbolToIsotropicExchange(secondJ, 3, 0);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList.TzSort()
        .Symmetrize(group::Group::S2, {{1, 0, 3, 2}})
        .Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

    EXPECT_THROW(
        runner::ConsistentModelOptimizationList(model, optimizationList),
        std::invalid_argument);
}

TEST(consistentModelOptimizationList_tests, throw_2222_gfactor_inconsistent_symmetry) {
    std::vector<int> mults = {2, 2, 2, 2};
    model::ModelInput model(mults);
    double J_value = 10;
    auto J = model.getSymbols().addSymbol("J", J_value);
    auto firstg = model.getSymbols().addSymbol("g1", 2.0);
    auto secondg = model.getSymbols().addSymbol("g2", 3.0);
    model.getSymbols()
        .assignSymbolToIsotropicExchange(J, 0, 1)
        .assignSymbolToIsotropicExchange(J, 1, 2)
        .assignSymbolToIsotropicExchange(J, 2, 3)
        .assignSymbolToIsotropicExchange(J, 3, 0)
        .assignSymbolToGFactor(firstg, 0)
        .assignSymbolToGFactor(firstg, 1)
        .assignSymbolToGFactor(firstg, 2)
        .assignSymbolToGFactor(secondg, 3);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList.TzSort()
        .Symmetrize(group::Group::S2, {{1, 0, 3, 2}})
        .Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

    EXPECT_THROW(
        runner::ConsistentModelOptimizationList(model, optimizationList),
        std::invalid_argument);
}

TEST(consistentModelOptimizationList_tests, throw_2222_gfactor_accidental_symmetry) {
    std::vector<int> mults = {2, 2, 2, 2};
    model::ModelInput model(mults);
    double J_value = 10;
    double g_value = 2.0;
    auto J = model.getSymbols().addSymbol("J", J_value);
    auto firstg = model.getSymbols().addSymbol("g1", g_value);
    auto secondg = model.getSymbols().addSymbol("g2", g_value);
    model.getSymbols()
        .assignSymbolToIsotropicExchange(J, 0, 1)
        .assignSymbolToIsotropicExchange(J, 1, 2)
        .assignSymbolToIsotropicExchange(J, 2, 3)
        .assignSymbolToIsotropicExchange(J, 3, 0)
        .assignSymbolToGFactor(firstg, 0)
        .assignSymbolToGFactor(firstg, 1)
        .assignSymbolToGFactor(firstg, 2)
        .assignSymbolToGFactor(secondg, 3);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList.TzSort()
        .Symmetrize(group::Group::S2, {{1, 0, 3, 2}})
        .Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

    EXPECT_THROW(
        runner::ConsistentModelOptimizationList(model, optimizationList),
        std::invalid_argument);
}
