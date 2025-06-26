#include <cmath>
#include "gtest/gtest.h"
#include "src/common/physical_optimization/OptimizationList.h"
#include "src/common/runner/Runner.h"

void expect_mu_squared_approximate_equality(runner::Runner& exact, runner::Runner& ftlm) {
    std::vector<double> temperatures;
    for (int power = -3; power <= 3; ++power) {
        for (int number = 1; number <= 9; ++number) {
            temperatures.push_back(number * std::pow(10, power));
        }
    }

    for (auto temp : temperatures) {
        auto exact_value = exact.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(temp);
        auto ftlm_value = ftlm.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(temp);
        EXPECT_NEAR(exact_value.mean(), ftlm_value.mean(), 2 * ftlm_value.stdev_total()) << "Temperature: " << temp;
    }
}

TEST(ftlm_integration_tests, 10x2_FM_ring_SameG) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    model::ModelInput model(mults);
    auto g = model.addSymbol("g", 2.0);
    auto J = model.addSymbol("J", +10.0);
    for (int center = 0; center < mults.size(); ++center) {
        model.assignSymbolToGFactor(g, center);
        model.assignSymbolToIsotropicExchange(J, center, (center + 1) % mults.size());
    }

    common::physical_optimization::OptimizationList exact_optimization_list;
    exact_optimization_list.TzSort().EliminatePositiveProjections();
    runner::Runner exact_runner(model, exact_optimization_list);

    common::physical_optimization::OptimizationList ftlm_optimization_list;
    ftlm_optimization_list.TzSort().EliminatePositiveProjections().FTLMApproximate({100, 128, 100});
    runner::Runner ftlm_runner(model, ftlm_optimization_list);
    
    expect_mu_squared_approximate_equality(exact_runner, ftlm_runner);
}

TEST(ftlm_integration_tests, 12x2_FM_ring_SameG) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    model::ModelInput model(mults);
    auto g = model.addSymbol("g", 2.0);
    auto J = model.addSymbol("J", +10.0);
    for (int center = 0; center < mults.size(); ++center) {
        model.assignSymbolToGFactor(g, center);
        model.assignSymbolToIsotropicExchange(J, center, (center + 1) % mults.size());
    }

    common::physical_optimization::OptimizationList exact_optimization_list;
    exact_optimization_list.TzSort().EliminatePositiveProjections();
    runner::Runner exact_runner(model, exact_optimization_list);

    common::physical_optimization::OptimizationList ftlm_optimization_list;
    ftlm_optimization_list.TzSort().EliminatePositiveProjections().FTLMApproximate({200, 256, 100});
    runner::Runner ftlm_runner(model, ftlm_optimization_list);
    
    expect_mu_squared_approximate_equality(exact_runner, ftlm_runner);
}

TEST(ftlm_integration_tests, 14x2_FM_ring_SameG) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    model::ModelInput model(mults);
    auto g = model.addSymbol("g", 2.0);
    auto J = model.addSymbol("J", +10.0);
    for (int center = 0; center < mults.size(); ++center) {
        model.assignSymbolToGFactor(g, center);
        model.assignSymbolToIsotropicExchange(J, center, (center + 1) % mults.size());
    }

    common::physical_optimization::OptimizationList exact_optimization_list;
    exact_optimization_list.TzSort().EliminatePositiveProjections();
    runner::Runner exact_runner(model, exact_optimization_list);

    common::physical_optimization::OptimizationList ftlm_optimization_list;
    ftlm_optimization_list.TzSort().EliminatePositiveProjections().FTLMApproximate({200, 512, 100});
    runner::Runner ftlm_runner(model, ftlm_optimization_list);
    
    expect_mu_squared_approximate_equality(exact_runner, ftlm_runner);
}

TEST(ftlm_integration_tests, 10x2_FM_ring_DifferentG) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    model::ModelInput model(mults);
    auto g_one = model.addSymbol("g1", 2.0);
    auto g_two = model.addSymbol("g2", 3.0);
    auto J = model.addSymbol("J", +10.0);
    for (int center = 0; center < mults.size(); ++center) {
        model.assignSymbolToIsotropicExchange(J, center, (center + 1) % mults.size());
    }
    for (int center = 0; center < mults.size(); center+=2) {
        model.assignSymbolToGFactor(g_one, center);
    }
    for (int center = 1; center < mults.size(); center+=2) {
        model.assignSymbolToGFactor(g_two, center);
    }

    common::physical_optimization::OptimizationList exact_optimization_list;
    exact_optimization_list.TzSort().EliminatePositiveProjections();
    runner::Runner exact_runner(model, exact_optimization_list);

    common::physical_optimization::OptimizationList ftlm_optimization_list;
    ftlm_optimization_list.TzSort().EliminatePositiveProjections().FTLMApproximate({100, 128, 100});
    runner::Runner ftlm_runner(model, ftlm_optimization_list);
    
    expect_mu_squared_approximate_equality(exact_runner, ftlm_runner);
}

TEST(ftlm_integration_tests, 12x2_FM_ring_DifferentG) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    model::ModelInput model(mults);
    auto g_one = model.addSymbol("g1", 2.0);
    auto g_two = model.addSymbol("g2", 3.0);
    auto J = model.addSymbol("J", +10.0);
    for (int center = 0; center < mults.size(); ++center) {
        model.assignSymbolToIsotropicExchange(J, center, (center + 1) % mults.size());
    }
    for (int center = 0; center < mults.size(); center+=2) {
        model.assignSymbolToGFactor(g_one, center);
    }
    for (int center = 1; center < mults.size(); center+=2) {
        model.assignSymbolToGFactor(g_two, center);
    }

    common::physical_optimization::OptimizationList exact_optimization_list;
    exact_optimization_list.TzSort().EliminatePositiveProjections();
    runner::Runner exact_runner(model, exact_optimization_list);

    common::physical_optimization::OptimizationList ftlm_optimization_list;
    ftlm_optimization_list.TzSort().EliminatePositiveProjections().FTLMApproximate({200, 256, 100});
    runner::Runner ftlm_runner(model, ftlm_optimization_list);
    
    expect_mu_squared_approximate_equality(exact_runner, ftlm_runner);
}

TEST(ftlm_integration_tests, 14x2_FM_ring_DifferentG) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    model::ModelInput model(mults);
    auto g_one = model.addSymbol("g1", 2.0);
    auto g_two = model.addSymbol("g2", 3.0);
    auto J = model.addSymbol("J", +10.0);
    for (int center = 0; center < mults.size(); ++center) {
        model.assignSymbolToIsotropicExchange(J, center, (center + 1) % mults.size());
    }
    for (int center = 0; center < mults.size(); center+=2) {
        model.assignSymbolToGFactor(g_one, center);
    }
    for (int center = 1; center < mults.size(); center+=2) {
        model.assignSymbolToGFactor(g_two, center);
    }

    common::physical_optimization::OptimizationList exact_optimization_list;
    exact_optimization_list.TzSort().EliminatePositiveProjections();
    runner::Runner exact_runner(model, exact_optimization_list);

    common::physical_optimization::OptimizationList ftlm_optimization_list;
    ftlm_optimization_list.TzSort().EliminatePositiveProjections().FTLMApproximate({200, 512, 100});
    runner::Runner ftlm_runner(model, ftlm_optimization_list);
    
    expect_mu_squared_approximate_equality(exact_runner, ftlm_runner);
}