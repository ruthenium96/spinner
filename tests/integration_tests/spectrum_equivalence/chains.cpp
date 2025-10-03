#include <cstddef>
#include <string>
#include "src/group/Group.h"
#include "src/spin_algebra/Multiplicity.h"
#include "tests/integration_tests/spectrum_equivalence/equivalence_check.h"
#include "tests/tools/OptimizationListsGenerator.h"

using namespace common::physical_optimization;

void initialize_chain_system(model::ModelInput& model, 
    const std::vector<double>& J_values) {
    const size_t number_of_centers = model.getMults().size();
    std::vector<model::symbols::SymbolName> J_symbols;
    for (size_t i = 0; i < J_values.size(); ++i) {
        auto J_symbol = model.addSymbol("J" + std::to_string(i), J_values[i]);
        J_symbols.emplace_back(J_symbol);
    }
    if (number_of_centers % 2 == 1) {
        for (size_t i = 0; i < (number_of_centers - 1) / 2; ++i) {
            model.assignSymbolToIsotropicExchange(J_symbols[i], i, i + 1);
            model.assignSymbolToIsotropicExchange(J_symbols[i], 
                number_of_centers - 1 - i, 
                number_of_centers - 2 - i);
        }
    }
    if (number_of_centers % 2 == 0) {
        for (size_t i = 0; i < number_of_centers / 2 - 1; ++i) {
            model.assignSymbolToIsotropicExchange(J_symbols[i], i, i + 1);
            model.assignSymbolToIsotropicExchange(J_symbols[i], 
                number_of_centers - 1 - i, 
                number_of_centers - 2 - i);
        }
        model.assignSymbolToIsotropicExchange(J_symbols.back(),
                number_of_centers / 2, 
                number_of_centers / 2 - 1);
    }
}

void initialize_same_g_factor_for_chain_system(
    model::ModelInput& model, 
    double g_value) {
    auto g = model.addSymbol("g", g_value);
    for (int i = 0; i < model.getMults().size(); ++i) {
        model.assignSymbolToGFactor(g, i);
    }
}

void initialize_different_g_factors_for_chain_system(
    model::ModelInput& model, const std::vector<double>& g_values) {
    const size_t number_of_centers = model.getMults().size();
    std::vector<model::symbols::SymbolName> g_symbols;
    for (size_t i = 0; i < g_values.size(); ++i) {
        auto J_symbol = model.addSymbol("g" + std::to_string(i), g_values[i]);
        g_symbols.emplace_back(J_symbol);
    }
    if (number_of_centers % 2 == 1) {
        for (size_t i = 0; i < (number_of_centers - 1) / 2; ++i) {
            model.assignSymbolToGFactor(g_symbols[i], i);
            model.assignSymbolToGFactor(g_symbols[i], 
                number_of_centers - 1 - i);
        }
        model.assignSymbolToGFactor(g_symbols.back(), (number_of_centers - 1) / 2);
    }
    if (number_of_centers % 2 == 0) {
        for (size_t i = 0; i < number_of_centers / 2; ++i) {
            model.assignSymbolToGFactor(g_symbols[i], i);
            model.assignSymbolToGFactor(g_symbols[i], 
                number_of_centers - 1 - i);
        }
    }
}

group::Group generate_S2_group_for_chain(size_t size) {
    group::Permutation generator(size);
    for (int i = 0; i < size; ++i) {
        generator[size - 1 - i] = i;
    }
    return group::Group(group::Group::S2, {generator});
}

class two_center_chain : public SpectrumFinalEquivalenceTest {};

TEST_P(two_center_chain, NoGFactors) {
    std::vector<double> js = {-50, -40.1, -20, -1.01, 1, 33, 46};
    for (auto J_value : js) {
        model::ModelInput model(std::get<1>(GetParam()));
        initialize_chain_system(model, {J_value});

        runner::Runner runner_simple(model);

        runner::Runner runner(model, std::get<0>(GetParam()));
        expect_final_vectors_equivalence(runner_simple, runner);                
    }
}

INSTANTIATE_TEST_SUITE_P(
    spectrum_final_equivalence,
    two_center_chain,
    ::testing::Combine(
    ::testing::ValuesIn( 
        generate_all_optimization_lists({generate_S2_group_for_chain(2)})
    ),
    ::testing::Values(
        std::vector<spin_algebra::Multiplicity>({2, 2}),
        std::vector<spin_algebra::Multiplicity>({3, 3}),
        std::vector<spin_algebra::Multiplicity>({4, 4}),
        std::vector<spin_algebra::Multiplicity>({5, 5}),
        std::vector<spin_algebra::Multiplicity>({6, 6}),
        std::vector<spin_algebra::Multiplicity>({7, 7}),
        std::vector<spin_algebra::Multiplicity>({8, 8})
    )
    ),
    spectrum_final_equivalence_test_name_generator
);

class three_center_chain : public SpectrumFinalEquivalenceTest {};

TEST_P(three_center_chain, NoGFactors) {
    std::vector<double> js = {-50, -40.1, -20, -1.01, 1, 33, 46};
    for (auto J_value : js) {
        model::ModelInput model(std::get<1>(GetParam()));
        initialize_chain_system(model, {J_value});

        runner::Runner runner_simple(model);

        runner::Runner runner(model, std::get<0>(GetParam()));
        expect_final_vectors_equivalence(runner_simple, runner);                
    }
}

TEST_P(three_center_chain, DifferentGFactors) {
    std::vector<double> js = {-50, -40.1, -20, -1.01, 1, 33, 46};
    std::vector<std::vector<double>> gss = {
        {2.0, 2.0}, {2.0, 2.1}, {3.1, 4.1}
    };    
    for (auto J_value : js) {
        for (const auto& gs : gss) {
            model::ModelInput model(std::get<1>(GetParam()));
            initialize_chain_system(model, {J_value});
            initialize_different_g_factors_for_chain_system(model, gs);

            runner::Runner runner_simple(model);

            runner::Runner runner(model, std::get<0>(GetParam()));
            expect_final_vectors_equivalence(runner_simple, runner);
        }
    }
}

TEST_P(three_center_chain, SameGFactors) {
    std::vector<double> js = {-50, -40.1, -20, -1.01, 1, 33, 46};
    std::vector<double> gs = {2.0, 2.2, 3.1, 5.1};
    for (auto J_value : js) {
        for (auto g_value : gs) {
            model::ModelInput model(std::get<1>(GetParam()));
            initialize_chain_system(model, {J_value});
            initialize_same_g_factor_for_chain_system(model, g_value);

            runner::Runner runner_simple(model);

            runner::Runner runner(model, std::get<0>(GetParam()));
            expect_final_vectors_equivalence(runner_simple, runner);
        }
    }
}

INSTANTIATE_TEST_SUITE_P(
    spectrum_final_equivalence,
    three_center_chain,
    ::testing::Combine(
    ::testing::ValuesIn( 
        generate_all_optimization_lists({generate_S2_group_for_chain(3)})
    ),
    ::testing::Values(
        std::vector<spin_algebra::Multiplicity>({2, 2, 2}),
        std::vector<spin_algebra::Multiplicity>({3, 3, 3}),
        std::vector<spin_algebra::Multiplicity>({4, 4, 4}),
        std::vector<spin_algebra::Multiplicity>({5, 5, 5}),
        std::vector<spin_algebra::Multiplicity>({2, 3, 2}),
        std::vector<spin_algebra::Multiplicity>({2, 4, 2}),        
        std::vector<spin_algebra::Multiplicity>({3, 2, 3}),        
        std::vector<spin_algebra::Multiplicity>({4, 3, 4})
    )
    ),
    spectrum_final_equivalence_test_name_generator
);

class four_center_chain : public SpectrumFinalEquivalenceTest {};

TEST_P(four_center_chain, NoGFactors) {
    std::vector<std::pair<double, double>> js =
        {{10, 15}, {-10, 15}, {10, -15}, {-10, -15}, {-10, -10}, {10, 10}};

    for (auto [Jfirst, Jsecond] : js) {
        model::ModelInput model(std::get<1>(GetParam()));
        initialize_chain_system(model, {Jfirst, Jsecond});

        runner::Runner runner_simple(model);
        runner::Runner runner(model, std::get<0>(GetParam()));
        expect_final_vectors_equivalence(runner_simple, runner);                
    }
}

TEST_P(four_center_chain, DifferentGFactors) {
    std::vector<std::pair<double, double>> js =
        {{10, 15}, {-10, 15}, {10, -15}, {-10, -15}};
    std::vector<std::vector<double>> gss = {
        {2.0, 2.0}, {2.0, 2.1}, {3.1, 4.1}
    };

    for (auto [Jfirst, Jsecond] : js) {
        for (const auto gs : gss) {
            model::ModelInput model(std::get<1>(GetParam()));
            initialize_chain_system(model, {Jfirst, Jsecond});
            initialize_different_g_factors_for_chain_system(model, gs);

            runner::Runner runner_simple(model);
            runner::Runner runner(model, std::get<0>(GetParam()));
            expect_final_vectors_equivalence(runner_simple, runner);
        }           
    }
}

TEST_P(four_center_chain, SameGFactors) {
    std::vector<std::pair<double, double>> js =
        {{10, 15}, {-10, 15}, {10, -15}, {-10, -15}};
    std::vector<double> gs = {2.0, 2.2, 3.1, 5.1};

    for (auto [Jfirst, Jsecond] : js) {
        for (const auto g : gs) {
            model::ModelInput model(std::get<1>(GetParam()));
            initialize_chain_system(model, {Jfirst, Jsecond});
            initialize_same_g_factor_for_chain_system(model, g);

            runner::Runner runner_simple(model);
            runner::Runner runner(model, std::get<0>(GetParam()));
            expect_final_vectors_equivalence(runner_simple, runner);
        }           
    }
}

INSTANTIATE_TEST_SUITE_P(
    spectrum_final_equivalence,
    four_center_chain,
    ::testing::Combine(
    ::testing::ValuesIn(
        generate_all_optimization_lists({generate_S2_group_for_chain(4)})
    ),
    ::testing::Values(
        std::vector<spin_algebra::Multiplicity>({2, 2, 2, 2}),
        std::vector<spin_algebra::Multiplicity>({2, 3, 3, 2}),
        std::vector<spin_algebra::Multiplicity>({3, 3, 3, 3}),
        std::vector<spin_algebra::Multiplicity>({3, 2, 2, 3}),
        std::vector<spin_algebra::Multiplicity>({4, 4, 4, 4})
    )
),
spectrum_final_equivalence_test_name_generator
);

class five_center_chain : public SpectrumFinalEquivalenceTest {};

TEST_P(five_center_chain, NoGFactors) {
    std::vector<std::pair<double, double>> js =
        {{10, 15}, {-10, 15}, {10, -15}, {-10, -15}, {-10, -10}, {10, 10}};

    for (auto [Jfirst, Jsecond] : js) {
        model::ModelInput model(std::get<1>(GetParam()));
        initialize_chain_system(model, {Jfirst, Jsecond});

        runner::Runner runner_simple(model);
        runner::Runner runner(model, std::get<0>(GetParam()));
        expect_final_vectors_equivalence(runner_simple, runner);                
    }
}

TEST_P(five_center_chain, DifferentGFactors) {
    std::vector<std::pair<double, double>> js =
        {{10, 15}, {-10, 15}, {10, -15}, {-10, -15}};
    std::vector<std::vector<double>> gss = {
        {2.0, 2.0, 2.0}, {2.0, 2.1, 2.2}, {3.1, 4.1, 5.1}
    };

    for (auto [Jfirst, Jsecond] : js) {
        for (const auto gs : gss) {
            model::ModelInput model(std::get<1>(GetParam()));
            initialize_chain_system(model, {Jfirst, Jsecond});
            initialize_different_g_factors_for_chain_system(model, gs);

            runner::Runner runner_simple(model);
            runner::Runner runner(model, std::get<0>(GetParam()));
            expect_final_vectors_equivalence(runner_simple, runner);
        }           
    }
}

TEST_P(five_center_chain, SameGFactors) {
    std::vector<std::pair<double, double>> js =
        {{10, 15}, {-10, 15}, {10, -15}, {-10, -15}};
    std::vector<double> gs = {2.0, 2.2, 3.1, 5.1};

    for (auto [Jfirst, Jsecond] : js) {
        for (const auto g : gs) {
            model::ModelInput model(std::get<1>(GetParam()));
            initialize_chain_system(model, {Jfirst, Jsecond});
            initialize_same_g_factor_for_chain_system(model, g);

            runner::Runner runner_simple(model);
            runner::Runner runner(model, std::get<0>(GetParam()));
            expect_final_vectors_equivalence(runner_simple, runner);
        }           
    }
}

INSTANTIATE_TEST_SUITE_P(
    spectrum_final_equivalence,
    five_center_chain,
    ::testing::Combine(
    ::testing::ValuesIn(
        generate_all_optimization_lists({generate_S2_group_for_chain(5)})
    ),
    ::testing::Values(
        std::vector<spin_algebra::Multiplicity>({2, 2, 2, 2, 2}),
        std::vector<spin_algebra::Multiplicity>({2, 2, 4, 2, 2}),
        std::vector<spin_algebra::Multiplicity>({3, 3, 2, 3, 3}),
        std::vector<spin_algebra::Multiplicity>({3, 3, 4, 3, 3}),
        std::vector<spin_algebra::Multiplicity>({4, 2, 4, 2, 4})
    )
),
spectrum_final_equivalence_test_name_generator
);
